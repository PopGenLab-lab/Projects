import time
import click
import csv
import os

# --- Settings for temporary files ---
PARENT_TMP_DIR = "ref_tmp"
CHILD_TMP_DIR = "test_tmp"

def remove_dir_recursive(path):
    """
    Removes all files and directories recursively from given path.
    """
    if not os.path.exists(path):
        return

    for root, dirs, files in os.walk(path, topdown=False):
        for file in files:
            os.remove(os.path.join(root, file))
        for dir in dirs:
            os.rmdir(os.path.join(root, dir))
    os.rmdir(path)

def passes_filter(row):
    """
    Returns True if the row meets these criteria:
      - Start equals End.
      - Neither Ref nor Alt equals '-'.
    Assumes CSV columns: Chr, Start, End, Ref, Alt, Count.
    """
    try:
        start = int(row[1])
        end = int(row[2])
    except ValueError:
        return False  # ignore header or malformed rows
    return (row[1] == row[2]) and (row[3] != '-') and (row[4] != '-')


def group_and_count(input_filepath, tmp_dir):
    """
    Process an input file (parents or children):
      - Filter rows based on passes_filter.
      - Add the global total_count.
      - Write each record into a temporary CSV file grouped by chromosome.

    Files are written in tmp_dir, one file per chromosome.
    Each temporary file will have columns: Start, Ref, Alt, Count.
    Returns total_count across all chromosomes.
    """
    # Make sure the temporary directory exists.
    os.makedirs(tmp_dir, exist_ok=True)

    total_count = 0
    tmp_writers = {}
    tmp_files = {}

    with open(input_filepath, "r") as infile:
        reader = csv.reader(infile, delimiter="\t")
        header = next(reader, None)  # skip header
        for row in reader:
            if not passes_filter(row):
                continue
            try:
                chr_val = row[0]
                start = int(row[1])
                # End is not needed here (and equals start after filtering)
                ref = row[3]
                alt = row[4]
                count = int(row[13])
            except Exception as e:
                print('Row error:', e, '\nAt: chr ', row[0], ' start ', row[1], ' ref ', row[3], ' alt ', row[4])
                continue  # skip problematic rows

            total_count += count

            # Open a temporary writer for this chromosome if not already.
            if chr_val not in tmp_writers:
                tmp_filename = os.path.join(tmp_dir, f"{chr_val}.csv")
                tmp_file = open(tmp_filename, "w", newline="")
                writer = csv.writer(tmp_file)
                # Write header for later clarity (could be omitted if desired)
                writer.writerow(["Start", "Ref", "Alt", "Count"])
                tmp_writers[chr_val] = writer
                tmp_files[chr_val] = tmp_file

            # Write record (only needed fields).
            tmp_writers[chr_val].writerow([start, ref, alt, count])

    # Close all temporary files.
    for f in tmp_files.values():
        f.close()

    return total_count


def merge_chr_files(chr_val, parent_tmp_filepath, child_tmp_filepath, total_parent, total_child, out_dir):
    """
    For a given chromosome, open the parent's and child's temporary files,
    read them row by row (they are assumed to be sorted by Start)
    and compute the relative fitness per record using:

      parent_ratio = parent_count / total_parent
      child_ratio = child_count / total_child
      relative_fitness = child_ratio / parent_ratio

    Then, output a CSV file (named <chr_val>.csv) in out_dir with columns:
      Start, Ref, Alt, Relative_Fitness
    """
    parent_file = open(parent_tmp_filepath, "r")
    child_file = open(child_tmp_filepath, "r")
    parent_reader = csv.reader(parent_file)
    child_reader = csv.reader(child_file)

    # Skip header rows
    next(parent_reader, None)
    next(child_reader, None)

    # Output file (one per chromosome)
    os.makedirs(out_dir, exist_ok=True)
    out_filepath = os.path.join(out_dir, f"{chr_val}.csv")
    with open(out_filepath, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Start", "Ref", "Alt", "Relative_Fitness"])

        # Using the assumption that each key (Start, Ref, Alt) exists in both files,
        # We simply iterate record by record.
        for p_row, c_row in zip(parent_reader, child_reader):
            try:
                start = int(p_row[0])
                ref = p_row[1]
                alt = p_row[2]
                parent_count = int(p_row[3])

                # For child's file.
                child_count = int(c_row[3])
            except Exception as e:
                continue  # skip if conversion fails

            # Compute normalized ratios.
            parent_freq = parent_count / total_parent
            child_freq = child_count / total_child

            parent_freq = 1e-8 if parent_freq == 0 else parent_freq

            # Avoid division by zero: if parent_ratio is zero, use denovo rate
            relative_fitness = child_freq**2 / ((parent_freq**2 + parent_freq**3) - (parent_freq * child_freq**2))

            writer.writerow([start, ref, alt, relative_fitness])

    parent_file.close()
    child_file.close()


def merge_all_chromosomes(parent_tmp_dir, child_tmp_dir, total_parent, total_child, out_dir="output"):
    """
    For every chromosome that exists in both temporary directories,
    merge the parent's and child's records to compute relative fitness.
    """
    parent_chrs = set(os.listdir(parent_tmp_dir))
    child_chrs = set(os.listdir(child_tmp_dir))

    common_chrs = parent_chrs.intersection(child_chrs)

    if not common_chrs:
        print("No common chromosomes found in both files.")
        return

    for tmp_filename in common_chrs:
        # tmp_filename is something like "chr1.csv". Extract the chromosome label.
        chr_val = os.path.splitext(tmp_filename)[0]
        parent_tmp_path = os.path.join(parent_tmp_dir, tmp_filename)
        child_tmp_path = os.path.join(child_tmp_dir, tmp_filename)
        print(f"Merging records for chromosome {chr_val}...")
        merge_chr_files(chr_val, parent_tmp_path, child_tmp_path, total_parent, total_child, out_dir)
    print(f"All output files are saved in the '{out_dir}' directory.")



@click.command()
@click.option('--input-test', help='Test generation input file.')
@click.option('--input-ref', help='Reference generation input file.')
@click.option('--out-dir', default='results', help='Directory for output files.')
@click.option('--temp-dir', default='tmp', help='Directory for temporary files.')
@click.option('--generate-graphics', default=False, help='Generates graphics for each chromosome.')
@click.option('--keep-temp', default=False, help='Do not delete temporary files.')

def pyrelfit(input_test, input_ref, out_dir, temp_dir,generate_graphics, keep_temp):
    """
    Tool used to calculate relative fitness of SNP mutations.\n
    note: data must be (Chr, Pos) sorted
    """

    # get current time
    start_time = time.time()

    parent_tmp_dir = temp_dir + '/' + PARENT_TMP_DIR
    child_tmp_dir = temp_dir + '/' + CHILD_TMP_DIR

    if generate_graphics:
        print('Graphics not implemented.')

    print("Processing parent file...")
    total_parent = group_and_count(input_ref, parent_tmp_dir)
    print(f"Total parent count: {total_parent}")

    print("Processing child file...")
    total_child = group_and_count(input_test, child_tmp_dir)
    print(f"Total child count: {total_child}")

    print("Merging chromosome-specific files and computing relative fitness...")
    merge_all_chromosomes(parent_tmp_dir, child_tmp_dir, total_parent, total_child, out_dir)

    if not keep_temp:
        print("Deleting temporary files...")
        remove_dir_recursive(temp_dir)

    duration = time.time() - start_time
    print("Done in ", duration, " seconds.")
