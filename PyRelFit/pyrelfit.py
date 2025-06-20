import glob
import sys
import time
import click
import csv
import os
import matplotlib.pyplot as plt
from cyvcf2 import VCF
from multiprocessing import Pool
import re

GRAPHICS_OPACITY = 0.8
GRAPHICS_POINT_SIZE = 4
GRAPHICS_ENABLED = False
OUTLIERS_ENABLED = False

def comma_separated(ctx, param, value):
    return value.split(',')

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

def list_to_dict_max(lst):
    max_dict = {}
    for value, key in lst:
        max_dict[key] = max(max_dict.get(key, -float('inf')), value)
    return max_dict

def parse_generations(samples, generation_string):
    if not generation_string.startswith("/") or not generation_string.endswith("/"):
        raise ValueError("Generations must be formatted as /<id>/<regex>/<id>/<regex>/...")

    parts = generation_string.strip("/").split("/")
    if len(parts) % 2 != 0:
        raise ValueError("Invalid generation format: must be /<id>/<regex>/<id>/<regex>/...")

    gen_sample_map = {}
    for i in range(0, len(parts), 2):
        gen_id = parts[i]
        regex = re.compile(parts[i + 1])
        matching_samples = [sample for sample in samples if regex.search(sample)]
        gen_sample_map[gen_id] = matching_samples

    return gen_sample_map


def compute_counts(variant):
    total = variant.num_called * 2
    alt_count = variant.num_het + variant.num_hom_alt * 2

    return alt_count, total

def filter_split_unit(args):
    vcf_path, chrom, gen, samples, temp_dir = args
    try:
        vcf = VCF(vcf_path, samples=samples)

        with open(f"{temp_dir}/tmp.{chrom}.{gen}.csv", "a") as f:
            for variant in vcf(chrom):
                if not variant.is_snp: continue
                if not variant.ALT: continue

                pos = variant.POS
                ref = variant.REF
                alt = variant.ALT[0]

                alt_count, total = compute_counts(variant)

                csv.writer(f).writerow([pos, ref, alt, alt_count, total])
    except Exception as e:
        print(e, file=sys.stderr)


def filter_and_split(vcf_path, generations, temp_dir, cores, chromosomes=None):
    vcf = VCF(vcf_path, threads=cores)

    if chromosomes:
        chrom_pattern = re.compile(chromosomes)
        chromosomes = [c for c in vcf.seqnames if chrom_pattern.search(c)]
    else:
        chromosomes = vcf.seqnames
    gens = parse_generations(vcf.samples, generations)

    tasks = [(vcf_path, c, gen, gens[gen], temp_dir) for c in chromosomes for gen in gens.keys()]
    with Pool(processes=cores) as pool:
        pool.map(filter_split_unit, tasks, chunksize=1)


def process_pair(args):
    chrom, pair, temp_dir, out_dir = args
    gen1, gen2 = pair.split("_")
    max_w = 0.0

    with open(f"{temp_dir}/tmp.{chrom}.{gen1}.csv") as f1:
        with open(f"{temp_dir}/tmp.{chrom}.{gen2}.csv") as f2:
            with open(f"{out_dir}/{chrom}.{pair}.csv", "w") as out:
                writer = csv.writer(out)
                writer.writerow(["Pos", "Ref", "Alt", "RF"])

                for (r1, r2) in zip(csv.reader(f1), csv.reader(f2)):
                    pos, ref, alt = r1[:3]
                    ac1, total1 = map(int, r1[3:])
                    ac2, total2 = map(int, r2[3:])

                    f1 = ac1 / total1 if total1 else 0
                    f2 = ac2 / total2 if total2 and total2 != 0 else 1e-8

                    w = (f2 ** 2) / (2 * f1 ** 2 - f1 * f2 ** 2) if f1 and f1 != 0.0 and f1 != 0 else 0
                    max_w = max(max_w, w)

                    writer.writerow([pos, ref, alt, w])

            return max_w, pair


def merge_and_compute(generation_pairs, cores, temp_dir, out_dir):
    chromosomes = {f.split(".")[1] for f in glob.glob(f"{temp_dir}/tmp.*.{generation_pairs[0].split('_')[0]}.csv")}
    tasks = [(c, p, temp_dir, out_dir) for c in chromosomes for p in generation_pairs]

    with Pool(processes=cores) as pool:
        max_vals = pool.map(process_pair, tasks, chunksize=1)
    # MAX should be per generation_pair
    return list_to_dict_max(max_vals)


def normalize_file(args):
    filename, global_max = args
    if global_max == 0 or global_max == 1: return

    with open(filename) as f_in:
        with open(f"{filename}.tmp", "w") as f_out:
            reader = csv.reader(f_in)
            writer = csv.writer(f_out)
            header = next(reader)
            writer.writerow(header)

            x_vals = []
            y_vals = []

            if OUTLIERS_ENABLED:
                outlier_groups = {
                    "gt_8": [],
                    "gt_6": [],
                    "gt_4": [],
                    "gt_2": []
                }

            for row in reader:
                row[-1] = float(row[-1]) / global_max
                if GRAPHICS_ENABLED:
                    x_vals.append(float(row[0]))
                    y_vals.append(row[-1])

                if OUTLIERS_ENABLED:
                    out_entry = row.copy()
                    if float(row[-1]) > 0.8:
                        outlier_groups["gt_8"].append(out_entry)
                    elif float(row[-1]) > 0.6:
                        outlier_groups["gt_6"].append(out_entry)
                    elif float(row[-1]) > 0.4:
                        outlier_groups["gt_4"].append(out_entry)
                    elif float(row[-1]) > 0.2:
                        outlier_groups["gt_2"].append(out_entry)

                writer.writerow(row)

            if GRAPHICS_ENABLED:
                png_filename = os.path.splitext(filename)[0]

                plt.figure(figsize=(16, 9))
                plt.scatter(x_vals, y_vals, alpha=GRAPHICS_OPACITY, s=GRAPHICS_POINT_SIZE)
                plt.title("Normalized Data Plot")
                plt.xlabel("Position in Chr")
                plt.ylabel("Relative fitness")
                plt.ylim(0, 1)
                plt.savefig(f"{png_filename}.png")
                plt.close()

        os.replace(f"{filename}.tmp", filename)

    if OUTLIERS_ENABLED:
        base_name = os.path.basename(filename)
        out_dir = os.path.dirname(filename)

        for group, entries in outlier_groups.items():
            group_file = os.path.join(out_dir, f"{group}.{base_name}")
            with open(group_file, "w", newline="") as gf:
                group_writer = csv.writer(gf)
                group_writer.writerow(header)
                group_writer.writerows(entries)


def normalise(scale_list, generation_pairs, cores, out_dir):
    files = []

    for pair in generation_pairs:
        files.extend([(s,pair) for s in glob.glob(f"{out_dir}/*.{pair}.csv")])

    with Pool(processes=cores) as pool:
        pool.map(normalize_file, [(f, scale_list.get(p)) for f, p in files], chunksize=1)
        pool.close()


@click.command()
@click.option('-i', '--input-file', help='Input VCF file.')
@click.option('-c', '--cores', default=1, help='Number of processes to spawn.')
@click.option('-C', '--chromosomes', help='Regex pattern to match chromosome names (e.g., "^chr[0-9]+$").')
@click.option('-g', '--generations', required=True, help='Sample generations in format: /<id>/<regex>/<id>/<regex>/...')
@click.option('-p', '--generation-pairs', callback=comma_separated, default='1_3,2_3', help='Sample generation pairs by id (e.g., "1_3,2_3".')
@click.option('-o', '--out-dir', default='results', help='Directory for output files.')
@click.option('-t', '--temp-dir', default='tmp', help='Directory for temporary files.')
@click.option('-O', '--outliers', default=False, help='Generate files with outlier relative fitness groups.')
@click.option('-G', '--generate-graphics', default=False, help='Generates graphics for each chromosome.')
@click.option('--keep-temp', default=False, help='Do not delete temporary files.')
def pyrelfit(input_file, generations, generation_pairs, cores, temp_dir, out_dir, keep_temp, generate_graphics, chromosomes, outliers):
    """
    Tool used to calculate relative fitness of SNP mutations.\n
    """
    global GRAPHICS_ENABLED, OUTLIERS_ENABLED
    GRAPHICS_ENABLED = generate_graphics
    OUTLIERS_ENABLED = outliers

    # get current time
    start_time = time.time()

    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    filter_and_split(input_file, generations, temp_dir, cores, chromosomes)
    print("Step 1 in: ", time.time() - start_time, " seconds.")
    scale_list = merge_and_compute(generation_pairs, cores, temp_dir, out_dir)
    print("Step 2 in: ", time.time() - start_time, " seconds.")
    normalise(scale_list, generation_pairs, cores, out_dir)

    if not keep_temp:
        print("Deleting temporary files...")
        remove_dir_recursive(temp_dir)

    print("Done in ", time.time() - start_time, " seconds.")


if __name__ == "__main__":
    sys.exit(pyrelfit())