import glob
import time
import click
import csv
import os
import matplotlib.pyplot as plt
from cyvcf2 import VCF
from multiprocessing import Pool
from collections import defaultdict

# --- Settings for temporary files ---
GRAPHICS_OPACITY = 0.6
GRAPHICS_POINT_SIZE = 4
GRAPHICS_ENABLED = False

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

def parse_generations(samples, generations):
    groups = defaultdict(list)
    for sample in samples:
        try:
            gen = sample.split("-")[-1]
            if gen in generations:
                groups[gen].append(sample)
        except IndexError:
            pass
    return dict(groups)


def compute_counts(variant):
    total = variant.num_called * 2
    alt_count = variant.num_het + variant.num_hom_alt * 2

    return alt_count, total


def filter_and_split(vcf_path, generations, temp_dir):
    vcf = VCF(vcf_path)
    gens = parse_generations(vcf.samples, generations)

    for gen in generations:
        vcf = VCF(vcf_path, samples=gens.get(gen))
        for variant in vcf:
            if not variant.is_snp: continue

            chrom = variant.CHROM
            pos = variant.POS
            ref = variant.REF
            alt = variant.ALT[0] if variant.ALT else "."

            alt_count, total = compute_counts(variant)
            with open(f"{temp_dir}/tmp.{chrom}.{gen}.csv", "a") as f:
                csv.writer(f).writerow([pos, ref, alt, alt_count, total])


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

                    w = (f2 ** 2) / (2 * f1 ** 2 - f1 * f2 ** 2) if f1 != 0 else 0
                    max_w = max(max_w, w)

                    writer.writerow([pos, ref, alt, w])

            return max_w, pair


def merge_and_compute(generation_pairs, cores, temp_dir, out_dir):
    chromosomes = {f.split(".")[1] for f in glob.glob(f"{temp_dir}/tmp.*.{generation_pairs[0].split('_')[0]}.csv")}
    tasks = [(c, p, temp_dir, out_dir) for c in chromosomes for p in generation_pairs]

    with Pool(processes=cores) as pool:
        max_vals = pool.map(process_pair, tasks)
    # MAX should be per generation_pair
    return list_to_dict_max(max_vals)


def normalize_file(filename, global_max):
    if global_max == 0 or global_max == 1: return

    with open(filename) as f_in:
        with open(f"{filename}.tmp", "w") as f_out:
            reader = csv.reader(f_in)
            writer = csv.writer(f_out)
            header = next(reader)
            writer.writerow(header)

            x_vals = []
            y_vals = []

            for row in reader:
                row[-1] = float(row[-1]) / global_max
                if GRAPHICS_ENABLED:
                    x_vals.append(float(row[0]))
                    y_vals.append(row[-1])

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


def normalise(scale_list, generation_pairs, cores, out_dir):
    files = []

    for pair in generation_pairs:
        files.extend([(s,pair) for s in glob.glob(f"{out_dir}/*.{pair}.csv")])

    with Pool(processes=cores) as pool:
        pool.starmap(normalize_file, [(f, scale_list.get(p)) for f, p in files])


@click.command()
@click.option('-i', '--input-file', help='Input VCF file.')
@click.option('-c', '--cores', default=1, help='Number of processes to spawn.')
@click.option('-g', '--generations', type=list, default=['1', '2', '3'], help='Sample generations.')
@click.option('-p', '--generation-pairs', type=list, default=['1_3', '2_3'], help='Sample generation pairs.')
@click.option('-o', '--out-dir', default='results', help='Directory for output files.')
@click.option('-t', '--temp-dir', default='tmp', help='Directory for temporary files.')
@click.option('-G', '--generate-graphics', default=False, help='Generates graphics for each chromosome.')
@click.option('--keep-temp', default=False, help='Do not delete temporary files.')
def pyrelfit(input_file, generations, generation_pairs, cores, temp_dir, out_dir, keep_temp, generate_graphics):
    """
    Tool used to calculate relative fitness of SNP mutations.\n
    """
    global GRAPHICS_ENABLED
    GRAPHICS_ENABLED = generate_graphics

    # get current time
    start_time = time.time()


    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    filter_and_split(input_file, generations, temp_dir)
    scale_list = merge_and_compute(generation_pairs, cores, temp_dir, out_dir)
    normalise(scale_list, generation_pairs, cores, out_dir)

    if not keep_temp:
        print("Deleting temporary files...")
        remove_dir_recursive(temp_dir)

    duration = time.time() - start_time
    print("Done in ", duration, " seconds.")
