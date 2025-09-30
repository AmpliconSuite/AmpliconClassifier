import argparse
import csv
import glob
from pathlib import Path
from pyliftover import LiftOver

def liftover_pos(lo, chrom, pos, file_type="tsv"):
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    res = lo.convert_coordinate(chrom, pos - 1)  #convert to 0-based
    if not res:
        return None, None
    #pick first mapping if multiple (could refine later)
    new_chrom, new_pos, strand, _ = res[0]
    if file_type == "bed":
        return new_chrom.replace("chr", ""), int(new_pos)  #keep 0-based for bed start
    else:
        return new_chrom.replace("chr", ""), int(new_pos) + 1  #back to 1-based

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, help="Top-level directory to search for *_amplicon*_SV_summary.tsv")
    parser.add_argument("--out-dir", required=True, help="Output directory (mirrors input structure)")
    parser.add_argument("--source", default="hg19", help="Source genome (default: hg19)")
    parser.add_argument("--target", default="hg38", help="Target genome (default: hg38)")
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.out_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    #LiftOver object
    lo = LiftOver(args.source, args.target)

    #find *_amplicon*_SV_summary.tsv files recursively
    tsv_files = list(input_dir.rglob("*_amplicon*_SV_summary.tsv"))
    if not tsv_files:
        raise SystemExit(f"No *_amplicon*_SV_summary.tsv files found under {input_dir}")

    print(f"Found {len(tsv_files)} files to process.")

    #find bedfiles recursively
    bed_files = list(input_dir.rglob("*.bed"))
    if not bed_files:
        raise SystemExit(f"No *.bed files found under {input_dir}")

    print(f"Found {len(bed_files)} files to process.")

    # process tsv files
    for tsv_path in tsv_files:
        rel_path = tsv_path.relative_to(input_dir)
        out_path = output_dir / rel_path
        out_path.parent.mkdir(parents=True, exist_ok=True)

        print(f"Processing {tsv_path}")

        #load SV summary tsv
        with open(tsv_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        if not rows:
            print(f"Warning: empty file {tsv_path}, skipping.")
            continue

        #check required columns
        required_cols = {"chrom1","pos1","chrom2","pos2","pos1_flanking_coordinate","pos2_flanking_coordinate"}
        if not required_cols.issubset(reader.fieldnames):
            print(f"Error: Missing required columns in {tsv_path}, found {reader.fieldnames}")
            continue

        with open(out_path, "w", newline="") as out_f:
            writer = csv.DictWriter(out_f, fieldnames=reader.fieldnames, delimiter="\t")
            writer.writeheader()

            for row in rows:
                #pos1
                new_chrom1, new_pos1 = liftover_pos(lo, row["chrom1"], int(row["pos1"]))
                if new_chrom1:
                    row["chrom1"] = new_chrom1
                    row["pos1"] = new_pos1
                else:
                    row["chrom1"] = "NA"
                    row["pos1"] = "NA"

                #pos2
                new_chrom2, new_pos2 = liftover_pos(lo, row["chrom2"], int(row["pos2"]))
                if new_chrom2:
                    row["chrom2"] = new_chrom2
                    row["pos2"] = new_pos2
                else:
                    row["chrom2"] = "NA"
                    row["pos2"] = "NA"

                #flanking positions
                new_chrom1f, new_pos1f = liftover_pos(lo, row["chrom1"], int(row["pos1_flanking_coordinate"]))
                if new_pos1f:
                    row["pos1_flanking_coordinate"] = new_pos1f
                else:
                    row["pos1_flanking_coordinate"] = "NA"

                new_chrom2f, new_pos2f = liftover_pos(lo, row["chrom2"], int(row["pos2_flanking_coordinate"]))
                if new_pos2f:
                    row["pos2_flanking_coordinate"] = new_pos2f
                else:
                    row["pos2_flanking_coordinate"] = "NA"
                if "NA" not in (row["chrom1"], row["pos1"], row["chrom2"], row["pos2"], row["pos1_flanking_coordinate"], row["pos2_flanking_coordinate"]):
                    writer.writerow(row)
                else:
                    print(f"Skipping row due to unmapped coordinates: {row}")

        print(f"Saved lifted tsvs file to {out_path}")

    
    #process bed files
    for bed_path in bed_files:
        rel_path = bed_path.relative_to(input_dir)
        out_path = output_dir / rel_path
        out_path.parent.mkdir(parents=True, exist_ok=True)

        print(f"Processing {bed_path}")

        with open(bed_path) as fh:
            lines = fh.readlines()
        if not lines: #skip empty files
            print(f"Warning: empty file {bed_path}, skipping.")
            continue

        with open(out_path, "w") as out_f:
            for line in lines:
                parts = line.strip().split("\t")
                if len(parts) < 3: #not enough columns for BED
                    print(f"Warning: malformed BED line in {bed_path}: {line.strip()}, skipping.")
                    continue
                chrom = parts[0]
                start = int(parts[1]) + 1  #BED is 0-based, convert to 1-based
                end = int(parts[2])         #end is exclusive in BED, keep as is for liftover

                new_chrom, new_start = liftover_pos(lo, chrom, start, "bed")
                new_chrom, new_end = liftover_pos(lo, chrom, end) #end is 1-based

                if new_start and new_end:  #if both positions mapped
                    output_fields = [new_chrom, str(new_start), str(new_end)]
                    
                    if len(parts) > 3: #preserve any additional columns
                        output_fields.extend(parts[3:])
                    
                    out_f.write('\t'.join(output_fields) + '\n')
                else:
                    print(f"Skipping BED line due to unmapped coordinates: {line.strip()}")

        print(f"Saved lifted bed file to {out_path}")

    print("All files processed successfully.")

if __name__ == "__main__":
    main()