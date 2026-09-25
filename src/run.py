#!/usr/bin/env python3
import argparse
import os
import sys
import time
import traceback
from pathlib import Path
from pipeline import annotate_sequence

VALID_FASTA = {".fasta", ".fa", ".fna"}
VALID_GB = {".gb", ".gbk", ".genbank"}
VALID_EXTENSIONS = VALID_FASTA | VALID_GB
ORGANELLE_NAMES = {"chloroplast", "plastid", "mitochondrion", "mitochondria"}

# NCBI translation tables (7, 8 and 17-20 are unassigned)
GENETIC_CODES = {
    "1", "2", "3", "4", "5", "6", "9", "10", "11", "12", "13", "14", "15", "16",
    *map(str, range(21, 34)),
}


def organelle_arg(value):
    """Checks --organelle: an organelle name or an NCBI genetic code number."""
    cleaned = value.strip().lower()
    if cleaned in ORGANELLE_NAMES or cleaned in GENETIC_CODES:
        return cleaned
    raise argparse.ArgumentTypeError(
        f"'{value}' is not an organelle name ({', '.join(sorted(ORGANELLE_NAMES))}) "
        f"or a valid NCBI genetic code (1-6, 9-16, 21-33)"
    )


def main():
    # Force line buffering so Docker prints output immediately in real-time
    try:
        sys.stdout.reconfigure(line_buffering=True)
        sys.stderr.reconfigure(line_buffering=True)
    except AttributeError:
        pass

    parser = argparse.ArgumentParser(
        description="Annotate seaweed organelle genomes using MFannot and Aragorn."
    )

    # 1. Positional Target
    parser.add_argument(
        "target",
        help="Input sequence file (.fasta/.gb) OR directory for batch mode",
    )

    # 2. Shared Options (apply to both single-file and batch mode)
    shared_group = parser.add_argument_group("shared options")
    shared_group.add_argument(
        "--organelle",
        type=organelle_arg,
        help=(
            "Organelle type or NCBI genetic code number: chloroplast/plastid use code 11, "
            "mitochondrion/mitochondria code 4, or give the code directly (e.g. 1). "
            "Required for single files; in batch mode, applies to all files "
            "(cannot be combined with --plastid_suffix or --mito_suffix)"
        ),
    )
    shared_group.add_argument(
        "--format",
        default="gb",
        choices=["gb", "gff3"],
        help="Output format (default: gb)",
    )
    shared_group.add_argument(
        "--directory",
        default="/data",
        help="Base directory containing the input file/folder (default: /data)",
    )
    shared_group.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        default=None,
        help="Directory where output files will be written (default: same as input)",
    )
    shared_group.add_argument(
        "--save-log",
        action="store_true",
        help="Save execution log to file even if all annotations succeed",
    )
    shared_group.add_argument(
        "--circular",
        dest="circular",
        action="store_true",
        default=None,
        help="Treat genome(s) as circular and attempt rotation to start at rns",
    )
    shared_group.add_argument(
        "--linear",
        dest="circular",
        action="store_false",
        help="Treat genome(s) as linear",
    )
    shared_group.add_argument(
        "--no-keep-old",
        dest="keep_old",
        action="store_false",
        help="Discard pre-rotation sequence file instead of keeping it",
    )
    shared_group.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Also write the MFannot table (.tbl) and Aragorn output (.txt)",
    )

    # 3. Batch Mode Suffixes
    batch_group = parser.add_argument_group(
        "batch mode options (used when target is a directory)"
    )
    batch_group.add_argument(
        "--plastid_suffix",
        default=None,
        help="Filename suffix for chloroplast (e.g. '_pt'); non-matching files become mitochondrion",
    )
    batch_group.add_argument(
        "--mito_suffix",
        default=None,
        help="Filename suffix for mitochondrion (e.g. '_mt'); non-matching files become chloroplast",
    )
    batch_group.add_argument(
        "--circular_suffix",
        default=None,
        help="Filename suffix for circular files (e.g. '_circ'); non-matching files default to unspecified",
    )
    batch_group.add_argument(
        "--linear_suffix",
        default=None,
        help="Filename suffix for linear files (e.g. '_lin'); non-matching files default to unspecified",
    )

    args = parser.parse_args()

    # Resolve target path (check directly or inside --directory)
    if os.path.exists(args.target):
        resolved_target = os.path.abspath(args.target)
    elif os.path.exists(os.path.join(args.directory, args.target)):
        resolved_target = os.path.abspath(os.path.join(args.directory, args.target))
    else:
        parser.error(
            f"Target '{args.target}' not found (checked current path and inside '{args.directory}')."
        )

    is_batch = os.path.isdir(resolved_target)

    # Validate suffix constraint: maximum 1 suffix allowed across the command
    suffix_args = {
        "--plastid_suffix": args.plastid_suffix,
        "--mito_suffix": args.mito_suffix,
        "--circular_suffix": args.circular_suffix,
        "--linear_suffix": args.linear_suffix,
    }
    active_suffixes = [k for k, v in suffix_args.items() if v is not None]
    if len(active_suffixes) > 1:
        parser.error(
            f"You cannot combine multiple suffix arguments together ({', '.join(active_suffixes)}). "
            f"Only 1 suffix argument is allowed."
        )

    # ==========================
    # BATCH MODE (Target is a folder)
    # ==========================
    if is_batch:
        batch_dir = resolved_target
        out_dir = os.path.abspath(args.output_dir) if args.output_dir else batch_dir
        os.makedirs(out_dir, exist_ok=True)

        if not args.organelle and not args.plastid_suffix and not args.mito_suffix:
            parser.error(
                "Batch mode requires either --organelle, --plastid_suffix, or --mito_suffix."
            )

        # Organelle suffixes decide the organelle themselves, so --organelle would be ignored
        if args.organelle and (args.plastid_suffix or args.mito_suffix):
            parser.error(
                "--organelle cannot be combined with --plastid_suffix or --mito_suffix."
            )

        # Collect matching sequence files
        files = [
            f
            for f in sorted(os.listdir(batch_dir))
            if os.path.isfile(os.path.join(batch_dir, f))
            and Path(f).suffix.lower() in VALID_EXTENSIONS
        ]

        if not files:
            print(f"[WARNING] No FASTA or GenBank sequence files found in {batch_dir}", flush=True)
            return

        batch_start_time = time.time()

        # Print shared parameters ONCE at the top
        batch_header = (
            f"\n{'=' * 60}\n"
            f"SOGA BATCH ANNOTATION\n"
            f"  Started:     {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            f"  Directory:   {batch_dir}\n"
            f"  Output Dir:  {out_dir}\n"
            f"  Settings:    format={args.format}, keep_old={args.keep_old}, intermediates={args.keep_intermediates}\n"
            f"  Total Files: {len(files)} sequence file(s) found\n"
            f"{'=' * 60}\n"
        )
        print(batch_header, flush=True)

        batch_log_buffer = [batch_header]
        succeeded = []
        failed = []

        for idx, fname in enumerate(files, 1):
            stem = Path(fname).stem

            # Determine organelle
            if args.plastid_suffix:
                is_plastid = stem.endswith(args.plastid_suffix) or fname.endswith(args.plastid_suffix)
                organelle = "chloroplast" if is_plastid else "mitochondrion"
            elif args.mito_suffix:
                is_mito = stem.endswith(args.mito_suffix) or fname.endswith(args.mito_suffix)
                organelle = "mitochondrion" if is_mito else "chloroplast"
            else:
                organelle = args.organelle

            # Determine circularity
            if args.circular_suffix:
                is_circ = stem.endswith(args.circular_suffix) or fname.endswith(args.circular_suffix)
                circular = True if is_circ else None
            elif args.linear_suffix:
                is_lin = stem.endswith(args.linear_suffix) or fname.endswith(args.linear_suffix)
                circular = False if is_lin else None
            else:
                circular = args.circular

            file_header = f"[{idx}/{len(files)}] Processing '{fname}' (organelle={organelle}, circular={circular})"
            divider = "-" * len(file_header)
            print(f"{file_header}\n{divider}", flush=True)
            batch_log_buffer.append(f"\n{file_header}\n{divider}")

            try:
                # Runs sequence and collects buffered log
                file_log = annotate_sequence(
                    fname,
                    format=args.format,
                    directory=batch_dir,
                    organelle=organelle,
                    circular=circular,
                    keep_old=args.keep_old,
                    keep_intermediates=args.keep_intermediates,
                    output_dir=out_dir,
                    save_log=False,
                    return_log=True,
                )
                if file_log:
                    batch_log_buffer.append(file_log)
                succeeded.append(fname)
                print("", flush=True)
            except Exception as e:
                err_msg = str(e)
                tb = traceback.format_exc()
                print(f"[ERROR] Failed to annotate '{fname}': {err_msg}\n", flush=True)
                batch_log_buffer.append(f"[ERROR] '{fname}' failed:\n{tb}")
                failed.append((fname, err_msg))
                continue

        # Total elapsed time
        total_time = time.time() - batch_start_time
        tm, ts = divmod(total_time, 60)

        # Write unified log file if requested or if errors occurred
        log_file_path = Path(out_dir) / "batch_annotation.log"
        has_errors = len(failed) > 0

        if args.save_log or has_errors:
            log_file_path.write_text("\n\n".join(batch_log_buffer) + "\n")

        # Print Batch Summary Block
        print("=" * 60, flush=True)
        print("BATCH SUMMARY", flush=True)
        print(f"Total files:   {len(files)}", flush=True)
        print(f"Succeeded:     {len(succeeded)}", flush=True)
        print(f"Failed:        {len(failed)}", flush=True)
        print(f"Total time:    {int(tm)} min {ts:.1f} s", flush=True)
        print("-" * 60, flush=True)

        if not has_errors:
            print("All files were successfully annotated!", flush=True)
            if args.save_log:
                print(f"Log file saved to: {log_file_path}", flush=True)
        else:
            print("Failed files:", flush=True)
            for failed_file, err in failed:
                print(f"  - {failed_file} ({err})", flush=True)
            print(f"\nA detailed error log has been saved to:\n{log_file_path}", flush=True)

        print("=" * 60 + "\n", flush=True)

        # Non-zero exit so scripts can tell the batch had failures
        if has_errors:
            sys.exit(1)

    # ==========================
    # SINGLE FILE MODE (Target is a file)
    # ==========================
    else:
        if active_suffixes:
            parser.error("Suffix arguments can only be used in batch mode on a directory.")

        if not args.organelle:
            parser.error("the following arguments are required: --organelle")

        file_dir = os.path.dirname(resolved_target)
        file_name = os.path.basename(resolved_target)
        out_dir = os.path.abspath(args.output_dir) if args.output_dir else file_dir

        annotate_sequence(
            file_name,
            format=args.format,
            directory=file_dir,
            organelle=args.organelle,
            circular=args.circular,
            keep_old=args.keep_old,
            keep_intermediates=args.keep_intermediates,
            output_dir=out_dir,
            save_log=args.save_log,
            return_log=False,
        )


if __name__ == "__main__":
    main()