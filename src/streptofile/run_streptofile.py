#!/usr/bin/env python3

from pathlib import Path
import polars as pl
import subprocess
import argparse
from importlib import resources
from streptofile import emm_typer
from streptofile import mlstyper
from streptofile import virulence_profiler


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input",
                        nargs="*",
                        help = "Input fasta file(s)",
                        type = Path)
    parser.add_argument("-o", "--output",
                        help = "Output directory.",
                        type=Path,
                        required = True)
    parser.add_argument("--analyses",
                        help = "Analyses to be performed in tab separated string. F.ex 'mlst,emm,virulence'. Default 'all.",
                        type=str,
                        default = "all")
    parser.add_argument("--emm_db",
                        help = "Path to fasta file with emm allele sequences",
                        type=Path,
                        default = resources.files("streptofile") / "db" / "emm_typing" / "alltrimmed.tfa")
    parser.add_argument("--mlst_db",
                        help="Path to folder with MLST data (alleles fasta file and profiles tsv)",
                        type=Path,
                        default=resources.files("streptofile") / "db" / "mlst")
    parser.add_argument("--virulence_db",
                        help = "Path to fasta file with emm allele sequences",
                        type=Path,
                        default = resources.files("streptofile") / "db" / "virulence_profiling")
    parser.add_argument("--full_path",
                        help = "Print full path to fasta input in output tsv rather than just file name",
                        action= "store_true",
                        default=False)
    return parser.parse_args()



def type_batch(assembly_files: list[Path],
               output_folder: Path,
               analyses_to_run: list,
               emm_allele_fasta: Path | None,
               mlst_database_dir: Path | None,
               virulence_database_dir: Path | None,
               full_path: bool = False,
               ) -> pl.DataFrame:
    result_dfs = []
    if "emm" in analyses_to_run:
        emm_results = emm_typer.type_batch(assembly_files=assembly_files,
                                           emm_allele_fasta=emm_allele_fasta,
                                           output_dir=output_folder)
        result_dfs.append(emm_results)
    if "mlst" in analyses_to_run:
        mlst_results = mlstyper.type_batch(assembly_files=assembly_files,
                                            output_dir=output_folder,
                                            database_dir=mlst_database_dir)
        result_dfs.append(mlst_results)
    if "virulence" in analyses_to_run:
        virulence_results, virulence_presence_absence = virulence_profiler.profile_batch(assembly_files=assembly_files,
                                                                                         database_dir=virulence_database_dir,
                                                                                         output_dir=output_folder)
        result_dfs.append(virulence_presence_absence)
    dfs_fixed = [result_dfs[0]] + [df.drop("sample") for df in result_dfs[1:]]
    combined_results = pl.concat(dfs_fixed, how="horizontal")
    return(combined_results)

def main():
    args = parse_args()
    print(f"Running streptofile on {len(args.input)} samples")
    args.analyses = args.analyses.lower()
    if args.analyses.lower == "all":
        args.analyses = "emm,mlst,virulence"
    print(f"Including analyses: {args.analyses}")
    analyses_to_run = args.analyses.split(",")
    results_df = type_batch(assembly_files=args.input,
                            output_folder=args.output,
                            analyses_to_run=analyses_to_run,
                            emm_allele_fasta=args.emm_db,
                            mlst_database_dir=args.mlst_db,
                            virulence_database_dir=args.virulence_db,
                            full_path=args.full_path)
    output_file = args.output / "results.tsv"
    print(f"Printing results to {output_file}")
    results_df.write_csv(file = output_file,
                         separator = "\t")
    print("Done")

if __name__ == "__main__":
    main()


