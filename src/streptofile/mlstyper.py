#!/usr/bin/env python3

from pathlib import Path
import polars as pl
import subprocess
import argparse
from importlib import resources


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        nargs="*",
        help="Input fasta file(s)",
        type=Path,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output directory.",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-d",
        "--database",
        help="Path to folder with MLST data (alleles fasta file and profiles tsv)",
        type=Path,
        default=resources.files("streptofile") / "db" / "mlst",
    )
    parser.add_argument(
        "--full_path",
        help="Print full path to fasta input in output tsv rather than just file name",
        action="store_true",
        default=False,
    )
    return parser.parse_args()


def run_mlst_blast(
    assembly_file: Path,
    mlst_allele_db: Path,
    output_file: Path,
    ) -> bool | None:
    output_file = Path(output_file)
    if not output_file.exists():
        cmd = (
            f'blastn -query {assembly_file} -subject {mlst_allele_db} '
            f'-out {output_file} -max_target_seqs 100000 -ungapped -perc_identity 80 -word_size 50 '
            f'-outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send sseq evalue bitscore slen"'
        )
        result = subprocess.run(cmd, shell=True)
        if result.returncode != 0:
            return False
        else:
            return True
    else:
        return True


def load_mlst_profiles(profiles_tsv: Path) -> pl.DataFrame:
    """
    Load MLST profiles TSV into a data frame.
    The gene columns are inferred as all columns except ST and clonal_complex.
    """
    try:
        profiles_df = pl.read_csv(
            profiles_tsv,
            separator="\t",
            has_header=True,
        )
    except pl.exceptions.NoDataError:
        raise Exception(f"Could not open MLST profiles tsv at {profiles_tsv}")

    if "ST" not in profiles_df.columns:
        raise Exception(f"MLST profiles tsv at {profiles_tsv} must contain an ST column")

    profiles_df = profiles_df.with_columns(
        [pl.col(col).cast(pl.Utf8) for col in profiles_df.columns]
    )
    return profiles_df


def extract_mlst_type(
    mlst_blast_tsv: Path,
    profiles_tsv: Path,
) -> pl.DataFrame:
    """
    Load MLST BLAST results and return a one-row dataframe with ST,
    top hit for each allele, and notes.

    Assumes BLAST was run with:
      - assembly as query
      - alleles.fasta as subject/database

    Therefore:
      - qseqid = assembly contig
      - sseqid = MLST allele name, e.g. gki_4
    """
    mlst_blast_tsv = Path(mlst_blast_tsv)

    profiles_df = load_mlst_profiles(profiles_tsv)
    mlst_genes = [col for col in profiles_df.columns if col not in ["ST", "clonal_complex"]]

    default_result = {
        "ST": ["ND"],
        **{gene: ["-"] for gene in mlst_genes},
        "mlst_notes": [""],
    }

    if not mlst_blast_tsv.exists():
        default_result["mlst_notes"] = ["No blast output found for MLST alleles"]
        return pl.DataFrame(default_result)

    try:
        blast_df = pl.read_csv(
            mlst_blast_tsv,
            separator="\t",
            has_header=False,
            new_columns=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "qlen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "sseq",
                "evalue",
                "bitscore",
                "slen",
            ],
            schema_overrides={
                "qseqid": pl.Utf8,
                "sseqid": pl.Utf8,
                "pident": pl.Float64,
                "length": pl.Int64,
                "qlen": pl.Int64,
                "qstart": pl.Int64,
                "qend": pl.Int64,
                "sstart": pl.Int64,
                "send": pl.Int64,
                "sseq": pl.Utf8,
                "evalue": pl.Float64,
                "bitscore": pl.Float64,
                "slen": pl.Int64,
            },
        )
    except pl.exceptions.NoDataError:
        raise Exception(f"Could not open MLST blast output tsv at {mlst_blast_tsv}")

    if blast_df.height == 0:
        default_result["mlst_notes"] = ["No blast hits found for MLST alleles"]
        return pl.DataFrame(default_result)

    blast_df = blast_df.with_columns(
        [
            pl.col("sseqid").str.extract(r"^(.+?)_([^_]+)$", group_index=1).alias("gene"),
            pl.col("sseqid").str.extract(r"^(.+?)_([^_]+)$", group_index=2).alias("allele"),
            (pl.col("length") / pl.col("slen") * 100).alias("percent_coverage"),
            pl.min_horizontal("sstart", "send").alias("subject_start"),
            pl.max_horizontal("sstart", "send").alias("subject_end"),
        ]
    ).with_columns(
        [
            (
                (pl.col("pident") == 100.0)
                & (pl.col("length") == pl.col("slen"))
                & (pl.col("subject_start") == 1)
                & (pl.col("subject_end") == pl.col("slen"))
            ).alias("exact_match")
        ]
    ).filter(
        pl.col("gene").is_in(mlst_genes)
    )

    if blast_df.height == 0:
        default_result["mlst_notes"] = ["No recognizable MLST allele names found in BLAST output"]
        return pl.DataFrame(default_result)

    best_hits_df = (
        blast_df
        .sort(
            by=["gene", "exact_match", "bitscore", "pident", "percent_coverage", "length"],
            descending=[False, True, True, True, True, True],
        )
        .group_by("gene")
        .first()
        .with_columns(
            [
                pl.col("pident").round(2).alias("pident"),
                pl.col("percent_coverage").round(2).alias("percent_coverage"),
                pl.when(pl.col("exact_match"))
                .then(pl.col("allele"))
                .otherwise(pl.concat_str([pl.col("allele"), pl.lit("*")], separator=""))
                .alias("allele_call"),
            ]
        )
        .select(
            [
                "gene",
                "allele",
                "allele_call",
                "pident",
                "length",
                "slen",
                "exact_match",
            ]
        )
    )

    allele_calls = {gene: "-" for gene in mlst_genes}
    exact_alleles = {gene: "-" for gene in mlst_genes}
    notes = []

    for row in best_hits_df.iter_rows(named=True):
        gene = row["gene"]
        allele_calls[gene] = row["allele_call"]
        exact_alleles[gene] = row["allele"]
        if not row["exact_match"]:
            notes.append(
                f"{gene} best hit is allele {row['allele']} with pident {row['pident']} and length {row['length']}/{row['slen']}"
            )

    missing_genes = [gene for gene in mlst_genes if allele_calls[gene] == "-"]
    if missing_genes:
        notes.append("Missing MLST genes: " + ", ".join(missing_genes))

    st = "ND"
    if not missing_genes:
        exact_match_df = profiles_df.filter(
            pl.all_horizontal([pl.col(gene) == exact_alleles[gene] for gene in mlst_genes])
        )

        if exact_match_df.height > 0:
            st = exact_match_df.get_column("ST").to_list()[0]
            if any(value.endswith("*") for value in allele_calls.values() if value != "-"):
                st = f"{st}*"
        else:
            profiles_with_distance = profiles_df.with_columns(
                pl.sum_horizontal(
                    [
                        pl.when(pl.col(gene) == exact_alleles[gene]).then(0).otherwise(1)
                        for gene in mlst_genes
                    ]
                ).alias("n_mismatches")
            )

            slv_df = (
                profiles_with_distance
                .filter(pl.col("n_mismatches") == 1)
                .select("ST")
                .unique()
                .sort("ST")
            )

            if slv_df.height > 0:
                slv_sts = slv_df.get_column("ST").to_list()
                notes.append("Single locus variant of ST(s): " + ", ".join(slv_sts))
            else:
                notes.append("No exact ST match found")

    if len(notes) == 0:
        notes_text = "All loci matched exactly"
    else:
        notes_text = "; ".join(notes)

    result_df = pl.DataFrame(
        {
            "ST": [st],
            **{gene: [allele_calls[gene]] for gene in mlst_genes},
            "mlst_notes": [notes_text],
        }
    )
    return result_df

def type_sample(
    assembly_file: Path,
    output_folder: Path,
    mlst_allele_db: Path,
    profiles_tsv: Path,
) -> pl.DataFrame:
    """
    Wrapper for single sample.
    Runs blastn and returns dataframe of MLST results.
    """
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    blast_output_tsv = output_folder.joinpath("mlst_blast.tsv")
    blast_complete = run_mlst_blast(
        assembly_file=assembly_file,
        mlst_allele_db=mlst_allele_db,
        output_file=blast_output_tsv,
    )
    if blast_complete:
        mlst_results_df = extract_mlst_type(
            mlst_blast_tsv=blast_output_tsv,
            profiles_tsv=profiles_tsv,
        )
        return mlst_results_df
    else:
        profiles_df = load_mlst_profiles(profiles_tsv)
        mlst_genes = [col for col in profiles_df.columns if col != "ST"]
        return pl.DataFrame(
            {
                "ST": ["ND"],
                **{gene: ["-"] for gene in mlst_genes},
                "mlst_notes": ["Failed to run blast"],
            }
        )


def type_batch(
    assembly_files: list[Path],
    database_dir: Path,
    output_dir: Path,
    full_path: bool = False,
) -> pl.DataFrame:
    mlst_results = []
    mlst_allele_db = database_dir / "alleles.fasta"
    profiles_tsv = database_dir / "profiles.tsv"

    for assembly_file in assembly_files:
        sample_output_dir = output_dir / assembly_file.stem
        mlst_results_df = type_sample(
            assembly_file=assembly_file,
            output_folder=sample_output_dir,
            mlst_allele_db=mlst_allele_db,
            profiles_tsv=profiles_tsv,
        )

        if full_path:
            sample_name = str(assembly_file)
        else:
            sample_name = assembly_file.name

        mlst_results_df = mlst_results_df.with_columns(
            pl.lit(sample_name).alias("sample")
        ).select(
            ["sample"] + mlst_results_df.columns
        )
        mlst_results.append(mlst_results_df)

    combined_results = pl.concat(mlst_results)
    return combined_results


def main():
    args = parse_args()
    print(f"Running MLST profiler on {len(args.input)} samples")
    combined_results = type_batch(
        assembly_files=args.input,
        database_dir=args.database,
        output_dir=args.output,
        full_path=args.full_path,
    )
    output_file = args.output / "results.tsv"
    print(f"Printing results to {output_file}")
    combined_results.write_csv(file=output_file, separator="\t")
    print("Done")


if __name__ == "__main__":
    main()