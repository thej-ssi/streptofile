#!/usr/bin/env python3

from pathlib import Path
import polars as pl
import subprocess
import argparse
from importlib import resources


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
    parser.add_argument("-d", "--database",
                        help = "Path to fasta file with emm allele sequences",
                        type=Path,
                        default = resources.files("streptofile") / "db" / "emm_typing" / "alltrimmed.tfa")
    parser.add_argument("--full_path",
                        help = "Print full path to fasta input in output tsv rather than just file name",
                        action= "store_true",
                        default=False)
    return parser.parse_args()



def run_emm_blast(assembly_file: Path, emm_allele_fasta: Path, output_file: Path) -> tuple[str] | None:
    output_file = Path(output_file)
    if not output_file.exists():
        cmd = f'blastn -query {emm_allele_fasta} -subject {assembly_file} -qcov_hsp_perc 90 -out {output_file} -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send sseq evalue bitscore"'
        result = subprocess.run(cmd, shell=True)
        if result.returncode != 0:
            return False
        else:
            return True
    else:
        return True


def extract_emm_type(emm_blast_tsv: Path) -> pl.DataFrame:
    emm_types_in_emm_plus_mrp_operons = []  # to update
    mrp_types_in_emm_plus_mrp_operons = [
        "134",
        "156",
        "159",
        "164",
        "174",
        "205",
    ]  # to update

    emm_blast_tsv = Path(emm_blast_tsv)

    def result_df(
        emm_type: str = "-",
        enn_type: str = "-",
        mrp_type: str = "-",
        genes_in_operon: str = "-",
        emm_typing_notes: str = "",
    ) -> pl.DataFrame:
        return pl.DataFrame(
            {
                "EMM_type": [emm_type],
                "ENN_type": [enn_type],
                "MRP_type": [mrp_type],
                "genes_in_operon": [genes_in_operon],
                "emm_typing_notes": [emm_typing_notes],
            }
        )

    if not emm_blast_tsv.exists():
        return result_df(emm_typing_notes="No blast output found for EMM genes")

    try:
        blast_df = pl.read_csv(
            emm_blast_tsv,
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
            ],
        )
    except pl.exceptions.NoDataError:
        return result_df(emm_typing_notes="Empty blast output, no EMM genes detected")

    notes = []

    blast_df = (
        blast_df.with_columns(
            [
                (pl.col("length") / pl.col("qlen") * 100).alias("plen"),
                pl.when(pl.col("sstart") < pl.col("send"))
                .then(((pl.col("sstart") - pl.col("qstart") + 1) / 100).floor())
                .otherwise(((pl.col("send") + pl.col("qend") - 180) / 100).floor())
                .cast(pl.Int64)
                .alias("extended_sstart"),
                pl.col("qseqid").str.replace(r"^EMM", "").alias("allele"),
                pl.when((pl.col("length") < pl.col("qlen")) | (pl.col("pident") < 100))
                .then(pl.col("qseqid").str.replace(r"^EMM", "") + pl.lit("*"))
                .otherwise(pl.col("qseqid").str.replace(r"^EMM", ""))
                .alias("typed_allele"),
                pl.col("qseqid")
                .str.replace(r"^EMM", "")
                .str.split(".")
                .list.get(0)
                .alias("main_type"),
            ]
        )
        .filter(pl.col("bitscore") > 280)
    )

    blast_df_unique = (
        blast_df.sort("bitscore", descending=True)
        .group_by("extended_sstart")
        .first()
    )

    if blast_df_unique.height == 0:
        notes.append("No blast hits found for EMM genes")
        return result_df(emm_typing_notes=", ".join(notes))

    genes_in_operon = ",".join(blast_df_unique.get_column("qseqid").to_list())

    if blast_df_unique["sseqid"].n_unique() != 1:
        emm_genes = blast_df_unique.get_column("typed_allele").to_list()
        notes.append(
            "Unable to determine EMM type because EMM and EMM-like genes found on multiple contigs. Alleles found: "
            + "/".join(emm_genes)
        )
        return result_df(
            genes_in_operon=genes_in_operon,
            emm_typing_notes=", ".join(notes) if notes else "All exact allele matches",
        )

    first_row = blast_df_unique.row(0, named=True)
    ascending = first_row["sstart"] < first_row["send"]
    blast_df_unique = blast_df_unique.sort("sstart", descending=not ascending)

    n = blast_df_unique.height

    if n == 1:
        row = blast_df_unique.row(0, named=True)
        emm_type = row["typed_allele"]

        if row["length"] < row["qlen"] or row["pident"] < 100:
            notes.append(
                f"EMM{row['allele']} with {round(row['pident'], 2)} and length {row['length']}/{row['qlen']}"
            )

        return result_df(
            emm_type=emm_type,
            genes_in_operon=genes_in_operon,
            emm_typing_notes=", ".join(notes) if notes else "All exact allele matches",
        )

    if n == 2:
        row0 = blast_df_unique.row(0, named=True)
        row1 = blast_df_unique.row(1, named=True)

        emm_type = row0["typed_allele"]
        enn_type = row1["typed_allele"]
        mrp_type = "-"

        if row0["length"] < row0["qlen"] or row0["pident"] < 100:
            notes.append(
                f"EMM{row0['allele']} with pident {round(row0['pident'], 2)} and length {row0['length']}/{row0['qlen']}"
            )

        if row1["length"] < row1["qlen"] or row1["pident"] < 100:
            notes.append(
                f"ENN{row1['allele']} with pident {round(row1['pident'], 2)} and length {row1['length']}/{row1['qlen']}"
            )

        emm_maintype = row0["main_type"]
        mrp_maintype = row1["main_type"]

        if (
            mrp_maintype in emm_types_in_emm_plus_mrp_operons
            or emm_maintype in mrp_types_in_emm_plus_mrp_operons
        ):
            mrp_type = emm_type
            emm_type = enn_type
            enn_type = "-"
            notes.append("EMM redesignated due to known MRP+EMM operon")

        return result_df(
            emm_type=emm_type,
            enn_type=enn_type,
            mrp_type=mrp_type,
            genes_in_operon=genes_in_operon,
            emm_typing_notes=", ".join(notes) if notes else "All exact allele matches",
        )

    if n == 3:
        row0 = blast_df_unique.row(0, named=True)
        row1 = blast_df_unique.row(1, named=True)
        row2 = blast_df_unique.row(2, named=True)

        mrp_type = row0["typed_allele"]
        emm_type = row1["typed_allele"]
        enn_type = row2["typed_allele"]

        if row0["length"] < row0["qlen"] or row0["pident"] < 100:
            notes.append(
                f"MRP{row0['allele']} with pident {round(row0['pident'], 2)} and length {row0['length']}/{row0['qlen']}"
            )

        if row1["length"] < row1["qlen"] or row1["pident"] < 100:
            notes.append(
                f"EMM{row1['allele']} with pident {round(row1['pident'], 2)} and length {row1['length']}/{row1['qlen']}"
            )

        if row2["length"] < row2["qlen"] or row2["pident"] < 100:
            notes.append(
                f"ENN{row2['allele']} with pident {round(row2['pident'], 2)} and length {row2['length']}/{row2['qlen']}"
            )

        return result_df(
            emm_type=emm_type,
            enn_type=enn_type,
            mrp_type=mrp_type,
            genes_in_operon=genes_in_operon,
            emm_typing_notes=", ".join(notes) if notes else "All exact allele matches",
        )

    notes.append("More than three EMM genes found")
    details = (
        blast_df_unique.select(
            pl.concat_str(
                [
                    pl.col("qseqid"),
                    pl.lit(" with pident "),
                    pl.col("pident").round(2).cast(pl.Utf8),
                    pl.lit(" and length "),
                    pl.col("length").cast(pl.Utf8),
                    pl.lit("/"),
                    pl.col("qlen").cast(pl.Utf8),
                ],
                separator="",
            ).alias("combined")
        )
        .get_column("combined")
        .to_list()
    )
    notes.extend(details)

    return result_df(
        genes_in_operon=genes_in_operon,
        emm_typing_notes=", ".join(notes) if notes else "All exact allele matches",
    )


def type_sample(
        assembly_file: Path,
        output_folder: Path,
        emm_allele_fasta: Path,
        ) -> pl.DataFrame:

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    blast_output_tsv = output_folder.joinpath("emm_blast.tsv")
    blast_complete = run_emm_blast(
        assembly_file=assembly_file,
        emm_allele_fasta=emm_allele_fasta,
        output_file=blast_output_tsv,
    )

    if blast_complete:
        return extract_emm_type(emm_blast_tsv=blast_output_tsv)

    return pl.DataFrame(
        {
            "EMM_type": ["-"],
            "ENN_type": ["-"],
            "MRP_type": ["-"],
            "genes_in_operon": ["-"],
            "emm_typing_notes": ["Failed to run blast"],
        }
    )

def type_batch(
    assembly_files: list[Path],
    emm_allele_fasta: Path,
    output_dir: Path,
    full_path: bool = False,
) -> pl.DataFrame:
    result_dfs = []

    for assembly_file in assembly_files:
        sample_output_dir = output_dir / assembly_file.stem

        sample_name = str(assembly_file) if full_path else assembly_file.name

        sample_result = (
            type_sample(
                assembly_file=assembly_file,
                output_folder=sample_output_dir,
                emm_allele_fasta=emm_allele_fasta,
            )
            .with_columns(pl.lit(sample_name).alias("sample"))
            .select(["sample", "EMM_type", "ENN_type", "MRP_type", "genes_in_operon", "emm_typing_notes"])
        )

        result_dfs.append(sample_result)

    if not result_dfs:
        return pl.DataFrame(
            schema={
                "sample": pl.Utf8,
                "EMM_type": pl.Utf8,
                "ENN_type": pl.Utf8,
                "MRP_type": pl.Utf8,
                "genes_in_operon": pl.Utf8,
                "emm_typing_notes": pl.Utf8,
            }
        )

    return pl.concat(result_dfs, how="vertical")


def main():
    args = parse_args()
    print(f"Running emm typer on {len(args.input)} samples")
    emm_results_df = type_batch(assembly_files=args.input,
                                   emm_allele_fasta=args.database,
                                   output_dir=args.output,
                                   full_path=args.full_path)
    output_file = args.output / "results.tsv"
    print(f"Printing results to {output_file}")
    emm_results_df.write_csv(file = output_file,
                             separator = "\t")
    print("Done")

if __name__ == "__main__":
    main()


