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
                        help = "Path to folder with virulence gene data (fasta file and tsv)",
                        type=Path,
                        default = resources.files("streptofile") / "db" / "virulence_profiling")
    parser.add_argument("--full_path",
                        help = "Print full path to fasta input in output tsv rather than just file name",
                        action= "store_true",
                        default=False)
    return parser.parse_args()

def run_virulence_gene_blast(assembly_file: Path, virulence_gene_fasta_file: Path, output_file: Path) -> tuple(str) | None:
    output_file = Path(output_file)
    if not output_file.exists():
        cmd = f'blastn -query {virulence_gene_fasta_file} -subject {assembly_file} -qcov_hsp_perc 60 -perc_identity 80 -out {output_file} -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send sseq evalue bitscore"'
        result = subprocess.run(cmd, shell=True)
        if result.returncode != 0:
            return False
        else:
            return True
    else:
        return True

def load_virulence_gene_table(virulence_gene_tsv: Path) -> dict:
    try:
        virulence_gene_df = pl.read_csv(
            virulence_gene_tsv,
            separator="\t",
            has_header=True,
        )
    except pl.exceptions.NoDataError:
        raise Exception(f"Could not open virulence gene database tsv at {virulence_gene_tsv}")
    #virulence_genes = {
    #row["Gene"]: {k: v for k, v in row.items() if k != "Gene"}
    #for row in virulence_gene_df.to_dicts()
    #}
    return(virulence_gene_df)



def extract_virulence_gene_presence(virulence_blast_tsv: Path,
                                    length_threshold = 60, 
                                    pident_threshold = 80) -> pl.DataFrame:
    
    virulence_blast_tsv = Path(virulence_blast_tsv)
    if not virulence_blast_tsv.exists():
        virulence_blast_tsv["virulence_typing_notes"] = "No blast output found for virulence genes"
        return None

    try:
        blast_df = pl.read_csv(
            virulence_blast_tsv,
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
        raise Exception(f"Could not open virulence gene database tsv at {virulence_blast_tsv}")
    
    blast_df = (
    blast_df
    .with_columns([
        (pl.col("length") / pl.col("qlen") * 100).alias("plen")
    ])
    .filter(pl.col("plen") > length_threshold, pl.col("pident") > pident_threshold)
    )
    results_df = (
        blast_df
        .sort("bitscore", descending=True)
        .group_by("qseqid")
        .first()
        .select(["qseqid","plen","pident","sseqid","sstart","send","sseq"])
        .with_columns(
        pl.col("plen").round(2).alias("plen"),
        pl.col("pident").round(2).alias("pident")
        )
        )
    print(results_df)
    return(results_df)

def profile_sample(assembly_file: Path,
                   output_folder: Path,
                   virulence_gene_fasta: Path,
                   ) -> pl.DataFrame:
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok= True)
    blast_output_tsv = output_folder.joinpath("virulence_blast.tsv")
    blast_complete = run_virulence_gene_blast(assembly_file = assembly_file, virulence_gene_fasta_file=virulence_gene_fasta,
                                              output_file=blast_output_tsv)
    if blast_complete:
        virulence_results_df = extract_virulence_gene_presence(virulence_blast_tsv=blast_output_tsv)
        return virulence_results_df
    else:
        return None

def profile_batch(assembly_files: list[Path],
                  database_dir: Path,
                  output_dir: Path,
                  full_path: bool = False,
                  length_threshold = 60,
                  pident_threshold = 80,
                  ) -> tuple[pl.DataFrame]:
    virulence_results = []
    virulence_gene_fasta = database_dir / "virulence_gene_references.fasta"
    virulence_gene_tsv = database_dir / "virulence_genes.tsv"
    sample_names = []
    for assembly_file in assembly_files:
        sample_output_dir = output_dir / assembly_file.stem
        blast_tsv_path = sample_output_dir / "virulence_gene_blast.tsv"
        virulence_results_df = profile_sample(assembly_file=assembly_file,
                                                  output_folder=sample_output_dir,
                                                  virulence_gene_fasta = virulence_gene_fasta,
                                                  )
        if full_path:
            sample_name = str(assembly_file)
        else:
            sample_name = assembly_file.name
        sample_names.append(sample_name)
        virulence_results_df = virulence_results_df.with_columns(
            pl.lit(sample_name).alias("sample")
        ).select(
            ["sample"] + virulence_results_df.columns
        )
        virulence_results.append(virulence_results_df)
    combined_results = pl.concat(virulence_results)
    virulence_gene_table = load_virulence_gene_table(virulence_gene_tsv)
    combined_results = combined_results.rename({"qseqid": "Gene",
                                                "plen": "percent_coverage",
                                                "pident": "percent_identity",
                                                "sseqid": "contig",
                                                "sstart": "contig_start",
                                                "send": "contig_end",
                                                "sseq": "sequence"})
    combined_results = combined_results.join(virulence_gene_table, on="Gene", how="left")
    gene_names = virulence_gene_table.get_column("Gene").to_list()
    full = pl.DataFrame({
    "Gene": gene_names
    }).join(
        pl.DataFrame({"sample": sample_names}),
        how="cross"
    )
    counts = combined_results.group_by(["Gene", "sample"]).len()

    presence_absence_matrix = (
        full
        .join(counts, on=["Gene", "sample"], how="left")
        .with_columns(pl.col("len").fill_null(0))
        .pivot(index="sample", on="Gene", values="len")
    )
    return(combined_results, presence_absence_matrix)

def print_results(combined_results: pl.DataFrame,
                  presence_absence_matrix: pl.DataFrame,
                  output_folder: Path) -> None:
    emm_results_df = pl.from_dicts(virulence_results_dict)
    emm_results_df.write_csv(file = output_file,
                             separator = "\t")
    return None


def main():
    args = parse_args()
    combined_results, presence_absence_matrix = profile_batch(
        assembly_files=args.input,
        database_dir=args.database,
        output_dir = args.output
    )
    output_file = args.output / "results.tsv"
    matrix_output_file = args.output / "results.matrix.tsv"
    combined_results.write_csv(file = output_file, separator= "\t")
    presence_absence_matrix.write_csv(file = matrix_output_file, separator = "\t")
    

if __name__ == "__main__":
    main()


