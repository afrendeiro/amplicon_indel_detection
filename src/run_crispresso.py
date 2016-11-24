import os
import pandas as pd


def reverse_complement(sequence):
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    return str(Seq(sequence, IUPAC.unambiguous_dna).reverse_complement())


def forward_or_reversecomplement(guide_rna, amplicon):
    try:
        amplicon.index(guide_rna)
        return guide_rna
    except ValueError:
        return reverse_complement(guide_rna)


def crispresso(fastq, amplicon, guide, name, output_dir, cpus=4, window=20):
    """
    """
    from pipelines import toolkit as tk
    import textwrap

    run_name = "_".join(["CRISPResso", name])
    log_file = os.path.join(output_dir, run_name + '.log.txt')
    job_file = os.path.join(output_dir, run_name + '.job.sh')

    cmd = tk.slurmHeader(run_name, log_file, cpusPerTask=4)

    cmd += """
        CRISPResso -w {} -r1 {} -a {} -g {} -n {} -o {} -p {}
    """.format(window, fastq, amplicon, guide, name, output_dir, cpus)
    cmd += tk.slurmFooter()

    with open(job_file, "w") as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurmSubmitJob(job_file)


# If template (amplicon sequence is known):
annotation = pd.read_csv(os.path.join("metadata", "annotation.csv"))

for sample in annotation.index:
    df = annotation.ix[sample]
    print(df["sample_name"])

    crispresso(
        df['sample_name'] + ".trimmed.fastq.gz",
        df["amplicon_sequence"],
        forward_or_reversecomplement(df["gRNA_seq"], df["amplicon_sequence"]),
        df["sample_name"], "/home/arendeiro/scratch/dropseq_runs/crispr_validations", 2)
