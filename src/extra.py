import os
import numpy as np
import pandas as pd


def get_amplicon_sequence(primer1, primer2):
    """
    Blat primers to get amplicon sequence.
    """
    raise NotImplementedError

    def blat(query, database="/data/groups/lab_bock/shared/resources/genomes/hg38/hg38.fa"):
        """
        """
        from subprocess import Popen, PIPE

        query_file = "tmp.fa"
        with open(query_file, "w") as handle:
            handle.writelines(">query\n" + query)
        output = "tmp.psl"

        process = Popen(['blat', query_file, database, output], stdout=PIPE, stderr=PIPE)
        result, stderr = process.communicate()

        return result

    # Blat primers
    s1 = blat(primer1)
    s2 = blat(primer2)

    # Get region in between
    region = (s1[0], s1[1], s2[2])
    amplicon = ""

    return {region: amplicon}


def generate_edits_reads(templates, weights=None, n_cells=1000, n_reads=10000, read_length=150):
    """
    """
    if weights is not None:
        assert len(templates) == len(weights)
    else:
        weights = [1.0 / len(templates)] * len(templates)

    # Additional params
    alphabet = ["A", "C", "T", "G"]
    geometric_p = 0.2

    reads = list()
    edits = list()

    # Generate fake edits in cells
    print("Generating edits in cells")
    for i in range(n_cells):
        # choose which template to use
        template = np.random.choice(templates, p=weights)
        # get random indels (higher prob for deletions)
        indel = np.random.geometric(geometric_p) * np.random.choice([-1, 1], p=[0.75, 0.25])
        # get a random position in the template
        position = np.random.choice(range(len(template)))
        # construct an edit
        if indel >= 0:
            insertion = "".join([np.random.choice(alphabet) for _ in range(indel)])
            edit = template[:position] + insertion + template[position:]
        else:
            edit = template[:position + indel] + template[position:]
        # append
        edits.append(edit)

    # Generate fake reads from the edits
    print("Generating reads from PCR")
    for i in range(n_reads):
        # choose an edit randomly
        edit = np.random.choice(edits)
        # choose a read start position (mimicking the enzyme cuts along the amplified region)
        position = np.random.choice(range(len(template)))
        # choose a strand
        direction = np.random.choice([-1, 1])
        # make a read
        if direction > 0:
            read = edit[position:position + read_length]
        else:
            read = edit[position - read_length:position]
        # skip small reads and append
        if len(read) > 10:
            reads.append(read)

    return (edits, reads)


def get_gtf():
    link = "ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz"
    cmd = "wget %s" % link
    os.system(cmd)

    os.system("gzip -d *gtf.gz")


# Read annotation in
annotation = pd.read_csv(os.path.join("metadata", "annotation.csv"))

# If only primers are given:
amplicons = get_amplicon_sequence()

# Generate fake reads from edits in cells
(edits, reads) = generate_edits_reads(annotation['sequence'], n_cells=2000, n_reads=10000)

# save edits (just for fun) and reads
with open("edits.txt", "w") as handle:
    handle.writelines("\n".join(edits))
with open("reads.txt", "w") as handle:
    handle.writelines("\n".join(reads))
