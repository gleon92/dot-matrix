import matplotlib.pyplot as plt
from Bio import SeqIO


def dot_matrix(file, window, threshold):
    in_handle = open(file)
    record_iterator = SeqIO.parse(in_handle, "fasta")
    rec_1 = next(record_iterator)
    rec_2 = next(record_iterator)
    chunks_1 = {}
    chunks_2 = {}
    for (seq, section_dict) in [
        (str(rec_1.seq).upper(), chunks_1),
        (str(rec_2.seq).upper(), chunks_2),
    ]:
        for i in range(len(seq) - window):
            section = seq[i: i + window]
            try:
                section_dict[section].append(i)
            except KeyError:
                section_dict[section] = [i]

    matches = set(chunks_1).intersection(chunks_2)
    x = []
    y = []
    for section in matches:
        for i in chunks_1[section]:
            for j in chunks_2[section]:
                x.append(i)
                y.append(j)
    plt.scatter(x, y)
    plt.xlabel("%s (length %i bp)" % (rec_1.id, len(rec_1)))
    plt.ylabel("%s (length %i bp)" % (rec_2.id, len(rec_2)))
    plt.show()


dot_matrix('flna.fasta', 5, 20)
