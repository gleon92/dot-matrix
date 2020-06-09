import matplotlib.pyplot as plt
from Bio import SeqIO


def identity(win_1, win_2, threshold):
    res = len(set(win_1) & set(win_2)) / max(len(win_1), len(win_2)) * 100
    # print(res)
    if res >= threshold:
        return True
    else:
        return False


def dot_matrix(file, window, threshold):
    in_handle = open(file)
    record_iterator = SeqIO.parse(in_handle, "fasta")
    rec_1 = next(record_iterator)
    rec_2 = next(record_iterator)
    chunks_1 = [str(rec_1.seq).upper()[i: i + window] for i in range(0, len(str(rec_1.seq).upper()), window)]
    chunks_2 = [str(rec_2.seq).upper()[i: i + window] for i in range(0, len(str(rec_2.seq).upper()), window)]
    x = []
    y = []

    for i in range(len(chunks_1)):
        for j in range(len(chunks_2)):
            if identity(chunks_1[i], chunks_2[j], threshold):
                x.extend([n for n in range(i*window, (i*window)+window)])
                y.extend([n for n in range(j*window, (j*window)+window)])

    plt.scatter(x, y)
    plt.xlabel("%s (length %i bp)" % (rec_1.id, len(rec_1)))
    plt.ylabel("%s (length %i bp)" % (rec_2.id, len(rec_2)))
    plt.title("window: %s | threshold: %s" % (window, threshold))
    plt.show()


dot_matrix('flna.fasta', 25, 50)
# print(identity('MSSSH', 'GSASC', 21))

