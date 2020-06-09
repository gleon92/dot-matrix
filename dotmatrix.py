import pylab
from Bio import SeqIO


def dot_matrix(file, window, threshold):
    in_handle = open(file)
    record_iterator = SeqIO.parse(in_handle, "fasta")
    rec_one = next(record_iterator)
    rec_two = next(record_iterator)
    dict_one = {}
    dict_two = {}
    for (seq, section_dict) in [
        (str(rec_one.seq).upper(), dict_one),
        (str(rec_two.seq).upper(), dict_two),
    ]:
        for i in range(len(seq) - window):
            section = seq[i: i + window]
            try:
                section_dict[section].append(i)
            except KeyError:
                section_dict[section] = [i]

    matches = set(dict_one).intersection(dict_two)
    x = []
    y = []
    for section in matches:
        for i in dict_one[section]:
            for j in dict_two[section]:
                x.append(i)
                y.append(j)

    pylab.cla()
    pylab.gray()
    pylab.scatter(x, y)
    pylab.xlim(0, len(rec_one) - window)
    pylab.ylim(0, len(rec_two) - window)
    pylab.xlabel("%s (length %i bp)" % (rec_one.id, len(rec_one)))
    pylab.ylabel("%s (length %i bp)" % (rec_two.id, len(rec_two)))
    pylab.title("DotMatrix FLNA  Window:%s)" % window)
    pylab.show()

    # for char1 in data_1:
    #     print(char1)


dot_matrix('flna.fasta', 5, 20)
