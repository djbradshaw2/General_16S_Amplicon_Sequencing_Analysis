from glob import glob
import gzip
import sys

from Bio import SeqIO



# dirname="exported-reads_trimmed.qza"
manifest_filename=sys.argv[1]


max=0
min=1000000
total = 0
count = 0

# read the manifest
filenames = []

with open(manifest_filename) as man_fh:
    man_fh.readline()
    for line in man_fh:
        (name, filename, direction) = line.split(",")
        filenames.append(filename)



for filename in filenames:
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as in_fh:
          for record in SeqIO.parse(in_fh,"fastq"):
            slen = len(record.seq)
            if slen<min:
              min = slen

            if slen>max:
              max = slen

            total += slen
            count += 1
    else:
        with open(filename, "rt") as in_fh:
          for record in SeqIO.parse(in_fh,"fastq"):
            slen = len(record.seq)
            if slen<min:
              min = slen

            if slen>max:
              max = slen

            total += slen
            count += 1



print("File: {} Min: {}, Max:{}, Mean: {}".format(manifest_filename, min, max, total/count))
