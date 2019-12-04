#!/usr/bin/env python3
''' Create a qiime2 manifest file from a directory that contains fastq files,
accepted extensions are .fastq, .fq, .fastq.gz, or .fq.gz.  Ignores files and
folders that begin with "."

Requires a string directory name, can provide an output filename as a second
parameter.  If no filename given with write to manifest.csv.
'''
__author__ = "Nick Dickens"
__copyright__ = "Copyright 2018, Nicholas J. Dickens"
__email__ = "dickensn@fau.edu"
__license__ = "MIT"

import sys
import os


def create_manifest(file_list):
    ''' Write the file list out as a manifest file.  Uses R2 or val_2 to determine
    forward reads. Takes a list of filepaths
    '''
    # current dir is self.working_dir
    with open("manifest.csv","w") as manifest_fh:
        # write the header
        manifest_fh.write("sample-id,absolute-filepath,direction\n")
        for filename in file_list:
            direction = "forward"
            if "_R2" in filename or "val_2" in filename:
                direction = "reverse"

            # take up to the first underscore as the sample name
            #sample_name = filename.partition("_")[0]
            sample_name = filename.split("-")[1]
            absolute_path = os.path.abspath(filename)

            manifest_fh.write("{},{},{}\n".format(sample_name,absolute_path,direction))

def search_path(directory_name):
    ''' find all fastq files in a directory, descends sub-directories
    '''
    file_list = []
    print("Searching {}".format(directory_name))
    for path, dirnames, filenames in os.walk(directory_name):
        #print(path, dirnames, filenames)
        for filename in filenames:
            if not filename.startswith("."):
                if filename.endswith(".fastq") or filename.endswith(".fq") or \
                filename.endswith(".fq.gz") or filename.endswith(".fastq.gz"):
                    file_list.append(os.path.abspath("{}/{}".format(path, filename)))

    print("Found {} fastq files.".format(len(file_list)))
    file_list.sort()
    return file_list

if __name__ == "__main__":
    try:
        dirname = sys.argv[1]
    except:
        sys.exit("ERROR: Did you specify a directory name?")

    outfile = "manifest.csv"
    if len(sys.argv) > 2:
        outfile = sys.argv[2]


    file_list = search_path(dirname)
    #print(file_list)
    create_manifest(file_list)
