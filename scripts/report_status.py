#!/usr/bin/env python
"""
Generate a report of the current impute status.
Can also be used on the output folder of an in progress run.

Example usage::
$ ./scripts/report_status.py -t ./scripts/report_status_ideogram.tpl \
    -i output/RUN_NAME/ -o output/RUN_NAME/status.html
"""

from __future__ import print_function
from string import Template
from glob import glob
import datetime
import time
import os.path
import argparse

def get_all_chunks(infolder):
    """ Return a dict of all chunk tuples, key is chrm """
    fin = open(os.path.join(infolder, "chunks.txt"))
    result = {}
    for line in fin:
        chrm, start, stop = line.strip().split(",")
        if chrm in result:
            result[chrm] = sorted(result[chrm] + [(chrm, start, stop)])
        else:
            result[chrm] = [(chrm, start, stop)]
    fin.close()
    return result

def get_wanted_chrms(infolder):
    """ peek inside the 'qc' folder and get the wanted chrms """
    paths = glob(os.path.join(infolder, "qc/*"))
    chrms = map(os.path.basename, paths)
    return chrms

def read_summary(summaryfile):
    """ read a single chunk's summary file """
    chrm = summaryfile.split("/")[-2]
    start, stop = os.path.basename(summaryfile).split(".")[-2].split("_")[-1].split("-")
    stats = []
    fin = open(summaryfile)
    for line in fin:
        if "Interval  #Genotypes" in line:
            break
    else:
        return chrm, int(start), int(stop), stats, 'empty'
    for line in fin:
        if not line.strip():
            return chrm, int(start), int(stop), stats, 'good'
        stats.append("\t".join(line.strip().split()))
    return chrm, int(start), int(stop), stats, 'error'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-i', '--infolder',
                        help='root output folder for this chipimputation run')
    parser.add_argument('-t', '--template',
                        help='template file for status.html',
                        default="scripts/report_status_ideogram.tpl")
    parser.add_argument('-o', '--outfile', help='output file for status.html')
    args = parser.parse_args()
    PATH = os.path.join(args.infolder,
                        "impute_results/impute/**/*imputed_summary")
    chunkdict = {}
    COLORS = {
        'good': '#6F6A',
        'empty': '#66F9',
        'error': '#F33F'
    }
    for filename in glob(PATH):
        print(filename)
        chrm, start, stop, stats, status = read_summary(filename)
        chunk = (chrm, chrm, start, stop, COLORS.get(status, "#000F"))
        if chrm in chunkdict:
            chunkdict[chrm] = sorted(chunkdict[chrm] + [chunk])
        else:
            chunkdict[chrm] = [chunk]
        if stats:
            print("\n".join(stats))
    # HACK: shrink the last chunk of every chromosome
    for key in sorted(chunkdict):
        print(key, len(chunkdict[key]))
        last_chunk = chunkdict[key][-1]
        # HACK: worse hack, chrm 19 has a whole chunk that doesn't map
        if key == '19':
            chunkdict[key] = chunkdict[key][:-1]
        else:
            chunkdict[key][-1] = (last_chunk[0], last_chunk[1],
                                  last_chunk[2], last_chunk[2]+11111,
                                  last_chunk[4])
        print("\n".join(map(repr,chunkdict[key])))
    chunklist = []
    for key in sorted(chunkdict):
        chunklist = chunklist + chunkdict[key]
    chunkstrings= ["""\n{ name: '%s', chr: '%s', start: %i, stop: %i, color:'%s'}"""%x for x in chunklist]
    chunkstring = "[" + ",".join(sorted(chunkstrings)) + "]"
    ts = time.time()
    context = {}
    context["timestamp"] = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    context["run_name"] = os.path.split(args.infolder.strip("/"))[-1]
    chromosomes = get_wanted_chrms(args.infolder)
    all_chunks = get_all_chunks(args.infolder)
    progress = {}
    for chrm in chromosomes:
        done = len(chunkdict[chrm])
        wanted = len(all_chunks[chrm])
        progress[chrm] = "%i%% (%i/%i)"%(100.0*done/wanted, done, wanted)
    context["progress"] = repr(progress)
    context["all_chunks"] = all_chunks
    context["chunkstring"] = chunkstring
    context["chromosomes"] = repr(chromosomes)
    TPL = Template(open(args.template).read())
    with open(args.outfile, "wt") as fout:
        fout.write(TPL.substitute(**context))
