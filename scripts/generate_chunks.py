#!/usr/bin/env python
"""
Reads map file, chunk size
Returns file with chromosome chunk_start chunk_end
"""
import sys

def chunk_split(map_file, output, chunk):
    '''
    Return: chunck files in the output folder
    '''
    data = [ [dat.split()[0], dat.split()[3]] for dat in open(map_file).readlines() ]
    datas = {}
    out = open(output, 'w')
    for dat in data:
        chrm = int(dat[0])
        myPos = int(dat[1])
        if chrm not in datas:
            datas[chrm] = []
        datas[chrm].append(myPos)
    data = {}
    for chrm in sorted(datas):
        if chrm not in data:
            data[chrm] = []
        chunk = int(chunk)
        max_ = max(datas[chrm])
        min_ = min(datas[chrm]) - (min(datas[chrm]) % 10) + 1
        myPos = list(range(min_, max_, chunk))
        for pos in myPos:
            start_ = pos
            end_ = start_ + chunk - 1
            out.writelines(','.join([str(chrm), str(start_), str(end_)])+'\n')
    out.close()

bimFile = sys.argv[1]
outputFile = sys.argv[2]
chunk_size = sys.argv[3]
chunk_split(bimFile, outputFile, chunk_size)