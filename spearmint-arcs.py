#!/usr/bin/env python
# Run ARCS for Spearmint

import csv
import os
import shlex
import subprocess

# Run ABySS
def main(job_id, params):
    c = params.get('c', 1)
    e = params.get('e', 50000)
    r = params.get('r', 0.05)
    a = params.get('a', 1.0)
    l = params.get('l', 1)
    tsv_filename = "psitchensismt.psitchensis.bx.atleast4.c%d_e%d_r%f.arcs.a%f_l%d.links.scaffolds.stats.tsv" % (c, e, r, a, l)
    command = "make c=%d e=%d r=%f a=%f l=%d %s" % (c, e, r, a, l, tsv_filename)
    print command
    subprocess.call(shlex.split(command))
    max_ng50 = 0
    with open(tsv_filename) as csvfile:
        csvfile.readline()
        csvreader = csv.reader(csvfile, delimiter="\t")
        for row in csvreader:
            ng50 = int(row[4])
            max_ng50 = max(max_ng50, ng50)
    print "c=%d\te=%d\tr=%f\ta=%f\tl=%d\tNG50=%d" % (c, e, r, a, l, ng50)
    return -max_ng50

if __name__ == "__main__":
    main(None, {})
