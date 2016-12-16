#!/usr/bin/env python
# Run ABySS for Spearmint

import csv
import os
import shlex
import subprocess

# Run ABySS
def main(job_id, params):
    k = params['k']
    kc = params['kc']
    directory = "abyss/2.0.1/k%d/kc%d" % (k, kc)
    subprocess.call(shlex.split("make k=%d kc=%d %s/psitchensis-scaffolds.stats.tsv" % (k, kc, directory)))
    subprocess.call(shlex.split("make k=%d kc=%d %s/psitchensis-scaftigs.stats.tsv" % (k, kc, directory)))
    subprocess.call(shlex.split("make k=%d kc=%d %s/organelles.psitchensis-scaffolds.samtobreak.tsv" % (k, kc, directory)))
    subprocess.call(shlex.split("make k=%d kc=%d %s/organelles.psitchensis-scaftigs.samtobreak.tsv" % (k, kc, directory)))
    max_ng50 = 0
    with open("%s/psitchensis-scaffolds.stats.tsv" % directory) as csvfile:
        csvfile.readline()
        csvreader = csv.reader(csvfile, delimiter="\t")
        for row in csvreader:
            ng50 = int(row[4])
            max_ng50 = max(max_ng50, ng50)
    print "k=%d\tkc=%d\tN50=%d" % (k, kc, ng50)
    return -max_ng50

if __name__ == "__main__":
    main(None, {"k": 84, "kc": 3})
