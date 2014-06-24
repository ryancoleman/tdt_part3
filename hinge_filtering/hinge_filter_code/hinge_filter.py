#!/usr/bin/env python

#filters out molecules that don't make contact in the so-called hinge region
#of some proteins. not general purpose, but works well.
#ryan coleman, joel karpiak, gpl v2.

import sys
import os
import mol2  # output of DOCK3.7 is mol2 files, use readdockmol2data method
import pdb  # reads pdb file containing filter region
import string

inputRegionFile = sys.argv[1]  # hardcoding these for now
inputLigandFile = sys.argv[2]  # input file, from dock
outputLigandFile = sys.argv[3]  # output file, just the good ones
hingeFilterDistance = 4.5  # angstroms

ligandData, ligandHeader = mol2.readDockMol2file(
    inputLigandFile, allheaders=True)
for count in xrange(len(ligandData)):
  for line in ligandHeader[count]:
    print string.strip(line)
  ligandData[count].writeMol2File(sys.stdout)
