#!/usr/bin/env python

#filters out molecules that don't make contact in the so-called hinge region
#of some proteins. not general purpose, but works well.
#ryan coleman, joel karpiak, gpl v2.

import sys
import os
import mol2  # output of DOCK3.7 is mol2 files, use readdockmol2data method
import pdb  # reads pdb file containing filter region
import geometry  # distance checkers
import string

inputRegionFile = sys.argv[1]  # hardcoding these for now
inputLigandFile = sys.argv[2]  # input file, from dock
outputLigandFile = sys.argv[3]  # output file, just the good ones
hingeFilterDistance = 4.5  # angstroms
hingeFilterDistSquared = hingeFilterDistance**2.0
ligandData, ligandHeader = mol2.readDockMol2file(
    inputLigandFile, allheaders=True)
filterPdb = pdb.pdbData(inputRegionFile)
outfile = open(outputLigandFile, 'w')
for count in xrange(len(ligandData)):
  ligXyz = ligandData[count].atomXyz[0]  # the 0 is because some mol2s have
  #multiple conformations but never dock3.7 written mol2s so just get the 0th
  passedAll = True
  for xyz in filterPdb.coords:
    passedOne = False
    for oneLig in ligXyz:
      if geometry.distL2Squared3(oneLig, xyz) < hingeFilterDistSquared:
        passedOne = True
        break
    if not passedOne:
      passedAll = False
      break
  if passedAll:  # means it passes the filters
    for line in ligandHeader[count]:
      outfile.write(line)       
    ligandData[count].writeMol2File(outfile)
