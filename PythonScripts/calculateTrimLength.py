#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python SubFunction to calculate trim length based on TF sequence length
Part of the ChIP-exo data analysis pipeline
@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
"""
import sys

seqLength=int(sys.argv[1])

#calculate weight of TF assuming an average amino acid weight of 110 daltons, taken from https://www.promega.com/~/media/files/resources/technical%20references/amino%20acid%20abbreviations%20and%20molecular%20weights.pdf
tfWeigth=(seqLength/3*110)
#calculate radius of TF assuming shape to be spherical, see https://doi.org/10.1007/s12575-009-9008-x
tfRadius=0.066*tfWeigth**(1/3)
#Translate sphere radius into footprint on DNA, assuming half overlapping monomers => 3 x radius x 3.03bp/nm
tfFootprint=3*tfRadius*3.03
print(round(tfFootprint))

