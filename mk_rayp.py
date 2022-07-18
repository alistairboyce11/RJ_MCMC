#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################### RUN IMPORTS ############################

# Code to make rayp files for Pforward code.

import numpy
import random
import sys

if len(sys.argv)!=3:
    print('>>    python mk_rayp.py <LENGTH> <FILENAME>')
    print('>>    python mk_rayp.py 100 rayp2.in')
    print('Length or rayp-filename missing....')
    sys.exit('Not enough arguments...')
else:
    length=int(sys.argv[1])
    file_name=str(sys.argv[2])

############################################################################

# file_name='rayp2.in'
# length = 100
max = 8.840
min = 4.642
ran_floats = numpy.random.rand(length) * (max-min) + min
ran_floats


file = open(file_name, 'w')

file.write("%i \n" % (int(length)))
for i in range(length):
    file.write("%f \n" % (float(ran_floats[i])))
file.close()
print('Produced: '+str(file_name)+'...')