# -*- coding: utf-8 -*-

'''
set tool path as global variables
Test samtools and bedtools
'''

import os
from subprocess import check_output

try:
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import minmax_scale
except:
    print "Error: required pyton modules are not found."
    print "Need pandas, matplotlib, sklearn.preprocessing"
    quit()

def init(samtools, bedtools):
    global samtoolspath
    samtoolspath = samtools
    global bedtoolspath
    bedtoolspath = bedtools

def Tools():
    '''
    Test whether bedtools or samtools works
    '''
    with open('mytools.version', 'w') as script:
        print>>script, '#!/bin/sh'
        print>>script, '{}bedtools --version'.format(bedtoolspath)
    try:
        output1 = check_output('sh mytools.version', shell=True)
        os.remove('mytools.version')
    except:
        print 'bedtools error.'
        os.remove('mytools.version')
        quit()
    
    with open('mytools.version', 'w') as script:
        print>>script, '#!/bin/sh'
        print>>script, '{}samtools --version'.format(samtoolspath)
    try:
        output2 = check_output('sh mytools.version', shell=True)
        os.remove('mytools.version')
    except:
        print 'samtools error.'
        os.remove('mytools.version')
        quit()
    return output1, output2
