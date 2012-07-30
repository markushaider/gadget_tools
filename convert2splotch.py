#!/usr/bin/python
import os

# start from snapshot number
startnum = 0
endnum = 175

inputpath = '/home/markus/sim2/'
outputpath= '/home/markus/sim2/splotch/'
snapname='snapshot_cut_'

for i in range(startnum,endnum+1):
    #print i
    inputName = inputpath+snapname+str(i).zfill(3)
    outputName = outputpath+'splotch_'+str(i).zfill(3)
    print ' Converting '+inputName+' to '+ outputName
    print ' ---------------------------------------------------'
    os.system('./convert2gadgetformat2 '+inputName+' '+outputName)

