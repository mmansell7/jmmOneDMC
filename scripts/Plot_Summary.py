#!/gpfs_share/santiso/SOFTWARE/miniconda2/envs/santimulators_1/bin/python

################################################################
#  Plot_Summary.py :
#    Plots Energy vs Temperature and Pressure on a 3D plot
#      using data from a user-supplied data file.
#
#  
#  Arguments:
#    1. summaryFile: this is the relative path to the file 
#          containing the data to be plotted.  The data in
#          this file must be arranged in columns with the
#          order: jobID, P, T, E
#
################################################################

import glob
import numpy as np
import sys
import matplotlib.pyplot as plt
import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

summaryFile = sys.argv[1]

#data = np.loadtxt(summaryFile)
data = np.genfromtxt(summaryFile,names=True)
#print(data)
#print('\n\n')
#print(data[:]['Energy'])
#print('Data[:,0]: ')
#print(data[:,0])
#print('Data[:,1]: ')
#print(data[:,1])
#print('Data[:,2]: ')
#print(data[:,2])
#print('Data[:,3]: ')
#print(data[:,3])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data[:]['P'],data[:]['T'],data[:]['Energy'])

ax.set_xlabel('P')
ax.set_ylabel('T')
ax.set_zlabel('E')

fig.savefig('tmp.png')









#def randrange(n, vmin, vmax):
#    '''
#    Helper function to make an array of random numbers having shape (n, )
#    with each number distributed Uniform(vmin, vmax).
#    '''
#    return (vmax - vmin)*np.random.rand(n) + vmin
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
#n = 100
#
## For each set of style and range settings, plot n random points in the box
## defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
#for c, m, zlow, zhigh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
#    xs = randrange(n, 23, 32)
#    ys = randrange(n, 0, 100)
#    zs = randrange(n, zlow, zhigh)
#    ax.scatter(xs, ys, zs, c=c, marker=m)
#
#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')
#
#fig.savefig('tmp.png')
#



