#!/
### Plot_AutoCorrelation.py ###
#####################################################################
# This is a module which allows one to quickly calculate and plot
# auto-correlation of data (over a single block) as a function of
# block size.
#
# Matt Mansell : NCSU : 2019/01/18
#
#####################################################################

import numpy as np
import matplotlib.pyplot as plt

def Plot_AutoCorrelation(B):
    # B should be a 1xN numpy array
    c = 1000*np.ones(B.shape[1]/2+1)
    c[0] = 1
    for blockSize in range(1,B.shape[1]/2+1):
        numBlocks = B.shape[1]/blockSize
        m = 1000*np.ones((1,numBlocks))
        print('blockSize,numBlocks= ' + str(blockSize) + ',' + str(numBlocks))
        for jj in range(0,numBlocks):
            m[0,jj] = np.mean(B[0,(jj*blockSize):(jj+1)*blockSize])
            #print('jj,m = ' + str(jj) + ',' + str(m))    
        print('m: ' + str(m))
        c[blockSize] = np.corrcoef(m[0,:-1],m[0,1:])[0,1]
    	print('c[' + str(jj) + ']: ' + str(c[jj]))
    lines = plt.plot(c)
    plt.xlabel('Block size')
    plt.ylabel('Correlation Coefficient')
    plt.ylim((-1,1.5))
    fig = plt.savefig('CorrCoef.png')
    return c,lines

def unMeaned(B):
    length = B.shape[1]-2
    c = 1000*np.ones(length)
    c[0] = 1
    for delta in range(1,length):
        print('delta = ' + str(delta))
        c[delta] = np.corrcoef(B[0,:-delta],B[0,delta:])[0,1]
    lines = plt.plot(c)
    plt.xlabel(r'$\Delta$t (Correlation ''time'')')
    plt.ylabel('Correlation Coefficient')
    plt.ylim((-1,1))
    plt.savefig('CorrCoef_unMeaned.png')
    return c,lines

def test_unMeaned(length=20000,offset=0,width=2.0):
    A = offset + width*( np.random.random_sample(length) - 0.5 )
    A = A.reshape(1,length)
    c,lines = unMeaned(A)
    return c,lines

