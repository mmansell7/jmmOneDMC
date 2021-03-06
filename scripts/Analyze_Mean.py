#!/gpfs_share/santiso/SOFTWARE/miniconda2/envs/santimulators_1/bin/python
#BSUB -W 03:00
#BSUB -o output.%J
#BSUB -e errors.%J
#BSUB -q gubbins
#BSUB -J An.Mean.

################################################################
#  Analyze_Mean.py :
#    Creates plots of the mean energy and mean length for each
#      simulation having a subdirectory contained in a given
#      directory (the first argument). Saves these plots in the
#      current diretory.
#      The overall mean energies and lengths are written to 
#      "Summary_tmp.txt".
#
#  
#  Arguments:
#    1. dir: the directory in which simulation subdirectories
#             are located.
#         
#  Supplements:
#    ${dir}/EqSteps.dat: this file must be used to list the
#          number of equilibration steps to be removed from
#          any of the data files. The format is P T EquilSteps
#          with any whitespace separating the values.
#
#
#  Workflow:
#    A good practice is to set each line in ${dir}/EqSteps.dat
#      to 0  0  0 and run
#      ./Analyze_Mean.py ${dir}.  Then fill in EqSteps.dat by
#      visually assessing convergence of each simulation in
#      the newly generated plots.  Then, re-run 
#      ./Analyze_Mean.py to get the equilibrium means.
#
#
################################################################


import glob
import numpy as np
import sys
import matplotlib.pyplot as plt
import re
#import wxversion
##wxversion.select('2.8')

blockSize = int(1000000)
printInterval = 10000

Nstr='2000'
potentialType='LJ'
if potentialType == 'LJ':
	summaryFile = open("LJ_Means.dat","a")
	# Use the first option for command-line input of the target directory
	dir = "../data/LJ/m-1/Longest/"
elif potentialType == 'HARMONIC':
	summaryFile = open("HARMONIC_Means.dat","a")
	# Use the first option for command-line input of the target directory
	dir = "../data/HARMONIC/m1/"
summaryFile.write('P\tT\tEnergy\tEnergy2\tLength\tLength2\tLengthEnergy\n')

ms = 4  # marker size for plots

# Use the first option for command-line input of the target directory
# dir = sys.argv[1]

#eqSteps = int(sys.argv[2])

#print("eqSteps: " + str(eqSteps))
print("printInterval: " + str(printInterval))
print("blockSize: " + str(blockSize))
#print("Analyzing runs in " + dir + " using " + str(eqSteps)  + " equilibration steps." + "\n")
print("Analyzing runs in " + dir + "\n")

eqStepsFile = glob.glob(dir + "EqSteps.dat")[0]
print("Found equilibration steps file: " + str(eqStepsFile))
eqData = np.genfromtxt(eqStepsFile)
#print("eqData: " + str(eqData))

kB    = 1.38064852E-23;  # J/K
sigma = 3.4E-10;         # m
epsilon = 125.0*kB;      # J


#thermoFiles = glob.glob(dir + "*/*thermo*")
#dir = "/gpfs_backup/gubbins_data/jmmansel/jmmOneDMC/data/LJ/m-1/N200/"
#Nstr='200'
thermoFiles = glob.glob(dir + "P0.1*T0.8*/*thermo*")
print("Found files: " + str(thermoFiles) + "\n")
kk = 0
for ii in thermoFiles:
	kk = kk + 1
	print("Analyzing file: " + ii + "\n")
	T = np.float(re.search('T(.*)_',ii).group(1))
	P = np.float(re.search('P(.*)_T',ii).group(1))
	#jobID = re.search('T.*_([0-9]*)\/',ii).group(1)
	#print('P: ' + str(P) + '  T: ' + str(T) + '  jobID: ' + str(jobID))
	print('P: ' + str(P) + '  T: ' + str(T) + '  N: ' + Nstr)
	
	thermo_data_raw = np.genfromtxt(ii,names=True)
	#thermo_data_raw = np.genfromtxt(ii,names=["Step","Energy","Energy2","l","l2","Virial","Virial2","lEnergy"],skip_header=1,filling_values=0,usecols=(0,1,2,3,4,5,6))
	condition = np.logical_and(eqData[:,0] == P,eqData[:,1] == T)
	#condition = np.logical_and(np.absolute(np.subtract(eqData[:,0],P)) < 0.0001  ,np.absolute(np.subtract(eqData[:,1] - T)) < 0.0001)
	#condition = np.logical_and(eqData[:,0] == P  ,eqData[:,1] == T)
	#print("Condition: " + str(condition))
	#eqSteps = np.compress(condition,eqData,axis=0)
	eqSteps = np.compress(condition,eqData,axis=0)[0,2].astype(int)
	LEstartSteps = np.compress(condition,eqData,axis=0)[0,4].astype(int)
        LEstartBlock = (LEstartSteps - eqSteps)/blockSize
        #eqSteps = np.compress(np.logical_and(eqData[:,0] == P,eqData[:,1] == T),eqData,axis=0)
	#print(eqSteps)
	print('Equilibration steps: ' + str(eqSteps))
        thermo_data_raw = thermo_data_raw[eqSteps/printInterval+1:]
	#for rawSampleNumber in range(np.size(thermo_data_raw)):
		#print("rawSampleNumber: " + str(rawSampleNumber) + " rawSample: " + str(thermo_data_raw[rawSampleNumber]))
	
	energy_tot_mean  = np.mean(thermo_data_raw[:]['Energy'])
	energy2_tot_mean = np.mean(thermo_data_raw[:]['Energy2'])
	l_tot_mean       = np.mean(thermo_data_raw[:]['l'])
	l2_tot_mean      = np.mean(thermo_data_raw[:]['l2'])
	lenergy_tot_mean = np.mean(thermo_data_raw[(LEstartSteps - eqSteps)/printInterval:]['lE'])
        
	print("Total mean energy  : " + str(energy_tot_mean));
	print("Total mean length  : " + str(l_tot_mean));
	print("Total mean energy^2: " + str(energy2_tot_mean));
	print("Total mean length^2: " + str(l2_tot_mean));
	print("Total mean l*E     : " + str(lenergy_tot_mean));
	summaryFile.write(str(P) + '\t' + str(T) + '\t' + str(energy_tot_mean) + '\t' + str(energy2_tot_mean) + '\t' 
                            + str(l_tot_mean) + '\t' + str(l2_tot_mean) + '\t' + str(lenergy_tot_mean) + '\n')
	#summaryFile.write(str(jobID) + '\t' + str(P) + '\t' + str(T) + '\n')
	
	for rawSampleNumber in range(np.size(thermo_data_raw) - 1):
		if(thermo_data_raw[rawSampleNumber + 1]['Step'] - thermo_data_raw[rawSampleNumber]['Step'] != printInterval):
			
                        print("########## ERROR: INCONSISTENT STEP SPACING IN THERMO FILE ##############\n")
			print("########## FAILURE OCCURED AT STEP NUMBER: " + str(rawSampleNumber) + " #############")
			sys.exit()
	numBlocks = np.floor( np.size(thermo_data_raw) / (blockSize / printInterval) ).astype(int)
	print("numBlocks: " + str(numBlocks))
	thermo_data_block = np.empty((0,8))
	for blockNum in range(numBlocks):
		step_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Step'].astype(np.uint64),axis=0)
		energy_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Energy'].astype(np.float),axis=0)
		esq_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Energy2'].astype(np.float),axis=0)
		length_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['l'].astype(np.float))
		lsq_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['l2'].astype(np.float),axis=0)
		virial_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Virial'].astype(np.float),axis=0)
		virsq_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Virial2'].astype(np.float),axis=0)
		le_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['lE'].astype(np.float),axis=0)
		#print("Mean Step: " + str(step_data_block) + "  Energy: " + str(energy_data_block) + "  Energy2: " + str(esq_data_block)
		#	 + "  Length: " + str(length_data_block) + "  Length2: " + str(lsq_data_block) + "  Virial: " + str(virial_data_block)
		#	 + "  Virial2: " + str(virsq_data_block))
		thermo_data_block_tmp = np.array([step_data_block,energy_data_block,esq_data_block,length_data_block,lsq_data_block,virial_data_block,virsq_data_block,le_data_block])
		thermo_data_block_tmp.shape = 1,8
		thermo_data_block = np.append(thermo_data_block,thermo_data_block_tmp,axis=0)
		
		#print(str(thermo_data_block))
		#print(str(thermo_data_block[:,0]))
		#print(str(thermo_data_block[:,1]))
	
        # "ea" = "ensemble average"
        eaE    = thermo_data_block[:,1]
        eaESq  = thermo_data_block[:,2]
        eaL    = thermo_data_block[:,3]
        eaLSq  = thermo_data_block[:,4]
        eaLE   = thermo_data_block[:,7]
        cp     = eaESq - eaE*eaE + 2.0*P*(eaLE -eaL*eaE) + P*P*(eaLSq - eaL*eaL)
        cp     = cp/T/T
        #betaT  = eaLSq - eaL*eaL
        betaT  = 1.0/eaL/T*(eaLSq - eaL*eaL)
        alphaP = 1/T/T/eaL*( (eaLE - eaL*eaE) + P*(eaLSq - eaL*eaL) )
        gammaV = alphaP/betaT
        betaS  = betaT - alphaP*alphaP*T*eaL/cp
        muJT   = eaL/cp*(alphaP*T - 1.0)

        print('Shapes (eaE,eaESq,eaL,eaLSq,eaLE,cp,muJT,thermo_data_block[:,0]: ' + str(eaE.shape) + ', ' + str(eaESq.shape) + ', ' 
               + str(eaL.shape) + ', ' + str(eaLSq.shape) + ', ' + str(eaLE.shape) + ', ' + str(cp.shape) + ', ' + str(muJT.shape) + ', ' + str(thermo_data_block[:,0].shape))
	#for num in range(np.size(eaL)):
		#print("SampleNumber: " + str(num) + " rawSample: " + str(thermo_data_block[num,0]) + "L,L2,betaT: " + str(eaL[num]) + "," + str(eaLSq[num]) + "," + str(betaT[num]) )
        
	#plt.scatter(thermo_data_block[:,0],thermo_data_block[:,1],s=ms)
	plt.scatter(thermo_data_block[LEstartBlock:,0],epsilon * thermo_data_block[LEstartBlock:,1],s=ms)
        plt.ylim((-4.2E-19,-3.5E-19))
	fig = plt.gcf()
	fig.set_dpi(200)
	#plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.title('N: {0:d}, P: {1:.4g} N, T: {2:.0f} K'.format(2000,P * epsilon/sigma, T * epsilon/kB),fontsize=18,loc='right')
	plt.xlabel(r'Step Number',fontsize=18)
	#plt.ylabel(r'$E_{conf}$',fontsize=18)
	plt.ylabel(r'$E_{conf}$ (J)',fontsize=18)
#	# plt.xticks(fontsize=18)
#	#plt.tick_params(axis='both', which='major', labelsize=18)
#	#ax = plt.gca()
#	#ax.yaxis.offsetText.set_fontsize(18)
#	#ax.xaxis.offsetText.set_fontsize(18)
#	#plt.tight_layout()
#	#ax.text(0.05, 0.17, 'Block size = ' + '{:.2E}'.format(blockSize), transform=ax.transAxes, fontsize=10,
#        #	verticalalignment='top')
#	#ax.text(0.05, 0.11, 'Mean = ' + '{:.2f}'.format(), transform=ax.transAxes, fontsize=10,
#       	#	verticalalignment='top')
#	#ax.text(0.05, 0.05, 'StdErr = ' + '{:.2f}'.format(std), transform=ax.transAxes, fontsize=10,
#        #	verticalalignment='top')
#	# plt.show()
	plt.savefig('MeanE_P' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure', bbox_inches="tight", pad_inches=0.5)
	plt.close()
	
	#plt.scatter(thermo_data_block[:,0],thermo_data_block[:,3],linewidth=0.8)
	plt.scatter(thermo_data_block[:,0],thermo_data_block[:,3],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'L',fontsize=18)
#	# plt.xticks(fontsize=18)
#	#plt.tick_params(axis='both', which='major', labelsize=18)
#	#ax = plt.gca()
#	#ax.yaxis.offsetText.set_fontsize(18)
#	#ax.xaxis.offsetText.set_fontsize(18)
#	#plt.tight_layout()
#	#ax.text(0.05, 0.17, 'Block size = ' + '{:.2E}'.format(blockSize), transform=ax.transAxes, fontsize=10,
#        #	verticalalignment='top')
#	#ax.text(0.05, 0.11, 'Mean = ' + '{:.2f}'.format(), transform=ax.transAxes, fontsize=10,
#       	#	verticalalignment='top')
#	#ax.text(0.05, 0.05, 'StdErr = ' + '{:.2f}'.format(std), transform=ax.transAxes, fontsize=10,
#        #	verticalalignment='top')
#	# plt.show()
	plt.savefig('MeanL' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure')
	plt.close()
        
	plt.scatter(thermo_data_block[LEstartBlock:,0],thermo_data_block[LEstartBlock:,2],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'$E_{conf}^2$',fontsize=18)
	plt.savefig('MeanESq' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure', bbox_inches="tight", pad_inches=0.5)
	plt.close()
        
	plt.scatter(thermo_data_block[:,0],thermo_data_block[:,4],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'$L^2$',fontsize=18)
	plt.savefig('MeanLSq' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure')
	plt.close()
        
	plt.scatter(thermo_data_block[LEstartBlock:,0],thermo_data_block[LEstartBlock:,7],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'$L\times E_{conf}$',fontsize=18)
	plt.savefig('MeanLE' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure')
	plt.close()
        
	plt.scatter(thermo_data_block[LEstartBlock:,0],cp[LEstartBlock:],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'$c_{P,conf}',fontsize=18)
	plt.savefig('MeanCp' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure')
	plt.close()
        
	plt.scatter(thermo_data_block[:,0],betaT[:],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'$\beta_{T,conf}$',fontsize=18)
	plt.savefig('MeanBetaT' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure')
	plt.close()
        
	plt.scatter(thermo_data_block[LEstartBlock:,0],betaS[LEstartBlock:],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'$\beta_{S,conf}$',fontsize=18)
	plt.savefig('MeanBetaS' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure')
	plt.close()
        
	plt.scatter(thermo_data_block[LEstartBlock:,0],alphaP[LEstartBlock:],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'$\alpha_{P,conf}$',fontsize=18)
	plt.savefig('MeanAlphaP' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure')
	plt.close()
        
	plt.scatter(thermo_data_block[LEstartBlock:,0],gammaV[LEstartBlock:],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'$\gamma_{V,conf}$',fontsize=18)
	plt.savefig('MeanGammaV' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure')
	plt.close()
        
	plt.scatter(thermo_data_block[LEstartBlock:,0],muJT[LEstartBlock:],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'$\mu_{JT}$',fontsize=18)
	plt.savefig('MeanMuJT' + str(P) + '_T' + str(T) + '_N' + Nstr + '.png',dpi='figure')
	plt.close()
        
#summaryFile.close()
#summaryFile.close()

