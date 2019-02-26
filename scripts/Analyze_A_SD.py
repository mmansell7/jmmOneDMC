#!/gpfs_share/santiso/SOFTWARE/miniconda2/envs/santimulators_1/bin/python

#BSUB -W 00:59
#BSUB -o output.%J
#BSUB -e errors.%J
#BSUB -q debug
#BSUB -J AnalyzeSD

################################################################
#  Analyze_A_SD.py :
#    Creates plots of the StDev of the mean energy and length
#      as functions of the block size for each
#      simulation having a subdirectory contained in a given
#      directory (the first argument). Saves these plots in the
#      current diretory.
#      Uses only production steps.
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
#      ./Analyze_mean.py ${dir}.  Then fill in EqSteps.dat by
#      visually assessing convergence of each simulation in
#      the newly generated plots.  Then, run 
#      ./Analyze_A_SD.py to get assess the appropriate block
#      sizes.
#
#
################################################################


import glob
import numpy as np
import sys
import matplotlib.pyplot as plt
import re
import os

blockSizeMin = 100000
blockSizeMax = 1000000000

printInterval = 10000

# Use the line below to give the target directory as a command line argument.
#dir = sys.argv[1]
# Use the line below to hard-code the target directory (for batch submission).
dir = "../data/LJ/m-1/Longest/"
#dir = "../data/HARMONIC/m1/"


#eqSteps = int(sys.argv[2])

#print("eqSteps: " + str(eqSteps))
print("printInterval: " + str(printInterval))
#print("Analyzing runs in " + dir + " using " + str(eqSteps)  + " equilibration steps." + "\n")
print("Analyzing runs in " + dir + "." + "\n")

#thermoFiles = glob.glob(dir + "*/*thermo*")
thermoFiles = glob.glob(dir + "*P0.1*/*thermo*")
print("Found files: ")
for tf in thermoFiles:
	print(tf)
kk = 0

eqStepsFile = glob.glob(dir + "EqSteps.dat")[0]
print("Found equilibrium steps file: " + str(eqStepsFile))
eqData = np.genfromtxt(eqStepsFile)
#print("eqData: " + str(eqData))

summaryFile = open("Summary_SD_tmp.txt","a")
summaryFile.write('P\tT\tblockSize\tnumBlocks\tSDM_Energy\tSDM_Length\tSEM_Energy\tSEM_Length\n')

### Use the following line to analyze a subset of the available thermo files ###
#thermoFiles = thermoFiles[0:1]
################################################################################

for tf in thermoFiles:
#for tf in thermoFiles[::-1]:
	thisDir = re.split('thermo.dat',tf)[0]
	dataBlockingResultsFile = thisDir + 'DataBlockingResults.dat'
	
	T = np.float(re.search('T(.*)_',tf).group(1))
	P = np.float(re.search('P(.*)_T',tf).group(1))
	#jobID = re.search('T.*_([0-9]*)\/',tf).group(1)
	#print('P: ' + str(P) + '  T: ' + str(T) + '  jobID: ' + str(jobID))
	print('P: ' + str(P) + '  T: ' + str(T))
	
	# Check for existence of DataBlockingResults.dat, skip building if it already exists
	if os.path.isfile(dataBlockingResultsFile):
		stdevs = np.genfromtxt(dataBlockingResultsFile,names=True)
		print('Shape of stdevs: ' + str(stdevs.shape))
		if stdevs.shape[0] > 0:
			print(dataBlockingResultsFile + ' already exists.  Using the existing file.\nTo build a new file, first delete the existing one.\n')
			fe = 1
			stdevs = stdevs.view(np.float64).reshape(stdevs.shape[0],-1)  # Reshape stdevs from a structured array (1D array of 1x3 tuples) to unstructured 2D array
		else:
			print(dataBlockingResultsFile + ' already exists but is empty...It will be replaced with a new file.\n')
			fe = 0
	else:
		print(dataBlockingResultsFile + ' does not exists.  It will be created.\n')
		fe = 0
	
	if fe == 0:
		kk = kk + 1
		print("Analyzing file: " + tf + "\n")
		
		thermo_data_raw = np.genfromtxt(tf,names=True)
		condition = np.logical_and(eqData[:,0] == P,eqData[:,1] == T)
		eqSteps = np.compress(condition,eqData,axis=0)[0,2].astype(int)
		uncorrelatedBlockSize = np.compress(condition,eqData,axis=0)[0,3].astype(int)
	
		thermo_data_raw = thermo_data_raw[eqSteps/printInterval+1:]
		#for rawSampleNumber in range(np.size(thermo_data_raw)):
			#print("rawSampleNumber: " + str(rawSampleNumber) + " rawSample: " + str(thermo_data_raw[rawSampleNumber]))
		
		energy_tot_mean = np.mean(thermo_data_raw[:]['Energy'])
		#print("Total mean energy: " + str(energy_tot_mean));
		#summaryFile.write(str(jobID) + '\t' + str(P) + '\t' + str(T) + '\t' + str(energy_tot_mean) + '\n')
		#summaryFile.write(str(jobID) + '\t' + str(P) + '\t' + str(T) + '\n')
		stdevs = np.empty((0,4))
		for rawSampleNumber in range(np.size(thermo_data_raw) - 1):
			if(thermo_data_raw[rawSampleNumber + 1]['Step'] - thermo_data_raw[rawSampleNumber]['Step'] != printInterval):
				print("########## ERROR: INCONSISTENT STEP SPACING IN THERMO FILE ##############\n")
				sys.exit()
		if (uncorrelatedBlockSize < 1):
			blockSizeCounter = 1
			blockSize = blockSizeMin
		else:
			blockSizeCounter = 0
			blockSize = uncorrelatedBlockSize
		blockSizeMax = np.size(thermo_data_raw) * printInterval / 2
		while blockSize <= blockSizeMax :
			numBlocks = np.floor( np.size(thermo_data_raw) / (blockSize / printInterval) ).astype(int)
			#print("BlockSize: " + str(blockSize))
			#print("numBlocks: " + str(numBlocks))
			thermo_data_block = np.empty((0,7))
			for blockNum in range(numBlocks):
				#print('blockNum: ' + str(blockNum))
				step_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Step'].astype(np.uint32),axis=0)
				energy_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Energy'].astype(np.float),axis=0)
				esq_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Energy2'].astype(np.float),axis=0)
				length_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['l'].astype(np.float))
				lsq_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['l2'].astype(np.float),axis=0)
				virial_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Virial'].astype(np.float),axis=0)
				virsq_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Virial2'].astype(np.float),axis=0)
				#print("Mean Step: " + str(step_data_block) + "  Energy: " + str(energy_data_block) + "  Energy2: " + str(esq_data_block)
				#	 + "  Length: " + str(length_data_block) + "  Length2: " + str(lsq_data_block) + "  Virial: " + str(virial_data_block)
				#	 + "  Virial2: " + str(virsq_data_block))
				thermo_data_block_tmp = np.array([step_data_block,energy_data_block,esq_data_block,length_data_block,lsq_data_block,virial_data_block,virsq_data_block])
				thermo_data_block_tmp.shape = 1,7
				thermo_data_block = np.append(thermo_data_block,thermo_data_block_tmp,axis=0)
					
				#print(str(thermo_data_block))
				#print(str(thermo_data_block[:,0]))
				#print(str(thermo_data_block[:,1]))
			#std1 = np.power(np.mean(np.power(thermo_data_block[:,1],2)) - np.power(np.mean(thermo_data_block[:,1]),2),0.5)
			stdE = np.std(thermo_data_block[:,1])/np.sqrt(numBlocks)
			stdL = np.std(thermo_data_block[:,3])/np.sqrt(numBlocks)
			#print("Energy StDev: " + str(std1) + '  ' + str(stdE))
			#print('stdevs shape: ' + str(stdevs.shape))
			#tmp = [[blockSize,stdE,stdL]]
			#print('tmp size: ' + str(np.array(tmp).shape))
			stdevs = np.append(stdevs,[[blockSize,numBlocks,stdE,stdL]],axis=0)
			#stdevs = np.append(stdevs,tmp,axis=0)
			if (blockSizeCounter > 0):
				# blockSize = printInterval*np.ceil(1.3*blockSize/printInterval).astype(int)
				blockSize = (blockSize+printInterval)
			else:
				summaryFile.write(str(P) + '\t' + str(T) + '\t' + str(blockSize) + '\t' + str(numBlocks) + '\t' + str(stdE) + '\t' + str(stdL) + '\t' + str(stdE/np.sqrt(numBlocks)) + '\t' + str(stdL/np.sqrt(numBlocks)) + '\n')
				blockSize = blockSizeMin
			blockSizeCounter = blockSizeCounter + 1
		#print('Stdevs: ' + str(stdevs))
		
		print('Saving data blocking results to ' + dataBlockingResultsFile)
		np.savetxt(dataBlockingResultsFile,stdevs,header='blockSize\tnumBlocks\tStdErrE\tStdErrL',comments='',fmt=['%.8g','%.12g','%.10g','%.10g'],delimiter='\t')
		print('Finished gathering data blocking data for this simulation.\nGenerating plots for this simulation.\n')
		
	plt.scatter(stdevs[:,0],stdevs[:,2],linewidth=0.5)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.xlabel(r'Block Size',fontsize=18)
	plt.ylabel(r'StdErrOfMean(E)',fontsize=18)
	plt.title('P: ' + str(P) + ' T: ' + str(T),fontsize=16)
	# Check for existence of data blocking plots. In each case, skip re-making if it already exists.
	for ii in [('1E5',1E5),('5E5',5E5),('1E6',1E6),('5E6',5E6),('1E7',1E7),('5E7',5E7),('1E8',1E8)]:
		plotPath = thisDir + 'StdevE_P' + str(P) + '_T' + str(T) + '_' + ii[0] + '.png'
		if os.path.isfile(plotPath):
			print(plotPath + ' already exists.  Using the existing plot.\nTo create a fresh plot, first delete the existing one.\n')
		else:
			plt.xlim((0,ii[1]))
			#fig = plt.gcf()
			#fig.set_dpi(200)
			#plt.scatter(thermo_data_block[:,0],thermo_data_block[:,1],linewidth=0.8)
			#plt.title(A[jj],fontsize=18)
			#plt.xlabel(r'Step Number',fontsize=18)
			#plt.ylabel(r'System ',fontsize=18)
			# plt.xticks(fontsize=18)
			#plt.tick_params(axis='both', which='major', labelsize=18)
			#ax = plt.gca()
			#ax.yaxis.offsetText.set_fontsize(18)
			#ax.xaxis.offsetText.set_fontsize(18)
			#plt.tight_layout()
			#ax.text(0.05, 0.17, 'Block size = ' + '{:.2E}'.format(blockSize), transform=ax.transAxes, fontsize=10,
		        #	verticalalignment='top')
			#ax.text(0.05, 0.11, 'Mean = ' + '{:.2f}'.format(), transform=ax.transAxes, fontsize=10,
		       	#	verticalalignment='top')
			#ax.text(0.05, 0.05, 'StdErr = ' + '{:.2f}'.format(std), transform=ax.transAxes, fontsize=10,
		        #	verticalalignment='top')
			# plt.show()
			plt.savefig(plotPath,dpi='figure')
	print('Done generating data blocking plots for energy.\n')
	plt.close()
	
	plt.scatter(stdevs[:,0],stdevs[:,3],linewidth=0.5)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.xlabel(r'Block Size',fontsize=18)
	plt.ylabel(r'StdDev(L)',fontsize=18)
	plt.title('P: ' + str(P) + ' T: ' + str(T),fontsize=16)
	for ii in [('1E5',1E5),('5E5',5E5),('1E6',1E6),('5E6',5E6),('1E7',1E7),('5E7',5E7),('1E8',1E8)]:
		plotPath = thisDir + 'StdevL_P' + str(P) + '_T' + str(T) + '_' + ii[0] + '.png'
		if os.path.isfile(plotPath):
			print(plotPath + ' already exists.  Using the existing plot.\nTo create a fresh plot, first delete the existing one.\n')
		else:
			plt.xlim((0,ii[1]))
			#fig = plt.gcf()
			#fig.set_dpi(200)
			#plt.scatter(thermo_data_block[:,0],thermo_data_block[:,1],linewidth=0.8)
			#plt.title(A[jj],fontsize=18)
			#plt.xlabel(r'Step Number',fontsize=18)
			#plt.ylabel(r'System ',fontsize=18)
			# plt.xticks(fontsize=18)
			#plt.tick_params(axis='both', which='major', labelsize=18)
			#ax = plt.gca()
			#ax.yaxis.offsetText.set_fontsize(18)
			#ax.xaxis.offsetText.set_fontsize(18)
			#plt.tight_layout()
			#ax.text(0.05, 0.17, 'Block size = ' + '{:.2E}'.format(blockSize), transform=ax.transAxes, fontsize=10,
        		#	verticalalignment='top')
			#ax.text(0.05, 0.11, 'Mean = ' + '{:.2f}'.format(), transform=ax.transAxes, fontsize=10,
       			#	verticalalignment='top')
			#ax.text(0.05, 0.05, 'StdErr = ' + '{:.2f}'.format(std), transform=ax.transAxes, fontsize=10,
        		#	verticalalignment='top')
			# plt.show()
			plt.savefig(plotPath,dpi='figure')
	print('Done generating data blocking plots for length.\n')
	plt.close()
	
	print('Done with this simulation.\n\n')
	
summaryFile.close()
#summaryFile.close()
