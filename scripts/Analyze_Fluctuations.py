#!/gpfs_share/santiso/SOFTWARE/miniconda2/envs/santimulators_1/bin/python

#BSUB -W 00:59
#BSUB -o output.%J
#BSUB -e errors.%J
#BSUB -q debug
#BSUB -J AnalyzeSD

################################################################
#  Analyze_Fluctuations.py :
#    Creates plots of the Fluctuations of the energy and length
#      as well as the energy-length co-fluctuation from a
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
#
#
#  Workflow:
#
#
################################################################


import glob
import numpy as np
import sys
import matplotlib.pyplot as plt
import re
import os

blockSize = int(1000000)
printInterval = 10000


pathToDir = '../data/LJ/m-1/Longest/'
#thermoFiles = glob.glob(pathToDir + "*/*thermo*")
thermoFiles = glob.glob(pathToDir + "P0.1_*/*thermo*")
print("Found files: " + str(thermoFiles) + "\n")
kk = 0
for ii in thermoFiles:
	kk = kk + 1
	print("Analyzing file: " + ii + "\n")
	T = np.float(re.search('T(.*)_',ii).group(1))
	P = np.float(re.search('P(.*)_T',ii).group(1))
	print('P: ' + str(P) + '  T: ' + str(T))
	
	thermo_data_raw = np.genfromtxt(ii,names=True)
	#thermo_data_raw = np.genfromtxt(ii,names=["Step","Energy","Energy2","l","l2","Virial","Virial2","lEnergy"],skip_header=1,filling_values=0,usecols=(0,1,2,3,4,5,6))
        
        #condition = np.logical_and(eqData[:,0] == P,eqData[:,1] == T)
	#condition = np.logical_and(np.absolute(np.subtract(eqData[:,0],P)) < 0.0001  ,np.absolute(np.subtract(eqData[:,1] - T)) < 0.0001)
	#condition = np.logical_and(eqData[:,0] == P  ,eqData[:,1] == T)
	#print("Condition: " + str(condition))
	#eqSteps = np.compress(condition,eqData,axis=0)
	#eqSteps = np.compress(condition,eqData,axis=0)[0,2].astype(int)
	#eqSteps = np.compress(np.logical_and(eqData[:,0] == P,eqData[:,1] == T),eqData,axis=0)
	eqSteps = 1E7
        print('Equlibration Steps: ' + str(eqSteps))
        thermo_data_raw = thermo_data_raw[eqSteps/printInterval+1:]
	#for rawSampleNumber in range(np.size(thermo_data_raw)):
		#print("rawSampleNumber: " + str(rawSampleNumber) + " rawSample: " + str(thermo_data_raw[rawSampleNumber]))
	
	energySq_tot_mean     = np.mean(thermo_data_raw[:]['Energy2'])
	lSq_tot_mean           = np.mean(thermo_data_raw[:]['l2'])
	lTimesEnergy_tot_mean = np.mean(thermo_data_raw[:]['lEnergy'])
	print("Total mean energy^2: " + str(energySq_tot_mean));
	print("Total mean length^2: " + str(lSq_tot_mean));
	print("Total mean length x energy: " + str(lTimesEnergy_tot_mean));
	# summaryFile.write(str(P) + '\t' + str(T) + '\t' + str(energySq_tot_mean) + '\t' + str(l_tot_mean) + '\n')
	#summaryFile.write(str(jobID) + '\t' + str(P) + '\t' + str(T) + '\n')
	
	for rawSampleNumber in range(np.size(thermo_data_raw) - 1):
		if(thermo_data_raw[rawSampleNumber + 1]['Step'] - thermo_data_raw[rawSampleNumber]['Step'] != printInterval):
			print("########## ERROR: INCONSISTENT STEP SPACING IN THERMO FILE ##############\n")
			sys.exit()
	numBlocks = np.floor( np.size(thermo_data_raw) / (blockSize / printInterval) ).astype(int)
	print("numBlocks: " + str(numBlocks))
	thermo_data_block = np.empty((0,8))
	for blockNum in range(numBlocks):
		step_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Step'].astype(np.uint32),axis=0)
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
	
	#plt.scatter(thermo_data_block[:,0],thermo_data_block[:,1],linewidth=0.8)
	plt.scatter(thermo_data_block[:,0],thermo_data_block[:,1],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'Energy',fontsize=18)
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
	plt.savefig('MeanE_P' + str(P) + '_T' + str(T) + '.png',dpi='figure')
	plt.close()
	
	#plt.scatter(thermo_data_block[:,0],thermo_data_block[:,3],linewidth=0.8)
	plt.scatter(thermo_data_block[:,0],thermo_data_block[:,3],s=ms)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
	plt.xlabel(r'Step Number',fontsize=18)
	plt.ylabel(r'Length',fontsize=18)
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
	plt.savefig('MeanL' + str(P) + '_T' + str(T) + '.png',dpi='figure')
	plt.close()
#summaryFile.close()
#summaryFile.close()






##################################################################################
########################### Copied from Analyze_Mean.py ##########################
##################################################################################
#   
#   blockSize = int(1000000)
#   printInterval = 10000
#   
#   potentialType='LJ'
#   if potentialType == 'LJ':
#   	summaryFile = open("LJ_Means.dat","w")
#   	# Use the first option for command-line input of the target directory
#   	dir = "../data/LJ/m-1/Longest/"
#   elif potentialType == 'HARMONIC':
#   	summaryFile = open("HARMONIC_Means.dat","w")
#   	# Use the first option for command-line input of the target directory
#   	dir = "../data/HARMONIC/m1/"
#   summaryFile.write('P\tT\tEnergy\tLength\n')
#   
#   ms = 4  # marker size for plots
#   
#   # Use the first option for command-line input of the target directory
#   # dir = sys.argv[1]
#   
#   #eqSteps = int(sys.argv[2])
#   
#   #print("eqSteps: " + str(eqSteps))
#   print("printInterval: " + str(printInterval))
#   print("blockSize: " + str(blockSize))
#   #print("Analyzing runs in " + dir + " using " + str(eqSteps)  + " equilibration steps." + "\n")
#   print("Analyzing runs in " + dir + "\n")
#   
#   eqStepsFile = glob.glob(dir + "EqSteps.dat")[0]
#   print("Found equilibrium steps file: " + str(eqStepsFile))
#   eqData = np.genfromtxt(eqStepsFile)
#   #print("eqData: " + str(eqData))
#   
#   #thermoFiles = glob.glob(dir + "*/*thermo*")
#   thermoFiles = glob.glob(dir + "P0.1_*/*thermo*")
#   print("Found files: " + str(thermoFiles) + "\n")
#   kk = 0
#   for ii in thermoFiles:
#   	kk = kk + 1
#   	print("Analyzing file: " + ii + "\n")
#   	T = np.float(re.search('T(.*)_',ii).group(1))
#   	P = np.float(re.search('P(.*)_T',ii).group(1))
#   	#jobID = re.search('T.*_([0-9]*)\/',ii).group(1)
#   	#print('P: ' + str(P) + '  T: ' + str(T) + '  jobID: ' + str(jobID))
#   	print('P: ' + str(P) + '  T: ' + str(T))
#   	
#   	thermo_data_raw = np.genfromtxt(ii,names=True)
#   	#thermo_data_raw = np.genfromtxt(ii,names=["Step","Energy","Energy2","l","l2","Virial","Virial2","lEnergy"],skip_header=1,filling_values=0,usecols=(0,1,2,3,4,5,6))
#   	condition = np.logical_and(eqData[:,0] == P,eqData[:,1] == T)
#   	#condition = np.logical_and(np.absolute(np.subtract(eqData[:,0],P)) < 0.0001  ,np.absolute(np.subtract(eqData[:,1] - T)) < 0.0001)
#   	#condition = np.logical_and(eqData[:,0] == P  ,eqData[:,1] == T)
#   	#print("Condition: " + str(condition))
#   	#eqSteps = np.compress(condition,eqData,axis=0)
#    	eqSteps = np.compress(condition,eqData,axis=0)[0,2].astype(int)
#    	#eqSteps = np.compress(np.logical_and(eqData[:,0] == P,eqData[:,1] == T),eqData,axis=0)
#    	#print(eqSteps)
#    	#print('Equilibrium steps: ' + str(eqSteps))
#            thermo_data_raw = thermo_data_raw[eqSteps/printInterval+1:]
#    	#for rawSampleNumber in range(np.size(thermo_data_raw)):
#    		#print("rawSampleNumber: " + str(rawSampleNumber) + " rawSample: " + str(thermo_data_raw[rawSampleNumber]))
#    	
#    	energy_tot_mean = np.mean(thermo_data_raw[:]['Energy'])
#    	l_tot_mean      = np.mean(thermo_data_raw[:]['l'])
#    	print("Total mean energy: " + str(energy_tot_mean));
#    	print("Total mean length: " + str(l_tot_mean));
#    	summaryFile.write(str(P) + '\t' + str(T) + '\t' + str(energy_tot_mean) + '\t' + str(l_tot_mean) + '\n')
#    	#summaryFile.write(str(jobID) + '\t' + str(P) + '\t' + str(T) + '\n')
#    	
#    	for rawSampleNumber in range(np.size(thermo_data_raw) - 1):
#    		if(thermo_data_raw[rawSampleNumber + 1]['Step'] - thermo_data_raw[rawSampleNumber]['Step'] != printInterval):
#    			print("########## ERROR: INCONSISTENT STEP SPACING IN THERMO FILE ##############\n")
#    			sys.exit()
#    	numBlocks = np.floor( np.size(thermo_data_raw) / (blockSize / printInterval) ).astype(int)
#    	print("numBlocks: " + str(numBlocks))
#    	thermo_data_block = np.empty((0,8))
#    	for blockNum in range(numBlocks):
#    		step_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Step'].astype(np.uint32),axis=0)
#    		energy_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Energy'].astype(np.float),axis=0)
#    		esq_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Energy2'].astype(np.float),axis=0)
#    		length_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['l'].astype(np.float))
#    		lsq_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['l2'].astype(np.float),axis=0)
#    		virial_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Virial'].astype(np.float),axis=0)
#    		virsq_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['Virial2'].astype(np.float),axis=0)
#    		le_data_block = np.mean(thermo_data_raw[ blockNum*blockSize/printInterval : ((blockNum+1)*blockSize)/printInterval ]['lE'].astype(np.float),axis=0)
#    		#print("Mean Step: " + str(step_data_block) + "  Energy: " + str(energy_data_block) + "  Energy2: " + str(esq_data_block)
#    		#	 + "  Length: " + str(length_data_block) + "  Length2: " + str(lsq_data_block) + "  Virial: " + str(virial_data_block)
#    		#	 + "  Virial2: " + str(virsq_data_block))
#    		thermo_data_block_tmp = np.array([step_data_block,energy_data_block,esq_data_block,length_data_block,lsq_data_block,virial_data_block,virsq_data_block,le_data_block])
#    		thermo_data_block_tmp.shape = 1,8
#    		thermo_data_block = np.append(thermo_data_block,thermo_data_block_tmp,axis=0)
#    		
#    		#print(str(thermo_data_block))
#    		#print(str(thermo_data_block[:,0]))
#    		#print(str(thermo_data_block[:,1]))
#    	
#    	#plt.scatter(thermo_data_block[:,0],thermo_data_block[:,1],linewidth=0.8)
#    	plt.scatter(thermo_data_block[:,0],thermo_data_block[:,1],s=ms)
#    	fig = plt.gcf()
#    	fig.set_dpi(200)
#    	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
#    	plt.xlabel(r'Step Number',fontsize=18)
#    	plt.ylabel(r'Energy',fontsize=18)
#    #	# plt.xticks(fontsize=18)
#    #	#plt.tick_params(axis='both', which='major', labelsize=18)
#    #	#ax = plt.gca()
#    #	#ax.yaxis.offsetText.set_fontsize(18)
#    #	#ax.xaxis.offsetText.set_fontsize(18)
#    #	#plt.tight_layout()
#    #	#ax.text(0.05, 0.17, 'Block size = ' + '{:.2E}'.format(blockSize), transform=ax.transAxes, fontsize=10,
#    #        #	verticalalignment='top')
#    #	#ax.text(0.05, 0.11, 'Mean = ' + '{:.2f}'.format(), transform=ax.transAxes, fontsize=10,
#    #       	#	verticalalignment='top')
#    #	#ax.text(0.05, 0.05, 'StdErr = ' + '{:.2f}'.format(std), transform=ax.transAxes, fontsize=10,
#    #        #	verticalalignment='top')
#    #	# plt.show()
#    	plt.savefig('MeanE_P' + str(P) + '_T' + str(T) + '.png',dpi='figure')
#    	plt.close()
#    	
#    	#plt.scatter(thermo_data_block[:,0],thermo_data_block[:,3],linewidth=0.8)
#    	plt.scatter(thermo_data_block[:,0],thermo_data_block[:,3],s=ms)
#    	fig = plt.gcf()
#    	fig.set_dpi(200)
#    	plt.title('P: ' + str(P) + '  T: ' + str(T),fontsize=18)
#    	plt.xlabel(r'Step Number',fontsize=18)
#    	plt.ylabel(r'Length',fontsize=18)
#    #	# plt.xticks(fontsize=18)
#    #	#plt.tick_params(axis='both', which='major', labelsize=18)
#    #	#ax = plt.gca()
#    #	#ax.yaxis.offsetText.set_fontsize(18)
#    #	#ax.xaxis.offsetText.set_fontsize(18)
#    #	#plt.tight_layout()
#    #	#ax.text(0.05, 0.17, 'Block size = ' + '{:.2E}'.format(blockSize), transform=ax.transAxes, fontsize=10,
#    #        #	verticalalignment='top')
#    #	#ax.text(0.05, 0.11, 'Mean = ' + '{:.2f}'.format(), transform=ax.transAxes, fontsize=10,
#    #       	#	verticalalignment='top')
#    #	#ax.text(0.05, 0.05, 'StdErr = ' + '{:.2f}'.format(std), transform=ax.transAxes, fontsize=10,
#    #        #	verticalalignment='top')
#    #	# plt.show()
#    	plt.savefig('MeanL' + str(P) + '_T' + str(T) + '.png',dpi='figure')
#    	plt.close()
#    #summaryFile.close()
#    #summaryFile.close()

