#!/gpfs_share/santiso/SOFTWARE/miniconda2/envs/santimulators_1/bin/python

import glob
import numpy as np
import sys
import matplotlib.pyplot as plt
import re

blockSizeMin = 100000
blockSizeMax = 1000000000

printInterval = 10000


dir = sys.argv[1]
eqSteps = int(sys.argv[2])

print("eqSteps: " + str(eqSteps))
print("printInterval: " + str(printInterval))
print("Analyzing runs in " + dir + " using " + str(eqSteps)  + " equilibration steps." + "\n")

thermoFiles = glob.glob(dir + "*/*thermo*")
print("Found files: " + str(thermoFiles) + "\n")
kk = 0

eqStepsFile = glob.glob(dir + "EqSteps.dat")[0]
print("Found equilibrium steps file: " + str(eqStepsFile))
eqData = np.genfromtxt(eqStepsFile)
#print("eqData: " + str(eqData))


### Use the following line to analyze a subset of the available thermo files ###
#thermoFiles = thermoFiles[0:1]
################################################################################

for tf in thermoFiles:
	kk = kk + 1
	print("Analyzing file: " + tf + "\n")
	T = np.float(re.search('T(.*)_',tf).group(1))
	P = np.float(re.search('P(.*)_T',tf).group(1))
	#jobID = re.search('T.*_([0-9]*)\/',tf).group(1)
	#print('P: ' + str(P) + '  T: ' + str(T) + '  jobID: ' + str(jobID))
	print('P: ' + str(P) + '  T: ' + str(T))
	
	thermo_data_raw = np.genfromtxt(tf,names=True)
	condition = np.logical_and(eqData[:,0] == P,eqData[:,1] == T)
	eqSteps = np.compress(condition,eqData,axis=0)[0,2].astype(int)
	
	thermo_data_raw = thermo_data_raw[eqSteps/printInterval+1:]
	#for rawSampleNumber in range(np.size(thermo_data_raw)):
		#print("rawSampleNumber: " + str(rawSampleNumber) + " rawSample: " + str(thermo_data_raw[rawSampleNumber]))
	
	energy_tot_mean = np.mean(thermo_data_raw[:]['Energy'])
	print("Total mean energy: " + str(energy_tot_mean));
	#summaryFile.write(str(jobID) + '\t' + str(P) + '\t' + str(T) + '\t' + str(energy_tot_mean) + '\n')
	#summaryFile.write(str(jobID) + '\t' + str(P) + '\t' + str(T) + '\n')
	stdevs = np.empty((0,3))
	for rawSampleNumber in range(np.size(thermo_data_raw) - 1):
		if(thermo_data_raw[rawSampleNumber + 1]['Step'] - thermo_data_raw[rawSampleNumber]['Step'] != printInterval):
			print("########## ERROR: INCONSISTENT STEP SPACING IN THERMO FILE ##############\n")
			sys.exit()
	blockSize = blockSizeMin
	blockSizeMax = np.size(thermo_data_raw) * printInterval / 2
	while blockSize <= blockSizeMax :
		numBlocks = np.floor( np.size(thermo_data_raw) / (blockSize / printInterval) ).astype(int)
		print("BlockSize: " + str(blockSize))
		print("numBlocks: " + str(numBlocks))
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
		stdE = np.std(thermo_data_block[:,1])
		stdL = np.std(thermo_data_block[:,3])
		#print("Energy StDev: " + str(std1) + '  ' + str(stdE))
		print('stdevs shape: ' + str(stdevs.shape))
		#tmp = [[blockSize,stdE,stdL]]
		#print('tmp size: ' + str(np.array(tmp).shape))
		stdevs = np.append(stdevs,[[blockSize,stdE,stdL]],axis=0)
		#stdevs = np.append(stdevs,tmp,axis=0)
		blockSize = printInterval*np.ceil(1.3*blockSize/printInterval).astype(int)
	print('Stdevs: ' + str(stdevs))
	
	
	plt.scatter(stdevs[:,0],stdevs[:,1],linewidth=0.8)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.xlabel(r'Block Size',fontsize=18)
	plt.ylabel(r'StdDev(E)',fontsize=18)
	#plt.xlim((0,5E6))
	plt.title('P: ' + str(P) + ' T: ' + str(T),fontsize=16)
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
	plt.savefig('StdevE_P' + str(P) + '_T' + str(T) + '.png',dpi='figure')
	plt.close()
	
	
	plt.scatter(stdevs[:,0],stdevs[:,2],linewidth=0.8)
	fig = plt.gcf()
	fig.set_dpi(200)
	plt.xlabel(r'Block Size',fontsize=18)
	plt.ylabel(r'StdDev(L)',fontsize=18)
	#plt.xlim((0,5E6))
	plt.title('P: ' + str(P) + ' T: ' + str(T),fontsize=16)
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
	plt.savefig('StdevL_P' + str(P) + '_T' + str(T) + '.png',dpi='figure')
	plt.close()
#summaryFile.close()
#summaryFile.close()
