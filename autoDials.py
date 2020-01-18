# free to use or modify
#  Authors: 
#  Qun Liu: qun.liu@gmail.com
#  Lina Takamaru:  llt45@cornell.edu
#  References:  J. Appl. Cryst. (2020). 53 https://doi.org/10.1107/S160057671901673X
import os, subprocess, argparse

def dImport(frame, images):
	if os.path.exists('datablock.json'):
		os.remove('datablock.json')
	#cmd = 'dials.import %s output.datablock=datablock.json scan.image_range=1,%i' %(images,frame)
	cmd = 'dials.import %s ' %(images)
	print cmd
	p=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()
	#out=p.communicate()
        #, stdout=subprocess.PIPE

def dFindspots(frame,dmin,minspot,sigb,sigs):
	#cmd = 'dials.find_spots datablock.json min_spot_size=%i xds.sigma_background=%i xds.sigma_strong=%i scan_range=1,%i filter.d_min=%f nproc=4' %(minspot,sigb,sigs,frame,dmin)
	cmd = 'dials.find_spots datablock.json scan_range=1,%i filter.d_min=%f nproc=6' %(frame,dmin)
	print cmd
	subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()

def dFindspots2(frame,dmin,minspot,sigb,sigs):
	if os.path.exists('strong.pickle'):
		os.remove('strong.pickle')
	cmd = 'dials.find_spots datablock.json min_spot_size=%i xds.sigma_background=%i xds.sigma_strong=%i scan_range=1,%i filter.d_min=%f nproc=6' %(minspot,sigb,sigs,frame,dmin)
	#cmd = 'dials.find_spots datablock.json scan_range=1,%i filter.d_min=%f nproc=4' %(frame,dmin)
	print cmd
	subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()

def dIndex(frame,spg,fix):
	if os.path.exists('experiments.json'):
		os.remove('experiments.json')
	if fix == True:
		detector = 'detector.fix=all'
	else:
		detector = ''
	cmd = 'dials.index datablock.json strong.pickle indexing.method=fft3d  refinement.parameterisation.auto_reduction.min_nref_per_parameter=4 scan_range=1,%i space_group=%s %s' %(frame,spg,detector)
	print cmd
	subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()
#indexing.method=real_space_grid_search  unit_cell="106,214,40,90,90,90"

def dIndex2(frame,spg, uc,fix):
	if os.path.exists('experiments.json'):
		os.remove('experiments.json')
	if fix == True:
		detector = 'detector.fix=all'
	else:
		detector = ''
	if uc == None:
		cmd = 'dials.index datablock.json strong.pickle'
		#print cmd
		#subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()
	else:
		cmd = 'dials.index datablock.json strong.pickle indexing.method=real_space_grid_search unit_cell="%s"   refinement.parameterisation.auto_reduction.min_nref_per_parameter=4 scan_range=1,%i space_group=%s %s' %(uc[0],frame,spg,detector)
	print cmd
	subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()

def dRefine():
	if os.path.exists('refined_experiments.json'):
		os.remove('refined_experiments.json')	
	cmd = 'dials.refine experiments.json indexed.pickle refinement.parameterisation.auto_reduction.min_nref_per_parameter=4 outlier.algorithm=tukey  output.log=refined.log'
	print cmd
	subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()
#scan_varying=True

def dIntegrate(dmin,frame):
	if os.path.exists('integrated.pickle'):
		os.remove('integrated.pickle')	
	cmd = 'dials.integrate refined_experiments.json refined.pickle '
	cmd += 'prediction.d_min=%f nproc=6  scan_range=1,%i' %(dmin,frame)
	#significance_filter.enable=True  significance_filter.isigi_cutoff=0.5
	print cmd
	subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()

def dExport():
	if os.path.exists('integrated.mtz'):
		os.remove('integrated.mtz')	
	cmd = 'dials.export integrated.pickle integrated_experiments.json mtz.hklout=integrated.mtz'
	subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()

def dPointless():
	if os.path.exists('sorted.mtz'):
		os.remove('sorted.mtz')	
	cmd = 'pointless hklin integrated.mtz hklout sorted.mtz <<eof > pointless.log\n'
	cmd += ' ALLOW OUTOFSEQUENCEFILES\n'
	cmd += ' copy\n'
	cmd += ' eof\n'
	print cmd
	subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()

def dPointless2(fls, folders, tolerance, dir_output):
	os.chdir(dir_output)
	check = [dir_output + "/sorted.mtz", dir_output + "/tmp2.txt"]
	for files in check:
		if os.path.exists(files):
			os.remove(files)

	tmp2 = open(dir_output+"/tmp2.txt", "a")
	pt ="pointless <<eof > pointless.log" + "\n"
	for x in range(len(folders)):
		folder = folders[x]
		files = fls[x]
		tmp2.write("HKLIN  " + folder + "/integrated_" + folder + "_" + str(files) + ".mtz\n")
		pt+= "HKLIN "  + folder + "/integrated_" + folder + "_" + str(files) + ".mtz\n"
	pt+=("\n" + "#########\nHKLOUT sorted.mtz\nNAME PROJECT tor CRYSTAL xtal DATASET data\nCOPY\nTOLERANCE " + str(tolerance) + "\neof")
	subprocess.Popen(pt,shell=True).wait()

def dAimless(res):
	check = ["scaled.mtz", "scala.mtz"]
	for files in check:
		if os.path.exists(files):
			os.remove(files)

	aim = "aimless hklin sorted.mtz hklout scaled.mtz <<eof > aimless.log\nbins 10\nanalysis groupbatch 1\nresolution 40 " + str(res) + "\nrefine parallel auto\nUSESDPARAMETER NO\nend\neof\n"
	#print aim
	subprocess.Popen(aim,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()


def dAimless20(res):
	check = ["scaled.mtz", "scala.mtz"]
	for files in check:
		if os.path.exists(files):
			os.remove(files)

	aim = "aimless hklin sorted.mtz hklout scaled.mtz <<eof > aimless.log\nbins 20\nanalysis groupbatch 1\nresolution 40 " + str(res) + "\nrefine parallel auto\nUSESDPARAMETER NO\nend\neof\n"
	#print aim
	subprocess.Popen(aim,shell=True).wait()
	
	
def dAimless2(res, dir_output):
	os.chdir(dir_output)
	aim = 'aimless hklin sorted.mtz hklout scaled.mtz <<eof > aimless2.log \n'
	with open(dir_output + "/tmpa.txt") as tmpa:
		lns = tmpa.readlines()
	for x in range(len(lns)):
		aim += str(lns[x][:-1])+'\n'
	aim += 'bins 20\nresolution 40\nanalysis groupbatch 1 ' + str(res)+'\nrefine parallel auto\nUSESDPARAMETER NO\nend\neof'
	#print aim
	subprocess.Popen(aim,shell=True).wait()	

### run truncate, mtzvarious, and freeR
def dProcesses(nres, dir_output):
	os.chdir(dir_output)
	if os.path.exists(dir_output + "/scala_unique.mtz"):
		os.remove(dir_output + "/scala_unique.mtz")

	ctrunc = " ctruncate -mtzin scaled.mtz -mtzout scala.mtz -nres " + str(nres) + " -colin '/*/*/[IMEAN,SIGIMEAN]' -colano '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]' > ctruncate.log"
	subprocess.Popen(ctrunc,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()
	
	mtzv = "mtz2various HKLIN scala.mtz HKLOUT scala.sca <<eof > mtz2various.log\nOUTPUT SCALEPACK\nlabin  I(+)=I(+) SIGI(+)=SIGI(+) I(-)=I(-) SIGI(-)=SIGI(-)\nend\neof"
	subprocess.Popen(mtzv,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()

	freerf = "freerflag HKLIN scala.mtz HKLOUT scala_unique.mtz <<eof > freerflag.log \nFREERFRAC 0.05 \nUNIQUE \neof"
	subprocess.Popen(freerf,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()

if __name__ == '__main__':
	dPointless()





