#!/bin/env python 
# free to use or modify
#  Authors: 
#  Qun Liu: qun.liu@gmail.com
#  Lina Takamaru:  llt45@cornell.edu
#  References:  J. Appl. Cryst. (2020). 53 https://doi.org/10.1107/S160057671901673X

import os, sys, glob, pprint, collections, re, subprocess, shutil, math, argparse
import picc_rejection, uc_rejection, helpers, autoDials
from multiprocessing import Pool
from numpy import arange

"""
main usage of pymda 
 1. examples of processing hdf5 format Eiger data by using dials  
 /share/apps/pymda/pymda  --run_dials  --hdf5_path /share/d1/nsls2/ping/uxtal/FMX_20190905/BNL_136_13_C1  --wedges 10   --spg p2221
/share/apps/pymda/pymda  --run_dials  --hdf5_path /share/d1/nsls2/ping/uxtal/FMX_20190920  --wedges 10  --spg p2221 --thread 4

 2. examples of running microcrystal data assembly with 1 or 2 classes and with xtal and frame rejections
 /share/apps/pymda/pymda  --run_mda --dataprefix BNL --ucr 10 
 /share/apps/pymda/pymda  --run_mda --dataprefix BNL --ucr 10 --data sm_sep  --rjxtal --xtal_steps 10  --reso 5.0
 /share/apps/pymda/pymda  --run_mda --dataprefix BNL --ucr 1 --rjxtal --xtal_steps 5 --rjframe  --reso 3.0
 /share/apps/pymda/pymda  --run_mda --dataprefix BNL --ucr 2 --rjxtal --xtal_steps 5 --rjframe  --decay "4.0 3.0 2.0 1.0"
 /share/apps/pymda/pymda  --run_mda --run_picc --dataprefix BNL --ucr  1  --rjxtal --xtal_steps 5 --rjframe 
 /share/apps/pymda/pymda  --run_mda --dataprefix  BNL --data smUxtl --ucr 20  --rjxtal --xtal_steps 10 --reso 3 --rjframe --decay "5.0 3.0 2.0 1.0"
"""
def setOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument("--hdf5_path", dest="hdf5_path", metavar="DIRECTORY",default="blank", 
			help="Absolute path of raw hdf5 data directory, required to run_dials")
	parser.add_argument("--dataprefix", dest="dataprefix", metavar="DIRECTORY",default="test", 
			help="The prefix of directories that contain integrated data from run_dials, required for assembly")
	parser.add_argument("--dir_output", dest="dir_output", metavar="DIRECTORY",default=os.getcwd(), 
			help="the absolution path for output data directory, default is the current directory.")
	parser.add_argument("--data", dest="data", metavar="STRING", default='output', 
			help="use a meaningful prefix for assembled data sets")
	parser.add_argument("--frames", dest="frames", metavar="VALUE", default=10, type=int, 
			help="number of frames in a data set")
	parser.add_argument("--reso", dest="reso", metavar="VALUE", default=3.0, 
			help="resolution limit in angstrom for index and integration")
	parser.add_argument("--reso_cchalf", dest="reso2", metavar="VALUE", default=4.0,
			help="resolution limit in angstroms for pointless and aimless for statistics and wedges selection, should be ~1 angstrom lower than reso ")
	parser.add_argument("--wedges", dest='wedges', metavar='VALUE', default=10,type=int, 
			help='number of wedges to split a single data set')
	parser.add_argument("--spg", dest="spg", metavar='VALUE', default='p1', type=str,
			help="space group if known")
	parser.add_argument("--unit_cell", dest="unit_cell", metavar="STRING", default = None, nargs = 1, 
			help="six unit cell parameters,seperated by a comma if known")
	parser.add_argument("--minspot", dest="minspot", metavar="VALUE", default=4,type=int)
	parser.add_argument("--sigb", dest="sigb", metavar="VALUE", default=5,type=int)
	parser.add_argument("--sigs", dest="sigs", metavar="VALUE", default=2,type=int)
	parser.add_argument("--cc_cutoff", dest="cccutoff", metavar="VALUE", default=0.4,type=float )

	parser.add_argument("--decay",dest="decay", metavar="STRING", default="4.0 2.0 1.0",  type=str,
			help="decay is a list of space-seprated values to determine the conditions of frame rejection.")
	parser.add_argument("--picc", dest="picc",  metavar="INT", default=2,  type=int, 
			help="number of classes for classification based on 1-picc ")
	parser.add_argument("--run_picc", dest="run_picc", default=False, action='store_true',  
			help="where to use 1-picc for classification")
	parser.add_argument("--ucr", dest="ucr", default=2, metavar="VALUE", type=int,  
			help="number of clusters to be used for unit cell rejection")
	parser.add_argument("--tolerance",dest="tolerance",  metavar="VALUE",type=int, default=16,  
			help="unit-cell tolerance used for pointless")
	parser.add_argument("--xtal_steps",dest="xtal_steps",  metavar="VALUE", type=int, default=2,  
			help="number of steps/datasets for xtal rejection")
	parser.add_argument("--frame_steps", dest="frame_steps", metavar="VALUE", type=int, default=5,  
			help="number of steps used for frame rejection")
	parser.add_argument("--nres", dest="nres", metavar="VALUE", default=400,type=int,  
			help="Estimated total number of molecules in a.u.")
	parser.add_argument("--thread", dest="thread",  metavar="INT",type=int, default=1, 
			help="number of cores to be used for dials")
	parser.add_argument("--opt", dest="opt", default=False,  action='store_true', 
			help="enable grid search for spotfind in dials")
	parser.add_argument("--run_dials", dest="run_dials",default=False, action='store_true', 
			help="run dials for data processing, require to know the --path")
	parser.add_argument("--fix_detector", dest="fix_detector",default=False, action='store_true', 
			help="run dials indexing and refine with a fixed detector distance")
	parser.add_argument("--run_mda", dest="run_mda",default=False, action='store_true', 
			help="run microdiffraction data assembly, require to known the directores that contain data processed by dials")
	parser.add_argument("--rjxtal", dest="rjxtal", default=False, action='store_true', 
			help="run xtal rejection")
	parser.add_argument("--rjframe", dest="rjframe", default=False, action='store_true', 
			help="run frame rejection")
	parser.add_argument("--single", dest="single", default=False, action='store_true', 
			help="use single or ward for unit cell classification")
	args = parser.parse_args()
	return args

def f(x): #for multiprocessing only for dials 
	options = setOptions()
	wedges=options.wedges
	minspot = options.minspot
	sigb = options.sigb
	sigs = options.sigs
	reso1 = float(options.reso)
	#print("res:", reso1)
	reso2 = float(options.reso2)
	spg = options.spg
	uc=options.unit_cell
	frames = options.frames
	fix_detector = options.fix_detector
	dir_output = options.dir_output
	os.chdir(dir_output)
	#print "here we are before: ", os.getcwd()  
	if options.opt == True:
		print("optimizing spot finding parameters ...")
		maxopt = dOptimize(frame, reso1, spg, minspot,sigb,sigs)
		minspot = maxopt[1]
		sigb = maxopt[2]
		sigs = maxopt[3]
	root = x.split()[0]
	dirfile = x.split()[1]
	dname = dirfile
	os.chdir(dirfile)
	ipt = root+'/'+dirfile+'_master.h5'
	#print "ipt", ipt
##  import raw data  
	autoDials.dImport(frames, ipt)
	imt = open('dials.import.log', 'r')
	lines = imt.readlines()
	imt.close()	
	frameinfo = {}
	for line in lines:
		line = line.strip()
		if line[0:12] == 'num images: ':
			findC = line.find(':')
			frameinfo['frames'] = line[findC+2:]
			frames = int(line[findC+2:])
			print("data set:",  dirfile, "Total number of frames is " + str(frames))

##  index & integration 
	steps = int(math.ceil(float(frames)/wedges))
	if os.path.exists('datablock.json'):
		autoDials.dFindspots2(frames,reso1,minspot,sigb,sigs)
		for frame2 in range(steps, frames+2, steps):
			if frame2 > frames:
				frame2 = frames
			if os.path.exists('strong.pickle'):
				if uc == None:
					autoDials.dIndex(frame2,spg,fix_detector)
				else:
					autoDials.dIndex2(frame2,spg,uc,fix_detector)
				if os.path.exists('experiments.json'):	
					autoDials.dRefine()
					if os.path.exists('refined_experiments.json'):
						autoDials.dIntegrate(reso1,frame2)
						if os.path.exists('integrated.pickle'):
							autoDials.dExport()
							if os.path.exists('integrated.mtz'):
								autoDials.dPointless()
								autoDials.dAimless(reso2)
						if os.path.exists('integrated.mtz'):
							os.rename('integrated.mtz','integrated_%s_%i.mtz' %(dname,frame2))
						if os.path.exists('dials.integrate.log'):
							os.rename('dials.integrate.log','dials.integrate_%s_%i.log' %(dname,frame2))
						if os.path.exists('scaled.mtz'):
							os.rename('scaled.mtz','scaled_%s_%i.mtz' %(dname,frame2))
						if os.path.exists('aimless.log'):
							os.rename('aimless.log','aimless_%s_%i.log' %(dname,frame2))
	print("spotfind parameters:", minspot, sigb, sigs)
	os.chdir(dir_output)
	print("completed data set:", dirfile)
	
### optimize minispot, sigmab and sigmas for dials spotfinding 
def dOptimize(frame,dmin,spg, minspot, sigb, sigs):
	grid_list=[[]]
	for minspot in range(2,6,2):
		for sigb in range(3,7,2):
			for sigs in range(2,4): 
				#print minspot, sigb, sigs
				dFindspots2(frame,dmin,minspot,sigb,sigs)
				dIndex(frame, spg)
				dRefine()
				if os.path.exists('dials.integrate.log'):
					os.remove('dials.integrate.log')	
				dIntegrate(dmin,frame)
				if os.path.exists('dials.integrate.log'):
					iovers=subprocess.Popen("cat dials.integrate.log |awk '/i\/sigi \(profile fitting/ {print $5}'", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
					#iovers=random.uniform(0, 1)
					#print("grid search: ", minspot, sigb, sigs, iovers)				
					grid_ldatablock.json
	grid_list.append([iovers, minspot, sigb, sigs])
	grid_list2 = filter(None, grid_list)
	if os.path.exists('grid_list'):
                os.remove('grid_list')
	op1=open('grid_list','w')
	op1.write(''.join(str(x) for x in grid_list2))
	return max(grid_list2, key=lambda x: x[0])

############################################
if __name__ == '__main__':
	options = setOptions()
	hdf5_path = options.hdf5_path
	dataprefix = options.dataprefix
	dir_output = options.dir_output
	frames = options.frames
	reso = float(options.reso)
	reso2 = float(options.reso2)
	spg = options.spg
	cc_cutoff = float(options.cccutoff)
	dname = options.data
	decay = options.decay.split()
	picc =  options.picc
	run_picc = options.run_picc
	ucr = options.ucr
	tolerance = options.tolerance
	xtal_steps = options.xtal_steps
	frame_steps = options.frame_steps
	rjxtal = options.rjxtal
	rjframe = options.rjframe
	nres = options.nres
	thread = options.thread
	run_dials = options.run_dials
	run_mda = options.run_mda
	single = options.single
### run dials to process hdf5 data  
	if run_dials == True:
		#where raw data is
		os.chdir(dir_output)
		ipt_list = []
		for root, dirs, files in os.walk(hdf5_path):
			for filename in files:
				if filename.endswith('_master.h5') & ("_Raster_" not in filename):
					#directories.append(root)
					ipt_list.append(root+' '+filename[:-10])
					#ipt_list.append([root,filename[:-10]])
					if not os.path.exists(dir_output+'/'+filename[:-10]):
						os.makedirs(dir_output+'/'+filename[:-10])
		#print ipt_list

		#print("ipt")
		os.chdir(dir_output)
		#sys.exit()
		try:
			p = Pool(thread)
			p.map(f, ipt_list)
		finally:
			p.close()
	#sys.exit()
	

#### run mda for microdiffraction data assembly 
	if run_mda == True:
	
########## prepare necessary files for assembly ###########################

		# return to the working directory 
		os.chdir(dir_output)
		
		cc12sorted = helpers.createFs(dir_output, dataprefix, cc_cutoff)
		#print("cc1/2 sorted: ", cc12sorted)
		#count number of entries in list cc12sorted
		n_data = len(cc12sorted)

		#put average unit cell of each folder into unit_cell.txt
		#unit_cell.txt2 contains only numbers, excluding 90 or 120 in unit_cell.txt
		if n_data > 1:
			helpers.createAUC1(dir_output,dataprefix, cc_cutoff)
		else:
			helpers.createAUC2(dir_output)
		
		if ucr > 1: 
			clslist = uc_rejection.ucr(dir_output+"/unit_cell.txt2", ucr, dir_output,single)
		else: 
			clslist = uc_rejection.single(n_data, ucr, dir_output)
		#sys.exit()

	###### assemble data with a resolution cutoff
		r = str(reso)
		for clr in range(1,ucr+1):
			#leave this function in local variables to be used in the dials.dialsPointless function below
			if os.path.exists(dir_output + "/include.txt"):
				os.remove(dir_output + "/include.txt")
			include = open(dir_output+"/include.txt", "a")
			trackFold = []
			trackFiles = []
			clustnum = []
			for x in range(len(clslist)):
				if clslist[x][0] == clr:
					include.write(str(clslist[x][1]) + '\n')
					trackFold.append(cc12sorted[x][0])
					trackFiles.append(cc12sorted[x][1][0])
					#print("cccc:", clslist[x][1], cc12sorted[x][0], cc12sorted[x][1][0])
			include.close()
			autoDials.dPointless2(trackFiles, trackFold, tolerance, dir_output) #run on all files
			autoDials.dAimless20(r)
			autoDials.dProcesses(nres, dir_output)

			#Rename files for back up
			rename = dname+'_'+str(clr)+'_'+str(len(trackFold))+'_all_'+ r
			if os.path.exists(dir_output + '/scala_unique.mtz'):
				shutil.copy2(dir_output + '/scala_unique.mtz', rename +'_scala_unique.mtz')
			if os.path.exists(dir_output + '/aimless.log'):
				shutil.copy2(dir_output + '/aimless.log', rename +'_aimless.log')
			if os.path.exists(dir_output + '/scala.sca'):
				shutil.copy2(dir_output + '/scala.sca', rename + '_scala.sca')
			if os.path.exists(dir_output + '/tmp2.txt'):
				shutil.copy2(dir_output + '/tmp2.txt', rename + '_included.txt')
			if os.path.exists(dir_output + '/pointless.log'):
				shutil.copy2(dir_output + '/pointless.log', rename +'_pointless.log')
			if os.path.exists(dir_output + '/sorted.mtz'):
				shutil.copy2(dir_output + '/sorted.mtz', rename +'_sorted.mtz')
			print("finished with process of clsuter #", clr)
			
			#call picc_rejection and ask  for the number of classes
			if run_picc == True:
				helpers.createPicc(dir_output)
				picc_rejection.pr(dir_output + '/1-picc.txt', picc, dir_output)
				t_class = input("How many classes to be split, followed by [ENTER]: ")
				picc_rejection.pr(dir_output+'/1-picc.txt', int(t_class), dir_output)
				n_cluster = input("Input class number followed by [ENTER]: ")
				if int(n_cluster)>int(t_class):
					print("Cluster number is too big, exit")

				#leave this function in - local variables to be used in the dials.dialsPointless function below
				if os.path.exists(dir_output + "/include.txt"):
					os.remove(dir_output + "/include.txt")

				include = open(dir_output+"/include.txt", "a")
				with open(dir_output + "/cluster.txt") as f:
					content = f.readlines()
					trackFiles = []
					trackFold = []
					for b in content:
						bl = b.find("b")
						space = b.find(" ")
						filenumber = int(b[bl+2:])
						if b[bl+1:space] == str(n_cluster):
							print b[space+1:]
							trackFold.append(cc12sorted[filenumber-1][0])
							trackFiles.append(cc12sorted[filenumber-1][1][0])
							include.write(b[space+1:])
					include.close()
		
				autoDials.dPointless2(trackFiles, trackFold, tolerance, dir_output)
				autoDials.dAimless20(r)
				autoDials.dProcesses(nres, dir_output)

				#aimless.log, sorted.mtz, pointless.log are used for xtal rejection, tmp2.txt to be used to generate tmp3.txt for further rejection
				rename = dname+'_'+str(clr)+'_'+str(n_cluster)+'_picc_'+str(len(trackFold))+'_'+ r
				if os.path.exists(dir_output+'/scala_unique.mtz'):
					shutil.copy2(dir_output+'/scala_unique.mtz', rename +'_scala_unique.mtz')
				if os.path.exists(dir_output+'/aimless.log',):
					shutil.copy2(dir_output+'/aimless.log', rename+'_aimless.log')
				if os.path.exists(dir_output+'/scala.sca'):
					shutil.copy2(dir_output+'/scala.sca', rename+'_scala.sca')
				if os.path.exists(dir_output +'/tmp2.txt'):
					shutil.copy2(dir_output +'/tmp2.txt', rename+'_included.txt')
				print("finished with PICC rejection")

			##### XTAL REJECTION #####
			if rjxtal == True:
				first_time = True
				a=0
				while a==0:
					check = [dir_output+"/extra.txt", dir_output+"/batch.txt", dir_output+"/SmRmerge.txt", dir_output+"/tmp3.txt"]
					for files in check:
						if os.path.exists(files):
							os.remove(files)

					helpers.createExtra(dir_output)
					helpers.createSmr(dir_output)
					helpers.createBatch(dir_output)

					shutil.copyfile(dir_output+"/tmp2.txt", dir_output+"/tmp3.txt")
					#data = sum(1 for line in open('/home/lina/uxtal/combine/tmp3.txt'))
					helpers.createTmp3a(dir_output)
					Tmp3b = helpers.createTmp3b(dir_output)
					total = len(Tmp3b)
					
					#if not the first iteration, subtract xtal_steps and run rejection processes
					if first_time == False:
						#helpers.createSmrTmp4(path, xtal_steps, dir_output)
						if os.path.exists(dir_output+"/tmp4.txt"):
							os.remove(dir_output+"/tmp4.txt")

						#smr = open(dir_output + "/smr_included.txt", "a")
						tmp4 = open(dir_output+"/tmp4.txt", "a")

						trackFiles = []
						trackFold = []
						#total = total - xtal_steps
						Tmp3b = Tmp3b[:total-xtal_steps]
						for x in range(len(Tmp3b)):
							#smr.write(Tmp3b[x][0] + '\n')
							num = int(Tmp3b[x][0]) - 1 ### list starts from zero 
							print("tt2:",x, num, cc12sorted[num][0], cc12sorted[num][1][0])
							trackFold.append(cc12sorted[num][0])
							trackFiles.append(cc12sorted[num][1][0])
							tmp4.write('HKLIN ' + str(cc12sorted[num][0])+ "/integrated_" + str(cc12sorted[num][0]) + "_" + str(cc12sorted[num][1][0]) + ".mtz" + "\n")
						#smr.close()
						tmp4.close()
						autoDials.dPointless2(trackFiles, trackFold, tolerance, dir_output)
						autoDials.dAimless20(r)
						autoDials.dProcesses(nres, dir_output)
						try:
							n_cluster
						except NameError:
							n_cluster = None
						if run_picc == True:
							rename = dname+'_'+str(clr)+ '_' + str(n_cluster)+'_'+str(len(trackFold))+'_' + r
						else:
							rename = dname+'_'+str(clr)+ '_' +str(len(trackFold))+'_' + r
						if os.path.exists(dir_output+'/tmp4.txt'):
							shutil.copy2(dir_output+'/tmp4.txt', rename+'_included.txt')
						if os.path.exists(dir_output+'/scala_unique.mtz'):
							shutil.copy2(dir_output+'/scala_unique.mtz', rename+'_scala_unique.mtz')
						if os.path.exists(dir_output+'/scala.sca'):
							shutil.copy2(dir_output+'/scala.sca', rename + '_scala.sca')
						if os.path.exists(dir_output+'/aimless.log'):
							shutil.copy2(dir_output+'/aimless.log', rename+'_aimless.log') #aimless.log to be used for next round of SmRmerge calculation
						if os.path.exists(dir_output+'/tmp3_smr.txt'):
							shutil.copy2(dir_output+'/tmp3_smr.txt', rename+'_smr.txt')
						if os.path.exists(dir_output+'/sorted.mtz'):
							shutil.copy2(dir_output+'/sorted.mtz', rename+'_sorted.mtz') #sorted.mtz to be used for frame rejection
						if os.path.exists(dir_output+'/tmp4.txt'):
							shutil.copy2(dir_output+'/tmp4.txt', 'tmp2.txt') #tmp2.txt to be used to produce tmp3.txt for iterative xtal rejection

					else: ## if this is the first time run, add everything to the file and loop back to while
						if os.path.exists(dir_output+"/tmp4.txt"):
							os.remove(dir_output+"/tmp4.txt")
						tmp4 = open(dir_output+"/tmp4.txt", "a")
						trackFiles = []
						trackFold = []
						#total = len(Tmp3b)
						#Tmp3b = Tmp3b[:total-xtal_steps]
						#print("Tmp3b", Tmp3b, len(Tmp3b))
						for x in range(len(Tmp3b)):
							#smr.write(Tmp3b[x][0] + '\n')
							num = int(Tmp3b[x][0]) - 1  ## list starts from zero 
							print("tt:",x, num, cc12sorted[num][0], cc12sorted[num][1][0])
							trackFold.append(cc12sorted[num][0])
							trackFiles.append(cc12sorted[num][1][0])
							tmp4.write('HKLIN ' + str(cc12sorted[num][0])+ "/integrated_" + str(cc12sorted[num][0]) + "_" + str(cc12sorted[num][1][0]) + ".mtz" + "\n")
						tmp4.close()
						first_time = False 

					##### FRAME REJECTION #####
					if rjframe == True:
						os.chdir(dir_output)
						#print("rjframe", rjframe) 
						if rjxtal == True :
							check = [dir_output + "/SmRmerge.txt",  dir_output + "/batch.txt", dir_output + "/extra.txt"]
							for files in check:
								if os.path.exists(files):
									os.remove(files)
							helpers.createExtra(dir_output)
							helpers.createSmr(dir_output)
							helpers.createBatch(dir_output)

						data = sum(1 for line in open(dir_output+'/tmp2.txt'))
						#print("rjframe", rjframe) 
						#for d in arange(decay, 0.1, -0.5):
						for d in list(decay):
							check = [dir_output+"/tmpa.txt", dir_output+"/smax.txt", dir_output+"/inc.txt"]
							for files in check:
								if os.path.exists(files):
									os.remove(files)
									#print(files, "is removed.")

							with open(dir_output+"/batch.txt") as bt:
								ln = bt.readlines()

								for batchSE in ln:
									inc = open(dir_output+"/inc.txt", "a")
									smaxfile = open(dir_output+"/smax.txt", "a")
									inc.seek(0)
									inc.truncate()
									x = re.sub(' +', ' ', batchSE)
									xnew = x.split(' ')

									batch_start = int(xnew[1])
									batch_end = int(xnew[2])
									#print(".....", batch_start, batch_end)
									with open(dir_output+"/SmRmerge.txt") as smr:
										ln2 = smr.readlines()
										values = []
										### find minimum smr for each data
										for y in ln2:
											y = re.sub(' +', ' ', y)
											ynew = y.split(' ')
											cut = int(ynew[2])
											sm = float(ynew[13])
											if cut <= batch_end and cut >= batch_start:
												values.append(sm) ### extract frames to get smr value 
										#print("values", values)
										smin = min(values) ### find minimum smr value
										smax = smin*(1+float(d))  ### get percentage 
										smaxfile.write(str(smin)+' '+str(smax)+'\n')
										for z in ln2:
											z = re.sub(' +', ' ', z)
											znew = z.split(' ')
											cut = int(znew[2])
											if cut <= batch_end and cut >= batch_start:
												val = float(znew[13])
												if val<smax:
													#inc = open(dir_output+"/inc.txt", "a")
													inc.write(znew[2]+'\n')
													#print("write out", znew[2])
									inc.write(str(batch_start)+'\n')
									inc.write(str(batch_end)+'\n')
									inc.close()
									smaxfile.close()
									#sys.exit()
									
									tmpa = open(dir_output+"/tmpa.txt", "a")
									with open(dir_output+"/inc.txt") as f:
										ln = f.readlines()
										sortThis = [line.strip() for line in ln]
										sortThis = sorted(map(int, sortThis))
						#				tmpa.seek(0)
						#				tmpa.truncate()
										#print("sortThis", sortThis)
										for i in range(len(sortThis)-1):
											if sortThis[i+1]-sortThis[i]>frame_steps:
												#tmpa.write("test")
												tmpa.write("exclude batch " + str(sortThis[i]+1) + " to " + str(sortThis[i+1]-1)+'\n')
										sortThis = map(str, sortThis)
										x = '\n'.join(sortThis)
									tmpa.close()

											
							autoDials.dAimless2(r, dir_output)
							autoDials.dProcesses(nres, dir_output)

							if run_picc == True:
								rename = dname+'_'+str(clr) + '_' + str(n_cluster)+'_'+str(len(trackFold))+'_'+str(d)+'_'+r
							else:
								rename = dname+'_'+str(clr) + '_' + str(len(trackFold))+'_'+str(d)+'_'+ r 
				
							if os.path.exists(dir_output+'/tmpa.txt'):
								shutil.copy2(dir_output+'/tmpa.txt', rename + '_tmpa.log')
							if os.path.exists(dir_output+'/smax.txt'):
								shutil.copy2(dir_output+'/smax.txt', rename + '_smax.log')
							if os.path.exists(dir_output+'/aimless2.log'):
								shutil.copy2(dir_output+'/aimless2.log', rename + '_aimless.log')
							if os.path.exists(dir_output+'/scala.sca'):
								shutil.copy2(dir_output+'/scala.sca', rename + '_scala.sca')
							if os.path.exists(dir_output+'/scala_unique.mtz'):
								shutil.copy2(dir_output+'/scala_unique.mtz', rename + '_scala_unique.mtz')

					#reset the value of a based on steps and files
					#print "total:", len(Tmp3b) 
					if len(Tmp3b) > xtal_steps:
						#total = total - xtal_steps
						a=0
					else:
						a=1
	## output sorted results
	cmd = 'fgrep "DelAnom" *_%s_aimless.log |sort -k 6 -n |tee DelAnom_sorted.log' %(reso)
	#subprocess.Popen(cmd,shell=True).wait()
	cmd = 'fgrep "Mn(I) half" *_%s_aimless.log |sort -k 7 -n |tee CC12_sorted.log' %(reso)
	#subprocess.Popen(cmd,shell=True).wait()
