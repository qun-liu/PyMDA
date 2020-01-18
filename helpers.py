# free to use or modify
#  Authors: 
#  Qun Liu: qun.liu@gmail.com
#  Lina Takamaru:  llt45@cornell.edu
#  References:  J. Appl. Cryst. (2020). 53 https://doi.org/10.1107/S160057671901673X
import math, re, os, sys, glob, pprint, argparse

## write for_sorting.tmp2 and data_with_frame_selected.txt files 
def createFs(dir_output, dataprefix='Puck', cc_cutoff=0.4):
	os.chdir(dir_output)
	#print "cctes", os.getcwd()
	check = [dir_output+"/for_sorting.tmp2", dir_output+"/data_with_frame_selected.txt"]
	for files in check:
		if os.path.exists(files):
			os.remove(files)

	#pFiles is a list of folders beginning with Puck
	#pFiles = glob.glob('../Puck*')
#	pFiles=[]
#	list_file = dataprefix.split(' ')
#	print("list_file", list_file)
#	for list_file in list(dataprefix):
#		pFiles.append(glob.glob(list_file+"*"))
	pFiles = glob.glob(dataprefix+"*")
	#print("directories", pFiles)
	wedges = []
	maxFiles = []
	output = open(dir_output+"/for_sorting.tmp2", "a")
	#glob returns a list, so make sure it is not empty before continuing
	if len(pFiles) != 0:
		#dictionary that contains the file name, wedges, max cc12
		cc12 = {}
		#dictionary that contains file names of max cc 
		filesC = {}
		for folder in pFiles:
			#items is a list of files in the current folder
			items = os.listdir(folder)
			maxcc = -1
			for files in items:
				if files.startswith('aimless'):
					with open(os.path.join(folder, files)) as f:
						content = f.read().strip()
					if 'Mn(I) half-set correlation CC(1/2)' in content:
						store = content.index('Mn(I) half-set correlation CC(1/2)')
						last = content.index('Completeness', store)
						newstr = content[store:last-1]
						lastp = newstr.rfind(')')
						newstr = newstr[lastp+1:].strip()
						space = newstr.find(' ')
						cc = float(newstr[:space+1].strip())
						if cc > maxcc:
							maxcc = cc
							dot = files.find('.')
							unders = files.rfind('_')
							wedge = int(files[unders+1:dot])
							wedges.append(wedge)
							filesC[folder] = [wedge, maxcc, files]
							cc12[folder]=[wedge, maxcc,files]
		x = pprint.pformat(cc12)
		output.write(x)
	output.close()
	##
	sortedfile = open(dir_output+"/data_with_frame_selected.txt", "a")
	#contains the files in addition to wedge and maxcc
	sorted2 = sorted(filesC.items(), key=lambda x: x[1][1], reverse=True)
	#cc12 sorted by cc1/2, then by file number
	cc12_sorted = sorted(cc12.items(), key=lambda x: (-x[1][1], x[0][-3:]))

	#cc12sorted only includes data with cc12 values greater than min cc cutoff
	cc12sorted_index=[d for d in range(len(cc12_sorted)) if float(cc12_sorted[d][1][1]) > float(cc_cutoff)]
	cc12sorted=[cc12_sorted[e] for e in cc12sorted_index]
	out = pprint.pformat(cc12sorted)
	sortedfile.write(out)
	sortedfile.close()
	return cc12sorted

## write unit_cell.txt and unit_cell.txt2 
def createAUC1(dir_output,dataprefix, cc_cutoff):
	cwd = os.getcwd()

	check=[dir_output+"/unit_cell.txt", dir_output+"/unit_cell.txt2",dir_output+"/check.txt"]
	for files in check:
		if os.path.exists(files):
			os.remove(files)
	
	auc = open(dir_output+"/unit_cell.txt", "a")
	aucRed = open(dir_output + "/unit_cell.txt2", "a")
	out = []

	x=0
	cc12sorted = createFs(dir_output,dataprefix, cc_cutoff)
	for items in cc12sorted:
		with open(cwd +'/' + cc12sorted[x][0] + "/" + cc12sorted[x][1][2]) as f:
			content = f.read()
			x = x+1
		if 'Average unit cell' in content:
			store = content.index('Average unit cell:')
			last = content.index('Selection',store)
			newstr = content[store:last-1]
			out.append(newstr[:-5])
	output = pprint.pformat(out)
	auc.write(output)
	result = []
	for data in out:
		new = data[21:]
		new = re.sub(' +', ' ', new)
		floats = map(float, new.split(' '))
		floats = filter(lambda a: a != 90.0, floats)
		floats = filter(lambda a: a != 120.0, floats)
		newstr = str(floats)
		result.append(newstr)
	numbers = '\n'.join(result)
	numbers = numbers.replace('[', '')
	numbers = numbers.replace(']', '')	
	aucRed.write(numbers)
	auc.close()
	aucRed.close()
	### combine unit cell with data
	combined = []
	combine = open(dir_output+"/check.txt", "a")
	#print list(cc12sorted), result
	for d in range(len(cc12sorted)):
		#print cc12sorted[d]
		combined.append([result[d],cc12sorted[d][0], cc12sorted[d][1][1]])
	#print combined
	combine.write(pprint.pformat(combined))
	combine.close()


## if only one data set, pass through 
def createAUC2(dir_output):
	cwd = os.getcwd()

	if os.path.exists(dir_output + '/cluster.txt'):
		os.remove(dir_output + '/cluster.txt')

	clust = open(dir_output+"/cluster.txt", "a")
	clust.write("cluster1q 1")
	clust.close()


##  write 1-picc for clustering
def createPicc(dir_output):
	if os.path.exists(dir_output+"/1-picc.txt"):
		os.remove(dir_output+"/1-picc.txt")

	#Create 1-picc.txt file to be used in picc rejection
	picc = open(dir_output+"/1-picc.txt", "a")
	with open(dir_output + "/aimless.log") as f:
		content = f.read()
		cbr = content.rfind("correlations by resolution")
		overall = content.find("Overall", cbr)
		zero = content.find("0", overall)
		equal = content.find("=", overall)
		newstr = content[zero:equal-2]
		picc.write(newstr)
		picc.close()

##### XTAL REJECTION #####
def createExtra(dir_output):
	cwd = os.getcwd()

	tmp = open(dir_output+"/extra.txt", "a")
	with open(cwd + "/aimless.log") as f:
		content = f.read()
		cm = content.find('Cumulative multiplicity')
		cc = content.find('Correlation coefficients',cm)
		fdollar = content.find('$$', cm)
		batch = content.find('Batch', fdollar)
		ln = content.rfind('$$', cm, cc)
		new = content[batch-4:ln]
		tmp.write(new)
		tmp.close()

def createSmr(dir_output):
	cwd = os.getcwd()
	smr = open(dir_output+"/SmRmerge.txt", "a")
	with open(dir_output + "/extra.txt") as f:
		ln = f.readlines()
		for x in ln[1:]:
			x = re.sub(' +', ' ', x)
			xnew = x.split(' ')
			xnew[15] = str(math.sqrt(float(xnew[15])*float(xnew[15])))
			x = " ".join(xnew)
			smr.write(str(x))
		smr.close()

def createBatch(dir_output):
	cwd = os.getcwd()
	batch = open(dir_output+"/batch.txt", "a")
	with open(cwd+"/pointless.log") as f:
		ln = f.read()
		s = ln.find(">*> Summary of test data read in:")
		e = ln.find("===", s)
		degrees = ln.rfind("degrees", s, e)
		ln = ln[s:degrees]
		x = re.sub(' +', ' ', ln)	
		x = x.split('\n')
		for item in x:
			if item.startswith(" Run number: "):
				start = item.index(":")
				c = item.index(" consists")
				first = item[start+2:c]
				cob = item.index("batches ")
				dash = item.index("- ")
				s = item[cob+8:dash]
				l = item[dash+2:]
				batch.write(first+' ')
				batch.write(s)
				batch.write(l+"\n")
		batch.close()

## dataset number, smr start and end frame
def createTmp3a(dir_output):
	cwd = os.getcwd()

	if os.path.exists(dir_output+"/tmp3_smr.txt"):
		os.remove(dir_output+"/tmp3_smr.txt")

	with open(dir_output+"/batch.txt") as f:
		content = f.readlines()
		se = {}
		for element in content:
			fspace = element.find(' ')
			sspace = element.rfind(' ')
			start = int(element[fspace+1:sspace])
			end = int(element[sspace+1:])
			se[start] = end
		se = sorted(se.iteritems())

	smr = open(dir_output+"/tmp3_smr.txt", "a")
	with open(dir_output+"/SmRmerge.txt") as f2:
		content2 = f2.readlines()
		y=1
		for x in range(len(se)):
			start = se[x][0]
			end = se[x][1]
			tmp = []
			for ln in content2:
				ln = re.sub(' +', ' ', ln)
				ln = ln.split(' ')
				tf = int(ln[2])>=start
				tf2 = int(ln[2]) <= end
				if tf&tf2 == True:
					tmp.append(float(ln[15]))
			if len(tmp)>0:
				avg = sum(tmp)/len(tmp)
				smr.write(str(y) + '    ' + str(avg) + '    ' + str(start) + '    ' +str(end)+'\n')
				y = y+1
		smr.close()

## sorted dataset number, smr start and end frame
def createTmp3b(dir_output):
	#sort SmRmerge numbers and add them to tmp3_smr2.txt
	if os.path.exists(dir_output+"/tmp3_smr2.txt"):
		os.remove(dir_output+"/tmp3_smr2.txt")

	smr = open(dir_output+"/tmp3_smr2.txt", "a")
	with open(dir_output+"/tmp3_smr.txt") as f:
		content = f.readlines()
		content2 = []
		for ln in content:
			ln = re.sub(' +', ' ', ln)
			ln = ln.split(' ')
			content2.append(ln)
		result = sorted(content2, key=lambda x: x[1])
		for s in result:
			x = ' '.join(s)
			smr.write(x)
		smr.close()
	return result
	
def createBdecay(dir_output):
	if os.path.exists(dir_output+"/Bdecay.txt"):
		os.remove(dir_output+"/Bdecay.txt")

	bd = open(dir_output+"/Bdecay.txt", "a")
	with open(dir_output+"/aimless.log") as aim:
		content = aim.read()
		rb = content.find("Relative Bfactor")
		ab = content.find("Agreement between batches")
		dr = content.rfind("$$", rb, ab)
		bf = content.find("    1", rb, dr)
		new = content[bf:dr]
		bd.write(new)
		bd.close()

def createSmrTmp4(path, xtal_steps, dir_output):	
	cc12sorted = createFs(path, dir_output)
	if os.path.exists(dir_output+"/tmp4.txt"):
		os.remove(dir_output+"/tmp4.txt")

	smr = open(dir_output + "/smr_included.txt", "a")
	tmp4 = open(dir_output+"/tmp4.txt", "a")
	with open(dir_output+"/tmp3_smr2.txt") as f:
		trackFiles = []
		trackFold = []
		content = f.readlines()
		total = len(content)
		content = content[:total-xtal_steps]
		for s in content:
			fspace = s.find(' ')
			smr.write(s[:fspace]+'\n')
			number = int(s[:fspace])-1
			trackFold.append(str(cc12sorted[number][0]))
			trackFiles.append(cc12sorted[number][1][2])
			tmp4.write('HKLIN ../' + str(cc12sorted[number][0])+ "/integrated" + cc12sorted[number][1][2][7:-4] + ".mtz" + "\n")
		smr.close()
		tmp4.close()
