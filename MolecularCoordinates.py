"""
A module for working with the three-dimensional molecular coordinates.

In particular, this allows:
1.representation of 3D-structures as sets of xyz coordinates
2.conversion of 3D coordinates to distance matrices
3.comparison of conformations of the same molecule
"""


################################################################################

import math

class MolecularGeometry:
	"""
	A class for representing 3D-molecular structures:
	
	========= ==============================================================
	Attribute Meaning
	========= ==============================================================
	`atomVec`: a vector of atomic numbers corresponding to the xyz coordinates in xyzCoor
	`xyzCoor`: a matrix of x (column 1), y (column 2) and z (column 3) coordinates for atoms within a molecule
	'atomLabels': (optional to specify) a vector of unique integer labels from 1 to natoms corresponding to the xyz coordinates in xyzCoor; in the absence of input, it will be assumed that the labels are 1 for the first atom, 2 for the second atom, etc.
	'atoms':   number of atoms in the molecule (based on the length of atomVec)
	'coordDict' and 'atomTypeDict': dictionaries mapping atomLabel to coordinates and atom types, respectively
	'atomTypeRange': the different possible values in atomVec
	========= ==============================================================
	"""

	#initialization with connectivity information
	def __init__(self, atomVec, xyzCoor, connectivity=None, atomLabels=None):
		self.xyzCoor=xyzCoor
		self.atomVec=atomVec
		self.connectivity=connectivity or False #connectivity will be false if it is not provided
		self.atoms = len(atomVec)
		self.atomLabels = atomLabels or range(1, self.atoms+1)
		#assert rows of xyzCoor = length of atomlabels = length of atomVec = atoms
		#assert columns of xyzCoor = 3
		#assert atomLabels contains all integers 1 through n

		#map the labels to xyzCoordinates and atom types
		self.coordDict= {}
		self.atomTypeDict = {}
		self.atomTypeRange = []
		for i in range(self.atoms):
		    self.coordDict[self.atomLabels[i]]=self.xyzCoor[i][0:3]
		    self.atomTypeDict[self.atomLabels[i]]=self.atomVec[i]

		#find atomTypeRange
		self.atomTypeRange = set(atomVec)

	def getDistanceMatrix(self):
		"""
		Constructs a (symmetric) matrix representing distances between atoms

		Atoms are indexed by atomLabels - 1
		"""
		natoms=self.atoms
		dist = [[0.0 for i in range(natoms)] for j in range(natoms)]
		for i in range(0,natoms):
		    #dist[i][i]=0.0 #diagonal elements are zero
		    for j in range(i+1, natoms):
			icoord=self.coordDict[i+1]
			jcoord=self.coordDict[j+1]
			xdiff=icoord[0]-jcoord[0]
			ydiff=icoord[1]-jcoord[1]
			zdiff=icoord[2]-jcoord[2]
			dist[i][j]=math.sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff)
			dist[j][i]=dist[i][j] #matrix is symmetric

		return dist

	def getDistanceMappings(self):
		"""
		Generates dictionaries mapping atom label pairs to associated distances

		Outputs are heterogeneous mappings (different atom types in tuple) and homogeneous mappings (same atom types) to distances and mappings to atomTypeTuples
		For heterogeneous mappings, the tuples will be ordered with lowest # atom type first
		For homogeneous mappings, the tuples will be ordered with lowest # atom label first
		Example:
		>>> import MolecularCoordinates
		>>> a = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
		>>> b = a.getDistanceMappings()
		>>> b
		({(3, 1): 1.0, (4, 1): 1.0, (2, 1): 1.0}, {(3, 4): 1.4142135623730951, (2, 3): 1.4142135623730951, (2, 4): 1.4142135623730951}, {(3, 1): (1, 6), (4, 1): (1, 6), (2, 1): (1, 6)}, {(3, 4): (1, 1), (2, 3): (1, 1), (2, 4): (1, 1)})
		"""
		dist = self.getDistanceMatrix()
		natoms = self.atoms
		hetMap = {}
		homMap = {}
		hetMapType={}
		homMapType={}
		for i in range(0,natoms):
		    atomType1 = self.atomTypeDict[i+1]
		    for j in range(i+1, natoms):
			atomType2 = self.atomTypeDict[j+1]
			if (atomType1 == atomType2): #if the atom types are the same, put the distance in homMap
			    homMap[(i+1,j+1)]=dist[i][j]
			    homMapType[(i+1,j+1)]=(atomType1,atomType2)
			elif (atomType1 < atomType2):#otherwise, put the distance in hetMap, with order of tuple determined by lowest # atom type
			    hetMap[(i+1,j+1)]=dist[i][j]
			    hetMapType[(i+1,j+1)]=(atomType1,atomType2)
			else:
			    hetMap[(j+1,i+1)]=dist[i][j]
			    hetMapType[(j+1,i+1)]=(atomType2,atomType1)

		#an alternative: dictionary where each mapType maps to a dictionary, which in turn, maps labels to distances

		return hetMap, homMap, hetMapType, homMapType

	def writeMM4File(self, filename, moleculename):
		"""
		Writes MM4 file

		The molecule must have connectivity information or this will not work properly; first 60 characters of moleculename will be used; file will be written to path specified by filename
		"""
		#determine attached vs. connected atoms
		b = len(self.connectivity)#the number of bonds
		counter = [0 for i in range(self.atoms)]
		for i in range(b):#count the number of times each atom appears
			counter[self.connectivity[i][0]-1]=counter[self.connectivity[i][0]-1]+1
			counter[self.connectivity[i][1]-1]=counter[self.connectivity[i][1]-1]+1
		#now, counter should be filled with positive integers and "attached" atoms will have only 1 in their corresponding location in counter list
		#count and identify attached cases
		attachedAtoms=[]
		for i in range(self.atoms):
			if counter[i]==1:
				attachedAtoms.append(i+1)
		nattach=len(attachedAtoms)
		ncon=b-nattach #the remainder are connected bonds, which will be written one at a time

		str1a = '%-60s3%4d 0  0 0  0%5s'%(moleculename[0:60], self.atoms, "2.0")#first portion of first line (left-aligned molecule name truncated to 60 characters, option 3 for automatic SCF calcs, number of atoms (right justified), IPRINT=0, MDERIV (blank), NRSTR=0, INIT=0, NCONST=0, TMAX=2.0 minutes
		str1b = moleculename + '\n' #print the full molecule name after the first 80 characters of the first line
		str1=str1a+str1b
		str2 = ' %4d                    %5d                                  1              0\n'%(ncon,nattach)#print second line: KFIXTYP omitted (assumed zero), a=0, NCON (number of connected atom lists), DEL, ISPEED (omitted), NATTACH (number of attached atom (each consists of a pair of integers), ISTYPE and LABEL and NDC and NCALC are skipped, HFORM=1 (heat of formation calculated), MVDW is skipped, NDRIVE=0 (for now; later this will probably be set at 1)
		#build the strings for the third and fourth blocks (connected atom connectivity and attached atom connectivity, respectively)
		str3=''
		str4=''
		attachedCounter=0
		for i in range(b):
			if(self.connectivity[i][0] not in attachedAtoms and self.connectivity[i][1] not in attachedAtoms):#connected atom case
				str3=str3+'%5d%5d\n'%(self.connectivity[i][0],self.connectivity[i][1])
			else:#attached atom case
				attachedCounter=attachedCounter+1
				if(self.connectivity[i][1] in attachedAtoms):#the first attached atom case
					str4=str4+'%5d%5d'%(self.connectivity[i][0],self.connectivity[i][1])
				else:#(self.connectivity[i][0] in attachedAtoms); the other attached case
					str4=str4+'%5d%5d'%(self.connectivity[i][1],self.connectivity[i][0])
				if(attachedCounter%8==0 or attachedCounter==nattach):#put a new line character when we are at the end of the line (a multiple of 8) or if we have reached the end of the attached atoms
					str4=str4+'\n'
		#build the string for the fifth block (coordinates and atom types)
		str5=''
		for i in range(self.atoms):
			label=self.atomLabels[i]
			coord=self.coordDict[label]
			atomtyp=self.atomTypeDict[label]
			str5=str5+'%10.5f%10.5f%10.5f%2s%3d(%3d)\n'%(coord[0],coord[1],coord[2],atomicNumberToSymbol(atomtyp),atomicNumberToMM4Type(atomtyp),label)
		
		#write the file
		f = open(filename, 'w')
		f.write(str1+str2+str3+str4+str5)
		f.close()

		return

	def writeMOLFile(self, filename, moleculename):
		"""
		Writes MOL file (without connectivity)

		will not necessarily be a valid MOL file (with RAD tags, etc.); file will be written to path specified by filename
		"""
		
		#write the file
		f = open(filename, 'w')
		#write the header
		f.write("\n")
		f.write(moleculename+"\n")
		f.write("\n")
		#write the first line after the header
		f.write(str(self.atoms)+ '  0  0  0  0  0  0  0  0  0  0 V2000\n')
		#write the coordinates
		for i in range(self.atoms):
			label=self.atomLabels[i]
			coord=self.coordDict[label]
			atomtyp=self.atomTypeDict[label]
			f.write('%9.4f %9.4f %9.4f %2s  0  0  0  0  0  0  0  0  0  0  0  0\n'%(coord[0],coord[1],coord[2],atomicNumberToSymbol(atomtyp)))
		f.write('M  END\n')
		f.close()

		return


################################################################################
def checkConformationalEquivalence(mg1, mg2, Atol=-1, Rtol=-1):
	"""Checks conformational equivalence between two molecules within a user specified tolerance

	mg1 and mg2 are the two conformers to be checked for equivalence
	note this this does not distinguish between mirror images
	At least one of Atol and Rtol must be specified (for now, we only treat the Atol case); Rtol support will be added later
	output: 1) a boolean indicating whether the two structures are identical within the desired tolerance
		2) a list containing distinct atom mappings giving equivalent conformations within the specified tolerance
	"""
	#assert mg1 and mg2 contain the same number of each type of atom (and obviously the same total number of atoms)
	#assert Atol > 0 or Rtol > 0 (and maybe also that Atol < 0 or Rtol < 0; i.e. only one of Atol and Rtol should be specified)

	#generate distance mappings
	(hetMap1, homMap1, hetMapType1, homMapType1)=mg1.getDistanceMappings()
	(hetMap2, homMap2, hetMapType2, homMapType2)=mg2.getDistanceMappings()

	atomMap = {} #a dictionary to keep track of the currently defined atom mappings as the algorithm progresses
	successfulAtomMapList = [] #a "pseudo-global" list variable to keep track of all successful atomMaps identified

	(matchQ, matches) = checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, mg1.atomVec, mg2.atomVec, atomMap, successfulAtomMapList, Atol=Atol, Rtol=Rtol)
	
	return matchQ, matches

def checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomType1, atomType2, atomMap, successfulAtomMapList, Atol=-1, Rtol=-1):
	"""Recursive function to assign mappings between molecule 1 and molecule 2 based on distances


	"""
	#initialize boolean variable to False
	successfulMatchQ = False

	#use the hetMap, if possible; if not, use homMap
	if(len(hetMap1)>0):
		#for the first iteration, when atomMap is empty, pull an arbitrary item from hetMap1; on subsequent iterations, make sure one (it should be exactly one by the design of the algorithm) of the elements is in atom map
		if(len(atomMap)==0):#if this is the first iteration:
			mapping = hetMap1.popitem()
			hetMap2TargetLabel = -1 #-1 corresponds to first iteration, where we don't have a particular target atom
		else:
			for i in range(1,len(atomType1)+1):#iterate over all the atom labels, looking for something new (i.e. not in the atomMap yet)
				if(i not in atomMap):
					label1new=i
					break
			for i in range(1,len(atomType1)+1):#iterate over all the atom labels, looking for something old (i.e. a mapping has already been found) and with a different atom type
				if(i in atomMap and (atomType1[label1new-1] != atomType1[i-1])):#look for an "old" atom with a different atom type
					label1old=i
					break
			if (atomType1[label1new-1] < atomType1[label1old-1]):#use canonical ordering
				labels1=(label1new,label1old)
			else:
				labels1=(label1old,label1new)
			hetMap2TargetLabel = atomMap[label1old]
			mapping =  (labels1, hetMap1.pop(labels1))
					
		mappingLabels = mapping[0]
		mappingType = hetMapType1[mappingLabels]
		for a in range(1,len(atomType1)+1):#search in hetMap2 for cases with the same mapping type and a pre-established atom correspondence
			if(hetMap2TargetLabel==-1):#for the first iteration, we should check all atom pairs
				for b in range(a+1,len(atomType1)+1):
					if (atomType2[a-1] < atomType2[b-1]):#use canonical ordering
						i=(a,b)
					else:
						i=(b,a)#determine canonical ordering
					if(i in hetMap2 and hetMapType2[i]==mappingType):#this line does the mapping type check
						if(distanceMatchQ(mapping[1],hetMap2[i], Atol=Atol, Rtol=Rtol)): #for each case where this is true, check that all other het and hom mappings involving already identified atoms also satisfy the constriant, removing them in the process
							successfulMatchQ = hetBlock(i, hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, hetMapType1, homMapType1, hetMapType2, homMapType2, atomMap, successfulAtomMapList, successfulMatchQ, Atol, Rtol)
			else:#on other iterations
				if (atomType2[a-1] < atomType2[hetMap2TargetLabel-1]):#use canonical ordering
					i=(a,hetMap2TargetLabel)
				else:
					i=(hetMap2TargetLabel,a)
				if(i in hetMap2 and hetMapType2[i]==mappingType):#this line does the mapping type check
					if(distanceMatchQ(mapping[1],hetMap2[i], Atol=Atol, Rtol=Rtol)): #for each case where this is true, check that all other het and hom mappings involving already identified atoms also satisfy the constriant, removing them in the process
						successfulMatchQ = hetBlock(i, hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, hetMapType1, homMapType1, hetMapType2, homMapType2, atomMap, successfulAtomMapList, successfulMatchQ, Atol, Rtol)
		#if(not successfulMatchQ):
		#	print atomMap
		#restore the initial popped out mapping
		hetMap1[mapping[0]]=mapping[1]
		return successfulMatchQ, successfulAtomMapList
	elif(len(homMap1)>0):
		#this should only be encountered for molecules like fullerene, graphene, hydrogen, etc. with only one type of atom
		#similar to hetMap case above(differences indicated by +++), except there are two possible mappings on first assignment
		#for the first iteration, when atomMap is empty, pull an arbitrary item from hetMap1; on subsequent iterations, make sure one (it should be exactly one by the design of the algorithm) of the elements is in atom map
		if(len(atomMap)==0):#if this is the first iteration:
			mapping = homMap1.popitem()
			homMap2TargetLabel = -1 #-1 corresponds to first iteration, where we don't have a particular target atom
		else:
			for i in range(1,len(atomType1)+1):#iterate over all the atom labels, looking for something new
				if(i not in atomMap):
					label1new=i
					break
			for i in range(1,len(atomType1)+1):#iterate over all the atom labels, looking for something old
				if(i in atomMap):#look for an "old" atom (by structure of algorithm, if we get to this block, all atoms will have the same atom type, so we don't need to check that it is the same)
					label1old=i
					break
			if (label1new < label1old):#use canonical ordering
				labels1=(label1new,label1old)
			else:
				labels1=(label1old,label1new)
			homMap2TargetLabel = atomMap[label1old]
			mapping =  (labels1, homMap1.pop(labels1))
		mappingLabels = mapping[0]
		mappingType = homMapType1[mappingLabels]
		for a in range(1,len(atomType1)+1):#search in homMap2 for cases with the appropriate target atoms
			if(homMap2TargetLabel==-1):#for the first iteration, we should check all atom pairs
				for b in range(a+1,len(atomType1)+1):
					if (a < b):#use canonical ordering
						i=(a,b)
					else:
						i=(b,a)#determine canonical ordering
					if(i in homMap2):#this line does the mapping type check
						if(distanceMatchQ(mapping[1],homMap2[i], Atol=Atol, Rtol=Rtol)): #for each case where this is true, check that all other het and hom mappings involving already identified atoms also satisfy the constriant, removing them in the process
							successfulMatchQ = homBlock(i, hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, hetMapType1, homMapType1, hetMapType2, homMapType2, atomMap, successfulAtomMapList, successfulMatchQ, Atol, Rtol)
			else:#on other iterations
				if (a < homMap2TargetLabel):#use canonical ordering
					i=(a,homMap2TargetLabel)
				else:
					i=(homMap2TargetLabel,a)
				if(i in homMap2):
					if(distanceMatchQ(mapping[1],homMap2[i], Atol=Atol, Rtol=Rtol)): #for each case where this is true, check that all other het and hom mappings involving already identified atoms also satisfy the constriant, removing them in the process
						successfulMatchQ = homBlock(i, hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, hetMapType1, homMapType1, hetMapType2, homMapType2, atomMap, successfulAtomMapList, successfulMatchQ, Atol, Rtol)
		#if(not successfulMatchQ):
		#	print atomMap
		#restore the initial popped out mapping
		homMap1[mapping[0]]=mapping[1]
		return successfulMatchQ, successfulAtomMapList
	else: #when both lists are empty, a successful mapping has been found and we return true
		successfulAtomMapList.append(atomMap.copy())
		return True, successfulAtomMapList

def hetBlock(i, hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, hetMapType1, homMapType1, hetMapType2, homMapType2, atomMap, successfulAtomMapList, successfulMatchQ, Atol, Rtol):
	successfulMatchLocal = successfulMatchQ
	#update the dictionaries
	#hetMap1 has already been popped; homMaps don't need to be popped yet
	hetMap2copy = (i, hetMap2.pop(i)) #pop out the assigned mapping, making a copy in the process
	#determine updates to atomMap; in particular: map atoms in molecule 1 to atoms in molecule 2; both the below if statements will hold on the first iteration (when atomMapC is empty), and exactly one should be called on subsequent iterations
	newMaps = {}
	if(mappingLabels[0] not in atomMap):
		newMaps[mappingLabels[0]]=i[0]
	if(mappingLabels[1] not in atomMap):
		newMaps[mappingLabels[1]]=i[1]
	(mappedDistanceMatch,addHetMap1,addHomMap1,addHetMap2,addHomMap2)= smartMappedDistanceMatchQ(hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, atomMap, newMaps, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
	if(mappedDistanceMatch):#if so,  call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
		if(checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomType1, atomType2, atomMap, successfulAtomMapList, Atol=Atol, Rtol=Rtol)[0]):#if they return true, set successfulMatch to true
			successfulMatchLocal = True
		for j in newMaps:
			del atomMap[j]
	#return dictionaries to their previous state
	hetMap2[hetMap2copy[0]]=hetMap2copy[1]
	for j in addHetMap1:
		hetMap1[j]=addHetMap1[j]
	for j in addHomMap1:
		homMap1[j]=addHomMap1[j]
	for j in addHetMap2:
		hetMap2[j]=addHetMap2[j]
	for j in addHomMap2:
		homMap2[j]=addHomMap2[j]
	return successfulMatchLocal

def homBlock(i, hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, hetMapType1, homMapType1, hetMapType2, homMapType2, atomMap, successfulAtomMapList, successfulMatchQ, Atol, Rtol):
	successfulMatchLocal = successfulMatchQ
	#update the dictionaries
	#homMap1 has already been popped; hetMaps are empty
	homMap2copy = (i, homMap2.pop(i)) #pop out the assigned mapping
	newMaps={}
	#look at atomMap to see if any mappings already exist (at most one mapping should exist); no mappings will be found on the first iteration leading to the need for considering two possible mappings
	if((mappingLabels[0] not in atomMap) and (mappingLabels[1] not in atomMap)):
		#assign the two possible mappings
		newMaps[mappingLabels[0]]=i[0]#mapping 1
		newMaps[mappingLabels[1]]=i[1]
		(mappedDistanceMatch1,addHetMap1,addHomMap1,addHetMap2,addHomMap2) = smartMappedDistanceMatchQ(hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, atomMap, newMaps, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
		if(mappedDistanceMatch1):#call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
			if(checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomType1, atomType2, atomMap, successfulAtomMapList, Atol=Atol, Rtol=Rtol)[0]):#if they return true, set successfulMatch to true
				successfulMatchLocal = True
			for j in newMaps:
				del atomMap[j]
		#return dictionaries to their previous state
		homMap2[homMap2copy[0]]=homMap2copy[1]
		for j in addHetMap1:
			hetMap1[j]=addHetMap1[j]
		for j in addHomMap1:
			homMap1[j]=addHomMap1[j]
		for j in addHetMap2:
			hetMap2[j]=addHetMap2[j]
		for j in addHomMap2:
			homMap2[j]=addHomMap2[j]
		newMaps={}#mapping 2
		newMaps[mappingLabels[0]]=i[1]
		newMaps[mappingLabels[1]]=i[0]
		(mappedDistanceMatch2,addHetMap1,addHomMap1,addHetMap2,addHomMap2) = smartMappedDistanceMatchQ(hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, atomMap, newMaps, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
		if(mappedDistanceMatch2):#call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
			if(checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomType1, atomType2, atomMap, successfulAtomMapList, Atol=Atol, Rtol=Rtol)[0]):#if they return true, set successfulMatch to true
				successfulMatchLocal = True
			for j in newMaps:
				del atomMap[j]
	elif(mappingLabels[0] not in atomMap):#label 0 is the new one; label 1 is already mapped
		if(atomMap[mappingLabels[1]]==i[1]):#we need to figure out which atom in the tuple is already mapped
			newMaps[mappingLabels[0]]=i[0]
		elif(atomMap[mappingLabels[1]]==i[0]):
			newMaps[mappingLabels[0]]=i[1]
		else:
			print "Algorithm error: unable to correctly match labels"
		(mappedDistanceMatch,addHetMap1,addHomMap1,addHetMap2,addHomMap2) = smartMappedDistanceMatchQ(hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, atomMap, newMaps, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
		if(mappedDistanceMatch):#if so,  call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
			if(checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomType1, atomType2, atomMap, successfulAtomMapList, Atol=Atol, Rtol=Rtol)[0]):#if they return true, set successfulMatch to true
				successfulMatchLocal = True
			for j in newMaps:
				del atomMap[j]
	elif(mappingLabels[1] not in atomMap): #label 1 is the new one; label 0 is already mapped
		if(atomMap[mappingLabels[0]]==i[1]):#we need to figure out which atom in the tuple is already mapped
			newMaps[mappingLabels[1]]=i[0]
		elif(atomMap[mappingLabels[0]]==i[0]):
			newMaps[mappingLabels[1]]=i[1]
		else:
			print "Algorithm error: unable to correctly match labels"
		(mappedDistanceMatch,addHetMap1,addHomMap1,addHetMap2,addHomMap2) = smartMappedDistanceMatchQ(hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels, atomMap, newMaps, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
		if(mappedDistanceMatch):#if so,  call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
			if(checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomType1, atomType2, atomMap, successfulAtomMapList, Atol=Atol, Rtol=Rtol)[0]):#if they return true, set successfulMatch to true
				successfulMatchLocal = True
			for j in newMaps:
				del atomMap[j]
	else:#this should not happen
		print "Algorithm error: two atoms already mapped when at least one should not be mapped"
	#return dictionaries to their previous state
	homMap2[homMap2copy[0]]=homMap2copy[1]
	for j in addHetMap1:
		hetMap1[j]=addHetMap1[j]
	for j in addHomMap1:
		homMap1[j]=addHomMap1[j]
	for j in addHetMap2:
		hetMap2[j]=addHetMap2[j]
	for j in addHomMap2:
		homMap2[j]=addHomMap2[j]
	return successfulMatchLocal

def distanceMatchQ(val1, val2, Atol=-1, Rtol=-1):
	"""Checks whether two values are within acceptable deviation

	Assumes exactly one of Rtol and Atol is positive; the one that is positive will be used; Rtol is based on the minimum of the values; in this way, the result should not depend on the order of val1 and val2
	"""
	if(Atol>0):#use Atol
		if(abs(val2-val1)< Atol):
		    return True
		else:
		    return False
	else:#Rtol > 0, use Rtol
		if(abs(val2-val1)/min(val1,val2)< Rtol):
		    return True
		else:
		    return False

def mappedDistanceMatchQ(hetMap1, homMap1, hetMap2, homMap2, atomMap, Atol=-1, Rtol=-1):
	"""Given atomMap, checks whether distances in hetMap2 and homMap2 are consistent (within specified tol) with hetMap1; mappings are removed from the dictionaries as they are found to match

	atomMap is not modified; hetMapType vectors are not used
	"""
	if(len(atomMap) < 3): #if there are only two mappings (as there will be on the first iteration), there cannot be anything else to check (the two mappings will already have been checked)
		return True
	else:
		hetMap1_icopy = hetMap1.copy()#make a copy of hetMap1 for the purposes of iteration (according to Python docs, iterating while adding or removing entries may cause RuntimeError)
		for i in hetMap1_icopy:#first go through hetMap1
			if((i[0] in atomMap) and (i[1] in atomMap)):
				targetLabels = (atomMap[i[0]], atomMap[i[1]])
				if(not distanceMatchQ(hetMap1[i], hetMap2[targetLabels], Atol=Atol, Rtol=Rtol)):
				    return False
				else:
				    del hetMap1[i]
				    del hetMap2[targetLabels]

		homMap1_icopy = homMap1.copy()#make a copy of homMap1 for the purposes of iteration (according to Python docs, iterating while adding or removing entries may cause RuntimeError)
		for i in homMap1_icopy:#next go through the homogeneous pairings; note that this, along with above line is basically a copy of the above so we should eventually make a function for it, and call it twice with two different arguments
			if((i[0] in atomMap) and (i[1] in atomMap)):
				if (atomMap[i[0]] < atomMap[i[1]]):#for homogeneous case, these must be sorted to the "canonical" order (1st is lowest label number, 2nd is highest label number)
					targetLabels = (atomMap[i[0]], atomMap[i[1]])
				else:
					targetLabels = (atomMap[i[1]], atomMap[i[0]])
				if(not distanceMatchQ(homMap1[i], homMap2[targetLabels], Atol=Atol, Rtol=Rtol)):
				    return False
				else:
				    del homMap1[i]
				    del homMap2[targetLabels]

		return True #return true if the threshhold has not been exceeded

def smartMappedDistanceMatchQ(hetMap1, homMap1, hetMap2, homMap2, atomType1, atomType2, mappingLabels1, atomMap, newMaps, Atol=-1, Rtol=-1):
	"""Given atomMap, checks whether distances in hetMap2 and homMap2 are consistent (within specified tol) with hetMap1; mappings are removed from the dictionaries as they are found to match

	atomMap IS modified if the function returns True...newMaps is added to atomMap; unlike the regular mappedDistanceMatchQ function, this uses information about the new mappings in an attempt to be faster (we don't iterate over all the elements of the maps and we don't have to search through long atom maps
	newMaps is expected to have length 2 at most (this will occur on the first iteration) and length 1 otherwise; if length 2, the mapping has already been checked earlier within the calling checkDistance function
	mappingLabels1 includes the recently added mapping; we must avoid checking the distance in this case as these cases have already been popped out of the maps for 1 and 2
	I don't think atomType2 is needed, but it is included for symmetry
	"""
	#initialize dictionaries used to store elements that have been deleted that must be later returned; the function will ultimately return these dictionaries
	addHetMap1 = {}
	addHomMap1 = {}
	addHetMap2 = {}
	addHomMap2 = {}
	if(len(atomMap) < 1): #if there are no old mappings (as there will be on the first iteration), there cannot be anything else to check (the two mappings will already have been checked)
		#copy newMaps elements into atomMaps
		for i in newMaps:
			atomMap[i] = newMaps[i]
		return True, addHetMap1, addHomMap1, addHetMap2, addHomMap2
	else:
		for i in newMaps:#iterate over new atom distances, corresponding to a combination of a new mapped atom with an old mapped atom
			i2 = newMaps[i]
			iT1 = atomType1[i-1]
			for j in atomMap:
				j2 = atomMap[j]
				jT1 = atomType1[j-1]
				if(iT1==jT1):#homogeneous case (the atom types are the same)
					if(i<j): #canonical origin ordering is (i,j)
						originLabels = (i, j)
					else:#canonical origin ordering (j,i)
						originLabels = (j, i)
					if(i2<j2):#this if-else block determines canonical target ordering; canonical target ordering = (i2, j2)
						targetLabels = (i2, j2)
					else: #canonical target ordering = (j2, i2)
						targetLabels = (j2, i2)
					if(originLabels==mappingLabels1):#every time this function is called there should be one case that has just been checked in the calling check distance function; we do not want to try to check this case as the elements have already been removed from the distance maps
						pass
					elif(not distanceMatchQ(homMap1[originLabels], homMap2[targetLabels], Atol=Atol, Rtol=Rtol)):
						return False, addHetMap1, addHomMap1, addHetMap2, addHomMap2
					else:
						addHomMap1[originLabels]=homMap1.pop(originLabels)
						addHomMap2[targetLabels]=homMap2.pop(targetLabels)
				else:#heterogeneous case
					if(iT1<jT1): #canonical ordering is (i,j)
						originLabels = (i, j)
						targetLabels = (i2, j2)
					else: #(atomType1[i-1]>atomType1[j-1]): canonical ordering is (j,i)
						originLabels = (j, i)
						targetLabels = (j2, i2)
					if(originLabels==mappingLabels1):#every time this function is called there should be one case that has just been checked in the calling check distance function; we do not want to try to check this case as the elements have already been removed from the distance maps
						pass
					elif(not distanceMatchQ(hetMap1[originLabels], hetMap2[targetLabels], Atol=Atol, Rtol=Rtol)):
						return False, addHetMap1, addHomMap1, addHetMap2, addHomMap2
					else:
						addHetMap1[originLabels]=hetMap1.pop(originLabels)
						addHetMap2[targetLabels]=hetMap2.pop(targetLabels)

		#copy newMaps elements into atomMaps
		for i in newMaps:
			atomMap[i] = newMaps[i]

		return True, addHetMap1, addHomMap1, addHetMap2, addHomMap2 #return true if the threshhold has not been exceeded

def calcDistanceDeviationsGivenMapping(mg1, mg2, atomMap):
	"""Calculate absolute and relative atom-to-atom distance deviations between two molecules, given an atom mapping

	mg1 and mg2 are the two conformers to be compared
	atomMap is a dictionary mapping mg1 atom labels to mg2 atom labels
	output: 1: a dictionary indicating the absolute distance deviations between mg1 and mg2 (more specifically d1-d2, where d1 and d2 are the atom-to-atom distances in mg1 and mg2, respectively) for each atom pair (keys correspond to mg1 labels)
		2: a dictionary indicating the relative distance deviations between mg1 and mg2 (more specifically (d1-d2)/min(d1,d2), where d1 and d2 are the atom-to-atom distances in mg1 and mg2, respectively) for each atom pair (keys correspond to mg1 labels)
	"""
	#assert mg1 and mg2 contain the same number of each type of atom (and obviously the same total number of atoms)

	#generate distance mappings
	(hetMap1, homMap1, hetMapType1, homMapType1)=mg1.getDistanceMappings()
	(hetMap2, homMap2, hetMapType2, homMapType2)=mg2.getDistanceMappings()

	distDevAbs = {}
	distDevRel = {}
	#first, go through the hetMaps
	for i in hetMap1:
		d1 = hetMap1[i]
		d2 = hetMap2[(atomMap[i[0]],atomMap[i[1]])]
		distDevAbs[i]=d1-d2
		distDevRel[i]=(d1-d2)/min(d1,d2)
	#second, go through the homMaps
	for i in homMap1:
		d1 = homMap1[i]
		j0 = atomMap[i[0]]
		j1 = atomMap[i[1]]
		if j0 < j1:
			d2 = homMap2[(j0,j1)]
		else:
			d2 = homMap2[(j1,j0)]
		distDevAbs[i]=d1-d2
		distDevRel[i]=(d1-d2)/min(d1,d2)

	return distDevAbs, distDevRel

def dictionaryMaxAbs(dict):
	"""Calculates the maximum absolute value among values in a dictionary"""
	#would it be faster to iterate through values, taking absolute value for each one and comparing to maximum found so far?
	v = dict.values()
	return max(abs(min(v)),abs(max(v)))


def readMOLFileWithConnectivity(filename):
	"""Given a MOL file, constructs a MolecularGeometry object with connectivity information

	Currently only supports C, H, and O atoms
	"""
	f = open(filename, 'r')
	#first three lines are irrelevant
	f.readline()
	f.readline()
	f.readline()

	line=f.readline()
	n = int(line[0:3])#read the number of atoms
	b = int(line[3:6])#read the number of bonds
	#initialize the atomTypes, connectivity, and atomCoords arrays
	#atomCoords = [[0.0 for j in range(3)] for i in range(n)]
	atomCoords = [[] for i in range(n)]
	connectivity = [[0,0] for i in range(b)]
	atomTypes = [0 for i in range(n)]

	#read info from the mole file into the arrays
	for i in range(n):
		splitLine = f.readline().split()
		atomCoords[i] = [float(splitLine[0]), float(splitLine[1]), float(splitLine[2])]
		atomTypes[i] = atomicSymbolToNumber(splitLine[3])

	#read connectivity info from mole file
	for i in range(b):
		line = f.readline()
		connectivity[i][0] = int(line[0:3])
		connectivity[i][1] = int(line[3:6])

	f.close() #close the file

	return MolecularGeometry(atomTypes,atomCoords,connectivity) #return the MolecularGeometry object

def readMOLFile(filename):
	"""Given a MOL file, constructs a MolecularGeometry object

	Currently only supports C, H, and O atoms
	"""
	f = open(filename, 'r')
	#first three lines are irrelevant
	f.readline()
	f.readline()
	f.readline()

	n = int(f.readline()[0:3])#read the number of atoms
	#initialize the atomTypes and atomCoords arrays
	#atomCoords = [[0.0 for j in range(3)] for i in range(n)]
	atomCoords = [[] for i in range(n)]
	atomTypes = [0 for i in range(n)]

	#read info from the mole file into the arrays
	for i in range(n):
		splitLine = f.readline().split()
		atomCoords[i] = [float(splitLine[0]), float(splitLine[1]), float(splitLine[2])]
		atomTypes[i] = atomicSymbolToNumber(splitLine[3])
	
	f.close() #close the file

	return MolecularGeometry(atomTypes,atomCoords) #return the MolecularGeometry object

def readMM4File(filename,startline):
	"""Given a MM4 file, constructs a MolecularGeometry object

	Currently only supports C, H, and O atoms; startline is an index (starting with 1) indicating where the block of interest starts; in most cases, except conformer output sets, this will be 1
	"""
	f = open(filename, 'r')
	#get to the block of interest
	for i in range (startline-1):
		f.readline()

	n = int(f.readline()[61:65])#read the number of atoms
	#initialize the atomTypes and atomCoords arrays
	atomCoords = [[] for i in range(n)]
	atomTypes = [0 for i in range(n)]

	#look for the beginning of the coordinates by checking for C H or O
	line = f.readline()
	while ( line[30:32].strip()!='C' and line[30:32].strip()!='H' and line[30:32].strip()!='O'):
		line = f.readline()
	#read coordinate info from the file into the arrays
	for i in range(n):
		#splitLine = line.split()
		atomCoords[i] = [float(line[0:10]), float(line[10:20]), float(line[20:30])]
		atomTypes[i] = atomicSymbolToNumber(line[30:32].strip())
		line = f.readline()

	f.close() #close the file

	return MolecularGeometry(atomTypes,atomCoords) #return the MolecularGeometry object


def atomicSymbolToNumber(symbol):
	"""Converts atomic symbol string to atomic number integer (e.g. 'C'->6)

	Currently only supports C, H, and O
	"""
	if(symbol=='C'):
		return 6
	elif(symbol=='O'):
		return 8
	else:#hydrogen
		return 1

def atomicNumberToSymbol(number):
	"""Converts atomic number integer to atomic symbol string  (e.g. 6->'C')

	Currently only supports C, H, and O
	"""
	if(number==6):
		return 'C'
	elif(number==8):
		return 'O'
	else:#hydrogen
		return 'H'

def atomicNumberToMM4Type(number):
	"""Converts atomic number integer to MM4 atom type integer  (e.g. 6 (carbon)->1)

	Currently only supports C, H, and O; relies on MM4 option KFIXTYP to refine atom types to be more specific/accurate
	"""
	if(number==6):
		return 1
	elif(number==8):
		return 6
	else:#hydrogen
		return 5

#if __name__ == '__main__':
#	a = readMOLFile('JP10A.mol')
#	b = readMOLFile('JP10B.mol')
#	(q, n) = checkConformationalEquivalence(b, a, Atol=0.10)
#	a = MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
#	b = MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
#	(q, n) = checkConformationalEquivalence(b, a, Atol=0.01)