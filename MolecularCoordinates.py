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
	
	def __init__(self, atomVec, xyzCoor, atomLabels=None):
		self.xyzCoor=xyzCoor
		self.atomVec=atomVec
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

	atomMap = {}
	successfulAtomMapList = []

	(matchQ, matches) = checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomMap, successfulAtomMapList, Atol=Atol, Rtol=Rtol)
	
	return matchQ, matches

def checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomMap, successfulAtomMapList, Atol=-1, Rtol=-1):
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
			for i in hetMap1:
				if(i[0] in atomMap):
					mapping =  (i, hetMap1.pop(i))
					hetMap2TargetLabel = atomMap[i[0]]
					break
				elif(i[1] in atomMap):
					mapping = (i, hetMap1.pop(i))
					hetMap2TargetLabel = atomMap[i[1]]
					break
		mappingLabels = mapping[0]
		mappingType = hetMapType1.pop(mappingLabels)
		hetMapType2_icopy = hetMapType2.copy()#make a copy of hetMapType2 for the purposes of iteration (according to Python docs, iterating while adding or removing entries may cause RuntimeError)
		for i in hetMapType2_icopy:#search in hetMapType2 for cases with the same mapping type and a pre-established atom correspondence
			if(hetMap2TargetLabel==i[0] or hetMap2TargetLabel==i[1] or hetMap2TargetLabel==-1): # this line checks whether the target atom label is consistent with previously established mappings;
				if(hetMapType2[i]==mappingType):#this line does the mapping type check
					if(distanceMatchQ(mapping[1],hetMap2[i], Atol=Atol, Rtol=Rtol)): #for each case where this is true, check that all other het and hom mappings involving already identified atoms also satisfy the constriant, removing them in the process
						#copy all the dictionaries
						hetMap1C = hetMap1.copy()
						homMap1C = homMap1.copy()
						hetMapType1C = hetMapType1.copy()
						homMapType1C = homMapType1.copy()
						hetMap2C = hetMap2.copy()
						homMap2C = homMap2.copy()
						hetMapType2C = hetMapType2.copy()
						homMapType2C = homMapType2.copy()
						atomMapC = atomMap.copy()
						#update the dictionaries
						#hetMap1 and hetMapType1 have already been popped; homMaps don't need to be popped yet
						del hetMap2C[i] #pop out the assigned mapping
						del hetMapType2C[i]
						#update atomMapC; in particular: map atoms in molecule 1 to atoms in molecule 2; both the below if statements will hold on the first iteration (when atomMapC is empty), and exactly one should be called on subsequent iterations
						if(mappingLabels[0] not in atomMapC):
							atomMapC[mappingLabels[0]]=i[0]
						if(mappingLabels[1] not in atomMapC):
							atomMapC[mappingLabels[1]]=i[1]
						mappedDistanceMatch = mappedDistanceMatchQ(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
						if(mappedDistanceMatch):#if so,  call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
							if(checkDistance(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, successfulAtomMapList, Atol=Atol, Rtol=Rtol)):#if they return true, set successfulMatch to true
								successfulMatchQ = True
		#if(not successfulMatchQ):
		#	print atomMap
		return successfulMatchQ, successfulAtomMapList
	elif(len(homMap1)>0):
		#this should only be encountered for molecules like fullerene, graphene, hydrogen, etc. with only one type of atom
		#similar to hetMap case above(differences indicated by +++), except there are two possible mappings on first assignment
		#for the first iteration, when atomMap is empty, pull an arbitrary item from hetMap1; on subsequent iterations, make sure one (it should be exactly one by the design of the algorithm) of the elements is in atom map
		if(len(atomMap)==0):#if this is the first iteration:
			mapping = homMap1.popitem()
			homMap2TargetLabel = -1 #-1 corresponds to first iteration, where we don't have a particular target atom
		else:
			for i in homMap1:
				if(i[0] in atomMap):
					mapping =  (i, homMap1.pop(i))
					homMap2TargetLabel = atomMap[i[0]]
					break
				elif(i[1] in atomMap):
					mapping = (i, homMap1.pop(i))
					homMap2TargetLabel = atomMap[i[1]]
					break
		mappingLabels = mapping[0]
		mappingType = homMapType1.pop(mappingLabels)
		homMapType2_icopy = homMapType2.copy()#make a copy of homMapType2 for the purposes of iteration (according to Python docs, iterating while adding or removing entries may cause RuntimeError)
		for i in homMapType2_icopy:#search in homMapType2 for cases with the same mapping type;
			if(homMap2TargetLabel==i[0] or homMap2TargetLabel==i[1] or homMap2TargetLabel==-1): # this line checks whether the target atom label is consistent with previously established mappings
				if(homMapType2[i]==mappingType):
					if(distanceMatchQ(mapping[1],homMap2[i], Atol=Atol, Rtol=Rtol)): #for each case where this returns true, check that all other het and hom mappings involving already identified atoms also satisfy the constriant, removing them in the process
						#copy all the dictionaries
						#+++ here, hetMaps are empty and don't need to be copied, but it is easier for consistent code to do so
						hetMap1C = hetMap1.copy()
						homMap1C = homMap1.copy()
						hetMapType1C = hetMapType1.copy()
						homMapType1C = homMapType1.copy()
						hetMap2C = hetMap2.copy()
						homMap2C = homMap2.copy()
						hetMapType2C = hetMapType2.copy()
						homMapType2C = homMapType2.copy()
						atomMapC = atomMap.copy()
						#update the dictionaries
						#homMap1 and homMapType1 have already been popped; hetMaps are empty
						del homMap2C[i] #pop out the assigned mapping
						del homMapType2C[i]
						#+++look at atomMapC to see if any mappings already exist (at most one mapping should exist); no mappings will be found on the first iteration leading to the need for considering two possible mappings
						if((mappingLabels[0] not in atomMapC) and (mappingLabels[1] not in atomMapC)):
							#in this case, we need to do two recursive calls, so we must make extra copies
							hetMap1D = hetMap1.copy()
							homMap1D = homMap1.copy()
							hetMapType1D = hetMapType1.copy()
							homMapType1D = homMapType1.copy()
							hetMap2D = hetMap2.copy()
							homMap2D = homMap2.copy()
							hetMapType2D = hetMapType2.copy()
							homMapType2D = homMapType2.copy()
							atomMapD = atomMap.copy()
							del homMap2D[i] #pop out the assigned mapping for the second copy
							del homMapType2D[i]
							#assign the two possible mappings
							atomMapC[mappingLabels[0]]=i[0]#mapping 1
							atomMapC[mappingLabels[1]]=i[1]
							atomMapD[mappingLabels[0]]=i[1]#mapping 2
							atomMapD[mappingLabels[1]]=i[0]
							mappedDistanceMatch1 = mappedDistanceMatchQ(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
							mappedDistanceMatch2 = mappedDistanceMatchQ(hetMap1D, homMap1D, hetMapType1D, homMapType1D, hetMap2D, homMap2D, hetMapType2D, homMapType2D, atomMapD, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
							if(mappedDistanceMatch1):#call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
								if(checkDistance(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, successfulAtomMapList, Atol=Atol, Rtol=Rtol)):#if they return true, set successfulMatch to true
									successfulMatchQ = True
							if(mappedDistanceMatch2):#call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
								if(checkDistance(hetMap1D, homMap1D, hetMapType1D, homMapType1D, hetMap2D, homMap2D, hetMapType2D, homMapType2D, atomMapD, successfulAtomMapList, Atol=Atol, Rtol=Rtol)):#if they return true, set successfulMatch to true
									successfulMatchQ = True
						elif(mappingLabels[0] not in atomMapC):#label 0 is the new one; label 1 is already mapped
							if(atomMapC[mappingLabels[1]]==i[1]):#we need to figure out which atom in the tuple is already mapped
								atomMapC[mappingLabels[0]]=i[0]
							elif(atomMapC[mappingLabels[1]]==i[0]):
								atomMapC[mappingLabels[0]]=i[1]
							else:
								print "Algorithm error: unable to correctly match labels"
							mappedDistanceMatch = mappedDistanceMatchQ(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
							if(mappedDistanceMatch):#if so,  call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
								if(checkDistance(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, successfulAtomMapList, Atol=Atol, Rtol=Rtol)):#if they return true, set successfulMatch to true
									successfulMatchQ = True
						elif(mappingLabels[1] not in atomMapC): #label 1 is the new one; label 0 is already mapped
							if(atomMapC[mappingLabels[0]]==i[1]):#we need to figure out which atom in the tuple is already mapped
								atomMapC[mappingLabels[1]]=i[0]
							elif(atomMapC[mappingLabels[0]]==i[0]):
								atomMapC[mappingLabels[1]]=i[1]
							else:
								print "Algorithm error: unable to correctly match labels"
							mappedDistanceMatch = mappedDistanceMatchQ(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
							if(mappedDistanceMatch):#if so,  call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
								if(checkDistance(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, successfulAtomMapList, Atol=Atol, Rtol=Rtol)):#if they return true, set successfulMatch to true
									successfulMatchQ = True

						else:#this should not happen
							print "Algorithm error: two atoms already mapped when at least one should not be mapped"
		#if(not successfulMatchQ):
		#	print atomMap
		return successfulMatchQ, successfulAtomMapList
	else: #when both lists are empty, a successful mapping has been found and we return true
		successfulAtomMapList.append(atomMap)
		return True, successfulAtomMapList

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

def mappedDistanceMatchQ(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomMap, Atol=-1, Rtol=-1):
	"""Given atomMap, checks whether distances in hetMap2 and homMap2 are consistent (within specified tol) with hetMap1; mappings are removed from the dictionaries as they are found to match

	atomMap is not modified
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
				    del hetMapType1[i]
				    del hetMap2[targetLabels]
				    del hetMapType2[targetLabels]

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
				    del homMapType1[i]
				    del homMap2[targetLabels]
				    del homMapType2[targetLabels]

		return True #return true if the threshhold has not been exceeded

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


def readMOLFile(filename):
	"""Given a MOL file, constructs a MolecularGeometry object

	Currently only supports C, H, and O atoms
	"""
	f = open(filename, 'r')
	#first three lines are irrelevant
	f.readline()
	f.readline()
	f.readline()

	n = int(f.readline().split()[0])#read the number of atoms
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

def atomicSymbolToNumber(symbol):
	"""Converts atomic symbol string to atomic number integer(e.g. 'C'->1)

	Currently only supports C, H, and O
	"""
	if(symbol=='C'):
		return 6
	elif(symbol=='O'):
		return 8
	else:#hydrogen
		return 1
