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
		2) an integer describing the number of distinct atom mappings giving equivalent conformations within the specified tolerance 
	"""
	#assert mg1 and mg2 contain the same number of each type of atom (and obviously the same total number of atoms)
	#assert Atol > 0 or Rtol > 0 (and maybe also that Atol < 0 or Rtol < 0; i.e. only one of Atol and Rtol should be specified)

	natoms = mg1.atoms
	#generate distance mappings
	(hetMap1, homMap1, hetMapType1, homMapType1)=mg1.getDistanceMappings()
	(hetMap2, homMap2, hetMapType2, homMapType2)=mg2.getDistanceMappings()

	atomMap = []
	matchQ = checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomMap, Atol=-1, Rtol=-1)
	#nmatches = size of number of mappings

	return matchQ, nmatches

def checkDistance(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomMap, Atol=-1, Rtol=-1):
	"""Recursive function to assign mappings between molecule 1 and molecule 2 based on distances


	"""
	#initialize boolean variable to false
	successfulMatchQ = false

	#use the hetMap, if possible; if not, use homMap
	if(len(hetMap1)>0):
		mapping = hetMap1.popitem()
		mappingLabels = mapping[0]
		mappingType = hetMapType1.pop(mappingLabels)
		hetMapType2_icopy = hetMapType2.copy()#make a copy of hetMapType2 for the purposes of iteration (according to Python docs, iterating while adding or removing entries may cause RuntimeError)
		for i in hetMapType2_icopy.iterkeys():#search in hetMapType2 for cases with the same mapping type;
			if(hetMapType2[i]==mappingType): #when they are encountered, perform a distanceMatchQ check
				if(distanceMatchQ(mapping[1],hetMap2[i], Atol=Atol, Rtol=Rtol)): #for each case where this returns true, check that all other het and hom mappings involving already identified atoms also satisfy the constriant, removing them in the process
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
					del hetMap2TypeC[i]
					#update atomMapC; in particular: map atoms in molecule 1 to atoms in molecule 2; both the below if statements will hold on the first iteration (when atomMapC is empty), and exactly one or two should be called on subsequent iterations (***we could modify so that exactly one atom is new on subsequent iterations...would this be better/faster?)
					if(!(mappingLabels[0] in atomMapC)):
						atomMapC[mappingLabels[0]]=i[0]
					if(!(mappingLabels[1] in atomMapC)):
						atomMapC[mappingLabels[1]]=i[1]
					mappedDistanceMatch = mappedDistanceMatchQ(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
					if(mappedDistanceMatch):#if so,  call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
						if(checkDistance(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)):#if they return true, set successfulMatch to true
							successfulMatchQ = true
		return successfulMatchQ
	elif(len(homMap1)>0):
		#this should only be encountered for molecules like fullerene, graphene, hydrogen, etc. with only one type of atom
		#similar to hetMap case above(differences indicated by +++), except there are two possible mappings on first assignment
		mapping = homMap1.popitem()
		mappingLabels = mapping[0]
		mappingType = homMapType1.pop(mappingLabels)
		homMapType2_icopy = homMapType2.copy()#make a copy of homMapType2 for the purposes of iteration (according to Python docs, iterating while adding or removing entries may cause RuntimeError)
		for i in homMapType2_icopy.iterkeys():#search in homMapType2 for cases with the same mapping type;
			if(homMapType2[i]==mappingType): #when they are encountered, perform a distanceMatchQ check
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
					del homMap2TypeC[i]
					#+++look at atomMapC to see if any mappings already exist (at most one mapping should exist); no mappings will be found on the first iteration, and possibly on subsequent iterations, leading to the need for considering two possible mappings; (***we could modify so that exactly one atom is new on subsequent iterations...this would PROBABLY be better/faster
					if(!(mappingLabels[0] in atomMapC)&&!(mappingLabels[1] in atomMapC)):
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
						del homMap2TypeD[i]
						#assign the two possible mappings
						atomMapC[mappingLabels[0]]=i[0]#mapping 1
						atomMapC[mappingLabels[1]]=i[1]
						atomMapD[mappingLabels[0]]=i[1]#mapping 2
						atomMapD[mappingLabels[1]]=i[0]
						mappedDistanceMatch1 = mappedDistanceMatchQ(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
						mappedDistanceMatch2 = mappedDistanceMatchQ(hetMap1D, homMap1D, hetMapType1D, homMapType1D, hetMap2D, homMap2D, hetMapType2D, homMapType2D, atomMapD, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
						if(mappedDistanceMatch1):#call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
							if(checkDistance(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)):#if they return true, set successfulMatch to true
								successfulMatchQ = true
						if(mappedDistanceMatch2):#call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
							if(checkDistance(hetMap1D, homMap1D, hetMapType1D, homMapType1D, hetMap2D, homMap2D, hetMapType2D, homMapType2D, atomMapD, Atol=Atol, Rtol=Rtol)):#if they return true, set successfulMatch to true
								successfulMatchQ = true
					elif(!(mappingLabels[0] in atomMapC)):#label 0 is the new one; label 1 is already mapped
						if(atomMapC[mappingLabels[1]]==i[1]):#we need to figure out which atom in the tuple is already mapped
							atomMapC[mappingLabels[0]]=i[0]
						elif(atomMapC[mappingLabels[1]]==i[0]):
							atomMapC[mappingLabels[0]]=i[1]
						else:
							print "Algorithm error: unable to correctly match labels"
						mappedDistanceMatch = mappedDistanceMatchQ(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
						if(mappedDistanceMatch):#if so,  call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
							if(checkDistance(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)):#if they return true, set successfulMatch to true
								successfulMatchQ = true
					elif(!(mappingLabels[1] in atomMapC)): #label 1 is the new one; label 0 is already mapped
						if(atomMapC[mappingLabels[0]]==i[1]):#we need to figure out which atom in the tuple is already mapped
							atomMapC[mappingLabels[1]]=i[0]
						elif(atomMapC[mappingLabels[0]]==i[0]):
							atomMapC[mappingLabels[1]]=i[1]
						else:
							print "Algorithm error: unable to correctly match labels"
						mappedDistanceMatch = mappedDistanceMatchQ(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)#note that, by design, this will modify (remove elements from) the dictionaries
						if(mappedDistanceMatch):#if so,  call a new instance of checkDistance with copies (dict.copy()) of variables with appropriately adjusted/popped values
							if(checkDistance(hetMap1C, homMap1C, hetMapType1C, homMapType1C, hetMap2C, homMap2C, hetMapType2C, homMapType2C, atomMapC, Atol=Atol, Rtol=Rtol)):#if they return true, set successfulMatch to true
								successfulMatchQ = true

					else:#this should not happen
						print "Algorithm error: two atoms already mapped when at least one should not be mapped"

		return successfulMatchQ
	else:
		#***add (now complete) mapping to "global variable"
		return true

def distanceMatchQ(val1, val2, Atol=-1, Rtol=-1):
	"""Checks whether two values are within acceptable deviation

	Currently, this is only written for Atol
	"""
	if(abs(val2-val1)< Atol):
	    return true
	else:
	    return false

def mappedDistanceMatchQ(hetMap1, homMap1, hetMapType1, homMapType1, hetMap2, homMap2, hetMapType2, homMapType2, atomMap, Atol=-1, Rtol=-1):
	"""Given atomMap, checks whether distances in hetMap2 and homMap2 are consistent (within specified tol) with hetMap1; mappings are removed from the dictionaries as they are found to match

	atomMap is not modified
	"""
	if(len(atomMap) < 3): #if there are only two mappings (as there will be on the first iteration), there cannot be anything else to check (the two mappings will already have been checked)
		return true
	else:
		hetMap1_icopy = hetMap1.copy()#make a copy of hetMap1 for the purposes of iteration (according to Python docs, iterating while adding or removing entries may cause RuntimeError)
		for i in hetMap1_icopy.iterkeys():#first go through hetMap1
			if((i[0] in atomMap) && (i[1] in atomMap)):
				targetLabels = (atomMap(i[0]), atomMap(i[1]))
				if(!distanceMatchQ(hetMap1[i], hetMap2[targetLabels], Atol=Atol, Rtol=Rtol)):
				    return false
				else:
				    del hetMap1[i]
				    del hetMapType1[i]
				    del hetMap2[targetLabels]
				    del hetMapType2[i]

		homMap1_icopy = homMap1.copy()#make a copy of homMap1 for the purposes of iteration (according to Python docs, iterating while adding or removing entries may cause RuntimeError)
		for i in homMap1_icopy.iterkeys():#next go through the homogeneous pairings; note that this, along with above line is basically a copy of the above so we should eventually make a function for it, and call it twice with two different arguments
			if((i[0] in atomMap) && (i[1] in atomMap)):
				targetLabels = (atomMap(i[0]), atomMap(i[1]))
				if(!distanceMatchQ(homMap1[i], homMap2[targetLabels], Atol=Atol, Rtol=Rtol)):
				    return false
				else:
				    del homMap1[i]
				    del homMapType1[i]
				    del homMap2[targetLabels]
				    del homMapType2[i]

		return true #return true if the threshhold has not been exceeded
	    