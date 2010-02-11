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
			dist[j][i]=dist[i][j]

		return dist

	def getDistanceMappings(self):
		"""
		Generates dictionaries mapping atom label pairs to associated distances

		Outputs are heterogeneous mappings (different atom types in tuple) and homogeneous mappings (same atom types) to distances and mappings to atomTypeTuples
		For heterogeneous mappings, the atom type tuple will be ordered with lowest # atom type first

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

	#generate distance mappings
	(hetMap, homMap, hetMapType, homMapType)=mg1.getDistanceMappings()
	(hetMap, homMap, hetMapType, homMapType)=mg2.getDistanceMappings()

	return matchQ, nmatches
