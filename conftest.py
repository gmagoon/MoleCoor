import unittest
import MolecularCoordinates
import conftest


class  ConfTestCase(unittest.TestCase):

	def testSimpleHetConf(self):
		""" A simple test (using the same molecule) of the heterogeneous case

		A pyramidal form for CH3
		"""
		(q, n) = conftest.SimpleHetConf()
		self.assertEqual(q, True)
		self.assertEqual(len(n), 6)

	def testSimpleHomConf(self):
		""" A simple test (using the same molecule) of the homogenous case

		A cube of hydrogens
		"""
		(q, n) = conftest.SimpleHomConf()
		self.assertEqual(q, True)
		self.assertEqual(len(n), 48) #I believe the answer is 48, though I am not positive off the top of my head, and it could be higher

	def testMirrorImageConf(self):
		""" A test of mirror images

		Mirror image conformations of gauche n-butane (each independently optimized using PM3 in Gaussian03 using default convergence criteria)
		"""
		(q, n) = conftest.MirrorImageConf()
		self.assertEqual(q, True)
		self.assertEqual(len(n), 2) #I believe the correct answer is 2

	def testBuckminsterfullerene(self):
		""" A test of Buckminsterfullerene with itself

		"""
		(q, n) = conftest.Buckminsterfullerene()
		self.assertEqual(q, True)
		self.assertEqual(len(n), 120) #A result of 120 seems reasonable (symmetry number is 60)

	def testDistinctJP10Conf(self):
		""" A test of distinct JP-10 conformations

		The conformations differ in the puckering of one of the rings; conformations come from CBS-QB3 calculations with opt=verytight and int=ultrafine (JP10_A.log and JP10_B.log in my records)
		"""
		(q, n) = conftest.DistinctJP10Conf()
		self.assertEqual(q, False)

	def testAtomTypeSwapConf(self):
		""" A test of molecules with swapped atom types

		The two "conformations" have atoms in the same position, but the types of atoms in each position are different, and thus should produce a value of False when compared
		"""
		(q, n) = conftest.AtomTypeSwapConf()
		self.assertEqual(q, False)

	def testLongLinearChainConf(self):
		""" A test of a long linear chain

		There should only be one unique atom mapping here by the construction of the distances (each atom-to-atom distance should be unique ; the function LongLinearChainConf is actually more useful for timing of scaling for large numbers of atoms, but it doesn't hurt to also check that the results are what we expect
		"""
		(q, n) = conftest.LongLinearChainConf(30)#use a chain of length 30 for this test
		self.assertEqual(q, True)
		self.assertEqual(len(n), 1)

#	def testOptConf(self):
#		""" A test of equivalent conformations of a large molecule obtained from optimization with different potential energy calculation methods
#
#		"""
#		#a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
#		#b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
#		#alternatively, we could play around with Atol/Rtol to see how tight they need to be
#		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.05)
#		print q
#
#	def testInitialConditionConf(self):
#		""" A test of equivalent conformations of a large molecule optimized with different initial conditions
#
#		The same potential energy calculation method is used in both cases. Default convergence criterion of Gaussian03 are used.
#		"""
#		#a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
#		#b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
#		#it may be appropriate to use Rtol here, particularly for a large molecule
#		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.05)
#		self.assertEqual(q, True)

	def testCalcDistDev(self):
		""" A test of distance deviations, given a mapping

		Uses mirror image conformations of gauche n-butane (each independently optimized using PM3 in Gaussian03 using default convergence criteria)
		"""
		ac =    [[-1.5846,   -0.5216,    0.1338],
			[-1.1512,   -1.5146,   -0.0477 ],
			[-1.7781,   -0.4440,    1.2120 ],
			[-2.5555,   -0.4907,   -0.3770 ],
			[-0.6751,    0.5855,   -0.3509 ],
			[-1.1705,    1.5633,   -0.1896 ],
			[-0.5293,    0.5028,   -1.4463 ],
			[ 0.6751,    0.5854,    0.3509 ],
			[ 1.1705,    1.5633,    0.1895 ],
			[ 0.5294,    0.5028,    1.4463 ],
			[ 1.5846,   -0.5216,   -0.1338 ],
			[ 2.5554,   -0.4908,    0.3772 ],
			[ 1.1510,   -1.5146,    0.0475 ],
			[ 1.7783,   -0.4438,   -1.2120]]
		aa = [6, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1]
		bc =  [[  -1.5845,   -0.5216,   -0.1338],
			   [-1.7784,   -0.4438,   -1.2120 ],
			   [-1.1509,   -1.5146,    0.0474 ],
			   [-2.5553,   -0.4910,    0.3773 ],
			   [-0.6751,    0.5855,    0.3508 ],
			   [-0.5293,    0.5028,    1.4463 ],
			   [-1.1705,    1.5634,    0.1895 ],
			   [ 0.6751,    0.5855,   -0.3508 ],
			   [ 1.1705,    1.5633,   -0.1893 ],
			   [ 0.5294,    0.5030,   -1.4463 ],
			   [ 1.5845,   -0.5216,    0.1338 ],
			   [ 1.1510,   -1.5146,   -0.0478 ],
			   [ 2.5555,   -0.4908,   -0.3770 ],
			   [ 1.7780,   -0.4441,    1.2121 ]]
		ba = [6, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1]
		a = MolecularCoordinates.MolecularGeometry(aa,ac)
		b = MolecularCoordinates.MolecularGeometry(ba,bc)
		#atomMap = {1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:7, 8:8, 9:9, 10:10, 11:11, 12:12, 13:13, 14:14} #a naive mapping
		atomMap = {1: 1, 2: 3, 3: 2, 4: 4, 5: 5, 6: 7, 7: 6, 8: 8, 9: 9, 10: 10, 11: 11, 12: 13, 13: 12, 14: 14} #a "correct" mapping (there is also one other "correct" mapping)
		(distDevAbs,distDevRel) = MolecularCoordinates.calcDistanceDeviationsGivenMapping(a, b, atomMap)
		distDevAbsMax = MolecularCoordinates.dictionaryMaxAbs(distDevAbs)
		distDevRelMax = MolecularCoordinates.dictionaryMaxAbs(distDevRel)
		self.assertAlmostEqual(distDevAbsMax, 0.000528777418128, 5)
		self.assertAlmostEqual(distDevRelMax, 0.000145546853437, 5)

def SimpleHetConf():
	a = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
	b = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
	return MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)


def SimpleHomConf():
	a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
	b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
	return MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)

def MirrorImageConf():
	ac =    [[-1.5846,   -0.5216,    0.1338],
			[-1.1512,   -1.5146,   -0.0477 ],
			[-1.7781,   -0.4440,    1.2120 ],
			[-2.5555,   -0.4907,   -0.3770 ],
			[-0.6751,    0.5855,   -0.3509 ],
			[-1.1705,    1.5633,   -0.1896 ],
			[-0.5293,    0.5028,   -1.4463 ],
			[ 0.6751,    0.5854,    0.3509 ],
			[ 1.1705,    1.5633,    0.1895 ],
			[ 0.5294,    0.5028,    1.4463 ],
			[ 1.5846,   -0.5216,   -0.1338 ],
			[ 2.5554,   -0.4908,    0.3772 ],
			[ 1.1510,   -1.5146,    0.0475 ],
			[ 1.7783,   -0.4438,   -1.2120]]
	aa = [6, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1]
	bc =  [[  -1.5845,   -0.5216,   -0.1338],
		   [-1.7784,   -0.4438,   -1.2120 ],
		   [-1.1509,   -1.5146,    0.0474 ],
		   [-2.5553,   -0.4910,    0.3773 ],
		   [-0.6751,    0.5855,    0.3508 ],
		   [-0.5293,    0.5028,    1.4463 ],
		   [-1.1705,    1.5634,    0.1895 ],
		   [ 0.6751,    0.5855,   -0.3508 ],
		   [ 1.1705,    1.5633,   -0.1893 ],
		   [ 0.5294,    0.5030,   -1.4463 ],
		   [ 1.5845,   -0.5216,    0.1338 ],
		   [ 1.1510,   -1.5146,   -0.0478 ],
		   [ 2.5555,   -0.4908,   -0.3770 ],
		   [ 1.7780,   -0.4441,    1.2121 ]]
	ba = [6, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1]
	a = MolecularCoordinates.MolecularGeometry(aa,ac)
	b = MolecularCoordinates.MolecularGeometry(ba,bc)
	return MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.001)

def Buckminsterfullerene():
	a = MolecularCoordinates.readMOLFile('c60_c.mol')
	b = MolecularCoordinates.readMOLFile('c60_c.mol')
	return MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)

def DistinctJP10Conf():
	a = MolecularCoordinates.readMOLFile('JP10A.mol')
	b = MolecularCoordinates.readMOLFile('JP10B.mol')
	return MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.10)

def AtomTypeSwapConf():
	a = MolecularCoordinates.MolecularGeometry([8,1,1],[[0,0,0],[1,0,0],[0,1,0]])
	b = MolecularCoordinates.MolecularGeometry([1,8,1],[[0,0,0],[1,0,0],[0,1,0]])
	return MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)

def LongLinearChainConf(n):
	#n represents the number of atoms in the chain
	#all the atom-to-atom distances should be unique by the "binary" construction of the coordinates
	atomtypes = [1 for i in range(0,n)]
	atomcoor = [[2**i,0,0] for i in range(0,n)]
	a = MolecularCoordinates.MolecularGeometry(atomtypes,atomcoor)
	b = MolecularCoordinates.MolecularGeometry(atomtypes,atomcoor)
	return MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.5)

if __name__ == '__main__':
	from timeit import Timer
	startup = """import MolecularCoordinates
import conftest
	"""
	test1 = "(q, n) = conftest.SimpleHetConf()"
	test2 = "(q, n) = conftest.SimpleHomConf()"
	test3 = "(q, n) = conftest.MirrorImageConf()"
	test4 = "(q, n) = conftest.Buckminsterfullerene()"
	test5 = "(q, n) = conftest.DistinctJP10Conf()"
	test6 = "(q, n) = conftest.LongLinearChainConf(75)"
	test7 = "(q, n) = conftest.LongLinearChainConf(150)"
	t = Timer(test1,startup)
	times = t.repeat(repeat=5,number=1000)
	print "test1 took %.3f seconds (%s)"%(min(times), times)
	t = Timer(test2, startup)
	times = t.repeat(repeat=5,number=10)
	print "test2 took %.3f seconds (%s)"%(min(times), times)
	t = Timer(test3, startup)
	times = t.repeat(repeat=5,number=100)
	print "test3 took %.3f seconds (%s)"%(min(times), times)
	t = Timer(test4, startup)
	times = t.repeat(repeat=1,number=1)
	print "test4 took %.3f seconds (%s)"%(min(times), times)
	t = Timer(test5, startup)
	times = t.repeat(repeat=5,number=100)
	print "test5 took %.3f seconds (%s)"%(min(times), times)
	t = Timer(test6, startup)
	times = t.repeat(repeat=3,number=1)
	print "test6 took %.3f seconds (%s)"%(min(times), times)
	t = Timer(test7, startup)
	times = t.repeat(repeat=3,number=1)
	print "test7 took %.3f seconds (%s)"%(min(times), times)
	for i in range(1,301):
		t = Timer("(q, n) = conftest.LongLinearChainConf(%s)"%(i), startup)
		times = t.repeat(repeat=1,number=1)
		print (min(times))
	print "\nContinuing with tests..."
	unittest.main(testRunner = unittest.TextTestRunner(verbosity=2))
