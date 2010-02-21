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
		(q, n) = conftest.SimpleHetConf()
		self.assertEqual(q, True)
		self.assertEqual(len(n), 48) #I belive the answer is 48, though I am not positive off the top of my head, and it could be higher

	def testMirrorImageConf(self):
		""" A test of mirror images

		Mirror image conformations of gauche n-butane (each independently optimized using PM3 in Gaussian03 using default convergence criteria)
		"""
		(q, n) = conftest.MirrorImageConf()
		self.assertEqual(q, True)
		self.assertEqual(len(n), 2) #I belive the correct answer is 2

#	def testBuckminsterfullerene(self):
#		""" A test of Buckminsterfullerene with itself
#
#
#		"""
#		#a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
#		#b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
#		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)
#		self.assertEqual(q, True)
#
#	def testDistinctJP10Conf(self):
#		""" A test of distinct JP-10 conformations
#
#		The conformations differ in the puckering of one of the rings
#		"""
#		#a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
#		#b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
#		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.05)
#		self.assertEqual(q, False)
#
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

if __name__ == '__main__':
	from timeit import Timer
	startup = """import MolecularCoordinates
import conftest
	"""
	test1 = "(q, n) = conftest.SimpleHetConf()"
	test2 = "(q, n) = conftest.SimpleHomConf()"
	test3 = "(q, n) = conftest.MirrorImageConf()"
	t = Timer(test1,startup)
	times = t.repeat(repeat=5,number=1000)
	print "test1 took %.3f seconds (%s)"%(min(times), times)
	t = Timer(test2, startup)
	times = t.repeat(repeat=5,number=10)
	print "test2 took %.3f seconds (%s)"%(min(times), times)
	t = Timer(test3, startup)
	times = t.repeat(repeat=5,number=100)
	print "test3 took %.3f seconds (%s)"%(min(times), times)
	print "\nContinuing with tests..."
	unittest.main(testRunner = unittest.TextTestRunner(verbosity=2))
