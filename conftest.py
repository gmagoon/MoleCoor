import unittest
import MolecularCoordinates


class  ConfTestCase(unittest.TestCase):

	def testSimpleHetConf(self):
		""" A simple test (using the same molecule) of the heterogeneous case

		A pyramidal form for CH3
		"""
		a = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
		b = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)
		self.assertEqual(q, True)
		self.assertEqual(n, 6)

	def testSimpleHomConf(self):
		""" A simple test (using the same molecule) of the homogenous case

		A cube of hydrogens
		"""
		a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)
		self.assertEqual(q, True)
		self.assertEqual(n, 48) #I belive the answer is 48, though I am not positive off the top of my head, and it could be higher

	def testMirrorImageConf(self):
		""" A test of mirror images

		Mirror image conformations of gauche n-butane
		"""
		#a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		#b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)
		self.assertEqual(q, True)

	def testBuckminsterfullerene(self):
		""" A test of Buckminsterfullerene with itself


		"""
		#a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		#b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)
		self.assertEqual(q, True)

	def testDistinctJP10Conf(self):
		""" A test of distinct JP-10 conformations

		The conformations differ in the puckering of one of the rings
		"""
		#a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		#b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.05)
		self.assertEqual(q, False)

	def testOptConf(self):
		""" A test of equivalent conformations obtained from optimization with different potential energy calculation methods

		"""
		#a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		#b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		#alternatively, we could play around with Atol/Rtol to see how tight they need to be
		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.05)
		print q

	def testInitialConditionConf(self):
		""" A test of equivalent conformations of a large molecule optimized with different initial conditions

		The same potential energy calculation method is used in both cases. Default convergence criterion of Gaussian03 are used.
		"""
		#a = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		#b = MolecularCoordinates.MolecularGeometry([1,1,1,1,1,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]])
		#it may be appropriate to use Rtol here, particularly for a large molecule
		(q, n) = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.05)
		self.assertEqual(q, True)


if __name__ == '__main__':
	unittest.main(testRunner = unittest.TextTestRunner(verbosity=2))

