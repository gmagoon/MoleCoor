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
		self.assertEqual(n, '24 or 48?')

if __name__ == '__main__':
	unittest.main()

