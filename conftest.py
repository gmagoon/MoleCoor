import unittest
import MolecularCoordinates


class  ConfTestCase(unittest.TestCase):

	def testSimpleConf(self):
		a = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
		b = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
		q = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)
		print q
		assert 1 != 2;
		#self.assertEqual(x, y, "Msg");
		#self.fail("TODO: Write test")

if __name__ == '__main__':
	a = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
	b = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
	q = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)
	print q
	unittest.main()

