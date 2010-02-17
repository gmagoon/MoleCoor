import unittest
import MolecularCoordinates


class  ConfTestCase(unittest.TestCase):

    def simpleConfTest(self):
	a = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
	b = MolecularCoordinates.MolecularGeometry([6,1,1,1],[[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
	q = MolecularCoordinates.checkConformationalEquivalence(b, a, Atol=0.01)
	print q
        #assert x != y;
        #self.assertEqual(x, y, "Msg");

if __name__ == '__main__':
    unittest.main()

