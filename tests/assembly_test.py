import unittest

from integer_polyomino import integer_polyomino as ip


class MainTest(unittest.TestCase):
	def test_assembly_dimer(self):
		self.assertEqual(ip.AssemblePlasticGenotype([0,0,0,1,0,0,0,2], **parameters), [(2, 0)])

	def test_assembly_dimer_frequency(self):
		self.assertEqual(ip.AssemblePlasticGenotypeFrequency([0,0,0,1,0,0,0,2]), {(2,0): parameters['phenotype_builds']})

	def test_assembly_tetramer(self):
		self.assertEqual(ip.AssemblePlasticGenotype([0,0,1,2], **parameters)


if __name__ == "__main__":
	unittest.main()
