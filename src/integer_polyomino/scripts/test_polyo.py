import sys, os

sys.path.append("/home/phd/Projects/Github/polyomino_integer/module")

import polyo

def test_tetramer():
    assert(polyo.AssemblePlasticGenotype([0,0,1,2]) == [(4,0)])

def test_dimer():
    assert(polyo.AssemblePlasticGenotype([0,0,0,1,0,0,0,2]) == [(2,0)])

def test_trimer():
    assert(polyo.AssemblePlasticGenotype([0,0,0,1,0,0,2,2]) == [(3,0)])

if __name__ == '__main__':
    test_dimer()
    test_trimer()
    test_tetramer()
