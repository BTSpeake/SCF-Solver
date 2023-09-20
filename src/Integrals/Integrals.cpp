#include "Integrals.h"

#include <iostream>

Integrals::Integrals(const Basis& bas, const Molecule& mol) {
	// Number of integrals of each type 
	_nS = bas.nFunctions() * bas.nFunctions();
	_nK = _nS;
	_nN = _nS * mol.getNumberOfAtoms();
	_nR = _nS * _nS;

	calculateShellPairData(bas);

	// Initialise the integral matrices
	_S = new std::complex<double>[_nS];
}

Integrals::~Integrals() {
	// Delete the heap allocated integral matrices
	for (auto shlp : _shellPairs) {
		delete shlp;
	}
	_shellPairs.clear();
	delete[] _S;
}

void Integrals::calculateShellPairData(const Basis& bas) {
	int nPairs = bas.nShells() * bas.nShells();
	_shellPairs.reserve(nPairs);
	for (int i = 0; i < bas.nShells(); i++) {
		for (int j = 0; j < bas.nShells(); j++) {
			_shellPairs.push_back(new ShellPair(bas.getShell(i), bas.getShell(j), bas.isSpherical()));
		}
	}
}

void Integrals::calculateOneElectron(const Basis& bas) {

}