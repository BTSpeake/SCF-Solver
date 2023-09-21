#pragma once 

#include <complex>
#include <vector>

#include "../Basis/Basis.h"
#include "../Molecule/Molecule.h"
#include "../Basis/ShellPair.h"

class Integrals {
public:
	Integrals(const Basis& bas, const Molecule& mol);
	~Integrals();
private:
	int _nBas = 0;
	int _nS = 0; 
	int _nK = 0;
	int _nN = 0;
	int _nR = 0; 
	std::vector<ShellPair*> _shellPairs;
	std::complex<double>* _S = 0;

	void calculateShellPairData(const Basis& bas);
	void calculateOneElectron(const Basis& bas);
	void calculateOverlapPair(const ShellPair* shlp);
};