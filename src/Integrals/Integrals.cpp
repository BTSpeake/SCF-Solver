#include "Integrals.h"

#include <iostream>

/*
Molecular integral algorithms over complex orbitals.

References:
https://doi.org/10.1021/acs.jctc.7b00540
https://doi.org/10.1063/1.2996525
https://doi.org/10.1063/1.459751
https://doi.org/10.1063/1.450106
https://doi.org/10.1063/1.455717

*/

#include <iostream>

Integrals::Integrals(const Basis& bas, const Molecule& mol) {
	// Number of integrals of each type 
	_nBas = bas.nFunctions();
	_nS = _nBas * _nBas;
	_nK = _nS;
	_nN = _nS * mol.getNumberOfAtoms();
	_nR = _nS * _nS;

	calculateShellPairData(bas);

	// Initialise the integral matrices
	_S = new std::complex<double>[_nS];

	calculateOneElectron(bas); 

	for (int i = 0; i < _nS; i++) {
		if (i % 4 == 0) {
			std::cout << std::endl;
		}
		std::cout << _S[i].real() << ",    ";
	}
	std::cout << std::endl;
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
	double mag[3] = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < bas.nShells(); i++) {
		for (int j = 0; j < bas.nShells(); j++) {
			calculateOverlapPair(_shellPairs[(i * bas.nShells()) + j]);
		}
	}
}

void Integrals::calculateOverlapPair(const ShellPair* shlp) {
	const int li = shlp->getli();
	const int lj = shlp->getlj();

	std::complex<double>* rpi = new std::complex<double>[3 * shlp->getk()];
	std::complex<double>* rpj = new std::complex<double>[3 * shlp->getk()];
	std::complex<double>* Srr = new std::complex<double>[(lj + 1) * (li + 1) * 3 * shlp->getk()];

	for (int k = 0; k < shlp->getk(); k++) {
		int ki = 3 * k;
		rpi[ki + 0] = shlp->getrp(ki + 0) - shlp->getri(0);
		rpi[ki + 1] = shlp->getrp(ki + 1) - shlp->getri(1);
		rpi[ki + 2] = shlp->getrp(ki + 2) - shlp->getri(2);

		rpj[ki + 0] = shlp->getrp(ki + 0) - shlp->getrj(0);
		rpj[ki + 1] = shlp->getrp(ki + 1) - shlp->getrj(1);
		rpj[ki + 2] = shlp->getrp(ki + 2) - shlp->getrj(2);

		Srr[k] = 1.0;
		Srr[shlp->getk() + k] = 1.0;
		Srr[(2 * shlp->getk()) + k] = shlp->getu(k);

	}


	for (int i = 0; i <= li; i++) {
		// i index aliases
		int ii	= (i	) * 3 * shlp->getk();
		int ip1 = (i + 1) * 3 * shlp->getk();
		int im1 = (i - 1) * 3 * shlp->getk();
		if (i < li) {
			for (int k = 0; k < shlp->getk(); k++) {
				// k index alias
				int ki = 3 * k;

				Srr[ip1 + k]							= rpi[ki]						* Srr[ii + k];
				Srr[ip1 + shlp->getk() + k]			= rpi[shlp->getk() + k]			* Srr[ii + shlp->getk() + k];
				Srr[ip1 + (2 * shlp->getk()) + k]	= rpi[(2 * shlp->getk()) + k]	* Srr[ii + (2 * shlp->getk()) + k];

				if (i > 0) {
					for (int a = 0; a < i; a++) {
						Srr[ip1 + k]							+= shlp->gets(k) * Srr[im1 + k];
						Srr[ip1 + shlp->getk() + k]			+= shlp->gets(k) * Srr[im1 + shlp->getk() + k];
						Srr[ip1 + (2 * shlp->getk()) + k]	+= shlp->gets(k) * Srr[im1 + (2 * shlp->getk()) + k];
					}
				}
			}
		}

		for (int j = 0; j < lj; j++) {
			// j index aliases 
			int jj = (j) * (li + 1) * 3 * shlp->getk();
			int jp1 = (j + 1) * (li + 1) * 3 * shlp->getk();
			int jm1 = (j - 1) * (li + 1) * 3 * shlp->getk();

			for (int k = 0; k < shlp->getk(); k++) {
				int ki = 3 * k;

				Srr[jp1 + ii + k]						= rpi[ki]						* Srr[jj + ii + k];
				Srr[jp1 + ii + shlp->getk() + k]			= rpi[shlp->getk() + k]			* Srr[jj + ii + shlp->getk() + k];
				Srr[jp1 + ii + (2 * shlp->getk()) + k]	= rpi[(2 * shlp->getk()) + k]	* Srr[jj + ii + (2 * shlp->getk()) + k];

				if (i > 0) {
					for (int a = 0; a < i; a++) {
						Srr[jp1 + ii + k]						+= shlp->gets(k) * Srr[jj + im1 + k];
						Srr[jp1 + ii + shlp->getk() + k]			+= shlp->gets(k) * Srr[jj + im1 + shlp->getk() + k];
						Srr[jp1 + ii + (2 * shlp->getk()) + k]	+= shlp->gets(k) * Srr[jj + im1 + (2 * shlp->getk()) + k];
					}
				}

				if (j > 0) {
					for (int a = 0; a < j; a++) {
						Srr[jp1 + ii + k]						+= shlp->gets(k) * Srr[jm1 + ii + k];
						Srr[jp1 + ii + shlp->getk() + k]			+= shlp->gets(k) * Srr[jm1 + ii + shlp->getk() + k];
						Srr[jp1 + ii + (2 * shlp->getk()) + k]	+= shlp->gets(k) * Srr[jm1 + ii + (2 * shlp->getk()) + k];
					}
				}

			}
		}
	}

	for (int i = 0; i < shlp->nCartesiani(); i++) {
		// i index aliases
		int ix = shlp->getlmni(i, 0) * 3 * shlp->getk();
		int iy = shlp->getlmni(i, 1) * 3 * shlp->getk();
		int iz = shlp->getlmni(i, 2) * 3 * shlp->getk();
		for (int j = 0; j < shlp->nCartesianj(); j++) {
			// j index aliases 
			int jx = shlp->getlmnj(j, 0) * (li + 1) * 3 * shlp->getk();
			int jy = shlp->getlmnj(j, 1) * (li + 1) * 3 * shlp->getk();
			int jz = shlp->getlmnj(j, 2) * (li + 1) * 3 * shlp->getk();

			std::complex<double> sum = 0.0;
			for (int k = 0; k < shlp->getk(); k++) {
				sum += (Srr[jx + ix + k] * Srr[jy + iy + shlp->getk() + k] * Srr[jz + iz + (2 * shlp->getk()) + k]);
			}

			// S matrix indices 
			int si = (shlp->getCartesianStarti() + i) * _nBas;
			int sj = shlp->getCartesianStartj() + j;
			_S[si + sj] = shlp->getNormali(i) * shlp->getNomralj(j) * sum;
		}
	}


	delete[] rpi;
	delete[] rpj;
	delete[] Srr;
}