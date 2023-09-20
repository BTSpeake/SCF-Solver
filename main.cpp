#include "src/Molecule/Molecule.h"
#include "src/Basis/Basis.h"
#include "src/Utils/Vector3.h"
#include "src/Integrals/Integrals.h"

#include <iostream>

void generateMolecule(Molecule& mol);

int main()
{
	std::cout << "...SCF Solver..." << std::endl;

	Molecule mol;
	
	int funcMap[4] = { 2, 1, 2, 1 };
	int atomMap[4] = { 1, 1, 2, 2 };
	int lmax = 0;
	int nShells = 4;
	double pExp[6] = {
		5.44717800E+00,  8.24547000E-01,  1.83192000E-01,  5.44717800E+00,  8.24547000E-01,
		1.83192000E-01
	};
	double cc[6] = {
		1.56284981E-01,  9.04690888E-01,  1.00000000E+00,  1.56284981E-01,  9.04690888E-01,
		1.00000000E+00
	};
	double r[12] = {
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 1.3, 0.0, 0.0, 1.3
	};
	
	bool contracted = false; 
	bool spherical = false;
	Basis bas(funcMap, lmax, atomMap, nShells, pExp, cc, r, contracted, spherical);
	// Basis input needs generalising to link with mol object and read in from a formatted basis file

	generateMolecule(mol);
	double nnr = mol.calculateNNRepulsion();

	Integrals intObj(bas, mol);


	std::cout << "NNRep = " << nnr << std::endl;
	std::cout << "nShells = " << bas.nShells() << std::endl;
	std::cout << "nBas = " << bas.nFunctions() << std::endl;


}


void generateMolecule(Molecule& mol) {
	// Create a H2 atom 
	mol.addAtom(1, 0.0, 0.0, 0.0);
	mol.addAtom(1, 0.0, 0.0, 1.3);
}
