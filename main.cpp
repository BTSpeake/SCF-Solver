#include "src/Molecule/Molecule.h"
#include "src/Utils/Vector3.h"

#include <iostream>

int main()
{
	std::cout << "...SCF Solver..." << std::endl;

	Molecule mol;
	Vector3 v1;
	Vector3 v2;

	mol.addAtom("HE", v1);
	mol.addAtom(8, v2);

	std::cout << "Molecule has " << mol.getNumberOfAtoms() << " atoms" << std::endl;

	std::cout << "Atom 1: " << mol.getAtomNumber(0) << std::endl;
	std::cout << "Position: " << mol.getAtomPosition(0)[0] << " " << mol.getAtomPosition(0)[1] << std::endl;

}
