#pragma once 

#include <fstream>

#include "../Molecule/Molecule.h"

class MoleculeParser {
public:
	MoleculeParser(const char* fname, Molecule& molObj);
	~MoleculeParser();
private:
	std::ifstream ifile;
	double AtB = 1.889725989;

	void parse(Molecule& molObj);
};