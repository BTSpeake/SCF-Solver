#include "MoleculeParser.h"

#include <string>

MoleculeParser::MoleculeParser(const char* fname, Molecule& molObj) {
	ifile.open(fname);
	if (ifile.is_open()) {
		parse(molObj);
	}
}

MoleculeParser::~MoleculeParser() {
	if (ifile.is_open()) {
		ifile.close();
	}
}

void MoleculeParser::parse(Molecule& molObj) {
	unsigned int lc = 0;
	std::string line;
	std::getline(ifile, line);
	std::getline(ifile, line);

	unsigned int atN = 0;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	//while (std::getline(ifile, line)) {
	//	std::istringstream iss(line);
	//}
}