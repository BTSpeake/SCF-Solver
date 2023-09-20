#pragma once 
#include "Atom.h"
#include <vector>

class Molecule {
public:
	Molecule();
	~Molecule();
	unsigned int getNumberOfAtoms() const;
	void addAtom(const char* symbol, Vector3& position);
	void addAtom(int atN, Vector3& position);
	void addAtom(Atom* atom);
	void addAtom(unsigned int atN, double x, double y, double z);
	const Vector3& getAtomPosition(unsigned int i) const;
	unsigned int getAtomNumber(unsigned int i) const;
	double calculateNNRepulsion() const;
private:
	std::vector<Atom*> _atoms;
	int charge = 0;
};