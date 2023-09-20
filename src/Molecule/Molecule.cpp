#include "Molecule.h"

Molecule::Molecule() {}

Molecule::~Molecule() {
	for (auto atom : _atoms) {
		delete atom;
	}
	_atoms.clear();
}

unsigned int Molecule::getNumberOfAtoms() const {
	return _atoms.size();
}

void Molecule::addAtom(const char* symbol, Vector3& position) {
	Atom* pAtom = new Atom(symbol, position);
	_atoms.push_back(pAtom);
}

void Molecule::addAtom(int atN, Vector3& position) {
	Atom* pAtom = new Atom(atN, position);
	_atoms.push_back(pAtom);
}

void Molecule::addAtom(Atom* atom) {
	_atoms.push_back(atom);
}

Vector3& Molecule::getAtomPosition(unsigned int i) const {
	return _atoms[i]->getPosition();
}

unsigned int Molecule::getAtomNumber(unsigned int i) const {
	return _atoms[i]->getAtomicNumber();
}