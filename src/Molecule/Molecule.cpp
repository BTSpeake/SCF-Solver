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

void Molecule::addAtom(unsigned int atN, double x, double y, double z) {
	Atom* pAtom = new Atom(atN, x, y, z);
	_atoms.push_back(pAtom);
}

const Vector3& Molecule::getAtomPosition(unsigned int i) const {
	return _atoms[i]->getPosition();
}

unsigned int Molecule::getAtomNumber(unsigned int i) const {
	return _atoms[i]->getAtomicNumber();
}

double Molecule::calculateNNRepulsion() const {
	double nnr = 0.0;
	for (int i = 0; i < _atoms.size(); i++) {
		for (int j = (i + 1); j < _atoms.size(); j++) {
			Vector3 Rij; 
			Rij = _atoms[i]->getPosition() - _atoms[j]->getPosition();
			double r = Rij.normal();
			nnr += _atoms[i]->getAtomicNumber() * _atoms[j]->getAtomicNumber() / r;
		}
	}
	return nnr;
}