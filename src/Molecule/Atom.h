#pragma once 
#include "../Utils/Vector3.h"

class Atom {
public:
	Atom(const char* symbol, Vector3& position);
	Atom(const unsigned int atN, Vector3& position);
	Atom(const unsigned int atN, double x, double y, double z);
	~Atom();
	const Vector3& getPosition() const;
	const char* getSymbol() const;
	unsigned int getAtomicNumber() const;
private:
	Vector3 _position;
	unsigned int _atN;

	void setAtomicNumberFromSymbol(const char* symbol);
	const char* numberToSymbol();
};