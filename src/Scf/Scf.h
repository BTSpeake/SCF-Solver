#pragma once

#include "../Utils/Vector3.h"
#include "../Basis/Basis.h"

class Scf {
public:

private:
	// Field vectors 
	Vector3 _B, _E, _G;
	// Basis object 
	Basis _basObj;
	// Energy values 
	double _e = 0.0; 
	// Integral matrices
	double* _ovi = 0; 
	double* _kei = 0;
	double* _nai = 0; 
	double* _erf = 0; 
	double* _eri = 0;
	double* _cmo_a = 0; 
	double* _cmo_b = 0;
	// Controls
	bool _restricted = true;
};