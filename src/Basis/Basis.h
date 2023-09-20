#pragma once 

#include <vector>

#include "shell.h"

class Basis {
public:
	Basis(int* funcMap, int lmax, int* atomMap, int nShells, double* pExp, double* cc, double* r, bool contracted, bool spherical);
	~Basis();

	Basis(const Basis&) = delete;
	Basis& operator=(const Basis&) = delete;

	int nShells() const; 
	Shell* getShell(int i) const; 
	int nSpherical() const; 
	int nCartesian() const; 
	int nFunctions() const; 
	bool isSpherical() const;

private:
	void createContractedShells(int* funcMap, int* atomMap, int& nShls, double* pExp, double* cc, double* r);
	void createPrimitiveShells(int* funcMap, int* atomMap, int& nShls, double* pExp, double* cc, double* r);
	void setC2S();
	double C2Sco(int& lx, int& ly, int& lz, int& m, double* fac);


	int _lmax = 0; 
	int _nSph = 0;
	int _nCar = 0; 
	int _c2shDim[2] = { 0, 0 };
	double* _c2sh = 0; 
	bool _contracted; 
	bool _spherical; 
	std::vector<Shell*> _orbShells;

};