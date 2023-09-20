#pragma once 

#include <complex>

#include "Shell.h"

class ShellPair {
public:
	ShellPair(Shell* shli, Shell* shlj, bool spherical);
	~ShellPair();

	Shell* _shli; 
	Shell* _shlj;

	int l; // la , lb
	int k;
	std::complex<double>* r = 0; // [3, kp]
	std::complex<double>* u = 0; // [kp]
	std::complex<double>* s = 0; // [kp] 
	// kp = shli->nExp * shlj->nExp
	
	std::complex<double> A[3];
	int np; //(lp+1)*(lp+2)*(lp+3)//6 - ioa   ioa = (la+0)*(la+1)*(la+2)//6

	std::complex<double>* sa = 0; 
	std::complex<double>* sb = 0; // [kp] --- is this needed or can we link to shli/shlj

	//vA - shli.kpw   vB - shlj.kpw 
	std::complex<double> k1[3]; // difference in v (vB - vA)
	// ra/rb = shli->r shlj->r
	bool occ = true; // true if rab^2 == 0 and k1^2 == 0

	// xa/xa is the 'slice' from the start function index to the final function index 

private:
	void calculatePrimitive();
	void calculateContracted();
};
