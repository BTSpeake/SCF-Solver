#pragma once 

#include <complex>

#include "Shell.h"


// Optimisation -> primitive basis fuinctions wouldn't require heap allocation of variables
// Optimisation? -> make arrays std::vectors ??

class ShellPair {
public:
	ShellPair(Shell* shli, Shell* shlj, bool spherical);
	~ShellPair();
	const int getli() const;
	const int getlj() const;
	int getk() const;

	std::complex<double> getrp(int i) const;
	double getri(int i) const;
	double getrj(int i) const;
	std::complex<double> getu(int i) const;
	std::complex<double> gets(int i) const;
	int nCartesiani() const;
	int nCartesianj() const;
	int getlmni(int k, int i) const;
	int getlmnj(int j, int i) const;
	const double& getNormali(int i) const;
	const double& getNomralj(int i) const;
	int getCartesianStarti() const;
	int getCartesianStartj() const;
	

private:
	void calculatePrimitive();
	void calculateContracted();
	void calculateHrrMatrix();

	Shell* _shli;
	Shell* _shlj;

	int _l; // la , lb
	int _k;
	std::complex<double>* _r = 0; // [3, kp]
	std::complex<double>* _u = 0; // [kp]
	std::complex<double>* _s = 0; // [kp] 
	// kp = shli->nExp * shlj->nExp

	std::complex<double> _A[3];
	int np; //(lp+1)*(lp+2)*(lp+3)//6 - ioa   ioa = (la+0)*(la+1)*(la+2)//6

	std::complex<double>* _sa = 0;
	std::complex<double>* _sb = 0; // [kp] --- is this needed or can we link to shli/shlj

	//vA - shli.kpw   vB - shlj.kpw 
	std::complex<double> _k1[3]; // difference in v (vB - vA)
	// ra/rb = shli->r shlj->r
	bool _occ = true; // true if rab^2 == 0 and k1^2 == 0

	// xa/xa is the 'slice' from the start function index to the final function index
};
