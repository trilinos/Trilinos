#include "Teuchos_Utils.hpp"
#include "Teuchos_Out.hpp"

namespace Teuchos
{

// These characters are needed for BLAS & LAPACK routines for compile time debugging
char ESideChar[] = {'L' , 'R' };
char ETranspChar[] = {'N' , 'T' , 'C' };
char EUploChar[] = {'U' , 'L' };
char EDiagChar[] = {'U' , 'N' };
char EFactChar[] = {'F', 'N' };
char ENormChar[] = {'O', 'I' };
char ECompQChar[] = {'N', 'I', 'V' };
char EJobChar[] = {'E', 'V', 'B' };
char EJobSChar[] = {'E', 'S' };
char EJobVSChar[] = {'V', 'N' };
char EHowmnyChar[] = {'A', 'S' };
char ECMachChar[] = {'E', 'S', 'B', 'P', 'N', 'R', 'M', 'U', 'L', 'O' };
char ESortChar[] = {'N', 'S'};

double Utils::chopVal_ = 1.0e-16;


double Utils::chop(const double& x) 
{
	if (fabs(x) < chopVal_) return 0;
	return x;
}

string Utils::toString(const int& x)
{
	char s[100];
	sprintf(s, "%d", x);
	return string(s);
}

string Utils::toString(const double& x)
{
	char s[100];
	sprintf(s, "%g", x);
	return string(s);
}


} // end namespace Teuchos
