#include "Teuchos_Utils.hpp"
#include "Teuchos_Out.hpp"

using namespace Teuchos;


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



