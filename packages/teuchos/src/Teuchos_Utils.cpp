// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_Utils.hpp"

namespace Teuchos {

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

string Utils::toString(const unsigned int& x)
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
