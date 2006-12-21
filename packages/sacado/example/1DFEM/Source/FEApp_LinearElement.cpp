// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "FEApp_LinearElement.hpp"

FEApp::LinearElement::
LinearElement() :
  xl(0.0),
  xr(0.0),
  left_GID(0),
  right_GID(0)
{
}

FEApp::LinearElement::
~LinearElement()
{
}

unsigned int
FEApp::LinearElement::
numNodes() const
{
  return 2;
}

void
FEApp::LinearElement::
createNodes(double x_left, double x_right, unsigned int first_node_gid)
{
  xl = x_left;
  xr = x_right;
  left_GID = first_node_gid;
  right_GID = first_node_gid+1;
}

unsigned int
FEApp::LinearElement::
nodeGID(unsigned int i) const 
{
  if (i == 0)
    return left_GID;
  else
    return right_GID;
}

void
FEApp::LinearElement::
evaluateShapes(const std::vector<double>& xi,
	       std::vector< std::vector<double> >& phi) const
{
  for (unsigned int i=0; i<xi.size(); i++) {
    if (phi[i].size() < 2)
      phi[i].resize(2);

    phi[i][0] = 0.5 * (1.0 - xi[i]);
    phi[i][1] = 0.5 * (1.0 + xi[i]);
  }
}

void
FEApp::LinearElement::
evaluateShapeDerivs(const std::vector<double>& xi,
		    std::vector< std::vector<double> >& dphidxi) const
{
  for (unsigned int i=0; i<xi.size(); i++) {
    if (dphidxi[i].size() < 2)
      dphidxi[i].resize(2);

    dphidxi[i][0] = -0.5;
    dphidxi[i][1] =  0.5;
  }
}

void
FEApp::LinearElement::
evaluateJacobian(const std::vector<double>& xi, std::vector<double>& jac) const
{
  double j = 0.5 * (xr-xl);
  for (unsigned int i=0; i<xi.size(); i++)
    jac[i] = j;
}
