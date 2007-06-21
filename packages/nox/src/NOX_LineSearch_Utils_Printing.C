// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_LineSearch_Utils_Printing.H"
#include "NOX_Utils.H"
#include "NOX_Common.H"

NOX::LineSearch::Utils::Printing::
Printing(const Teuchos::RCP<NOX::Utils>& u) :
  NOX::Utils(*u)
{

}

NOX::LineSearch::Utils::Printing::~Printing()
{

}

void NOX::LineSearch::Utils::Printing::
reset(const Teuchos::RCP<NOX::Utils>& u)
{
  NOX::Utils* tmp = this;
  tmp = u.get();
}

void NOX::LineSearch::Utils::Printing::
printOpeningRemarks(const string& lineSearchName) const
{
  if (this->isPrintType(NOX::Utils::InnerIteration)) 
    {
      this->out() << "\n" << NOX::Utils::fill(72) << "\n" 
		  << "-- " << lineSearchName << " -- \n";
    }
}


void NOX::LineSearch::Utils::Printing::
printStep(int n, double step, double oldf, double newf, const string s, 
	  bool unscaleF) const
{
  if (isPrintType(NOX::Utils::InnerIteration)) 
  {
    this->out() << setw(3) << n << ":";
    this->out() << NOX::Utils::fill(1,' ') << "step = " << sciformat(step);
    if (unscaleF == true) {
      this->out() << NOX::Utils::fill(1,' ') << "old f = " 
		  << sciformat(sqrt(2. * oldf));
      this->out() << NOX::Utils::fill(1,' ') << "new f = " 
		  << sciformat(sqrt(2. * newf));
    }
    else {
      this->out() << NOX::Utils::fill(1,' ') << "old f = " << sciformat(oldf);
      this->out() << NOX::Utils::fill(1,' ') << "new f = " << sciformat(newf);
    }
    if (!s.empty()) 
    {
      this->out() << " " << s << "\n";
      this->out() << NOX::Utils::fill(72);
    }
    this->out() << endl;
  }
}
