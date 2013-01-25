// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
printOpeningRemarks(const std::string& lineSearchName) const
{
  if (this->isPrintType(NOX::Utils::InnerIteration)) 
    {
      this->out() << "\n" << NOX::Utils::fill(72) << "\n" 
		  << "-- " << lineSearchName << " -- \n";
    }
}


void NOX::LineSearch::Utils::Printing::
printStep(int n, double step, double oldf, double newf, const std::string s, 
	  bool unscaleF) const
{
  if (isPrintType(NOX::Utils::InnerIteration)) 
  {
    this->out() << std::setw(3) << n << ":";
    this->out() << NOX::Utils::fill(1,' ') << "step = " << sciformat(step);
    if (unscaleF == true) {
      this->out() << NOX::Utils::fill(1,' ') << "old f = " 
		  << sciformat(std::sqrt(2. * oldf));
      this->out() << NOX::Utils::fill(1,' ') << "new f = " 
		  << sciformat(std::sqrt(2. * newf));
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
    this->out() << std::endl;
  }
}
