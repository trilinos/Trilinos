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

#ifndef NOX_LINESEARCH_UTILS_PRINTING_H
#define NOX_LINESEARCH_UTILS_PRINTING_H

#include "NOX_LineSearch_Generic.H"
#include "NOX_Utils.H"
#include "NOX_Common.H"

// Forward declaration
namespace NOX {
  class Utils;
}

namespace NOX {

namespace LineSearch {

namespace Utils {

//! %Common line search utilites for printing line search information to the screen.
/*!  
  All line searches should print output results in a similar format.
  This utility provides common output routines.
*/

class Printing : public NOX::Utils {

public:

  //! Default constructor
  Printing(const Teuchos::RCP<NOX::Utils>& u);

  //! Destructor
  virtual ~Printing();

  void reset(const Teuchos::RCP<NOX::Utils>& u);

  //! Prints the opening information.
  void printOpeningRemarks(const std::string& lineSearchName) const;

  //! Print out step information for the inner iterations of a line search algorithm.
  /*!
    Example of output from the inner iterations of a Polynomial line search:

\verbatim
************************************************************************
-- Polynomial Line Search -- 
  1: step = 1.000e+00 oldf = 2.403e+00 newf = 1.076e+03
  2: step = 1.000e-01 oldf = 2.403e+00 newf = 4.440e+00
  3: step = 1.000e-02 oldf = 2.403e+00 newf = 2.394e+00 (STEP ACCEPTED!)
************************************************************************
\endverbatim
    
    \param unscaleF - If this is true (the default), than the values
    printed are \f$ \sqrt{2 * {\rm oldf}} \f$ and \f$ \sqrt{2 * {\rm
    newf}} \f$. This is to accomodate the standard merit function, \f$
    \phi(x) = \frac{1}{2} \|F(x)\|^2 \f$.

  */
  void printStep(int n, double step, double oldf, double newf, 
		 const std::string s = "", bool unscaleF = true) const;

};
} // namespace Utils
} // namespace LineSearch
} // namespace NOX

#endif
