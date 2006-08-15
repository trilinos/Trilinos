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
  Printing(const Teuchos::RefCountPtr<NOX::Utils>& u);

  //! Destructor
  virtual ~Printing();

  void reset(const Teuchos::RefCountPtr<NOX::Utils>& u);

  //! Prints the opening information.
  void printOpeningRemarks(const string& lineSearchName) const;

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
		 const string s = "", bool unscaleF = true) const;

};
} // namespace Utils
} // namespace LineSearch
} // namespace NOX

#endif
