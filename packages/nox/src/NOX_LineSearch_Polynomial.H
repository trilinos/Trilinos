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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_LINESEARCH_POLYNOMIAL_H
#define NOX_LINESEARCH_POLYNOMIAL_H

#include "NOX_LineSearch_Generic.H" // base class

namespace NOX {
namespace LineSearch {

//! Use quad/cubic polynomial fit from Dennis & Schnabel
/*!

  \todo Russ needs to add a citation for the source of this line
  search. Is it Alg 6.3.1mod from D&S? We also need to add a basic
  description of the search.
  
  \todo Does this work with non-Newton directions?

  This line search can be called via NOX::LineSearch::Manager.
 
  The following parameters can be specified for this line search:
 
  - "Default Step" - starting step length (defaults to 1.0)

  - "Minimum Step" - minimum acceptable step length (defaults to 1.0e-12)

  - "Recovery Step" - step to take when the line search fails
    (defaults to value for "Default Step")

  - "Max Iters" - maximum number of iterations (i.e., RHS computations)

*/

class Polynomial : public Generic {

public:

  //! Constructor
  Polynomial(Parameter::List& params);

  //! Destructor
  ~Polynomial();
  
  // derived
  bool reset(Parameter::List& params);

  // derived
  bool compute(Abstract::Group& newgrp, double& step, 
		  const Abstract::Group& oldgrp, 
		  const Abstract::Vector& dir) ;
  
private:

   //! Minimum step length (i.e., when we give up)
  double minstep;
 
 //! Default step
  double defaultstep;

  //! Default step for linesearch failure
  double recoverystep;

  //! Maximum itersions
  int maxiters;

};
} // namespace LineSearch
} // namespace NOX
#endif
