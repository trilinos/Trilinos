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

#ifndef NOX_MORETHUENTE_H
#define NOX_MORETHUENTE_H

#include "NOX_Linesearch_Generic.H"

namespace NOX {
namespace Linesearch {

//! Use quad/cubic polynomial fit from Dennis & Schnabel

class MoreThuente : public Generic {

public:

  //! Constructor
  MoreThuente(const Parameter::List& params);

  //! Destructor
  ~MoreThuente();
  
  // derived
  void reset(const Parameter::List& params);

  // derived
  bool operator()(Abstract::Group& newgrp, double& step, 
		  const Abstract::Group& oldgrp, 
		  const Abstract::Vector& dir) ;

  //! \todo This needs to be documented
  int MTStep(double*, double*, double*, double*, double*, double*, double*,
           double*, double*, double*, double*, int*);
//   int MTStep();
  
  //! \todo This needs to be documented
  double mymin(double, double);

  //! \todo This needs to be documented
  double mymax(double, double);


private:

   //! Minimum step length (i.e., when we give up)
  double minstep;
 
 //! Default step
  double defaultstep;

  //! Default step for linesearch failure
  double recoverystep;

};
} // namespace Linesearch
} // namespace NOX
#endif
