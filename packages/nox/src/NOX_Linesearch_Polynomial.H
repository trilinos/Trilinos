// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_LINESEARCH_POLYNOMIAL_H
#define _NOX_LINESEARCH_POLYNOMIAL_H

#include "NOX_Linesearch_Generic.H"

namespace NOX {
namespace Linesearch {

//! Use quad/cubic polynomial fit from Dennis & Schnabel

class Polynomial : public Generic {

public:

  //! Constructor
  Polynomial(const Parameter::List& params);

  //! Destructor
  ~Polynomial();
  
  // derived
  void reset(const Parameter::List& params);

  // derived
  bool operator()(Abstract::Group& newgrp, double& step, 
		  const Abstract::Group& oldgrp, 
		  const Abstract::Vector& dir) ;
  
private:

   //! Minimum step length (i.e., when we give up)
  double minstep;
 
 //! Default step
  double defaultstep;

  //! Default step for linesearch failure
  double recoverystep;

};
}
}
#endif
