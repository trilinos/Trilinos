// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_LINESEARCH_FULLSTEP_H
#define _NOX_LINESEARCH_FULLSTEP_H

#include "NOX_Linesearch_Generic.H"

namespace NOX {
namespace Linesearch {

//! Simplest line search - always take the full step.
class FullStep : public Generic {

public:

  //! Constructor
  FullStep(const Parameter::List& params);

  //! Destructor
  ~FullStep();
  
  // derived
  void reset(const Parameter::List& params);

  // derived
  bool operator()(Abstract::Group& newgrp, double& step, 
		  const Abstract::Group& oldgrp, 
		  const Abstract::Vector& dir) const;
private:

  //! Full step length
  double defaultstep;

};
}
}

#endif
