// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_LINESEARCH_HALVING_H
#define _NOX_LINESEARCH_HALVING_H

#include "NOX_Linesearch_Generic.H"

namespace NOX {
namespace Linesearch {

//! Halve step until reduction in the RHS

class Halving : public Generic {

public:

  //! Constructor
  Halving(const Parameter::List& params);

  //! Destructor
  ~Halving();
  
  // derived
  void reset(const Parameter::List& params);

  // derived
  bool operator()(Abstract::Group& newgrp, double& step, 
		  const Abstract::Group& oldgrp, 
		  const Abstract::Vector& dir) const;
  
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
