// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_MORETHUNTE_H
#define _NOX_MORETHUNTE_H

#include "NOX_Linesearch_Generic.H"

namespace NOX {
namespace Linesearch {

//! Use quad/cubic polynomial fit from Dennis & Schnabel

class MoreThunte : public Generic {

public:

  //! Constructor
  MoreThunte(const Parameter::List& params);

  //! Destructor
  ~MoreThunte();
  
  // derived
  void reset(const Parameter::List& params);

  // derived
  bool operator()(Abstract::Group& newgrp, double& step, 
		  const Abstract::Group& oldgrp, 
		  const Abstract::Vector& dir) ;

  int MTStep(double*, double*, double*, double*, double*, double*, double*,
           double*, double*, double*, double*, int*);
//   int MTStep();
  
  double mymin(double, double);

  double mymax(double, double);


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
