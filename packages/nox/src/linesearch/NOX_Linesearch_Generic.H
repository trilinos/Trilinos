// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_LINESEARCH_GENERIC_H
#define _NOX_LINESEARCH_GENERIC_H

#include "NOX_Parameter_List.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"

namespace NOX {
namespace Linesearch {

//! Generic line search interface
class Generic {

public:

  //! Default constructor
  Generic() {};

  //! Destructor
  virtual ~Generic() {};

  //! Reset parameters
  virtual void reset(const Parameter::List& params) = 0;
  
  //! Perform a line search from oldgrp along dirction dir. newgrp = oldgrp + step * dir.
  virtual bool operator()(Abstract::Group& newgrp, double& step,
			  const Abstract::Group& oldgrp, 
			  const Abstract::Vector& dir) const = 0;

};
};
};
#endif
