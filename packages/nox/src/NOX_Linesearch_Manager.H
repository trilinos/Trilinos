// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_LINESEARCH_MANAGER_H
#define _NOX_LINESEARCH_MANAGER_H

#include "NOX_Linesearch_Generic.H"

namespace NOX {
namespace Linesearch {

//! Manager for all different line searches
class Manager : public Generic {

public:

  //! Constructor
  Manager(const Parameter::List& params);

  //! Destructor
  ~Manager();

  // derived
  void reset(const Parameter::List& params);

  // derived
  bool operator()(Abstract::Group& newgrp, double& step, 
		  const Abstract::Group& oldgrp, const Abstract::Vector& dir) ;

private:

  //! Update method and return true if a new method name is specified.
  bool isNewMethod(const Parameter::List& params);

  //! Allocate a method of the type specified by method
  void newPtr(const Parameter::List& params);

private:
  //! Method name
  string method;

  //! Pointer to method
  Generic* ptr;

};
};
};
#endif
