// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_LINESEARCHMANAGER_H
#define _NLS_LINESEARCHMANAGER_H

#include "NLS_LineSearch.H"

//! Manager for all different line searches
class NLS_LineSearchManager : public NLS_LineSearch {

public:

  //! Constructor
  NLS_LineSearchManager(const NLS_ParameterList& params);

  //! Destructor
  ~NLS_LineSearchManager();

  // derived
  void reset(const NLS_ParameterList& params);

  // derived
  bool operator()(NLS_Group& newgrp, double& step, 
		  const NLS_Group& oldgrp, const NLS_Vector& dir) const;

private:

  //! Update method and return true if a new method name is specified.
  bool isNewMethod(const NLS_ParameterList& params);

  //! Allocate a method of the type specified by method
  void newPtr(const NLS_ParameterList& params);

private:
  //! Method name
  string method;

  //! Pointer to method
  NLS_LineSearch* ptr;

};

#endif
