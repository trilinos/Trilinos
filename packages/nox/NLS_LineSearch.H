// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_LINESEARCH_H
#define _NLS_LINESEARCH_H

#include "NLS_ParameterList.H"
#include "NLS_Vector.H"
#include "NLS_Group.H"

//! Abstract line search method
class NLS_LineSearch {

public:

  //! Default constructor
  NLS_LineSearch();

  //! Recommended constructor
  NLS_LineSearch(const NLS_ParameterList& params);

  //! Destructor
  virtual ~NLS_LineSearch() {};

  //! Reset parameters
  virtual void reset(const NLS_ParameterList& params) = 0;
  
  //! Perform a line search from oldgrp along dirction dir. The solution should be stored in newgrp.
  virtual bool search(const NLS_Group& oldgrp, const NLS_Vector& dir, NLS_Group& newgrp) const = 0;

};

#endif
