// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_INTERVALHALVING_H
#define _NLS_INTERVALHALVING_H

#include "NLS_LineSearch.H"

//! Abstract line search method
class NLS_IntervalHalving : public NLS_LineSearch {

public:

  //! Constructor
  NLS_IntervalHalving(const NLS_ParameterList& params);

  //! Destructor
  ~NLS_IntervalHalving();
  
  // derived
  void reset(const NLS_ParameterList& params);

  // derived
  bool search(const NLS_Group& oldgrp, const NLS_Vector& dir, NLS_Group& newgrp) const;

private:

  //! Full step length
  double maxstep;
  //! Minimum step length (i.e., when we give up)
  double minstep;
};

#endif
