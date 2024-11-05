/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_QR_FIX_HPP
#define ML_QR_FIX_HPP

#include "ml_include.h"
#include <vector>

struct ML_qr_fix_struct {

  int                 numAggsWithDeadDofs;
  // Records which coarse DOFs in each aggregate should be recorded as "dead".
  // coarseDOFState[i][j] gives the state of coarse DOF j unknown in aggregate i
  std::vector< std::vector<bool> > coarseDOFState;
  // Records whether each aggregate is too small to support entire nullspace.
  // true = aggregate size is smaller than nullspace dimension
  // false = aggregate size is greater than or equal to nullspace dimension
  std::vector<bool> aggTooSmall;
};

typedef struct ML_qr_fix_struct ML_qr_fix;

#endif //ifndef ML_QR_FIX_HPP

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

