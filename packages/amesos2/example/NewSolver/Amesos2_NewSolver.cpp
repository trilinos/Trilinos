// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Amesos2_NewSolver_decl.hpp"

#ifdef HAVE_AMESOS2_EXPLICIT_INSTANTIATION

#  include "Amesos2_NewSolver_def.hpp"
#  include "Amesos2_ExplicitInstantiationHelpers.hpp"

namespace Amesos2 {
#ifdef HAVE_AMESOS2_EPETRA
  AMESOS2_SOLVER_EPETRA_INST(NewSolver);
#endif

  // Remove explicit instantiations on the a Tpetra scalar type if it
  // is not supported by the TPL.

#ifdef HAVE_TPETRA_INST_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(NewSolver,float,int,int);
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(NewSolver,double,int,int);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  AMESOS2_SOLVER_TPETRA_INST(NewSolver,std::complex<float>,int,int);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  AMESOS2_SOLVER_TPETRA_INST(NewSolver,std::complex<double>,int,int);
#endif
}

#endif  // HAVE_AMESOS2_EXPLICIT_INSTANTIATION
