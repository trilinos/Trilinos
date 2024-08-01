// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_ConfigDefs.hpp"
#if defined(HAVE_IFPACK2_HYPRE) && defined(HAVE_IFPACK2_MPI)

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE.h"
#include <map>
#include "Ifpack2_Hypre_FunctionParameters.hpp"

namespace Ifpack2 {

  #include "Ifpack2_HypreParameterMap.hpp"

  void IFPACK2_CHK_ERRV(int code) {
    if(code<0) {
      std::ostringstream ofs;
      ofs << "Ifpack2::Hypre: Error with code "<<code<<std::endl;
      throw std::runtime_error(ofs.str());
    }
  }

  void IFPACK2_CHK_ERR(int code) {
    if(code<0) {
      std::ostringstream ofs;
      ofs << "Ifpack2::Hypre: Error with code "<<code<<std::endl;
      throw std::runtime_error(ofs.str());
    }
  }

}

#endif // HAVE_IFPACK2_HYPRE && HAVE_IFPACK2_MPI
