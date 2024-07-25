// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Amesos2_Factory.hpp"

namespace Amesos2 {

  /**********************
   *   QUERY function   *
   **********************/

  bool query (const char* solverName) {
    std::string solver = solverName;
    return query (solver);
  }

  bool query (const std::string& solver_name) {
    // Make the input name canonical by converting to lower case.
    const std::string solverName = tolower (solver_name);
    //
    // Go through all implemented solvers.  If the name matches one of
    // them, return true; else return false.
    //
#ifdef HAVE_AMESOS2_BASKER
    if (solverName == "amesos2_basker" ||
        solverName == "basker") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_SHYLU_NODEBASKER
    if (solverName == "amesos2_shylubasker" ||
        //solverName == "ShyLUBasker" || // unnecessary - tolower called on solver name prior
        solverName == "shylubasker") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_KLU2
    if (solverName == "amesos2_klu2" ||
        solverName == "klu2" ||
        solverName == "amesos2_klu" ||
        solverName == "klu") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUDIST
    if (solverName == "amesos2_superludist" ||
        solverName == "superludist" ||
        solverName == "amesos2_superlu_dist" ||
        solverName == "superlu_dist") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUMT
    if (solverName == "amesos2_superlumt" ||
        solverName == "superlumt" ||
        solverName == "amesos2_superlu_mt" ||
        solverName == "superlu_mt") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_UMFPACK
    if (solverName == "amesos2_umfpack" || solverName == "umfpack") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_SHYLU_NODETACHO
    if (solverName == "amesos2_tacho" || solverName == "tacho") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLU
    if (solverName == "amesos2_superlu" || solverName == "superlu") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_PARDISO_MKL
    if (solverName == "amesos2_pardiso_mkl" ||
        solverName == "pardiso_mkl" ||
        solverName == "amesos2_pardisomkl"  ||
        solverName == "pardisomkl") {
      return true;
    }
#endif
#ifdef HAVE_AMESOS2_CSS_MKL
    if (solverName == "amesos2_css_mkl" ||
        solverName == "css_mkl" ||
        solverName == "amesos2_cssmkl"  ||
        solverName == "cssmkl") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_LAPACK
    if (solverName == "amesos2_lapack" || solverName == "lapack") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_MUMPS
    if(solverName  == "MUMPS" || solverName == "mumps" ||
       solverName == "amesos2_MUMPS" || solverName == "amesos2_MUMPS")
      {
        return true;
      }
#endif

#ifdef HAVE_AMESOS2_STRUMPACK
    if(solverName == "strumpack" || solverName == "amesos2_strumpack")
      {
        return true;
      }
#endif

#if defined (HAVE_AMESOS2_CHOLMOD) && defined (HAVE_AMESOS2_EXPERIMENTAL)
    if (solverName == "amesos2_cholmod" || solverName == "cholmod") {
      return true;
    }
#endif

#ifdef HAVE_AMESOS2_CUSOLVER
    if (solverName == "amesos2_cusolver" || solverName == "cusolver") {
      return true;
    }
#endif

    // Otherwise, the solver is not available
    return false;
  }

  std::string tolower (const std::string& s)
  {
    std::string rtn = s;
    const size_t len = rtn.length ();
    for (size_t i = 0; i < len; ++i) {
      rtn[i] = ::tolower (rtn[i]);
    }
    return rtn;
  }
} // end namespace Amesos2
