/**
 * \file   Amesos2_TypeDecl.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Mon Jul 18 11:57:14 2011
 * 
 * \brief  Enum and other types declarations for Amesos2
 * 
 * 
 */

#ifndef AMESOS2_TYPEDECL_HPP
#define AMESOS2_TYPEDECL_HPP

namespace Amesos2 {

  typedef enum {
    CLEAN = 0,
    PREORDERING,
    SYMBFACT,
    NUMFACT,
    SOLVE
  } EPhase;

}

#endif	// AMESOS2_TYPEDECL_HPP
