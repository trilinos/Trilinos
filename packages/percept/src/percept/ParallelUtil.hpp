// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_ParallelUtil_hpp
#define percept_ParallelUtil_hpp

#include <utility>

#include <percept/Util.hpp>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif


  namespace percept { 


    /** Copied/Modeled after MPIX.[hc]pp in ~/code/stk/src/rsmesh and is an extension to the Sierra framework's
     *    Fmwk_global_min capability.
     */

    //========================================================================================================================
    // external interface
    namespace {

      template<class T>
      inline
      void percept_global_lex_min(stk::ParallelMachine comm,  int n , T local_min[] , T global_min[] );
    }

  } // percept

//========================================================================================================================
// implementation

#include <percept/ParallelUtilDef.hpp>

#endif
