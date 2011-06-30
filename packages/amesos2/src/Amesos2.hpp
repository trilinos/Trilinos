// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include "Amesos2_config.h"
#include "Amesos2_Factory.hpp"

// Define some groups to be used in our doxygen documentation

/**
 * \defgroup amesos2_adapters Amesos2 Linear Algebra Object Adapters
 *
 * Amesos2 derives much of its flexibility from its adapters.  Amesos2
 * solver instances are templated on matrix and multivector types.  As
 * long as an Amesos2 adapter exists for those types, then Amesos2 can
 * interact with those objects.
 *
 * Amesos2 has two types of adapters:
 * - \ref amesos2_matrix_adapters "Matrix Adapters", and
 * - \ref amesos2_multivec_adapters "MultiVector Adapters"
 *
 * The adapters provide a unifying interface for Amesos2 to use,
 * regardless of the actual type of the object.  In this way, a solver
 * interface itself does not need to be changed in order to work with
 * a new linear algebra object.
 * 
 * New adapters can be created whenever the need arises.  Just contact
 * the Amesos2 developers if you have a linear algebra object that you
 * want Amesos2 to be able to interact with.
 */

/**
 * \defgroup amesos2_matrix_adapters Amesos2 Matrix Adapters
 * \ingroup amesos2_adapters
 * 
 */

/**
 * \defgroup amesos2_multivec_adapters Amesos2 MultiVector Adapters
 * \ingroup amesos2_adapters
 * 
 */

///////////////////////////////////////////////////////////////////////////

/**
 * \defgroup amesos2_solvers Amesos2 Solvers
 *
 * Perhaps the most interesting part of Amesos2 from a user's
 * perspective, but the Amesos2 solver system was designed to be
 * useful for both users and developers.  The system can be split into
 * two distinct but inter-related parts: The \ref
 * amesos2_solver_framework "solver framework", and the \ref
 * amesos2_solver_interfaces "solver interfaces".
 */

/**
 * \defgroup amesos2_solver_framework Amesos2 Solver Framework
 * \ingroup amesos2_solvers
 *
 * The Amesos2 Solver Framework provides an infrastructure that \ref
 * amesos2_solver_interfaces "concrete solver interfaces" depend upon.
 * At the same time this framework provides a sort of
 * fill-in-the-blank framework for solving a system of linear
 * equations that depends upon the solver-interfaces to fill in those
 * blanks.
 *
 * Its purpose is to abstract all the solver-independent features in
 * order to reduce the burden on those developing new solver
 * interfaces.  In doing this it reduces the amount of maintained code
 * and focuses a developers concerns.  A developer writing a new
 * solver interface does not have to be bothered to recreate
 * infrastructure but can instead rely on the framework provided by
 * Amesos2.
 *
 * \note In a way, the \ref amesos2_adapters "Amesos2 adapters" could
 * be viewed as part of the solver framework, but they are interesting
 * enough in their own right to be listed as a separate module.
 */

/**
 * \defgroup amesos2_solver_interfaces Amesos2 Solver Interfaces
 * \ingroup amesos2_solvers
 */


//////////////////////////////// Examples /////////////////////////////////

/*
 * We put all the example blocks in a single place, since doxygen doesn't do
 * an incredible job with putting links to the examples where you would like
 * them.
 */

/**
 * \example SimpleSolve.cpp
 *
 * Shows how to create an Amesos2 solver using the Amesos::create()
 * factory method interface, followed by solving a small linear
 * system.
 */

/**
 * \example SimpleSolve_File.cpp
 *
 * Shows how you could use Amesos2 with Tpetra's MatrixMarket
 * functionality.
 */

/**
 * \example SimpleSolve_WithParameters.cpp
 *
 * An example of how to set solver parameters with the Superlu interface.
 */

/**
 * \example MultipleSolves_File.cpp
 *
 * This example gives an example of how to re-solve a system after
 * having changed it, either changing the matrix itself, or changing
 * part of the RHS B.
 */
