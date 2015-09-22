/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file mrtr_solver.H
 *
 * \class MOERTEL::Solver
 *
 * \brief A class to solve mortar constraint problems
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifndef MOERTEL_SOLVER_H
#define MOERTEL_SOLVER_H

#include <ctime>
#include <iostream>
#include <vector>

// Trilinos includes
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Amesos.h"
#include "Amesos_BaseSolver.h"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"

#include "AztecOO.h"

//#include "mrtr_manager.H"
#include "mrtr_utils.H"
#include "mrtr_ml_preconditioner.H"

/*!
\brief MOERTEL: namespace of the Moertel package

The Moertel package depends on \ref Epetra, \ref EpetraExt, \ref Teuchos,
\ref Amesos, \ref ML and \ref AztecOO:<br>
Use at least the following lines in the configure of Trilinos:<br>
\code
--enable-moertel 
--enable-epetra 
--enable-epetraext
--enable-teuchos 
--enable-ml
--enable-aztecoo --enable-aztecoo-teuchos 
--enable-amesos
\endcode

*/


namespace MOERTEL
{
// forward declarations
class Manager;

/*!
\class Solver

\brief <b>The solver interface to Amesos, ML and Aztec</b>


\author Glen Hansen (gahanse@sandia.gov)

*/
class Solver 
{
public:

  // @{ \name Constructors and destructors
  
  /*!
  \brief 
  
  */
  explicit Solver(Epetra_Comm& comm, int outlevel);
  
  /*!
  \brief Destroys an instance of this class
  
  Destructor
  
  */
  virtual ~Solver();
  
  //@}
  
  // @{ \name Query methods
  
  /*!
  \brief Returns the Epetra_Comm object associated with this class
  */
  Epetra_Comm& Comm() const {return comm_;}

  /*!
  \brief Returns the Level of output (0 - 10) the user specified in the constructor
  */
  int OutLevel() { return outlevel_; }

  //@}

  // @{ \name Solution methods
  
  /*!
  \brief Set solver parameters
  */
  void SetParameters(Teuchos::ParameterList *params) { params_ = params; return;}

  /*!
  \brief Set linear system
  */
  void SetSystem(Teuchos::RCP<Epetra_CrsMatrix> matrix,
                 Teuchos::RCP<Epetra_Vector> x,
                 Teuchos::RCP<Epetra_Vector> b);

  /*!
  \brief Solve system
  */
  bool Solve(Teuchos::RCP<Teuchos::ParameterList> params,
             Teuchos::RCP<Epetra_CrsMatrix> matrix,
             Teuchos::RCP<Epetra_Vector> x,
             Teuchos::RCP<Epetra_Vector> b,
             MOERTEL::Manager& manager);
  
  //@}

private:  

  bool Solve();
  bool Solve_Amesos(Teuchos::ParameterList& amesosparams);
  bool Solve_MLAztec(Teuchos::ParameterList& mlparams, Teuchos::ParameterList& aztecparams);

  // don't want = operator and copy-ctor
  Solver operator = (const Solver& old);
  Solver(MOERTEL::Solver& old);

private:

  int                                               outlevel_;      // output level (0-10)
  Epetra_Comm&                                      comm_;          // communicator (global, contains ALL procs)
  
  Teuchos::ParameterList*                           params_;        // solution parameters

  Teuchos::RCP<Epetra_CrsMatrix>                     matrix_;        // linear system
  bool                                              matrixisnew_;   // flag indicating whether system has been solved before
  Teuchos::RCP<Epetra_Vector>                        x_;             // solution vector
  Teuchos::RCP<Epetra_Vector>                        b_;             // rhs vector
  
  Teuchos::RCP<Epetra_LinearProblem>                 linearproblem_; // the linear system
  Teuchos::RCP<Amesos_BaseSolver>                    amesossolver_;  // the amesos solver
  Teuchos::RCP<MOERTEL::Mortar_ML_Preconditioner>    mlprec_;        // the ml preconditioner
  //Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner>  mlprec_;        // the ml preconditioner
  Teuchos::RCP<AztecOO>                              aztecsolver_;   // the aztec solver

  Teuchos::RCP<Epetra_CrsMatrix>                     origmatrix_;
  Teuchos::RCP<Epetra_CrsMatrix>                     WT_;
  Teuchos::RCP<Epetra_CrsMatrix>                     B_;
  Teuchos::RCP<Epetra_CrsMatrix>                     I_;  
  Teuchos::RCP<Epetra_Map>                           Annmap_;
  Teuchos::RCP<std::map<int,int> >                        lm_to_dof_;
};

} // namespace MOERTEL

#endif // MOERTEL_SOLVER_H
