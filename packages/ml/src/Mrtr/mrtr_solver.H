/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
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
 * \date Last update do Doxygen: 19-Dec-05
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
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"




using namespace std;
using namespace Teuchos;

/*!
\brief MOERTEL: namespace of the Moertel package

The Moertel package depends on \ref Epetra, \ref EpetraExt, \ref Teuchos,
\ref Amesos and \ref ML

It makes use of the Epetra define EPETRA_MPI

*/
namespace MOERTEL
{

/*!
\class Solver

\brief <b>  </b>


\author Michael W. Gee (mwgee@sandia.gov)

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
  void SetParameters(RefCountPtr<Teuchos::ParameterList> params) { params_ = params; return;}

  /*!
  \brief Set linear system
  */
  void SetSystem(RefCountPtr<Epetra_CrsMatrix> matrix,
                 RefCountPtr<Epetra_Vector> x,
                 RefCountPtr<Epetra_Vector> b);

  /*!
  \brief Solve system
  */
  bool Solve(RefCountPtr<Teuchos::ParameterList> params,
             RefCountPtr<Epetra_CrsMatrix> matrix,
             RefCountPtr<Epetra_Vector> x,
             RefCountPtr<Epetra_Vector> b);
  
  /*!
  \brief Solve system
  */
  bool Solve();
  
  //@}

private:  

  bool Solve_Amesos(ParameterList& amesosparams);

  // don't want = operator and copy-ctor
  Solver operator = (const Solver& old);
  Solver(MOERTEL::Solver& old);

private:

  int                                 outlevel_;      // output level (0-10)
  Epetra_Comm&                        comm_;          // communicator (global, contains ALL procs)
  
  RefCountPtr<Teuchos::ParameterList> params_;        // solution parameters

  RefCountPtr<Epetra_CrsMatrix>       matrix_;        // linear system
  bool                                matrixisnew_;   // flag indicating whether system has been solved before
  
  RefCountPtr<Epetra_Vector>          x_;             // solution vector
  RefCountPtr<Epetra_Vector>          b_;             // rhs vector
  
  RefCountPtr<Epetra_LinearProblem>   linearproblem_;  // the linear system

  RefCountPtr<Amesos_BaseSolver>      amesossolver_;  // the amesos solver
};

} // namespace MOERTEL

#endif // MOERTEL_SOLVER_H
