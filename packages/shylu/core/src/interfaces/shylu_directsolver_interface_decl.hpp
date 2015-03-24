//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// ************************************************************************
//@HEADER

/** \file shylu_directsolver_interface_decl.hpp

    \brief Eperta/Tpetra templated interface for call Amesos and Amesos2

    \author Joshua Dennis Booth
*/
#ifndef SHYLU_DIRECTSOLVER_INTERFACE_DECL_HPP
#define SHYLU_DIRECTSOLVER_INTERFACE_DECL_HPP

#include "ShyLUCore_config.h"

//ShyLU
#include <shylu.h>
#include <shylu_util.h>

//Epetra
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_LinearProblem.h>

//Isorropia
#include <Isorropia_config.h>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraProber.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>

//Tperta
#ifdef HAVE_SHYLUCORE_TPETRA
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#endif

//#include <Zoltan2_config.h>
#ifdef HAVE_SHYLUCORE_ZOLTAN2
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#endif

#include <Teuchos_XMLParameterListHelpers.hpp>

//Amesos
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

#ifdef HAVE_SHYLUCORE_AMESOS2
#include <Amesos2.hpp>
#endif

namespace ShyLU{

  /** \brief DirectSolverInterface class templated on Epetra/Tpetra Matrix and Vector
   *
   * This class acts as an interface that will allow Shylu to call either Amesos/Amesos2 
   * without having to address if the matrix/submatrix is either of Epetra or Tpetra form.
   * Currently:  Only supporting limited solvers in both Amesos and Amesos2, will be updated
   */
template <class Matrix, class Vector>
class DirectSolverInterface
{
public:
  ~DirectSolverInterface();

  /** \brief Main constructor of class
   *
   * This constructor requires a Teuchos ParameterList that provides information on solver.
   * It assumes that if Tpetra matrix is given then Amesos2 ParameterList must be given.
   * Likewise, if Epetra matrix is given then Amesos ParameterList must be given.
   */
  DirectSolverInterface(Matrix *inA, Teuchos::ParameterList* pList);
  int solve(Vector* b, Vector* x);
private:
  int solveAmesos(Vector* b, Vector *x );
  Teuchos::ParameterList *pList;
  Matrix *A;
  int maxproc;
#ifdef HAVE_SHYLUCORE_AMESOS2
  int solveAmesos2(Vector* b, Vector *x);
#endif
}; //end class

}// end namespace
#endif //end header if
