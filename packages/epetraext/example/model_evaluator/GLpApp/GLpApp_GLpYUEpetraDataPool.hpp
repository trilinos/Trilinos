/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
//@HEADER
*/

#ifndef GLPAPP_GLPYUEPETRADATAPOOL_H
#define GLPAPP_GLPYUEPETRADATAPOOL_H

//#include "Epetra_config.h"

#include <iostream>

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_LAPACK.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "GenSQP_DataPool.hpp"
#include "GenSQP_YUEpetraVector.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

namespace GLpApp {

class GLpYUEpetraDataPool : public GenSQP::DataPool
{
public:

  GLpYUEpetraDataPool(
    Teuchos::RCP<const Epetra_Comm>    const& commptr
    ,const double                              beta
    ,const double                              len_x     // Ignored if myfile is *not* empty
    ,const double                              len_y     // Ignored if myfile is *not* empty
    ,const int                                 local_nx  // Ignored if myfile is *not* empty
    ,const int                                 local_ny  // Ignored if myfile is *not* empty
    ,const char                                myfile[]
    ,const bool                                trace
    );

  /** \brief Calls functions to compute nonlinear quantities and the augmented system matrix.
             These computations are performed after every update of the SQP iterate.
  */
  void computeAll( const GenSQP::Vector &x );

  /** \brief Solves augmented system.*/
  int  solveAugsys( const Teuchos::RCP<const Epetra_MultiVector> & rhsy,
                    const Teuchos::RCP<const Epetra_MultiVector> & rhsu,
                    const Teuchos::RCP<const Epetra_MultiVector> & rhsp,
                    const Teuchos::RCP<Epetra_MultiVector> & y,
                    const Teuchos::RCP<Epetra_MultiVector> & u,
                    const Teuchos::RCP<Epetra_MultiVector> & p,
                    double tol );

  Teuchos::RCP<const Epetra_Comm> getCommPtr();

  Teuchos::RCP<Epetra_FECrsMatrix> getA();
  Teuchos::RCP<Epetra_FECrsMatrix> getB();
  Teuchos::RCP<Epetra_FECrsMatrix> getH();
  Teuchos::RCP<Epetra_FECrsMatrix> getR();
  Teuchos::RCP<Epetra_CrsMatrix> getAugmat();
  Teuchos::RCP<Epetra_FECrsMatrix> getNpy();

  Teuchos::RCP<Epetra_FEVector> getb();
  Teuchos::RCP<Epetra_FEVector> getq();
  Teuchos::RCP<Epetra_FEVector> getNy();

  /** \brief Calls the function that computes the nonlinear term.*/
  void computeNy(const Teuchos::RCP<const Epetra_MultiVector> & y);

  /** \brief Calls the function that computes the Jacobian of the nonlinear term.*/
  void computeNpy(const Teuchos::RCP<const Epetra_MultiVector> & y);

  /** \brief Assembles the augmented system (KKT-type) matrix.*/
  void computeAugmat();

  Teuchos::RCP<const Epetra_SerialDenseMatrix> getipcoords();
  Teuchos::RCP<const Epetra_IntSerialDenseVector> getipindx();
  Teuchos::RCP<const Epetra_SerialDenseMatrix> getpcoords();
  Teuchos::RCP<const Epetra_IntSerialDenseVector> getpindx();
  Teuchos::RCP<const Epetra_IntSerialDenseMatrix> gett();
  Teuchos::RCP<const Epetra_IntSerialDenseMatrix> gete();

  double getbeta();

  /** \brief Outputs the solution vector to files.*/
  void PrintVec( const Teuchos::RCP<const Epetra_Vector> & x );

private:

  Teuchos::RCP<const Epetra_Comm> commptr_;

  /** \brief Coordinates of nodes that are unique to this subdomain.*/
  Teuchos::RCP<Epetra_SerialDenseMatrix> ipcoords_;
  /** \brief Global nodes (interior, nonoverlapping) in this subdomain.*/
  Teuchos::RCP<Epetra_IntSerialDenseVector> ipindx_;
  /** \brief Coordinates of all nodes in this subdomain.*/
  Teuchos::RCP<Epetra_SerialDenseMatrix> pcoords_;
  /** \brief Global nodes (interior + shared, overlapping) in this subdomain.*/
  Teuchos::RCP<Epetra_IntSerialDenseVector> pindx_;
  /** \brief Elements (this includes all overlapping nodes).*/
  Teuchos::RCP<Epetra_IntSerialDenseMatrix> t_;
  /** \brief Edges.*/
  Teuchos::RCP<Epetra_IntSerialDenseMatrix> e_;

  /** \brief Volume stiffness matrix.*/
  Teuchos::RCP<Epetra_FECrsMatrix> A_;
  /** \brief Control/state mass matrix.*/
  Teuchos::RCP<Epetra_FECrsMatrix> B_;
  /** \brief Volume mass matrix.*/
  Teuchos::RCP<Epetra_FECrsMatrix> H_;
  /** \brief Edge mass matrix.*/
  Teuchos::RCP<Epetra_FECrsMatrix> R_;

  /** \brief Basis matrix for p_bar=B*p.*/
  Teuchos::RCP<Epetra_MultiVector> B_bar_;

  /** \brief Augmented system matrix:
   [ I  Jac* ]
   [Jac  0   ]
  */
  Teuchos::RCP<Epetra_CrsMatrix> Augmat_;

  /** \brief Jacobian of the nonlinear term.*/
  Teuchos::RCP<Epetra_FECrsMatrix> Npy_;

  /** \brief Right-hand side of the PDE.*/
  Teuchos::RCP<Epetra_FEVector> b_;
  /** \brief The desired state.*/
  Teuchos::RCP<Epetra_FEVector> q_;

  Teuchos::RCP<Epetra_FEVector> Ny_;

  /** \brief Regularization parameter.*/
  double beta_;

};

class Usr_Par {
public:

  Usr_Par();

  Epetra_SerialDenseMatrix Nodes;
  Epetra_SerialDenseVector Weights;

  Epetra_SerialDenseMatrix N;

  Epetra_SerialDenseMatrix Nx1;

  Epetra_SerialDenseMatrix Nx2;

  Epetra_SerialDenseMatrix S1;
  Epetra_SerialDenseMatrix S2;
  Epetra_SerialDenseMatrix S3;

  Epetra_SerialDenseVector Nw;

  Epetra_SerialDenseMatrix NNw;

  Epetra_SerialDenseMatrix * NNNw;

  Epetra_SerialDenseMatrix * NdNdx1Nw;

  Epetra_SerialDenseMatrix * NdNdx2Nw;

  ~Usr_Par() {
    delete [] NNNw;
    delete [] NdNdx1Nw;
    delete [] NdNdx2Nw;
  }

  void Print(std::ostream& os) const;
};

} // namespace GLpApp

#endif // GLPAPP_GLPYUEPETRADATAPOOL_H
