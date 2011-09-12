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
    Teuchos::RefCountPtr<const Epetra_Comm>    const& commptr
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
  int  solveAugsys( const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsy,
                    const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsu,
                    const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsp,
                    const Teuchos::RefCountPtr<Epetra_MultiVector> & y,
                    const Teuchos::RefCountPtr<Epetra_MultiVector> & u,
                    const Teuchos::RefCountPtr<Epetra_MultiVector> & p,
                    double tol );

  Teuchos::RefCountPtr<const Epetra_Comm> getCommPtr();

  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getA();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getB();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getH();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getR();
  Teuchos::RefCountPtr<Epetra_CrsMatrix> getAugmat();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getNpy();

  Teuchos::RefCountPtr<Epetra_FEVector> getb();
  Teuchos::RefCountPtr<Epetra_FEVector> getq();
  Teuchos::RefCountPtr<Epetra_FEVector> getNy();

  /** \brief Calls the function that computes the nonlinear term.*/
  void computeNy(const Teuchos::RefCountPtr<const Epetra_MultiVector> & y);

  /** \brief Calls the function that computes the Jacobian of the nonlinear term.*/
  void computeNpy(const Teuchos::RefCountPtr<const Epetra_MultiVector> & y);

  /** \brief Assembles the augmented system (KKT-type) matrix.*/
  void computeAugmat();
  
  Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> getipcoords();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> getipindx();
  Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> getpcoords();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> getpindx();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> gett();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> gete();

  double getbeta();
  
  /** \brief Outputs the solution vector to files.*/
  void PrintVec( const Teuchos::RefCountPtr<const Epetra_Vector> & x );

private:

  Teuchos::RefCountPtr<const Epetra_Comm> commptr_;
          
  /** \brief Coordinates of nodes that are unique to this subdomain.*/
  Teuchos::RefCountPtr<Epetra_SerialDenseMatrix> ipcoords_;
  /** \brief Global nodes (interior, nonoverlapping) in this subdomain.*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseVector> ipindx_;
  /** \brief Coordinates of all nodes in this subdomain.*/
  Teuchos::RefCountPtr<Epetra_SerialDenseMatrix> pcoords_;
  /** \brief Global nodes (interior + shared, overlapping) in this subdomain.*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseVector> pindx_;
  /** \brief Elements (this includes all overlapping nodes).*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseMatrix> t_;
  /** \brief Edges.*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseMatrix> e_;

  /** \brief Volume stiffness matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> A_;
  /** \brief Control/state mass matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> B_;
  /** \brief Volume mass matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> H_;
  /** \brief Edge mass matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> R_;

  /** \brief Basis matrix for p_bar=B*p.*/
  Teuchos::RefCountPtr<Epetra_MultiVector> B_bar_;

  /** \brief Augmented system matrix:
   [ I  Jac* ]
   [Jac  0   ]
  */
  Teuchos::RefCountPtr<Epetra_CrsMatrix> Augmat_;

  /** \brief Jacobian of the nonlinear term.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> Npy_;

  /** \brief Right-hand side of the PDE.*/
  Teuchos::RefCountPtr<Epetra_FEVector> b_;
  /** \brief The desired state.*/
  Teuchos::RefCountPtr<Epetra_FEVector> q_;

  Teuchos::RefCountPtr<Epetra_FEVector> Ny_;

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
  
  void Print(ostream& os) const;
};

} // namespace GLpApp

#endif // GLPAPP_GLPYUEPETRADATAPOOL_H
