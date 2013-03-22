// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_EpetraMultiVector.hpp>

#include "MueLu_EpetraOperator.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

int EpetraOperator::SetUseTranspose(bool UseTransposeBool) { return -1; }

int EpetraOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const { return -1; }

int EpetraOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  try {

    Epetra_MultiVector & temp_x = const_cast<Epetra_MultiVector &>(X);

    const Xpetra::EpetraMultiVector tX(rcpFromRef(temp_x));
    Xpetra::EpetraMultiVector       tY(rcpFromRef(Y));

    // check if X and Y points to the same memory
    if(X.Values() == Y.Values()) {
      // For AztecOO X and Y point to the same memory
      // use work vectors

      // reserve memory for deep copy vectors
      RCP<Xpetra::EpetraMultiVector> epX = Teuchos::rcp(new Xpetra::EpetraMultiVector(tX.getMap(), tX.getNumVectors())); // oops, we don't have a copy constructor?
      RCP<Xpetra::EpetraMultiVector> epY = Teuchos::rcp(new Xpetra::EpetraMultiVector(tY.getMap(), tY.getNumVectors()));

      // deep copy of RHS vector
      epX->update(1.0,tX,0.0);

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      epY->putScalar(0.0);

      // apply one V/W-cycle as preconditioner
      Hierarchy_->Iterate(*epX, 1, *epY, true);

      // deep copy solution from MueLu to AztecOO
      tY.update(1.0,*epY,0.0);
    } else {
      // X and Y point to different memory
      // avoid additonal working vectors, just pass through the vectors

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      tY.putScalar(0.0);
      Hierarchy_->Iterate(tX, 1, tY, true);
    }
  } catch(std::exception& e) {
    //TODO: error msg directly on std::cerr?
    std::cerr << "Caught an exception in MueLu::EpetraOperator::ApplyInverse():" << std::endl
        << e.what() << std::endl;
    return -1;
  }

  return 0;
}

double EpetraOperator::NormInf() const { return 0; }

const char * EpetraOperator::Label() const { return "MueLu::Hierarchy"; }

bool EpetraOperator::UseTranspose() const { return false; }

bool EpetraOperator::HasNormInf() const { return false; }

const Epetra_Comm & EpetraOperator::Comm() const {
  RCP<Matrix> A = Hierarchy_->GetLevel(0)->Get<RCP<Matrix> >("A");

  //TODO: This code is not pretty
  RCP<Xpetra::BlockedCrsMatrix<double, int, int> > epbA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<double, int, int> >(A);
  if(epbA != Teuchos::null) {
    RCP<const Xpetra::EpetraCrsMatrix> tmp_ECrsMtx = rcp_dynamic_cast<Xpetra::EpetraCrsMatrix >(epbA->getMatrix(0,0));
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
    RCP<Epetra_CrsMatrix> epA = tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
    return epA->Comm();
  }
  //

  RCP<Epetra_CrsMatrix>epA = Utils::Op2NonConstEpetraCrs(A);
  return epA->Comm();
}

const Epetra_Map & EpetraOperator::OperatorDomainMap() const {
  RCP<Matrix> A = Hierarchy_->GetLevel(0)->Get<RCP<Matrix> >("A");

  RCP<Xpetra::BlockedCrsMatrix<double, int, int> > epbA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<double, int, int> >(A);
  if(epbA != Teuchos::null)
    return Xpetra::toEpetra(epbA->getDomainMap());

  RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A);
  return epA->DomainMap();
}

const Epetra_Map & EpetraOperator::OperatorRangeMap() const {
  RCP<Matrix> A = Hierarchy_->GetLevel(0)->Get<RCP<Matrix> >("A");

  RCP<Xpetra::BlockedCrsMatrix<double, int, int> > epbA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<double, int, int> >(A);
  if(epbA != Teuchos::null)
    return Xpetra::toEpetra(epbA->getRangeMap());

  RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A);
  return epA->RangeMap();
}

RCP<MueLu::Hierarchy<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<double,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> > EpetraOperator::GetHierarchy() const { return Hierarchy_; }

} // namespace
