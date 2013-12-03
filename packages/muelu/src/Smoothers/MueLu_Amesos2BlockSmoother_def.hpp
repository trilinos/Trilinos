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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_AMESOS2BLOCKSMOOTHER_DEF_HPP
#define MUELU_AMESOS2BLOCKSMOOTHER_DEF_HPP

#include <algorithm>

#include "MueLu_ConfigDefs.hpp"
#if defined (HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
#include <Xpetra_Matrix.hpp>

#include <Amesos2_config.h>
#include <Amesos2.hpp>

#include "MueLu_Amesos2BlockSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Amesos2BlockSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Amesos2BlockSmoother(const std::string& type, const Teuchos::ParameterList& paramList)
    : type_(type) {
    this->SetParameterList(paramList);

    if (!type_.empty()) {
      // Transform string to "Abcde" notation
      std::transform(type_.begin(),   type_.end(),   type_.begin(), ::tolower);
      std::transform(type_.begin(), ++type_.begin(), type_.begin(), ::toupper);
    }
    if (type_ == "Superlu_dist")
      type_ = "Superludist";

    // Try to come up with something availble
    // Order corresponds to our preference
    // TODO: It would be great is Amesos2 provides directly this kind of logic for us
    if (type_ == "" || Amesos2::query(type_) == false) {
      std::string oldtype = type_;
#if defined(HAVE_AMESOS2_SUPERLU)
      type_ = "Superlu";
#elif defined(HAVE_AMESOS2_KLU2)
      type_ = "Klu";
#elif defined(HAVE_AMESOS2_SUPERLUDIST)
      type_ = "Superludist";
#else
      throw Exceptions::RuntimeError("Amesos2 have been compiled without SuperLU_DIST, SuperLU or Klu. By default, MueLu tries"
                                     "to use one of these libraries. Amesos2 must be compiled with one of these solvers or"
                                     "a valid Amesos2 solver have to be specified explicitly.");
#endif
      if (oldtype != "")
        this->GetOStream(Warnings0, 0) << "Warning: MueLu::Amesos2BlockSmoother: \"" << oldtype << "\" is not available. Using \"" << type_ << "\" instead" << std::endl;
      else
        this->GetOStream(Warnings0, 0) << "MueLu::Amesos2BlockSmoother: using \"" << type_ << "\"" << std::endl;
    }

    // Check the validity of the solver type parameter
    TEUCHOS_TEST_FOR_EXCEPTION(Amesos2::query(type_) == false, Exceptions::RuntimeError, "The Amesos2 library reported that the solver '" << type_ << "' is not available. "
                               "Amesos2 have been compiled without the support of this solver or the solver name is misspelled.");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Amesos2BlockSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~Amesos2BlockSmoother() { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Amesos2BlockSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level& currentLevel) const {
    this->Input(currentLevel, "A");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Amesos2BlockSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level& currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
      this->GetOStream(Warnings0, 0) << "Warning: MueLu::Amesos2BlockSmoother::Setup(): Setup() has already been called" << std::endl;

    RCP<Matrix> A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    RCP<Tpetra_CrsMatrix> tA = Utils::Op2NonConstTpetraCrs(A_);

    // make local matrix

    // local communicator
    Teuchos::RCP<const Teuchos::Comm<int> > localComm;
#ifdef HAVE_MPI
    localComm = Teuchos::rcp(new Teuchos::MpiComm<int> (MPI_COMM_SELF));
#else
    localComm = Teuchos::rcp(new Teuchos::SerialComm<int> ());
#endif

    // get row map and setup local matrix
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > globalRowMap = tA->getRowMap();
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > globalColMap = tA->getColMap();
    LocalOrdinal numRows = globalRowMap -> getNodeNumElements();
    localRowMap_ = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(numRows, 0, localComm,
										 Tpetra::GloballyDistributed, tA->getNode()));
    Teuchos::RCP<Tpetra_CrsMatrix> localA = Teuchos::rcp(new Tpetra_CrsMatrix(localRowMap_,localRowMap_,100));
    
    // extract rows
    for(LocalOrdinal i = 0; i < numRows; i++) {
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> values;
      tA -> getLocalRowView(i,indices,values);
      std::vector<LocalOrdinal> indices_vec;
      std::vector<Scalar> values_vec;
      indices_vec.resize(0);
      values_vec.resize(0);
      for(unsigned int j=0; j < indices.size(); j++ ) {
	LocalOrdinal local_col = indices[j];
	GlobalOrdinal global_col = globalColMap->getGlobalElement(local_col);
	if(globalRowMap->isNodeGlobalElement(global_col)) {
	  local_col = globalRowMap->getLocalElement(global_col);
	  indices_vec.push_back(local_col);
	  values_vec.push_back(values[j]);
	}
      }
      localA->insertLocalValues(i,
				Teuchos::ArrayView<LocalOrdinal>(indices_vec),
				Teuchos::ArrayView<Scalar>(values_vec));
    }
    localA->fillComplete();
    
    prec_ = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>(type_, localA);
    prec_ -> preOrdering();
    prec_ -> symbolicFactorization();
    prec_ -> numericFactorization();
    TEUCHOS_TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "Amesos2::create returns Teuchos::null");

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Amesos2BlockSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::Amesos2BlockSmoother::Apply(): Setup() has not been called");

    RCP<Tpetra_MultiVector> tX = Utils::MV2NonConstTpetraMV2(X);
    MultiVector & BNonC = const_cast<MultiVector&>(B);
    RCP<Tpetra_MultiVector> tB = Utils::MV2NonConstTpetraMV2(BNonC);

    // construct local vectors    
    size_t numvecs = tX->getNumVectors();
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > globalRowMap = tX->getMap();
    Teuchos::ArrayView<const GlobalOrdinal> globalRows = globalRowMap -> getNodeElementList();
    LocalOrdinal numRows = globalRowMap->getNodeNumElements();
    Teuchos::RCP<Tpetra_MultiVector> localX = Teuchos::rcp( new Tpetra_MultiVector(localRowMap_,numvecs) );
    Teuchos::RCP<Tpetra_MultiVector> localB = Teuchos::rcp( new Tpetra_MultiVector(localRowMap_,numvecs) );
    // extract values
    for(size_t j=0; j<numvecs; j++) {
      Teuchos::ArrayRCP<const Scalar> vecj = tB -> getData(j);
      for(LocalOrdinal i = 0; i < numRows; i++) {
	localB->replaceLocalValue(i,j,vecj[i]);
      }
    }

    // solve
    prec_->setX(localX);
    prec_->setB(localB);
    prec_->solve();

    if(InitialGuessIsZero==true) {
      // extract to global vector
      for(size_t j=0; j<numvecs; j++) {
	Teuchos::ArrayRCP<const Scalar> localview = localX->getData(j);
	for(LocalOrdinal i=0; i<numRows; i++) {
	  tX->replaceLocalValue(i,j,localview[i]);
	}
      }
    }
    else {
      // create temporary vector
      Teuchos::RCP<Tpetra_MultiVector> Xtemp = Teuchos::rcp( new Tpetra_MultiVector(globalRowMap,numvecs) );
      for(size_t j=0; j<numvecs; j++) {
	Teuchos::ArrayRCP<const Scalar> localview = localX->getData(j);
	for(LocalOrdinal i=0; i<numRows; i++) {
	  Xtemp->replaceLocalValue(i,j,localview[i]);
	}
      }
      // update X
      tX->update((Scalar)1.0, *Xtemp, (Scalar)1.0);
    }

    prec_->setX(Teuchos::null);
    prec_->setB(Teuchos::null);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Amesos2BlockSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp(new Amesos2BlockSmoother(*this));
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string Amesos2BlockSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;

    if (SmootherPrototype::IsSetup() == true) {
      out << prec_->description();

    } else {
      out << SmootherPrototype::description();
      out << "{type = " << type_ << "}";
    }
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Amesos2BlockSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl;
      Teuchos::OSTab tab2(out);
      out << this->GetParameterList();
    }

    if ((verbLevel & External) && prec_ != Teuchos::null) {
      Teuchos::OSTab tab2(out);
      out << *prec_ << std::endl;
    }

    if (verbLevel & Debug)
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
  }

} // namespace MueLu

#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_AMESOS2
#endif // MUELU_AMESOS2BLOCKSMOOTHER_DEF_HPP
