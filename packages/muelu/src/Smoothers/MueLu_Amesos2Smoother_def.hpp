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
#ifndef MUELU_AMESOS2SMOOTHER_DEF_HPP
#define MUELU_AMESOS2SMOOTHER_DEF_HPP

#include <algorithm>

#include "MueLu_ConfigDefs.hpp"
#if defined (HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
#include <Xpetra_Matrix.hpp>

#include <Amesos2_config.h>
#include <Amesos2.hpp>

#include "MueLu_Amesos2Smoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Amesos2Smoother(const std::string& type, const Teuchos::ParameterList& paramList)
    : type_(type), useTransformation_(false) {
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
#elif defined(HAVE_AMESOS2_BASKER)
      type_ = "Basker";
#else
      throw Exceptions::RuntimeError("Amesos2 has been compiled without SuperLU_DIST, SuperLU, Klu, or Basker. By default, MueLu tries"
                                     "to use one of these libraries. Amesos2 must be compiled with one of these solvers, "
                                     "or a valid Amesos2 solver has to be specified explicitly.");
#endif
      if (oldtype != "")
        this->GetOStream(Warnings0) << "MueLu::Amesos2Smoother: \"" << oldtype << "\" is not available. Using \"" << type_ << "\" instead" << std::endl;
      else
        this->GetOStream(Runtime1) << "MueLu::Amesos2Smoother: using \"" << type_ << "\"" << std::endl;
    }

    // Check the validity of the solver type parameter
    TEUCHOS_TEST_FOR_EXCEPTION(Amesos2::query(type_) == false, Exceptions::RuntimeError, "The Amesos2 library reported that the solver '" << type_ << "' is not available. "
                               "Amesos2 has been compiled without the support of this solver, or the solver name is misspelled.");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Amesos2Smoother() { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    this->Input(currentLevel, "A");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::Amesos2Smoother::Setup(): Setup() has already been called" << std::endl;

    RCP<Matrix> A = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    // Do a quick check if we need to modify the matrix
    RCP<const Map> rowMap = A->getRowMap();
    if (rowMap->getGlobalNumElements() != as<size_t>((rowMap->getMaxAllGlobalIndex() - rowMap->getMinAllGlobalIndex())+1)) {
      // If our system is non-conventional, here is the place where we try to fix it
      // One example is: if our maps contain a gap in them, for instance GIDs
      // are [0,..., 100, 10000, ..., 10010], then Amesos2 breaks down.
      //
      // The approach we take is to construct a second system with maps
      // replaced by their continuous versions.
      //
      // FIXME: for the moment, this functionality works for selected limited scenarios
      this->GetOStream(Runtime1) << "MueLu::Amesos2Smoother::Setup(): using system transformation" << std::endl;

      TEUCHOS_TEST_FOR_EXCEPTION(rowMap->getComm()->getSize() > 1, Exceptions::RuntimeError,
        "MueLu::Amesos2Smoother::Setup Fixing coarse matrix for Amesos2 for multiple processors has not been implemented yet.");
      TEUCHOS_TEST_FOR_EXCEPTION(!rowMap->isSameAs(*A->getColMap()), Exceptions::RuntimeError,
        "MueLu::Amesos2Smoother::Setup Fixing coarse matrix for Amesos2 when row map is different from column map has not been implemented yet.");

      RCP<CrsMatrixWrap> Acrs = rcp_dynamic_cast<CrsMatrixWrap>(A);
      TEUCHOS_TEST_FOR_EXCEPTION(Acrs.is_null(), Exceptions::RuntimeError,
        "MueLu::Amesos2Smoother::Setup Fixing coarse matrix for Amesos2 when matrix is not a Crs matrix has not been implemented yet.");

      useTransformation_ = true;

      ArrayRCP<const size_t> rowPointers;
      ArrayRCP<const LO>     colIndices;
      ArrayRCP<const SC>     values;
      Acrs->getCrsMatrix()->getAllValues(rowPointers, colIndices, values);

      // Create new map
      RCP<Map>       map     = MapFactory::Build(rowMap->lib(), rowMap->getGlobalNumElements(), 0, rowMap->getComm());
      RCP<Matrix>    newA    = rcp(new CrsMatrixWrap(map, map, 0, Xpetra::StaticProfile));
      RCP<CrsMatrix> newAcrs = rcp_dynamic_cast<CrsMatrixWrap>(newA)->getCrsMatrix();

      using Teuchos::arcp_const_cast;
      newAcrs->setAllValues(arcp_const_cast<size_t>(rowPointers), arcp_const_cast<LO>(colIndices), arcp_const_cast<SC>(values));
      newAcrs->expertStaticFillComplete(map, map);

      A.swap(newA);

      X_ = MultiVectorFactory::Build(map, 1);
      B_ = MultiVectorFactory::Build(map, 1);
    }

    RCP<Tpetra_CrsMatrix> tA = Utilities::Op2NonConstTpetraCrs(A);

    prec_ = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>(type_, tA);
    TEUCHOS_TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "Amesos2::create returns Teuchos::null");

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool /* InitialGuessIsZero */) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::Amesos2Smoother::Apply(): Setup() has not been called");

    RCP<Tpetra_MultiVector> tX, tB;
    if (!useTransformation_) {
      tX = Utilities::MV2NonConstTpetraMV2(X);
      tB = Utilities::MV2NonConstTpetraMV2(const_cast<MultiVector&>(B));
    } else {
      // Copy data of the original vectors into the transformed ones
      size_t numVectors = X.getNumVectors();
      size_t length     = X.getLocalLength();

      TEUCHOS_TEST_FOR_EXCEPTION(numVectors > 1, Exceptions::RuntimeError,
        "MueLu::Amesos2Smoother::Apply: Fixing coarse matrix for Amesos2 for multivectors has not been implemented yet.");
      ArrayRCP<const SC> Xdata  = X.  getData(0),         Bdata  = B.  getData(0);
      ArrayRCP<SC>       X_data = X_->getDataNonConst(0), B_data = B_->getDataNonConst(0);

      for (size_t i = 0; i < length; i++) {
        X_data[i] = Xdata[i];
        B_data[i] = Bdata[i];
      }

      tX = Utilities::MV2NonConstTpetraMV2(*X_);
      tB = Utilities::MV2NonConstTpetraMV2(*B_);
    }

    prec_->setX(tX);
    prec_->setB(tB);

    prec_->solve();

    prec_->setX(Teuchos::null);
    prec_->setB(Teuchos::null);

    if (useTransformation_) {
      // Copy data from the transformed vectors into the original ones
      size_t length     = X.getLocalLength();

      ArrayRCP<SC>       Xdata  = X.  getDataNonConst(0);
      ArrayRCP<const SC> X_data = X_->getData(0);

      for (size_t i = 0; i < length; i++)
        Xdata[i] = X_data[i];
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Copy() const
  {
    return rcp (new Amesos2Smoother (*this));
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;

    if (SmootherPrototype::IsSetup() == true) {
      out << prec_->description();

    } else {
      out << SmootherPrototype::description();
      out << "{type = " << type_ << "}";
    }
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
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

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
    if(!prec_.is_null())
      return prec_->getStatus().getNnzLU();
    else
      return 0.0;

  }
} // namespace MueLu

#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_AMESOS2
#endif // MUELU_AMESOS2SMOOTHER_DEF_HPP
