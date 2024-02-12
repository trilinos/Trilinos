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
#if defined(HAVE_MUELU_AMESOS2)
#include <Xpetra_Matrix.hpp>
#include <Xpetra_IO.hpp>

#include <Amesos2_config.h>
#include <Amesos2.hpp>

#include "MueLu_Amesos2Smoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Projection<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Projection(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Nullspace) {
  localMap_ = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Nullspace->getMap()->lib(),
                                                                           Nullspace->getNumVectors(),
                                                                           Nullspace->getMap()->getIndexBase(),
                                                                           Nullspace->getMap()->getComm(),
                                                                           Xpetra::LocallyReplicated);

  Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tempMV = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(localMap_, Nullspace->getNumVectors());
  const Scalar ONE                                                                     = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar ZERO                                                                    = Teuchos::ScalarTraits<Scalar>::zero();
  tempMV->multiply(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, ONE, *Nullspace, *Nullspace, ZERO);

  Kokkos::View<Scalar**, Kokkos::LayoutLeft, Kokkos::HostSpace> Q("Q", Nullspace->getNumVectors(), Nullspace->getNumVectors());
  int LDQ;
  {
    auto dots = tempMV->getHostLocalView(Xpetra::Access::ReadOnly);
    Kokkos::deep_copy(Q, dots);
    int strides[2];
    Q.stride(strides);
    LDQ = strides[1];
  }

  Teuchos::LAPACK<LocalOrdinal, Scalar> lapack;
  int info = 0;
  lapack.POTRF('L', Nullspace->getNumVectors(), Q.data(), LDQ, &info);
  TEUCHOS_ASSERT(info == 0);
  lapack.TRTRI('L', 'N', Nullspace->getNumVectors(), Q.data(), LDQ, &info);
  TEUCHOS_ASSERT(info == 0);

  Nullspace_ = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Nullspace->getMap(), Nullspace->getNumVectors());

  for (size_t i = 0; i < Nullspace->getNumVectors(); i++) {
    for (size_t j = 0; j <= i; j++) {
      Nullspace_->getVectorNonConst(i)->update(Q(i, j), *Nullspace->getVector(j), ONE);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Projection<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    projectOut(Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X) {
  const Scalar ONE  = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();

  // Project X onto orthonormal nullspace
  // Nullspace_ ^T * X
  Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tempMV = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(localMap_, X.getNumVectors());
  tempMV->multiply(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, ONE, *Nullspace_, X, ZERO);
  auto dots      = tempMV->getHostLocalView(Xpetra::Access::ReadOnly);
  bool doProject = true;
  for (size_t i = 0; i < X.getNumVectors(); i++) {
    for (size_t j = 0; j < Nullspace_->getNumVectors(); j++) {
      doProject = doProject || (Teuchos::ScalarTraits<Scalar>::magnitude(dots(j, i)) > 100 * Teuchos::ScalarTraits<Scalar>::eps());
    }
  }
  if (doProject) {
    for (size_t i = 0; i < X.getNumVectors(); i++) {
      for (size_t j = 0; j < Nullspace_->getNumVectors(); j++) {
        X.getVectorNonConst(i)->update(-dots(j, i), *Nullspace_->getVector(j), ONE);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Amesos2Smoother(const std::string& type, const Teuchos::ParameterList& paramList)
  : type_(type)
  , useTransformation_(false) {
  this->SetParameterList(paramList);

  if (!type_.empty()) {
    // Transform string to "Abcde" notation
    std::transform(type_.begin(), type_.end(), type_.begin(), ::tolower);
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
    this->declareConstructionOutcome(true, std::string("Amesos2 has been compiled without SuperLU_DIST, SuperLU, Klu, or Basker. By default, MueLu tries") +
                                               "to use one of these libraries. Amesos2 must be compiled with one of these solvers, " +
                                               "or a valid Amesos2 solver has to be specified explicitly.");
    return;
#endif
    if (oldtype != "")
      this->GetOStream(Warnings0) << "MueLu::Amesos2Smoother: \"" << oldtype << "\" is not available. Using \"" << type_ << "\" instead" << std::endl;
    else
      this->GetOStream(Runtime1) << "MueLu::Amesos2Smoother: using \"" << type_ << "\"" << std::endl;
  }

  // Check the validity of the solver type parameter
  this->declareConstructionOutcome(Amesos2::query(type_) == false, "The Amesos2 library reported that the solver '" + type_ + "' is not available. " +
                                                                       "Amesos2 has been compiled without the support of this solver, or the solver name is misspelled.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Amesos2Smoother() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<RCP<const FactoryBase> >("A", null, "Factory of the coarse matrix");
  validParamList->set<RCP<const FactoryBase> >("Nullspace", null, "Factory of the nullspace");
  validParamList->set<bool>("fix nullspace", false, "Remove zero eigenvalue by adding rank one correction.");
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("Amesos2", norecurse, "Parameters that are passed to Amesos2");
  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  ParameterList pL = this->GetParameterList();

  this->Input(currentLevel, "A");
  if (pL.get<bool>("fix nullspace"))
    this->Input(currentLevel, "Nullspace");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
  FactoryMonitor m(*this, "Setup Smoother", currentLevel);

  if (SmootherPrototype::IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::Amesos2Smoother::Setup(): Setup() has already been called" << std::endl;

  RCP<Matrix> A = Factory::Get<RCP<Matrix> >(currentLevel, "A");

  // Do a quick check if we need to modify the matrix
  RCP<const Map> rowMap = A->getRowMap();
  RCP<Matrix> factorA;
  Teuchos::ParameterList pL = this->GetParameterList();

  if (pL.get<bool>("fix nullspace")) {
    this->GetOStream(Runtime1) << "MueLu::Amesos2Smoother::Setup(): fixing nullspace" << std::endl;

    rowMap            = A->getRowMap();
    size_t gblNumCols = rowMap->getGlobalNumElements();

    RCP<MultiVector> NullspaceOrig = Factory::Get<RCP<MultiVector> >(currentLevel, "Nullspace");

    projection_                = rcp(new Projection<Scalar, LocalOrdinal, GlobalOrdinal, Node>(NullspaceOrig));
    RCP<MultiVector> Nullspace = projection_->Nullspace_;

    RCP<MultiVector> ghostedNullspace;
    RCP<const Map> colMap;
    RCP<const Import> importer;
    if (rowMap->getComm()->getSize() > 1) {
      this->GetOStream(Warnings0) << "MueLu::Amesos2Smoother::Setup(): Applying nullspace fix on distributed matrix. Try rebalancing to single rank!" << std::endl;
      ArrayRCP<GO> elements_RCP;
      elements_RCP.resize(gblNumCols);
      ArrayView<GO> elements = elements_RCP();
      for (size_t k = 0; k < gblNumCols; k++)
        elements[k] = Teuchos::as<GO>(k);
      colMap           = MapFactory::Build(rowMap->lib(), gblNumCols * rowMap->getComm()->getSize(), elements, Teuchos::ScalarTraits<GO>::zero(), rowMap->getComm());
      importer         = ImportFactory::Build(rowMap, colMap);
      ghostedNullspace = MultiVectorFactory::Build(colMap, Nullspace->getNumVectors());
      ghostedNullspace->doImport(*Nullspace, *importer, Xpetra::INSERT);
    } else {
      ghostedNullspace = Nullspace;
      colMap           = rowMap;
    }

    using ATS         = Kokkos::ArithTraits<SC>;
    using impl_Scalar = typename ATS::val_type;
    using impl_ATS    = Kokkos::ArithTraits<impl_Scalar>;
    using range_type  = Kokkos::RangePolicy<LO, typename NO::execution_space>;

    typedef typename Matrix::local_matrix_type KCRS;
    typedef typename KCRS::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename KCRS::values_type::non_const_type scalar_view_t;

    const impl_Scalar impl_SC_ZERO = impl_ATS::zero();

    size_t lclNumRows       = rowMap->getLocalNumElements();
    LocalOrdinal lclNumCols = Teuchos::as<LocalOrdinal>(gblNumCols);
    lno_view_t newRowPointers("newRowPointers", lclNumRows + 1);
    lno_nnz_view_t newColIndices("newColIndices", lclNumRows * gblNumCols);
    scalar_view_t newValues("newValues", lclNumRows * gblNumCols);

    impl_Scalar shift;
    {
      RCP<Vector> diag = VectorFactory::Build(A->getRowMap());
      A->getLocalDiagCopy(*diag);
      shift = diag->normInf();
    }

    // form normalization * nullspace * nullspace^T
    {
      auto lclNullspace        = Nullspace->getDeviceLocalView(Xpetra::Access::ReadOnly);
      auto lclGhostedNullspace = ghostedNullspace->getDeviceLocalView(Xpetra::Access::ReadOnly);
      Kokkos::parallel_for(
          "MueLu:Amesos2Smoother::fixNullspace_1", range_type(0, lclNumRows + 1),
          KOKKOS_LAMBDA(const size_t i) {
            if (i < lclNumRows) {
              newRowPointers(i) = i * gblNumCols;
              for (LocalOrdinal j = 0; j < lclNumCols; j++) {
                newColIndices(i * gblNumCols + j) = j;
                newValues(i * gblNumCols + j)     = impl_SC_ZERO;
                for (size_t I = 0; I < lclNullspace.extent(1); I++)
                  for (size_t J = 0; J < lclGhostedNullspace.extent(1); J++)
                    newValues(i * gblNumCols + j) += shift * lclNullspace(i, I) * impl_ATS::conjugate(lclGhostedNullspace(j, J));
              }
            } else
              newRowPointers(lclNumRows) = lclNumRows * gblNumCols;
          });
    }

    // add A
    if (colMap->lib() == Xpetra::UseTpetra) {
      auto lclA          = A->getLocalMatrixDevice();
      auto lclColMapA    = A->getColMap()->getLocalMap();
      auto lclColMapANew = colMap->getLocalMap();
      Kokkos::parallel_for(
          "MueLu:Amesos2Smoother::fixNullspace_2", range_type(0, lclNumRows),
          KOKKOS_LAMBDA(const size_t i) {
            for (size_t jj = lclA.graph.row_map(i); jj < lclA.graph.row_map(i + 1); jj++) {
              LO j          = lclColMapANew.getLocalElement(lclColMapA.getGlobalElement(lclA.graph.entries(jj)));
              impl_Scalar v = lclA.values(jj);
              newValues(i * gblNumCols + j) += v;
            }
          });
    } else {
      auto lclA = A->getLocalMatrixHost();
      for (size_t i = 0; i < lclNumRows; i++) {
        for (size_t jj = lclA.graph.row_map(i); jj < lclA.graph.row_map(i + 1); jj++) {
          LO j = colMap->getLocalElement(A->getColMap()->getGlobalElement(lclA.graph.entries(jj)));
          SC v = lclA.values(jj);
          newValues(i * gblNumCols + j) += v;
        }
      }
    }

    RCP<Matrix> newA       = rcp(new CrsMatrixWrap(rowMap, colMap, 0));
    RCP<CrsMatrix> newAcrs = rcp_dynamic_cast<CrsMatrixWrap>(newA)->getCrsMatrix();
    newAcrs->setAllValues(newRowPointers, newColIndices, newValues);
    newAcrs->expertStaticFillComplete(A->getDomainMap(), A->getRangeMap(),
                                      importer, A->getCrsGraph()->getExporter());

    factorA = newA;
    rowMap  = factorA->getRowMap();
  } else {
    factorA = A;
  }

  RCP<Tpetra_CrsMatrix> tA = Utilities::Op2NonConstTpetraCrs(factorA);

  prec_ = Amesos2::create<Tpetra_CrsMatrix, Tpetra_MultiVector>(type_, tA);
  TEUCHOS_TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "Amesos2::create returns Teuchos::null");
  RCP<Teuchos::ParameterList> amesos2_params = Teuchos::rcpFromRef(pL.sublist("Amesos2"));
  amesos2_params->setName("Amesos2");
  if ((rowMap->getGlobalNumElements() != as<size_t>((rowMap->getMaxAllGlobalIndex() - rowMap->getMinAllGlobalIndex()) + 1)) ||
      (!rowMap->isContiguous() && (rowMap->getComm()->getSize() == 1))) {
    if (!(amesos2_params->sublist(prec_->name()).template isType<bool>("IsContiguous")))
      amesos2_params->sublist(prec_->name()).set("IsContiguous", false, "Are GIDs Contiguous");
  }
  prec_->setParameters(amesos2_params);

  prec_->numericFactorization();

  SmootherPrototype::IsSetup(true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
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
    ArrayRCP<const SC> Xdata = X.getData(0), Bdata = B.getData(0);
    ArrayRCP<SC> X_data = X_->getDataNonConst(0), B_data = B_->getDataNonConst(0);

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
    size_t length = X.getLocalLength();

    ArrayRCP<SC> Xdata        = X.getDataNonConst(0);
    ArrayRCP<const SC> X_data = X_->getData(0);

    for (size_t i = 0; i < length; i++)
      Xdata[i] = X_data[i];
  }

  {
    Teuchos::ParameterList pL = this->GetParameterList();
    if (pL.get<bool>("fix nullspace")) {
      projection_->projectOut(X);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Copy() const {
  return rcp(new Amesos2Smoother(*this));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
  if (!prec_.is_null())
    return prec_->getStatus().getNnzLU();
  else
    return 0.0;
}
}  // namespace MueLu

#endif  // HAVE_MUELU_AMESOS2
#endif  // MUELU_AMESOS2SMOOTHER_DEF_HPP
