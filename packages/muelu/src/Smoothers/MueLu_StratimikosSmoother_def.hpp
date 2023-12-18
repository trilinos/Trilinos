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
#ifndef MUELU_STRATIMIKOSSMOOTHER_DEF_HPP
#define MUELU_STRATIMIKOSSMOOTHER_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_StratimikosSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#ifdef MUELU_RECURMG
#include "Stratimikos_MueLuHelpers.hpp"
#endif

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include "Teuchos_AbstractFactoryStd.hpp"
#include <Teuchos_ParameterList.hpp>
#include <unordered_map>

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::StratimikosSmoother(const std::string type, const Teuchos::ParameterList& paramList)
  : type_(type) {
  std::transform(type_.begin(), type_.end(), type_.begin(), ::toupper);
  ParameterList& pList = const_cast<ParameterList&>(paramList);

  if (pList.isParameter("smoother: recurMgOnFilteredA")) {
    recurMgOnFilteredA_ = true;
    pList.remove("smoother: recurMgOnFilteredA");
  }
  bool isSupported = type_ == "STRATIMIKOS";
  this->declareConstructionOutcome(!isSupported, "Stratimikos does not provide the smoother '" + type_ + "'.");
  if (isSupported)
    SetParameterList(paramList);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
  Factory::SetParameterList(paramList);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  this->Input(currentLevel, "A");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
  FactoryMonitor m(*this, "Setup Smoother", currentLevel);

  A_ = Factory::Get<RCP<Matrix> >(currentLevel, "A");
  SetupStratimikos(currentLevel);
  SmootherPrototype::IsSetup(true);
  this->GetOStream(Statistics1) << description() << std::endl;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::SetupStratimikos(Level& currentLevel) {
  RCP<const Thyra::LinearOpBase<Scalar> > thyraA;
  if (recurMgOnFilteredA_) {
    RCP<Matrix> filteredA;
    ExperimentalDropVertConnections(filteredA, currentLevel);
    thyraA = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(filteredA)->getCrsMatrix());
  } else
    thyraA = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A_)->getCrsMatrix());

  // Build Stratimikos solver
  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
  if (recurMgOnFilteredA_) {
#ifdef MUELU_RECURMG
    Stratimikos::enableMueLu<LocalOrdinal, GlobalOrdinal, Node>(linearSolverBuilder);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::StratimikosSmoother:: must compile with MUELU_RECURMG defined. Unfortunately, cmake does not always produce a proper link.txt file (which sometimes requires libmuelu.a before and after libmuelu-interface.a). After configuring, run script muelu/utils/misc/patchLinkForRecurMG to change link.txt files manually. If you want to create test example, add -DMUELU_RECURMG=ON to cmake arguments.");
#endif
  }

  linearSolverBuilder.setParameterList(rcpFromRef(const_cast<ParameterList&>(this->GetParameterList())));

  // Build a new "solver factory" according to the previously specified parameter list.
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
  solver_                                                         = Thyra::linearOpWithSolve(*solverFactory, thyraA);
#ifdef dumpOutRecurMGDebug
  char mystring[120];
  sprintf(mystring, "for i in A_[0123456789].m P_[0123456789].m;  do T=Xecho $i | sed Xs/.m$/%d.m/XX; mv $i $T; done", (int)currentLevel.GetLevelID());
  fflush(stdout);
  mystring[50] = '`';
  mystring[65] = '"';
  mystring[76] = '"';
  mystring[77] = '`';
  system(mystring);
#endif
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::ExperimentalDropVertConnections(RCP<Matrix>& filteredA, Level& currentLevel) {
  // strip out the veritcal connections.
  //
  // There is some code, which is currently turned off (via sumDropped). That
  // makes things more complicated. I want to keep it for now, and so it is
  // explained here. The basic idea is to try and maintain the character
  // of the horizontal stencil by summing dropped entries appropriately.
  // However, this does not correspond to plane relexation, and so I am
  // not sure it is really justified. Anyway, the picture below
  // gives a situation with a 15-pt stencil
  //
  //                              a
  //                             /
  //                            /
  //                     b ----c----- d
  //                          /
  //                         /
  //                        e
  //                              f         dropped a & l summed into f
  //                             /          dropped b & m summed into g
  //                            /           dropped c & n summed into i
  //                     g ----i----- j     dropped d & o summed into j
  //                          /             dropped e & p summed into k
  //                         /
  //                        k
  //                              l
  //                             /
  //                            /
  //                     m ----n----- o
  //                          /
  //                         /
  //                        p
  // To do this, we use umap to record locations within the middle layer associated
  // with each line ID (e.g. g corresponds to the 13th line). Then, when dropping (in a 2nd pass) we
  // use lineId and umap to find where dropped entries should be summed (e.g., b corresponds to the
  // 13th line and umap[13] points at location for g).
  // using TST = typename Teuchos::ScalarTraits<SC>;

  bool sumDropped = false;

  LO dofsPerNode = A_->GetFixedBlockSize();

  RCP<ParameterList> fillCompleteParams(new ParameterList);  // This code taken from Build method
  fillCompleteParams->set("No Nonlocal Changes", true);      // within MueLu_FilteredAFactory_def
  filteredA = MatrixFactory::Build(A_->getCrsGraph());
  filteredA->resumeFill();

  ArrayView<const LocalOrdinal> inds;
  ArrayView<const Scalar> valsA;
  ArrayView<Scalar> vals;
  Teuchos::ArrayRCP<LocalOrdinal> TVertLineId = Factory::Get<Teuchos::ArrayRCP<LocalOrdinal> >(currentLevel, "LineDetection_VertLineIds");
  Teuchos::ArrayRCP<LocalOrdinal> TLayerId    = Factory::Get<Teuchos::ArrayRCP<LocalOrdinal> >(currentLevel, "LineDetection_Layers");
  LocalOrdinal* VertLineId                    = TVertLineId.getRawPtr();
  LocalOrdinal* LayerId                       = TLayerId.getRawPtr();
  TEUCHOS_TEST_FOR_EXCEPTION((LayerId == NULL) || (VertLineId == NULL), Exceptions::RuntimeError, "MueLu::StratimikosSmoother:: no line information found on this level. Cannot use recurMgOnFilteredA on this level.");

  Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  for (size_t i = 0; i < A_->getRowMap()->getLocalNumElements(); i++) {
    A_->getLocalRowView(i, inds, valsA);
    size_t nnz = inds.size();
    ArrayView<const Scalar> vals1;
    filteredA->getLocalRowView(i, inds, vals1);
    vals = ArrayView<Scalar>(const_cast<Scalar*>(vals1.getRawPtr()), nnz);
    memcpy(vals.getRawPtr(), valsA.getRawPtr(), nnz * sizeof(Scalar));
    size_t inode, jdof, jnode, jdof_offset;
    inode = i / dofsPerNode;

    std::unordered_map<LocalOrdinal, LocalOrdinal> umap;  // umap[k] indicates where the dropped entry
                                                          // corresponding to kth line should be added
                                                          // within the row. See comments above.
    if (sumDropped) {
      for (size_t j = 0; j < nnz; j++) {
        jdof        = inds[j];
        jnode       = jdof / dofsPerNode;
        jdof_offset = jdof - jnode * dofsPerNode;
        if (LayerId[jnode] == LayerId[inode]) umap[dofsPerNode * VertLineId[jnode] + jdof_offset] = j;
      }
    }

    // drop non-middle layer entries. When sumDropped=true,  sum dropped entries to corresponding mid-layer entry
    for (size_t j = 0; j < nnz; j++) {
      jdof        = inds[j];
      jnode       = jdof / dofsPerNode;
      jdof_offset = jdof - jnode * dofsPerNode;
      if (LayerId[jnode] != LayerId[inode]) {
        if (sumDropped) {
          if (umap.find(dofsPerNode * VertLineId[jnode + jdof_offset]) != umap.end())
            vals[umap[dofsPerNode * VertLineId[jnode + jdof_offset]]] += vals[j];
        }
        vals[j] = ZERO;
      }
    }
  }
  filteredA->fillComplete(fillCompleteParams);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::StratimikosSmoother::Apply(): Setup() has not been called");

  // Apply
  if (InitialGuessIsZero) {
    X.putScalar(0.0);
    RCP<Thyra::MultiVectorBase<Scalar> > thyraX       = Teuchos::rcp_const_cast<Thyra::MultiVectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyraMultiVector(rcpFromRef(X)));
    RCP<const Thyra::MultiVectorBase<Scalar> > thyraB = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyraMultiVector(rcpFromRef(B));
    Thyra::SolveStatus<Scalar> status                 = Thyra::solve<Scalar>(*solver_, Thyra::NOTRANS, *thyraB, thyraX.ptr());
    RCP<MultiVector> thyXpX                           = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toXpetra(thyraX, X.getMap()->getComm());
    X                                                 = *thyXpX;

  } else {
    typedef Teuchos::ScalarTraits<Scalar> TST;
    RCP<MultiVector> Residual = Utilities::Residual(*A_, X, B);

    RCP<MultiVector> Cor                                = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(X.getMap(), X.getNumVectors(), true);
    RCP<Thyra::MultiVectorBase<Scalar> > thyraCor       = Teuchos::rcp_const_cast<Thyra::MultiVectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyraMultiVector(Cor));
    RCP<const Thyra::MultiVectorBase<Scalar> > thyraRes = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyraMultiVector(Residual);
    Thyra::SolveStatus<Scalar> status                   = Thyra::solve<Scalar>(*solver_, Thyra::NOTRANS, *thyraRes, thyraCor.ptr());
    RCP<MultiVector> thyXpCor                           = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toXpetra(thyraCor, X.getMap()->getComm());
    X.update(TST::one(), *thyXpCor, TST::one());
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<double, LocalOrdinal, GlobalOrdinal, Node> > StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  RCP<StratimikosSmoother> smoother = rcp(new StratimikosSmoother(*this));
  smoother->SetParameterList(this->GetParameterList());
  return smoother;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::string StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  if (SmootherPrototype::IsSetup()) {
    out << solver_->description();
  } else {
    out << "STRATIMIKOS {type = " << type_ << "}";
  }
  return out.str();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters1) {
    out0 << "Parameter list: " << std::endl;
    Teuchos::OSTab tab2(out);
    out << this->GetParameterList();
  }

  if (verbLevel & External)
    if (solver_ != Teuchos::null) {
      Teuchos::OSTab tab2(out);
      out << *solver_ << std::endl;
    }

  if (verbLevel & Debug) {
    out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
         << "-" << std::endl
         << "RCP<solver_>: " << solver_ << std::endl;
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t StratimikosSmoother<double, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

#endif  // HAVE_MUELU_STRATIMIKOS
#endif  // MUELU_STRATIMIKOSSMOOTHER_DEF_HPP
