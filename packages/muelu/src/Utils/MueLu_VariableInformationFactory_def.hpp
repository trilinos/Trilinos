/*
 * MueLu_VariableInformationFactory_def.hpp
 *
 *  Created on: 20.02.2012
 *      Author: tobias
 */

#ifndef MUELU_VARIABLEINFORMATIONFACTORY_DEF_HPP_
#define MUELU_VARIABLEINFORMATIONFACTORY_DEF_HPP_

#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsOperator.hpp>

#include "MueLu_VariableInformationFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  VariableInformationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::VariableInformationFactory(const std::string& ename, const FactoryBase* fac, bool bCoarseLevelInfo)
    : varName_(ename), factory_(fac), coarseLevel_(bCoarseLevelInfo)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  VariableInformationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~VariableInformationFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void VariableInformationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    if(coarseLevel_)
      coarseLevel.DeclareInput(varName_,factory_,this);
    else
      fineLevel.DeclareInput(varName_,factory_,this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void VariableInformationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {
    typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OOperator; //TODO
    typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOOperator; //TODO
    typedef MueLu::Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> UUtils;
    typedef MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> UUUtils;

    // check if variable is available
    if(coarseLevel_) {
      if(coarseLevel.IsAvailable(varName_,factory_) == false) {
        coarseLevel.Release(varName_,factory_);
        return;
      }
    } else {
      if(fineLevel.IsAvailable(varName_,factory_) == false) {
        fineLevel.Release(varName_,factory_);
        return;
      }
    }

    // extract variable from Level
    RCP<OOperator> Ain = Teuchos::null;
    if(coarseLevel_) Ain = coarseLevel.Get< RCP<OOperator> >(varName_, factory_);
    else             Ain = fineLevel.Get< RCP<OOperator> >(varName_, factory_);

    // check if input variable is a Xpetra::Operator
    if(Ain==Teuchos::null) {
      GetOStream(Warnings0, 0) << varName_ << " (" << factory_ << ") is not a Operator" << std::endl;
      return;
    }

    GetOStream(Statistics0,0) << varName_ << "(" << factory_ << "): " << Ain->getRowMap()->getGlobalNumElements() << "x" << Ain->getDomainMap()->getGlobalNumElements() << " matrix, " << Ain->getGlobalNumEntries() << " entries, ||" << varName_ << "||_F=" << Ain->getFrobeniusNorm() << std::endl;

    if(Ain->getRangeMap()->getGlobalNumElements() == Ain->getDomainMap()->getGlobalNumElements()) {
      // calculate At = A-A^T
      RCP<OOperator> At = UUtils::Transpose(Ain,true);
      RCP<OOperator> AtmA = rcp(new CrsOOperator(Ain->getRowMap(),Ain->getGlobalMaxNumRowEntries()));
      UUtils::TwoMatrixAdd(At, false, 1.0, Ain, false, -1.0, AtmA);
      AtmA->fillComplete();

      GetOStream(Statistics0,0) << varName_ << "(" << factory_ << "): ||" << varName_ << "-" << varName_ << "^T||_F = " << AtmA->getFrobeniusNorm() << ", entries: " << AtmA->getGlobalNumEntries() << std::endl;

      Gershgorin(Ain,At); // TODO: do the same for the transposed
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void VariableInformationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Gershgorin(const RCP<Operator>& A, const RCP<Operator>& At) const {

    Scalar minAbsEigValue = 1e15; // rough approximation of min/max absolute EW on current proc
    Scalar maxAbsEigValue = 0.0;

    // loop over local rows
    for(size_t row=0; row<A->getNodeNumRows(); row++) {
      size_t nnz = A->getNumEntriesInLocalRow(row);

      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      A->getLocalRowView(row, indices, vals);

      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::VariableInformationFactory::Gershgorin: number of nonzeros not equal to number of indices? Error.");

      Scalar gershgorin_radius = 0.0;
      Scalar diagonal_entry = 0.0;

      for(size_t i=0; i<(size_t)indices.size(); i++) {
        if(indices[i]==(LocalOrdinal)row) { // entry on matrix diagonal (row,row)
          diagonal_entry = vals[i];
        } else {  // off diagonal entry
          gershgorin_radius += abs(vals[i]);
        }
      }

      Scalar a = diagonal_entry - gershgorin_radius;
      Scalar b = diagonal_entry + gershgorin_radius;
      if(a*b < 0.0) { // a and b have different sign
        minAbsEigValue = 0.0; // zero is possible EW
        if(std::max(abs(a),abs(b)) > maxAbsEigValue) maxAbsEigValue = std::max(abs(a),abs(b));
      } else { // a and b have same sign
        if(std::min(abs(a),abs(b)) < minAbsEigValue) minAbsEigValue = std::min(abs(a),abs(b));
        if(std::max(abs(a),abs(b)) > maxAbsEigValue) maxAbsEigValue = std::max(abs(a),abs(b));
      }
    }

    if(At != Teuchos::null) {
      // loop over local rows
      for(size_t row=0; row<A->getNodeNumRows(); row++) {
        size_t nnz = At->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        At->getLocalRowView(row, indices, vals);

        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::VariableInformationFactory::Gershgorin: number of nonzeros not equal to number of indices? Error.");

        Scalar gershgorin_radius = 0.0;
        Scalar diagonal_entry = 0.0;

        for(size_t i=0; i<(size_t)indices.size(); i++) {
          if(indices[i]==(LocalOrdinal)row) { // entry on matrix diagonal (row,row)
            diagonal_entry = vals[i];
          } else {  // off diagonal entry
            gershgorin_radius += abs(vals[i]);
          }
        }

        Scalar a = diagonal_entry - gershgorin_radius;
        Scalar b = diagonal_entry + gershgorin_radius;
        if(a*b < 0.0) { // a and b have different sign
          minAbsEigValue = 0.0; // zero is possible EW
          if(std::max(abs(a),abs(b)) > maxAbsEigValue) maxAbsEigValue = std::max(abs(a),abs(b));
        } else { // a and b have same sign
          if(std::min(abs(a),abs(b)) < minAbsEigValue) minAbsEigValue = std::min(abs(a),abs(b));
          if(std::max(abs(a),abs(b)) > maxAbsEigValue) maxAbsEigValue = std::max(abs(a),abs(b));
        }
      }
    }

    Scalar globalMinAbsEWApprox = 0.0;
    Scalar globalMaxAbsEWApprox = 0.0;

    minAll(A->getRowMap()->getComm(), minAbsEigValue, globalMinAbsEWApprox);
    maxAll(A->getRowMap()->getComm(), maxAbsEigValue, globalMaxAbsEWApprox);

    GetOStream(Statistics0,0) << varName_ << "(" << factory_ << "): min |EW(" << varName_ << ")|=" << globalMinAbsEWApprox << ",  max |EW(" << varName_ << ")|=" << globalMaxAbsEWApprox << std::endl;

  }

} // namespace MueLu


#endif /* MUELU_VARIABLEINFORMATIONFACTORY_DEF_HPP_ */
