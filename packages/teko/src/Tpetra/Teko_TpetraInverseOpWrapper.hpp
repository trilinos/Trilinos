// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_TpetraInverseOpWrapper_hpp__
#define __Teko_TpetraInverseOpWrapper_hpp__

#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_ConfigDefs.hpp"
#include "Tpetra_MultiVector.hpp"

namespace Teko {
namespace TpetraHelpers {

class TpetraInverseOpWrapper : public TpetraOperatorWrapper {
 public:
  TpetraInverseOpWrapper(const RCP<const MappingStrategy>& forwardMaps)
      : TpetraOperatorWrapper(forwardMaps) {}

  TpetraInverseOpWrapper(const RCP<const Thyra::LinearOpBase<ST> >& thyraOp)
      : TpetraOperatorWrapper(thyraOp) {}

  /** */
  virtual void apply(const Tpetra::MultiVector<ST, LO, GO, NT>& X,
                     Tpetra::MultiVector<ST, LO, GO, NT>& Y,
                     Teuchos::ETransp mode = Teuchos::NO_TRANS,
                     ST alpha              = Teuchos::ScalarTraits<ST>::one(),
                     ST beta               = Teuchos::ScalarTraits<ST>::zero()) const;

  /** */
  virtual void applyInverse(const Tpetra::MultiVector<ST, LO, GO, NT>& X,
                            Tpetra::MultiVector<ST, LO, GO, NT>& Y,
                            Teuchos::ETransp mode = Teuchos::NO_TRANS,
                            ST alpha              = Teuchos::ScalarTraits<ST>::one(),
                            ST beta               = Teuchos::ScalarTraits<ST>::zero()) const;

  //    /** */
  //    virtual const Epetra_Map& OperatorDomainMap() const;
  //
  //    /** */
  //    virtual const Epetra_Map& OperatorRangeMap() const;
 protected:
  TpetraInverseOpWrapper() {}
};

}  // namespace TpetraHelpers
}  // namespace Teko

#endif
