// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_EpetraInverseOpWrapper_hpp__
#define __Teko_EpetraInverseOpWrapper_hpp__

#include "Teko_EpetraOperatorWrapper.hpp"

namespace Teko {
namespace Epetra {

class EpetraInverseOpWrapper : public EpetraOperatorWrapper {
 public:
  EpetraInverseOpWrapper(const RCP<const MappingStrategy>& forwardMaps)
      : EpetraOperatorWrapper(forwardMaps) {}

  EpetraInverseOpWrapper(const RCP<const Thyra::LinearOpBase<double> >& thyraOp)
      : EpetraOperatorWrapper(thyraOp) {}

  /** */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  /** */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //    /** */
  //    virtual const Epetra_Map& OperatorDomainMap() const;
  //
  //    /** */
  //    virtual const Epetra_Map& OperatorRangeMap() const;
 protected:
  EpetraInverseOpWrapper() {}
};

}  // namespace Epetra
}  // namespace Teko

#endif
