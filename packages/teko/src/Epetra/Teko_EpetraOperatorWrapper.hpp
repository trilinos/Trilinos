// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_EpetraOperatorWrapper_hpp__
#define __Teko_EpetraOperatorWrapper_hpp__

#include "Thyra_LinearOpBase.hpp"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"

#include <string>

namespace Teko {
namespace Epetra {
using Teuchos::RCP;

class EpetraOperatorWrapper;

/// Abstract Mapping strategy for an EpetraOperatorWrapper
class MappingStrategy {
 public:
  virtual ~MappingStrategy() {}

  /** \brief Copy an Epetra_MultiVector into a Thyra::MultiVectorBase
   *
   * Copy an Epetra_MultiVector into a Thyra::MultiVectorBase. The exact
   * method for copying is specified by the concrete implementations.
   *
   * \param[in]     epetraX Vector to be copied into the Thyra object
   * \param[in,out] thyraX  Destination Thyra object
   */
  virtual void copyEpetraIntoThyra(
      const Epetra_MultiVector& epetraX,
      const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& thyraX) const = 0;
  // const EpetraOperatorWrapper & eow) const = 0;

  /** \brief Copy an Thyra::MultiVectorBase into a Epetra_MultiVector
   *
   * Copy an Thyra::MultiVectorBase into an Epetra_MultiVector. The exact
   * method for copying is specified by the concrete implementations.
   *
   * \param[in]     thyraX  Source Thyra object
   * \param[in,out] epetraX Destination Epetra object
   */
  virtual void copyThyraIntoEpetra(const RCP<const Thyra::MultiVectorBase<double> >& thyraX,
                                   Epetra_MultiVector& epetraX) const = 0;
  // const EpetraOperatorWrapper & eow) const = 0;

  /** \brief Domain map for this strategy */
  virtual const RCP<const Epetra_Map> domainMap() const = 0;

  /** \brief Range map for this strategy */
  virtual const RCP<const Epetra_Map> rangeMap() const = 0;

  /** \brief Identifier string */
  virtual std::string toString() const = 0;
};

/// Flip a mapping strategy object around to give the "inverse" mapping strategy.
class InverseMappingStrategy : public MappingStrategy {
 public:
  /** \brief Constructor to build a inverse MappingStrategy from
   * a forward map.
   */
  InverseMappingStrategy(const RCP<const MappingStrategy>& forward) : forwardStrategy_(forward) {}

  virtual ~InverseMappingStrategy() {}

  virtual void copyEpetraIntoThyra(
      const Epetra_MultiVector& epetraX,
      const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& thyraX) const
  // const EpetraOperatorWrapper & eow) const
  {
    forwardStrategy_->copyEpetraIntoThyra(epetraX, thyraX);
  }

  virtual void copyThyraIntoEpetra(const RCP<const Thyra::MultiVectorBase<double> >& thyraX,
                                   Epetra_MultiVector& epetraX) const
  // const EpetraOperatorWrapper & eow) const
  {
    forwardStrategy_->copyThyraIntoEpetra(thyraX, epetraX);
  }

  /** \brief Domain map for this strategy */
  virtual const RCP<const Epetra_Map> domainMap() const { return forwardStrategy_->rangeMap(); }

  /** \brief Range map for this strategy */
  virtual const RCP<const Epetra_Map> rangeMap() const { return forwardStrategy_->domainMap(); }

  /** \brief Identifier string */
  virtual std::string toString() const {
    return std::string("InverseMapping(") + forwardStrategy_->toString() + std::string(")");
  }

 protected:
  /** \brief Forward mapping strategy object */
  const RCP<const MappingStrategy> forwardStrategy_;

 private:
  InverseMappingStrategy();
  InverseMappingStrategy(const InverseMappingStrategy&);
};

/// default mapping strategy for the basic EpetraOperatorWrapper
class DefaultMappingStrategy : public MappingStrategy {
 public:
  /** */
  DefaultMappingStrategy(const RCP<const Thyra::LinearOpBase<double> >& thyraOp,
                         const Epetra_Comm& comm);

  virtual ~DefaultMappingStrategy() {}

  /** \brief Copy an Epetra_MultiVector into a Thyra::MultiVectorBase
   *
   * Copy an Epetra_MultiVector into a Thyra::MultiVectorBase. The exact
   * method for copying is specified by the concrete implementations.
   *
   * \param[in]     epetraX Vector to be copied into the Thyra object
   * \param[in,out] thyraX  Destination Thyra object
   */
  virtual void copyEpetraIntoThyra(
      const Epetra_MultiVector& epetraX,
      const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& thyraX) const;
  // const EpetraOperatorWrapper & eow) const;

  /** \brief Copy an Thyra::MultiVectorBase into a Epetra_MultiVector
   *
   * Copy an Thyra::MultiVectorBase into an Epetra_MultiVector. The exact
   * method for copying is specified by the concrete implementations.
   *
   * \param[in]     thyraX  Source Thyra object
   * \param[in,out] epetraX Destination Epetra object
   */
  virtual void copyThyraIntoEpetra(const RCP<const Thyra::MultiVectorBase<double> >& thyraX,
                                   Epetra_MultiVector& epetraX) const;
  // const EpetraOperatorWrapper & eow) const;

  /** \brief Domain map for this strategy */
  virtual const RCP<const Epetra_Map> domainMap() const { return domainMap_; }

  /** \brief Range map for this strategy */
  virtual const RCP<const Epetra_Map> rangeMap() const { return rangeMap_; }

  /** \brief Identifier string */
  virtual std::string toString() const { return std::string("DefaultMappingStrategy"); }

 protected:
  RCP<const Thyra::VectorSpaceBase<double> > domainSpace_;  ///< Domain space object
  RCP<const Thyra::VectorSpaceBase<double> > rangeSpace_;   ///< Range space object

  RCP<const Epetra_Map> domainMap_;  ///< Pointer to the constructed domain map
  RCP<const Epetra_Map> rangeMap_;   ///< Pointer to the constructed range map
};

/** \brief
 * Implements the Epetra_Operator interface with a Thyra LinearOperator. This
 * enables the use of absrtact Thyra operators in AztecOO as preconditioners and
 * operators, without being rendered into concrete Epetra matrices. This is my own
 * modified version that was originally in Thyra.
 */
class EpetraOperatorWrapper : public Epetra_Operator {
 public:
  /** */
  EpetraOperatorWrapper(const RCP<const Thyra::LinearOpBase<double> >& thyraOp);
  EpetraOperatorWrapper(const RCP<const Thyra::LinearOpBase<double> >& thyraOp,
                        const RCP<const MappingStrategy>& mapStrategy);
  EpetraOperatorWrapper(const RCP<const MappingStrategy>& mapStrategy);

  /** */
  virtual ~EpetraOperatorWrapper() { ; }

  /** */
  int SetUseTranspose(bool useTranspose) {
    useTranspose_ = useTranspose;
    return 0;
  }

  /** */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  /** */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  /** */
  double NormInf() const;

  /** */
  const char* Label() const { return label_.c_str(); }

  /** */
  bool UseTranspose() const { return useTranspose_; }

  /** */
  bool HasNormInf() const { return false; }

  /** */
  const Epetra_Comm& Comm() const { return *comm_; }

  /** */
  const Epetra_Map& OperatorDomainMap() const { return *mapStrategy_->domainMap(); }

  /** */
  const Epetra_Map& OperatorRangeMap() const { return *mapStrategy_->rangeMap(); }

  //! Return the thyra operator associated with this wrapper
  const RCP<const Thyra::LinearOpBase<double> > getThyraOp() const { return thyraOp_; }

  //! Get the mapping strategy for this wrapper (translate between Thyra and Epetra)
  const RCP<const MappingStrategy> getMapStrategy() const { return mapStrategy_; }

  //! Get the number of block rows in this operator
  virtual int GetBlockRowCount();

  //! Get the number of block columns in this operator
  virtual int GetBlockColCount();

  //! Grab the i,j block
  Teuchos::RCP<const Epetra_Operator> GetBlock(int i, int j) const;

 protected:
  /** */
  EpetraOperatorWrapper();

  /** */
  RCP<const Epetra_Comm> getEpetraComm(const Thyra::LinearOpBase<double>& inOp) const;

  /** */
  void SetOperator(const RCP<const Thyra::LinearOpBase<double> >& thyraOp, bool buildMap = true);

  /** */
  void SetMapStrategy(const RCP<const MappingStrategy>& mapStrategy) { mapStrategy_ = mapStrategy; }

  /** */
  RCP<const MappingStrategy> mapStrategy_;

  /** */
  RCP<const Thyra::LinearOpBase<double> > thyraOp_;

  /** */
  bool useTranspose_;

  /** */
  RCP<const Epetra_Comm> comm_;

  /** */
  std::string label_;
};
}  // end namespace Epetra
}  // end namespace Teko

#endif
