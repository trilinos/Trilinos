// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_TpetraOperatorWrapper_hpp__
#define __Teko_TpetraOperatorWrapper_hpp__

#include "Thyra_LinearOpBase.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#include "Teko_ConfigDefs.hpp"

#include <string>

namespace Teko {
namespace TpetraHelpers {
using Teuchos::RCP;

class TpetraOperatorWrapper;

/// Abstract Mapping strategy for an TpetraOperatorWrapper
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
  virtual void copyTpetraIntoThyra(
      const Tpetra::MultiVector<ST, LO, GO, NT>& tpetraX,
      const Teuchos::Ptr<Thyra::MultiVectorBase<ST> >& thyraX) const = 0;
  // const TpetraOperatorWrapper & eow) const = 0;

  /** \brief Copy an Thyra::MultiVectorBase into a Epetra_MultiVector
   *
   * Copy an Thyra::MultiVectorBase into an Epetra_MultiVector. The exact
   * method for copying is specified by the concrete implementations.
   *
   * \param[in]     thyraX  Source Thyra object
   * \param[in,out] epetraX Destination Epetra object
   */
  virtual void copyThyraIntoTpetra(const RCP<const Thyra::MultiVectorBase<ST> >& thyraX,
                                   Tpetra::MultiVector<ST, LO, GO, NT>& tpetraX) const = 0;
  // const TpetraOperatorWrapper & eow) const = 0;

  /** \brief Domain map for this strategy */
  virtual const RCP<const Tpetra::Map<LO, GO, NT> > domainMap() const = 0;

  /** \brief Range map for this strategy */
  virtual const RCP<const Tpetra::Map<LO, GO, NT> > rangeMap() const = 0;

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

  virtual void copyTpetraIntoThyra(const Tpetra::MultiVector<ST, LO, GO, NT>& tpetraX,
                                   const Teuchos::Ptr<Thyra::MultiVectorBase<ST> >& thyraX) const
  // const TpetraOperatorWrapper & eow) const
  {
    forwardStrategy_->copyTpetraIntoThyra(tpetraX, thyraX);
  }

  virtual void copyThyraIntoTpetra(const RCP<const Thyra::MultiVectorBase<ST> >& thyraX,
                                   Tpetra::MultiVector<ST, LO, GO, NT>& tpetraX) const
  // const TpetraOperatorWrapper & eow) const
  {
    forwardStrategy_->copyThyraIntoTpetra(thyraX, tpetraX);
  }

  /** \brief Domain map for this strategy */
  virtual const RCP<const Tpetra::Map<LO, GO, NT> > domainMap() const {
    return forwardStrategy_->rangeMap();
  }

  /** \brief Range map for this strategy */
  virtual const RCP<const Tpetra::Map<LO, GO, NT> > rangeMap() const {
    return forwardStrategy_->domainMap();
  }

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

/// default mapping strategy for the basic TpetraOperatorWrapper
class DefaultMappingStrategy : public MappingStrategy {
 public:
  /** */
  DefaultMappingStrategy(const RCP<const Thyra::LinearOpBase<ST> >& thyraOp,
                         const Teuchos::Comm<Thyra::Ordinal>& comm);

  virtual ~DefaultMappingStrategy() {}

  /** \brief Copy an Epetra_MultiVector into a Thyra::MultiVectorBase
   *
   * Copy an Epetra_MultiVector into a Thyra::MultiVectorBase. The exact
   * method for copying is specified by the concrete implementations.
   *
   * \param[in]     epetraX Vector to be copied into the Thyra object
   * \param[in,out] thyraX  Destination Thyra object
   */
  virtual void copyTpetraIntoThyra(const Tpetra::MultiVector<ST, LO, GO, NT>& tpetraX,
                                   const Teuchos::Ptr<Thyra::MultiVectorBase<ST> >& thyraX) const;
  // const TpetraOperatorWrapper & eow) const;

  /** \brief Copy an Thyra::MultiVectorBase into a Epetra_MultiVector
   *
   * Copy an Thyra::MultiVectorBase into an Epetra_MultiVector. The exact
   * method for copying is specified by the concrete implementations.
   *
   * \param[in]     thyraX  Source Thyra object
   * \param[in,out] epetraX Destination Epetra object
   */
  virtual void copyThyraIntoTpetra(const RCP<const Thyra::MultiVectorBase<ST> >& thyraX,
                                   Tpetra::MultiVector<ST, LO, GO, NT>& tpetraX) const;
  // const TpetraOperatorWrapper & eow) const;

  /** \brief Domain map for this strategy */
  virtual const RCP<const Tpetra::Map<LO, GO, NT> > domainMap() const { return domainMap_; }

  /** \brief Range map for this strategy */
  virtual const RCP<const Tpetra::Map<LO, GO, NT> > rangeMap() const { return rangeMap_; }

  /** \brief Identifier string */
  virtual std::string toString() const { return std::string("DefaultMappingStrategy"); }

 protected:
  RCP<const Thyra::VectorSpaceBase<ST> > domainSpace_;  ///< Domain space object
  RCP<const Thyra::VectorSpaceBase<ST> > rangeSpace_;   ///< Range space object

  RCP<const Tpetra::Map<LO, GO, NT> > domainMap_;  ///< Pointer to the constructed domain map
  RCP<const Tpetra::Map<LO, GO, NT> > rangeMap_;   ///< Pointer to the constructed range map
};

/** \brief
 * Implements the Epetra_Operator interface with a Thyra LinearOperator. This
 * enables the use of absrtact Thyra operators in AztecOO as preconditioners and
 * operators, without being rendered into concrete Epetra matrices. This is my own
 * modified version that was originally in Thyra.
 */
class TpetraOperatorWrapper : public Tpetra::Operator<ST, LO, GO, NT> {
 public:
  /** */
  TpetraOperatorWrapper(const RCP<const Thyra::LinearOpBase<ST> >& thyraOp);
  TpetraOperatorWrapper(const RCP<const Thyra::LinearOpBase<ST> >& thyraOp,
                        const RCP<const MappingStrategy>& mapStrategy);
  TpetraOperatorWrapper(const RCP<const MappingStrategy>& mapStrategy);

  /** */
  virtual ~TpetraOperatorWrapper() { ; }

  /** */
  int SetUseTranspose(bool useTranspose) {
    useTranspose_ = useTranspose;
    return 0;
  }

  /** */
  void apply(const Tpetra::MultiVector<ST, LO, GO, NT>& X, Tpetra::MultiVector<ST, LO, GO, NT>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS, ST alpha = Teuchos::ScalarTraits<ST>::one(),
             ST beta = Teuchos::ScalarTraits<ST>::zero()) const;

  /** */
  void applyInverse(const Tpetra::MultiVector<ST, LO, GO, NT>& X,
                    Tpetra::MultiVector<ST, LO, GO, NT>& Y,
                    Teuchos::ETransp mode = Teuchos::NO_TRANS,
                    ST alpha              = Teuchos::ScalarTraits<ST>::one(),
                    ST beta               = Teuchos::ScalarTraits<ST>::zero()) const;

  /** */
  double NormInf() const;

  /** */
  const char* Label() const { return label_.c_str(); }

  /** */
  bool UseTranspose() const { return useTranspose_; }

  /** */
  bool HasNormInf() const { return false; }

  /** */
  const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> >& Comm() const { return comm_; }

  /** */
  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > getDomainMap() const {
    return mapStrategy_->domainMap();
  }

  /** */
  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > getRangeMap() const {
    return mapStrategy_->rangeMap();
  }

  //! Return the thyra operator associated with this wrapper
  const RCP<const Thyra::LinearOpBase<ST> > getThyraOp() const { return thyraOp_; }

  //! Get the mapping strategy for this wrapper (translate between Thyra and Epetra)
  const RCP<const MappingStrategy> getMapStrategy() const { return mapStrategy_; }

  //! Get the number of block rows in this operator
  virtual int GetBlockRowCount();

  //! Get the number of block columns in this operator
  virtual int GetBlockColCount();

  //! Grab the i,j block
  Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > GetBlock(int i, int j) const;

 protected:
  /** */
  TpetraOperatorWrapper();

  /** */
  RCP<const Teuchos::Comm<Thyra::Ordinal> > getThyraComm(const Thyra::LinearOpBase<ST>& inOp) const;

  /** */
  void SetOperator(const RCP<const Thyra::LinearOpBase<ST> >& thyraOp, bool buildMap = true);

  /** */
  void SetMapStrategy(const RCP<const MappingStrategy>& mapStrategy) { mapStrategy_ = mapStrategy; }

  /** */
  RCP<const MappingStrategy> mapStrategy_;

  /** */
  RCP<const Thyra::LinearOpBase<ST> > thyraOp_;

  /** */
  bool useTranspose_;

  /** */
  RCP<const Teuchos::Comm<Thyra::Ordinal> > comm_;

  /** */
  std::string label_;
};
}  // namespace TpetraHelpers
}  // end namespace Teko

#endif
