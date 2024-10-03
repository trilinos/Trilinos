// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_EPETRA_OPERATOR_WRAPPER_HPP
#define THYRA_EPETRA_OPERATOR_WRAPPER_HPP

#include "Thyra_LinearOpBase.hpp"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"


namespace Thyra {


/** \brief Implements the Epetra_Operator interface with a Thyra LinearOperator.
 *
 * This enables the use of absrtact Thyra operators in AztecOO as
 * preconditioners and operators, without being rendered into concrete Epetra
 * matrices.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class EpetraOperatorWrapper : public Epetra_Operator
{
public:

  /** \name Constructor, utilties. */
  //@{

  /** \brief . */
  EpetraOperatorWrapper(const RCP<const LinearOpBase<double> > &thyraOp);

  /** \brief . */
  RCP<const LinearOpBase<double> > getThyraOp() const
  {
    return thyraOp_;
  }

  /** \brief . */
  void copyEpetraIntoThyra(const Epetra_MultiVector &x,
    const Ptr<VectorBase<double> > &thyraVec) const;

  /** \brief . */
  void copyThyraIntoEpetra(const VectorBase<double> &thyraVec,
    Epetra_MultiVector &x) const;

  //@}

  /** \name Overridden from Epetra_Operator */
  //@{

  /** \brief . */
  int SetUseTranspose(bool UseTranspose_in)
    {
      useTranspose_ = UseTranspose_in;
      return 0;
    }

  /** \brief . */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

  /** \brief . */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

  /** \brief . */
  double NormInf() const ;

  /** \brief . */
  const char* Label() const {return label_.c_str();}

  /** \brief . */
  bool UseTranspose() const {return useTranspose_;}

  /** \brief . */
  bool HasNormInf() const {return false;}
    
  /** \brief . */
  const Epetra_Comm& Comm() const {return *comm_;}

  /** \brief . */
  const Epetra_Map& OperatorDomainMap() const {return *domainMap_;}

  /** \brief . */
  const Epetra_Map& OperatorRangeMap() const {return *rangeMap_;}

  //@}

private:

  bool useTranspose_;
  RCP<const LinearOpBase<double> > thyraOp_;
  RCP<const VectorSpaceBase<double> > range_;
  RCP<const VectorSpaceBase<double> > domain_;
  RCP<const Epetra_Comm> comm_;
  RCP<const Epetra_Map> rangeMap_;
  RCP<const Epetra_Map> domainMap_;

  std::string label_;

  static RCP<const Epetra_Comm> getEpetraComm(const LinearOpBase<double>& thyraOp);

};


/** \brief Wrap a Thyra operator in the Epetra_Operator interface, and then
 * wrap it again in a Thyra operator interface.
 *
 * This lets an arbitrary Thyra operator be given to the Thyra AztecOO
 * adapters.
 *
 * \relates EpetraOperatorWrapper
 */
RCP<const LinearOpBase<double> > 
makeEpetraWrapper(const RCP<const LinearOpBase<double> > &thyraOp);


}  // namespace Thyra


#endif // THYRA_EPETRA_OPERATOR_WRAPPER_HPP

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraAdapters package is deprecated"
#endif
#endif

