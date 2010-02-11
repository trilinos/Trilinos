// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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
  void copyEpetraIntoThyra(const Epetra_MultiVector &x,
    const Ptr<VectorBase<double> > &thyraVec) const;

  /** \brief . */
  void copyThyraIntoEpetra(const VectorBase<double> &thyraVec,
    Epetra_MultiVector &x) const;

  //@}

  /** \name Overridden from Epetra_Operator */
  //@{

  /** \brief . */
  int SetUseTranspose(bool UseTranspose)
    {
      useTranspose_ = UseTranspose;
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
