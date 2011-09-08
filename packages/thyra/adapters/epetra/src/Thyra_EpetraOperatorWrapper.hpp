// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
