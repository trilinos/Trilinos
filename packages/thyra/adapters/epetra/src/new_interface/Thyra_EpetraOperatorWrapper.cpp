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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Thyra_EpetraOperatorWrapper.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraMultiVector.hpp"
#include "Thyra_EpetraLinearOp.hpp"

namespace Thyra { 

// Constructor, utilties


EpetraOperatorWrapper::
EpetraOperatorWrapper(const RCP<const LinearOpBase<double> > &thyraOp)
  : useTranspose_(false),
    thyraOp_(thyraOp)
{
  TEUCHOS_TEST_FOR_EXCEPTION (thyraOp_.is_null(), std::logic_error,
                              "Error! EpetraOperatorWrapper requires a nonnull LinearOpBase ptr as input.\n");
  range_     = thyraOp->range();
  domain_    = thyraOp->domain();
  label_     = thyraOp->description();
  domainMap_ = EpetraOperatorVectorExtraction::getOrCreateEpetraMap(domain_);
  rangeMap_  = EpetraOperatorVectorExtraction::getOrCreateEpetraMap(range_);
}


// Overridden from Epetra_Operator


int EpetraOperatorWrapper::
Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y) const
{
  // Wrap inputs into EpetraMultiVector's
  const RCP<const MultiVectorBase<double>> eX = createConstMultiVector(Teuchos::rcpFromRef(X));
  const RCP<MultiVectorBase<double>> eY = createMultiVector(Teuchos::rcpFromRef(Y));

  Thyra::apply<double>( *thyraOp_, !useTranspose_ ? NOTRANS : CONJTRANS, *eX, eY.ptr());

  return 0;
}


int EpetraOperatorWrapper::ApplyInverse(const Epetra_MultiVector& X, 
  Epetra_MultiVector& Y) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "EpetraOperatorWrapper::ApplyInverse not implemented");
  TEUCHOS_UNREACHABLE_RETURN(1);
}


double EpetraOperatorWrapper::NormInf() const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "EpetraOperatorWrapper::NormInf not implemated");
  TEUCHOS_UNREACHABLE_RETURN(1.0);
}


// Free function
Teuchos::RCP<const Thyra::LinearOpBase<double> > 
makeEpetraWrapper(const RCP<const LinearOpBase<double> > &thyraOp)
{
  return constEpetraLinearOp(Teuchos::rcp(new EpetraOperatorWrapper(thyraOp)));
}

} // namespace Thyra
