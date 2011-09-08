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

#include "Thyra_EpetraOperatorWrapper.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"

#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "Teuchos_DefaultSerialComm.hpp"


namespace Thyra { 


// Constructor, utilties


EpetraOperatorWrapper::EpetraOperatorWrapper(
  const RCP<const LinearOpBase<double> > &thyraOp
  )
  : useTranspose_(false),
    thyraOp_(thyraOp),
    range_(thyraOp->range()),
    domain_(thyraOp->domain()),
    comm_(getEpetraComm(*thyraOp)),
    rangeMap_(get_Epetra_Map(*range_, comm_)),
    domainMap_(get_Epetra_Map(*domain_, comm_)),
    label_(thyraOp->description())
{;}


void EpetraOperatorWrapper::copyEpetraIntoThyra(const Epetra_MultiVector& x,
  const Ptr<VectorBase<double> > &thyraVec) const
{

  using Teuchos::rcpFromPtr;
  using Teuchos::rcp_dynamic_cast;

  const int numVecs = x.NumVectors();

  TEST_FOR_EXCEPTION(numVecs != 1, std::runtime_error,
    "epetraToThyra does not work with MV dimension != 1");

  const RCP<ProductVectorBase<double> > prodThyraVec =
    castOrCreateNonconstProductVectorBase(rcpFromPtr(thyraVec));

  const ArrayView<const double> epetraData(x[0], x.Map().NumMyElements());
  // NOTE: I tried using Epetra_MultiVector::operator()(int) to return an
  // Epetra_Vector object but it has a defect when Reset(...) is called which
  // results in a memory access error (see bug 4700).

  int offset = 0;
  const int numBlocks = prodThyraVec->productSpace()->numBlocks();
  for (int b = 0; b < numBlocks; ++b) {
    const RCP<VectorBase<double> > vec_b = prodThyraVec->getNonconstVectorBlock(b);
    const RCP<const SpmdVectorSpaceBase<double> > spmd_vs_b =
      rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(vec_b->space(), true);
    DetachedSpmdVectorView<double> view(vec_b);
    const ArrayRCP<double> thyraData = view.sv().values();
    const int localNumElems = spmd_vs_b->localSubDim();
    for (int i=0; i < localNumElems; ++i) {
      thyraData[i] = epetraData[i+offset];
    }
    offset += localNumElems;
  }

}


void EpetraOperatorWrapper::copyThyraIntoEpetra(
  const VectorBase<double>& thyraVec, Epetra_MultiVector& x) const 
{

  using Teuchos::rcpFromRef;
  using Teuchos::rcp_dynamic_cast;

  const int numVecs = x.NumVectors();

  TEST_FOR_EXCEPTION(numVecs != 1, std::runtime_error,
    "epetraToThyra does not work with MV dimension != 1");

  const RCP<const ProductVectorBase<double> > prodThyraVec =
    castOrCreateProductVectorBase(rcpFromRef(thyraVec));

  const ArrayView<double> epetraData(x[0], x.Map().NumMyElements());
  // NOTE: See above!

  int offset = 0;
  const int numBlocks = prodThyraVec->productSpace()->numBlocks();
  for (int b = 0; b < numBlocks; ++b) {
    const RCP<const VectorBase<double> > vec_b = prodThyraVec->getVectorBlock(b);
    const RCP<const SpmdVectorSpaceBase<double> > spmd_vs_b =
      rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(vec_b->space(), true);
    ConstDetachedSpmdVectorView<double> view(vec_b);
    const ArrayRCP<const double> thyraData = view.sv().values();
    const int localNumElems = spmd_vs_b->localSubDim();
    for (int i=0; i < localNumElems; ++i) {
      epetraData[i+offset] = thyraData[i];
    }
    offset += localNumElems;
  }

}


// Overridden from Epetra_Operator


int EpetraOperatorWrapper::Apply(const Epetra_MultiVector& X,
  Epetra_MultiVector& Y) const
{

  const RCP<const VectorSpaceBase<double> >
    opRange = ( !useTranspose_ ? range_ : domain_ ),
    opDomain = ( !useTranspose_ ? domain_ : range_ );

  const RCP<VectorBase<double> > tx = createMember(opDomain);
  copyEpetraIntoThyra(X, tx.ptr());

  const RCP<VectorBase<double> > ty = createMember(opRange);

  Thyra::apply<double>( *thyraOp_, !useTranspose_ ? NOTRANS : CONJTRANS, *tx, ty.ptr());

  copyThyraIntoEpetra(*ty, Y);

  return 0;

}


int EpetraOperatorWrapper::ApplyInverse(const Epetra_MultiVector& X, 
  Epetra_MultiVector& Y) const
{
  TEST_FOR_EXCEPTION(true, std::runtime_error,
    "EpetraOperatorWrapper::ApplyInverse not implemented");
  return 1;
}


double EpetraOperatorWrapper::NormInf() const 
{
  TEST_FOR_EXCEPTION(true, std::runtime_error,
    "EpetraOperatorWrapper::NormInf not implemated");
  return 1.0;
}


// private


RCP<const Epetra_Comm> 
EpetraOperatorWrapper::getEpetraComm(const LinearOpBase<double>& thyraOp)
{

  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::SerialComm;
#ifdef HAVE_MPI
  using Teuchos::MpiComm;
#endif

  RCP<const VectorSpaceBase<double> > vs = thyraOp.range();
  
  const RCP<const ProductVectorSpaceBase<double> > prod_vs =
    rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(vs);

  if (nonnull(prod_vs))
    vs = prod_vs->getBlock(0);

  const RCP<const SpmdVectorSpaceBase<double> > spmd_vs =
    rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(vs, true);

  return get_Epetra_Comm(*spmd_vs->getComm());

}


} // namespace Thyra


Teuchos::RCP<const Thyra::LinearOpBase<double> > 
Thyra::makeEpetraWrapper(const RCP<const LinearOpBase<double> > &thyraOp)
{
  return epetraLinearOp(
    Teuchos::rcp(new EpetraOperatorWrapper(thyraOp))
    );
}
