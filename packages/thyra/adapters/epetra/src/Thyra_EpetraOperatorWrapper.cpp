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
