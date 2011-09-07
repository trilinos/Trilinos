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

#ifndef THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_FACTORY_HPP


#include "Thyra_DefaultBlockedTriangularLinearOpWithSolveFactory_decl.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_PhysicallyBlockedLinearOpWithSolveBase.hpp" // Interface
#include "Thyra_DefaultBlockedTriangularLinearOpWithSolve.hpp" // Implementation
#include "Thyra_DefaultLinearOpSource.hpp"


namespace Thyra {


// Overridden from Constructors/Initializers/Accessors

  
template<class Scalar>
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::DefaultBlockedTriangularLinearOpWithSolveFactory(
  const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(lowsf));
#endif
  lowsf_.initialize(lowsf);
}

  
template<class Scalar>
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::DefaultBlockedTriangularLinearOpWithSolveFactory(
  const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(lowsf));
#endif
  lowsf_.initialize(lowsf);
}


template<class Scalar>
RCP<LinearOpWithSolveFactoryBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF()
{
  return lowsf_.getNonconstObj();
}


template<class Scalar>
RCP<const LinearOpWithSolveFactoryBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF() const
{
  return lowsf_.getConstObj();
}


// Overridden from Teuchos::Describable


template<class Scalar>
std::string
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::description() const
{
  std::ostringstream oss;
  oss << this->Teuchos::Describable::description()
      << "{"
      << "lowsf=";
  if (!is_null(lowsf_.getConstObj()))
    oss << lowsf_.getConstObj()->description();
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


// Overridden from ParameterListAcceptor


template<class Scalar>
void
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  lowsf_.getNonconstObj()->setParameterList(paramList);
}


template<class Scalar>
RCP<ParameterList>
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::getNonconstParameterList()
{
  return lowsf_.getNonconstObj()->getNonconstParameterList();
}


template<class Scalar>
RCP<ParameterList> 
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  return lowsf_.getNonconstObj()->unsetParameterList();
}


template<class Scalar>
RCP<const ParameterList>
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return lowsf_.getConstObj()->getParameterList();
}


template<class Scalar>
RCP<const ParameterList>
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::getValidParameters() const
{
  return lowsf_.getConstObj()->getValidParameters();
}


// Overridden from LinearOpWithSolveFactoyBase


template<class Scalar>
bool
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::acceptsPreconditionerFactory() const
{
  return false;
}


template<class Scalar>
void
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::setPreconditionerFactory(
  const RCP<PreconditionerFactoryBase<Scalar> > &precFactory,
  const std::string &precFactoryName
  )
{
  TEST_FOR_EXCEPTION(true,std::logic_error,
    "Error, we don't support a preconditioner factory!");
}


template<class Scalar>
RCP<PreconditionerFactoryBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::getPreconditionerFactory() const
{
  return Teuchos::null;
}


template<class Scalar>
void DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::unsetPreconditionerFactory(
  RCP<PreconditionerFactoryBase<Scalar> > *precFactory,
  std::string *precFactoryName
  )
{
  TEST_FOR_EXCEPTION(true,std::logic_error,
    "Error, we don't support a preconditioner factory!");
}


template<class Scalar>
bool
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::isCompatible(
  const LinearOpSourceBase<Scalar> &fwdOpSrc
  ) const
{
  TEST_FOR_EXCEPT(true);
  return false;
}


template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::createOp() const
{
  return defaultBlockedTriangularLinearOpWithSolve<Scalar>();
}


template<class Scalar>
void
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::initializeOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{

  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==Op);
#endif

  // Set the verbosity settings for the wrapped LOWSF object!
  lowsf_.getConstObj()->setOStream(this->getOStream());
  lowsf_.getConstObj()->setVerbLevel(this->getVerbLevel());

  // Get the block interface to get at the blocks
  typedef PhysicallyBlockedLinearOpBase<Scalar> PBLOB;
  const RCP<const PBLOB> blo =
    rcp_dynamic_cast<const PBLOB>(fwdOpSrc->getOp().assert_not_null());
  
  // Dynamic cast to get the DefaultBlockedTriangularLinearOpWithSolveBase
  // interface that we will fill.

  typedef DefaultBlockedTriangularLinearOpWithSolve<Scalar> DBTLOWS;
  DBTLOWS &btlows = dyn_cast<DBTLOWS>(*Op);

  // Determine if this is the first time through or if we have already
  // initialized before.  This will be needed to allow efficient reuse of the
  // LOWSB objects for the diagonal blocks.
  const bool firstTime = is_null(btlows.range());

  // If this is the first time through, we need to fill and create the block
  // structure
  if (firstTime)
    btlows.beginBlockFill(blo->productRange(),blo->productDomain());

  const int N = blo->productRange()->numBlocks();
  for ( int k = 0; k < N; ++k ) {
    const RCP<const LinearOpBase<Scalar> > fwdOp_k =
      blo->getBlock(k,k).assert_not_null();
    if (firstTime) {
      // This is the first time through so reate and initialize a new LOWSB
      // object for each block
      btlows.setNonconstLOWSBlock( k, k,
        linearOpWithSolve<Scalar>(*lowsf_.getConstObj(),fwdOp_k)
        );
    }
    else {
      // This is not the first time through so we need to just reinitiallize
      // the object that is already created.  This allows us to efficiently
      // reuse precreated structure and storage.
      RCP<LinearOpWithSolveBase<Scalar> >
        invOp_k = btlows.getNonconstLOWSBlock(k,k).assert_not_null();
      Thyra::initializeOp<Scalar>(*lowsf_.getConstObj(), fwdOp_k, invOp_k.ptr());
    }
  }

  // If this is the first time through, then we need to finalize the block
  // structure.
  if (firstTime)
    btlows.endBlockFill();

  // After the block structure has been setup, set the off-diagonal blocks.
  // Note that this also sets the diagonal blocks but these are ignored since
  // the LOWSB blocks created above override these.
  btlows.setBlocks(blo);

  // Set the verbosity settings
  btlows.setOStream(this->getOStream());
  btlows.setVerbLevel(this->getVerbLevel());

}


template<class Scalar>
void
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::initializeAndReuseOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op
  ) const
{
  TEST_FOR_EXCEPT(true);
}


template<class Scalar>
void
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::uninitializeOp(
  LinearOpWithSolveBase<Scalar> *Op,
  RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
  RCP<const PreconditionerBase<Scalar> > *prec,
  RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
  ESupportSolveUse *supportSolveUse
  ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcp_dynamic_cast;
  typedef DefaultBlockedTriangularLinearOpWithSolve<Scalar> DBTLOWS;
  TEST_FOR_EXCEPT(0==Op);
  DBTLOWS &btlowsOp = dyn_cast<DBTLOWS>(*Op);
  if (fwdOpSrc) {
    const RCP<const LinearOpBase<Scalar> > fwdOp = btlowsOp.getBlocks();
    if (!is_null(fwdOp))
      *fwdOpSrc = defaultLinearOpSource<Scalar>(fwdOp);
    else
      *fwdOpSrc = Teuchos::null;
  }
  if (prec) *prec = Teuchos::null;
  if (approxFwdOpSrc) *approxFwdOpSrc = Teuchos::null;
}


template<class Scalar>
bool
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::supportsPreconditionerInputType(
  const EPreconditionerInputType precOpType
  ) const
{
  // We don't support any external preconditioners!
  return false;
  // 20071006: rabartl: Note: We could support external preconditioners but it
  // will take some work.  We would have to extract out the individual
  // preconditioners from each block.  This would be pretty easy to do but I
  // am not going to do this until we have to.
}


template<class Scalar>
void
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const PreconditionerBase<Scalar> > &prec,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  TEST_FOR_EXCEPTION(true,std::logic_error,
    "Error, we don't support an external preconditioner!");
}


template<class Scalar>
void
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::initializeApproxPreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  TEST_FOR_EXCEPTION(true,std::logic_error,
    "Error, we don't support an external preconditioner!");
}


// protected


template<class Scalar>
void
DefaultBlockedTriangularLinearOpWithSolveFactory<Scalar>::informUpdatedVerbosityState() const
{
  lowsf_.getConstObj()->setVerbLevel(this->getVerbLevel());
  lowsf_.getConstObj()->setOStream(this->getOStream());
}


} // namespace Thyra


#endif // THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
