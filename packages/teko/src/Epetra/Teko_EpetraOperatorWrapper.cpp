// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER


#include "Teko_EpetraOperatorWrapper.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "Teuchos_DefaultSerialComm.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Thyra_EpetraThyraWrappers.hpp"

// #include "Thyra_LinearOperator.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Teko_EpetraThyraConverter.hpp"
#include "Teuchos_Ptr.hpp"


namespace Teko { 
namespace Epetra { 


using namespace Teuchos;
using namespace Thyra;

DefaultMappingStrategy::DefaultMappingStrategy(const RCP<const Thyra::LinearOpBase<double> > & thyraOp,const Epetra_Comm & comm)
{
   RCP<Epetra_Comm> newComm = rcp(comm.Clone());

   // extract vector spaces from linear operator
   domainSpace_ = thyraOp->domain();
   rangeSpace_ = thyraOp->range();

   domainMap_ = Teko::Epetra::thyraVSToEpetraMap(*domainSpace_,newComm);
   rangeMap_ = Teko::Epetra::thyraVSToEpetraMap(*rangeSpace_,newComm);
}

void DefaultMappingStrategy::copyEpetraIntoThyra(const Epetra_MultiVector& x, const Ptr<Thyra::MultiVectorBase<double> > & thyraVec) const
{
   Teko::Epetra::blockEpetraToThyra(x,thyraVec);
}

void DefaultMappingStrategy::copyThyraIntoEpetra(const RCP<const Thyra::MultiVectorBase<double> > & thyraVec, Epetra_MultiVector& v) const
{
   Teko::Epetra::blockThyraToEpetra(thyraVec,v);
}

EpetraOperatorWrapper::EpetraOperatorWrapper()
{
   useTranspose_ = false;
   mapStrategy_ = Teuchos::null;
   thyraOp_ = Teuchos::null;
   comm_ = Teuchos::null;
   label_ = Teuchos::null;
}

EpetraOperatorWrapper::EpetraOperatorWrapper(const RCP<const Thyra::LinearOpBase<double> > & thyraOp)
{ 
   SetOperator(thyraOp); 
}

EpetraOperatorWrapper::EpetraOperatorWrapper(const RCP<const Thyra::LinearOpBase<double> > & thyraOp,const RCP<const MappingStrategy> & mapStrategy)
   : mapStrategy_(mapStrategy)
{
   SetOperator(thyraOp);
}

EpetraOperatorWrapper::EpetraOperatorWrapper(const RCP<const MappingStrategy> & mapStrategy)
   : mapStrategy_(mapStrategy)
{
   useTranspose_ = false;
   thyraOp_ = Teuchos::null;
   comm_ = Teuchos::null;
   label_ = Teuchos::null;
}

void EpetraOperatorWrapper::SetOperator(const RCP<const LinearOpBase<double> > & thyraOp,bool buildMap)
{
   useTranspose_ = false;
   thyraOp_ = thyraOp;
   comm_ = getEpetraComm(*thyraOp);
   label_ = thyraOp_->description();
   if(mapStrategy_==Teuchos::null && buildMap)
      mapStrategy_ = Teuchos::rcp(new DefaultMappingStrategy(thyraOp,*comm_));
}

double EpetraOperatorWrapper::NormInf() const 
{
  TEST_FOR_EXCEPTION(true, std::runtime_error,
                     "EpetraOperatorWrapper::NormInf not implemated");
  return 1.0;
}

int EpetraOperatorWrapper::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
   if (!useTranspose_)
   {
       // allocate space for each vector
       RCP<Thyra::MultiVectorBase<double> > tX;
       RCP<Thyra::MultiVectorBase<double> > tY; 

       tX = Thyra::createMembers(thyraOp_->domain(),X.NumVectors()); 
       tY = Thyra::createMembers(thyraOp_->range(),X.NumVectors());

       Thyra::assign(tX.ptr(),0.0);
       Thyra::assign(tY.ptr(),0.0);

       // copy epetra X into thyra X
       mapStrategy_->copyEpetraIntoThyra(X, tX.ptr());
       mapStrategy_->copyEpetraIntoThyra(Y, tY.ptr()); // if this matrix isn't block square, this probably won't work!

       // perform matrix vector multiplication
       thyraOp_->apply(Thyra::NOTRANS,*tX,tY.ptr(),1.0,0.0);

       // copy thyra Y into epetra Y
       mapStrategy_->copyThyraIntoEpetra(tY, Y);
   }
   else
   {
       TEUCHOS_ASSERT(false);
   }
 
   return 0;
}


int EpetraOperatorWrapper::ApplyInverse(const Epetra_MultiVector& X, 
                                      Epetra_MultiVector& Y) const
{
  TEST_FOR_EXCEPTION(true, std::runtime_error,
                     "EpetraOperatorWrapper::ApplyInverse not implemented");
  return 1;
}


RCP<const Epetra_Comm> 
EpetraOperatorWrapper::getEpetraComm(const Thyra::LinearOpBase<double>& inOp) const
{
  RCP<const VectorSpaceBase<double> > vs = inOp.domain();

  RCP<const SpmdVectorSpaceBase<double> > spmd;
  RCP<const VectorSpaceBase<double> > current = vs;
  while(current!=Teuchos::null) {
     // try to cast to a product vector space first
     RCP<const ProductVectorSpaceBase<double> > prod
           = rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(current);

     // figure out what type it is
     if(prod==Teuchos::null) {
        // hopfully this is a SPMD vector space
        spmd = rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(current);

        break;
     }
     else // get first convenient vector space
        current = prod->getBlock(0);
  }

  TEST_FOR_EXCEPTION(spmd==Teuchos::null, std::runtime_error, 
                     "EpetraOperatorWrapper requires std::vector space "
                     "blocks to be SPMD std::vector spaces");

  return Thyra::get_Epetra_Comm(*spmd->getComm());
/*
  const Thyra::ConstLinearOperator<double> thyraOp = rcpFromRef(inOp); 

  RCP<Epetra_Comm> rtn;
  // VectorSpace<double> vs = thyraOp.domain().getBlock(0);
  RCP<const VectorSpaceBase<double> > vs = thyraOp.domain().getBlock(0).constPtr();

  // search for an SpmdVectorSpaceBase object
  RCP<const SpmdVectorSpaceBase<double> > spmd;
  RCP<const VectorSpaceBase<double> > current = vs;
  while(current!=Teuchos::null) {
     // try to cast to a product vector space first
     RCP<const ProductVectorSpaceBase<double> > prod
           = rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(current);

     // figure out what type it is
     if(prod==Teuchos::null) {
        // hopfully this is a SPMD vector space
        spmd = rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(current);

        break;
     }
     else {
        // get first convenient vector space
        current = prod->getBlock(0);
     }
  }

  TEST_FOR_EXCEPTION(spmd==Teuchos::null, std::runtime_error, 
                     "EpetraOperatorWrapper requires std::vector space "
                     "blocks to be SPMD std::vector spaces");

  const SerialComm<Thyra::Ordinal>* serialComm 
    = dynamic_cast<const SerialComm<Thyra::Ordinal>*>(spmd->getComm().get());

#ifdef HAVE_MPI
  const MpiComm<Thyra::Ordinal>* mpiComm 
    = dynamic_cast<const MpiComm<Thyra::Ordinal>*>(spmd->getComm().get());

  TEST_FOR_EXCEPTION(mpiComm==0 && serialComm==0, std::runtime_error, 
                     "SPMD std::vector space has a communicator that is "
                     "neither a serial comm nor an MPI comm");

  if (mpiComm != 0)
    {
      rtn = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
    }
  else
    {
      rtn = rcp(new Epetra_SerialComm());
    }
#else
  TEST_FOR_EXCEPTION(serialComm==0, std::runtime_error, 
                     "SPMD std::vector space has a communicator that is "
                     "neither a serial comm nor an MPI comm");
  rtn = rcp(new Epetra_SerialComm());
  
#endif

  TEST_FOR_EXCEPTION(rtn.get()==0, std::runtime_error, "null communicator created");
  return rtn;
*/
}

int EpetraOperatorWrapper::GetBlockRowCount()
{
   const RCP<const Thyra::BlockedLinearOpBase<double> > blkOp 
         = Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<double> >(getThyraOp());

   return blkOp->productRange()->numBlocks();
}

int EpetraOperatorWrapper::GetBlockColCount()
{
   const RCP<const Thyra::BlockedLinearOpBase<double> > blkOp 
         = Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<double> >(getThyraOp());

   return blkOp->productDomain()->numBlocks();
}

Teuchos::RCP<const Epetra_Operator> EpetraOperatorWrapper::GetBlock(int i,int j) const
{
   const RCP<const Thyra::BlockedLinearOpBase<double> > blkOp 
         = Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<double> >(getThyraOp());

   return Thyra::get_Epetra_Operator(*blkOp->getBlock(i,j));
}

} // namespace Epetra
} // namespace Teko
