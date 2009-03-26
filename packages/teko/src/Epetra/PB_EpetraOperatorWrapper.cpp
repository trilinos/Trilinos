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

#include "PB_EpetraOperatorWrapper.hpp"
#include "Thyra_SpmdVectorBase.hpp"
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

#include "PB_EpetraHelpers.hpp"
#include "PB_EpetraThyraConverter.hpp"
#include "Teuchos_Ptr.hpp"


namespace PB { 
namespace Epetra { 


using namespace Teuchos;
using namespace Thyra;

DefaultMappingStrategy::DefaultMappingStrategy(const RCP<const Thyra::LinearOpBase<double> > & thyraOp, Epetra_Comm & comm)
{
   // extract vector spaces from linear operator
   domainSpace_ = thyraOp->domain();
   rangeSpace_ = thyraOp->range();

   domainMap_ = PB::Epetra::thyraVSToEpetraMap(*domainSpace_,Teuchos::rcpFromRef(comm));
   rangeMap_ = PB::Epetra::thyraVSToEpetraMap(*rangeSpace_,Teuchos::rcpFromRef(comm));
}

void DefaultMappingStrategy::copyEpetraIntoThyra(const Epetra_MultiVector& x, const Ptr<Thyra::MultiVectorBase<double> > & thyraVec,
                                      const EpetraOperatorWrapper & eow) const
{
   // epetraToThyra(x,Teuchos::ptr_dynamic_cast<Thyra::VectorBase<double> >(thyraVec));
   PB::Epetra::blockEpetraToThyra(x,thyraVec);
}

void DefaultMappingStrategy::copyThyraIntoEpetra(const RCP<const Thyra::MultiVectorBase<double> > & thyraVec, Epetra_MultiVector& v,
                                      const EpetraOperatorWrapper & eow) const
{
   // thyraToEpetra(rcp_dynamic_cast<const Thyra::VectorBase<double> >(thyraVec), v);
   PB::Epetra::blockThyraToEpetra(thyraVec,v);
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
   comm_ = getEpetraComm(thyraOp);
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
       mapStrategy_->copyEpetraIntoThyra(X, tX.ptr(),*this);

       // perform matrix vector multiplication
       thyraOp_->apply(Thyra::NONCONJ_ELE,*tX,&*tY);

       // copy thyra Y into epetra Y
       mapStrategy_->copyThyraIntoEpetra(tY, Y,*this);
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


RCP<Epetra_Comm> 
EpetraOperatorWrapper::getEpetraComm(const ConstLinearOperator<double>& thyraOp) const
{
  RCP<Epetra_Comm> rtn;
  VectorSpace<double> vs = thyraOp.domain().getBlock(0);

  RCP<const SpmdVectorSpaceBase<double> > spmd 
    = rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(vs.constPtr());

  TEST_FOR_EXCEPTION(!isSPMD(vs), std::runtime_error, 
                     "EpetraOperatorWrapper requires std::vector space "
                     "blocks to be SPMD std::vector spaces");


  const SerialComm<int>* serialComm 
    = dynamic_cast<const SerialComm<int>*>(spmd->getComm().get());

#ifdef HAVE_MPI
  const MpiComm<int>* mpiComm 
    = dynamic_cast<const MpiComm<int>*>(spmd->getComm().get());

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
}


} // namespace Epetra
} // namespace PB
