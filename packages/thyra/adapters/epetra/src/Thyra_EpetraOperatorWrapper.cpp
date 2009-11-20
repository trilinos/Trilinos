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


namespace Thyra { 


using namespace Teuchos;


EpetraOperatorWrapper::EpetraOperatorWrapper(const ConstLinearOperator<double>& thyraOp)
  : useTranspose_(false),
    thyraOp_(thyraOp),
    domain_(thyraOp_.domain()),
    range_(thyraOp_.range()),
    comm_(getEpetraComm(thyraOp)),
    domainMap_(thyraVSToEpetraMap(thyraOp_.domain(), comm_)),
    rangeMap_(thyraVSToEpetraMap(thyraOp_.range(), comm_)),
    label_(thyraOp_.description())
{;}


double EpetraOperatorWrapper::NormInf() const 
{
  TEST_FOR_EXCEPTION(true, std::runtime_error,
                     "EpetraOperatorWrapper::NormInf not implemated");
  return 1.0;
}


RCP<Epetra_Map> 
EpetraOperatorWrapper
::thyraVSToEpetraMap(const VectorSpace<double>& vs,
                     const RCP<Epetra_Comm>& comm) const 
{
  int globalDim = vs.dim();
  int myLocalElements = 0;
  
  /* find the number of local elements, summed over all blocks */
  if (isSPMD(vs))
    {
      myLocalElements = numLocalElements(vs);
    }
  else
    {
      for (int b=0; b<vs.numBlocks(); b++)
        {
          TEST_FOR_EXCEPTION(!isSPMD(vs.getBlock(b)), std::runtime_error, 
                             "EpetraOperatorWrapper requires std::vector space "
                             "blocks to be SPMD vectors");
          myLocalElements += numLocalElements(vs.getBlock(b));
        }
    }

  
  /* find the GIDs owned by this processor, taken from all blocks */
  Array<int> myGIDs(myLocalElements);
  
  int count=0;
  int blockOffset = 0;
  for (int b=0; b<vs.numBlocks(); b++)
    {
      int lowGIDInBlock = lowestLocallyOwnedIndex(vs.getBlock(b));
      int numLocalElementsInBlock = numLocalElements(vs.getBlock(b));
      for (int i=0; i<numLocalElementsInBlock; i++, count++)
        {
          myGIDs[count] = blockOffset+lowGIDInBlock+i;
        }

      blockOffset += vs.getBlock(b).dim();
    }

  /* create the std::map */
  RCP<Epetra_Map> rtn 
    = rcp(new Epetra_Map(globalDim, myLocalElements, &(myGIDs[0]), 0, *comm));

  return rtn;
}


void EpetraOperatorWrapper::copyEpetraIntoThyra(const Epetra_MultiVector& x,
                                                Vector<double> thyraVec) const
{
  int numVecs = x.NumVectors();
  TEST_FOR_EXCEPTION(numVecs != 1, std::runtime_error,
                     "epetraToThyra does not work with MV dimension != 1");

  double* const epetraData = x[0];

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  bool verbose = false;
  if (verbose)
    {
      *out << "in ETO::copyEpetraIntoThyra()" << std::endl;
      double normIn=0;
      x.Norm2(&normIn);
      *out << "before: norm(Epetra vec) = " << normIn << std::endl;
    }

  int offset = 0;
  for (int b=0; b<thyraVec.numBlocks(); b++)
    {
      Vector<double> block = thyraVec.getBlock(b);
      TEST_FOR_EXCEPTION(!isSPMD(Thyra::space(block)), std::runtime_error, "block is not SPMD");
      int localOffset = lowestLocallyOwnedIndex(Thyra::space(block));
      int localNumElems = numLocalElements(Thyra::space(block));

      RCP<SpmdVectorBase<double> > spmd 
        = rcp_dynamic_cast<SpmdVectorBase<double> >(block.ptr());

      /** get a non-const pointer to the data in the thyra std::vector */ 
      {
        /* limit the scope so that it gets destroyed, and thus committed,
         * before using the std::vector */
        Teuchos::RCP<DetachedVectorView<double> > view 
          = rcp(new DetachedVectorView<double>(block.ptr(), 
                                               Thyra::Range1D(localOffset,
                                                       localOffset+localNumElems-1), 
                                               true));
        double* thyraData = view->values();
        for (int i=0; i<localNumElems; i++)
          {
            thyraData[i] = epetraData[i+offset];
          }
      }
      offset += localNumElems;
    }
  if (verbose)
    {
      *out << "after: norm(Thyra vec) = " << norm2(thyraVec) << std::endl;
    }
}


void EpetraOperatorWrapper::
copyThyraIntoEpetra(const ConstVector<double>& thyraVec,
                    Epetra_MultiVector& v) const 
{
  int numVecs = v.NumVectors();
  TEST_FOR_EXCEPTION(numVecs != 1, std::runtime_error,
                     "copyThyraIntoEpetra does not work with MV dimension != 1");

  /* get a writable pointer to the contents of the Epetra std::vector */
  double* epetraData = v[0];

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  bool verbose = false;
  if (verbose)
    {
      *out << "in ETO::copyThyraIntoEpetra()" << std::endl;
      *out << "before: norm(Thyra vec) = " << norm2(thyraVec) << std::endl;
    }

  int offset = 0;
  for (int b=0; b<thyraVec.numBlocks(); b++)
    {
      ConstVector<double> block = thyraVec.getBlock(b);
      TEST_FOR_EXCEPTION(!isSPMD(Thyra::space(block)), std::runtime_error, 
                         "block is not SPMD");
      int localOffset = lowestLocallyOwnedIndex(Thyra::space(block));
      int localNumElems = numLocalElements(Thyra::space(block));
      RCP<const SpmdVectorBase<double> > spmd 
        = rcp_dynamic_cast<const SpmdVectorBase<double> >(block.constPtr());

      /** get a const pointer to the data in the thyra std::vector */ 
      {
        /* limit the scope so that it gets destroyed, and thus committed,
         * before using the std::vector */
        Teuchos::RCP<const ConstDetachedVectorView<double> > view 
          = rcp(new ConstDetachedVectorView<double>(block.constPtr(), 
                                               Thyra::Range1D(localOffset,
                                                              localOffset+localNumElems-1), 
                                               true));
        const double* thyraData = view->values();
        for (int i=0; i<localNumElems; i++)
          {
            epetraData[i+offset] = thyraData[i];
          }
      }
      offset += localNumElems;
    }
  
  if (verbose)
    {
      double normOut=0;
      v.Norm2(&normOut);
      *out << "after: norm(Epetra vec) = " << normOut << std::endl;
    }
}


int EpetraOperatorWrapper::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  bool verbose = false;
  if (verbose)
    {
      *out << "in ETO::Apply()" << std::endl;
      double normX=0;
      double normY=0;
      X.Norm2(&normX);
      Y.Norm2(&normY);
      *out << "before: norm(in) = " << normX << std::endl;
      *out << "before: norm(out) = " << normY << std::endl;
    }

  if (!useTranspose_)
    {
      Vector<double> tX = domain_.createMember();
      copyEpetraIntoThyra(X, tX);
      Vector<double> tY = range_.createMember();
      zeroOut(tY);
      thyraOp_.apply(tX, tY);
      copyThyraIntoEpetra(tY, Y);
    }
  else
    {
      Vector<double> tX = range_.createMember();
      copyEpetraIntoThyra(X, tX);
      Vector<double> tY = domain_.createMember();
      zeroOut(tY);
      thyraOp_.applyTranspose(tX, tY);
      copyThyraIntoEpetra(tY, Y);
    }

  if (verbose)
    {
      double normX=0;
      double normY=0;
      X.Norm2(&normX);
      Y.Norm2(&normY);
      *out << "after: norm(in) = " << normX << std::endl;
      *out << "after: norm(out) = " << normY << std::endl;
      *out << "leaving ETO::Apply()" << std::endl;
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


  const SerialComm<Ordinal>* serialComm 
    = dynamic_cast<const SerialComm<Ordinal>*>(spmd->getComm().get());

#ifdef HAVE_MPI
  const MpiComm<Ordinal>* mpiComm 
    = dynamic_cast<const MpiComm<Ordinal>*>(spmd->getComm().get());

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


} // namespace Thyra


Teuchos::RCP<const Thyra::LinearOpBase<double> > 
Thyra::makeEpetraWrapper(const ConstLinearOperator<double>& thyraOp)
{
  RCP<const LinearOpBase<double> > rtn =
    epetraLinearOp(
      Teuchos::rcp(new EpetraOperatorWrapper(thyraOp))
      );
  return rtn;
}
