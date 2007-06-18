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
  TEST_FOR_EXCEPTION(true, runtime_error,
                     "EpetraOperatorWrapper::NormInf not implemated");
  return 1.0;
}


RefCountPtr<Epetra_Map> 
EpetraOperatorWrapper
::thyraVSToEpetraMap(const VectorSpace<double>& vs,
                     const RefCountPtr<Epetra_Comm>& comm) const 
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
          TEST_FOR_EXCEPTION(!isSPMD(vs.getBlock(b)), runtime_error, 
                             "EpetraOperatorWrapper requires vector space "
                             "blocks to be SPMD vectors");
          myLocalElements += numLocalElements(vs.getBlock(b));
        }
    }

  
  /* find the GIDs owned by this processor, taken from all blocks */
  Array<int> myGIDs(myLocalElements);
  
  int count=0;
  for (int b=0; b<vs.numBlocks(); b++)
    {
      int lowGIDInBlock = lowestLocallyOwnedIndex(vs.getBlock(b));
      int numLocalElementsInBlock = numLocalElements(vs.getBlock(b));
      for (int i=0; i<numLocalElementsInBlock; i++, count++)
        {
          myGIDs[count] = lowGIDInBlock+i;
        }
    }

  /* create the map */
  RefCountPtr<Epetra_Map> rtn 
    = rcp(new Epetra_Map(globalDim, myLocalElements, &(myGIDs[0]), 0, *comm));

  return rtn;
}


void EpetraOperatorWrapper::copyEpetraIntoThyra(const Epetra_MultiVector& x,
                                                Vector<double> thyraVec) const
{
  int numVecs = x.NumVectors();
  TEST_FOR_EXCEPTION(numVecs != 1, runtime_error,
                     "epetraToThyra does not work with MV dimension != 1");

  double* const epetraData = x[0];

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  bool verbose = false;
  if (verbose)
    {
      *out << "in ETO::copyEpetraIntoThyra()" << endl;
      double normIn=0;
      x.Norm2(&normIn);
      *out << "before: norm(Epetra vec) = " << normIn << endl;
    }

  int offset = 0;
  for (int b=0; b<thyraVec.numBlocks(); b++)
    {
      Vector<double> block = thyraVec.getBlock(b);
      TEST_FOR_EXCEPTION(!isSPMD(Thyra::space(block)), runtime_error, "block is not SPMD");
      int localOffset = lowestLocallyOwnedIndex(Thyra::space(block));
      int localNumElems = numLocalElements(Thyra::space(block));

      RefCountPtr<SpmdVectorBase<double> > spmd 
        = rcp_dynamic_cast<SpmdVectorBase<double> >(block.ptr());

      /** get a non-const pointer to the data in the thyra vector */ 
      {
        /* limit the scope so that it gets destroyed, and thus committed,
         * before using the vector */
        Teuchos::RefCountPtr<DetachedVectorView<double> > view 
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
      *out << "after: norm(Thyra vec) = " << norm2(thyraVec) << endl;
    }
}


void EpetraOperatorWrapper::
copyThyraIntoEpetra(const ConstVector<double>& thyraVec,
                    Epetra_MultiVector& v) const 
{
  int numVecs = v.NumVectors();
  TEST_FOR_EXCEPTION(numVecs != 1, runtime_error,
                     "copyThyraIntoEpetra does not work with MV dimension != 1");

  /* get a writable pointer to the contents of the Epetra vector */
  double* epetraData = v[0];

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  bool verbose = false;
  if (verbose)
    {
      *out << "in ETO::copyThyraIntoEpetra()" << endl;
      *out << "before: norm(Thyra vec) = " << norm2(thyraVec) << endl;
    }

  int offset = 0;
  for (int b=0; b<thyraVec.numBlocks(); b++)
    {
      ConstVector<double> block = thyraVec.getBlock(b);
      TEST_FOR_EXCEPTION(!isSPMD(Thyra::space(block)), runtime_error, 
                         "block is not SPMD");
      int localOffset = lowestLocallyOwnedIndex(Thyra::space(block));
      int localNumElems = numLocalElements(Thyra::space(block));
      RefCountPtr<const SpmdVectorBase<double> > spmd 
        = rcp_dynamic_cast<const SpmdVectorBase<double> >(block.constPtr());

      /** get a const pointer to the data in the thyra vector */ 
      {
        /* limit the scope so that it gets destroyed, and thus committed,
         * before using the vector */
        Teuchos::RefCountPtr<const ConstDetachedVectorView<double> > view 
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
      *out << "after: norm(Epetra vec) = " << normOut << endl;
    }
}


int EpetraOperatorWrapper::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  bool verbose = false;
  if (verbose)
    {
      *out << "in ETO::Apply()" << endl;
      double normX=0;
      double normY=0;
      X.Norm2(&normX);
      Y.Norm2(&normY);
      *out << "before: norm(in) = " << normX << endl;
      *out << "before: norm(out) = " << normY << endl;
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
      *out << "after: norm(in) = " << normX << endl;
      *out << "after: norm(out) = " << normY << endl;
      *out << "leaving ETO::Apply()" << endl;
    }
  return 0;
}


int EpetraOperatorWrapper::ApplyInverse(const Epetra_MultiVector& X, 
                                      Epetra_MultiVector& Y) const
{
  TEST_FOR_EXCEPTION(true, runtime_error,
                     "EpetraOperatorWrapper::ApplyInverse not implemented");
  return 1;
}


RefCountPtr<Epetra_Comm> 
EpetraOperatorWrapper::getEpetraComm(const ConstLinearOperator<double>& thyraOp) const
{
  RefCountPtr<Epetra_Comm> rtn;
  VectorSpace<double> vs = thyraOp.domain().getBlock(0);

  RefCountPtr<const SpmdVectorSpaceBase<double> > spmd 
    = rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(vs.constPtr());

  TEST_FOR_EXCEPTION(!isSPMD(vs), runtime_error, 
                     "EpetraOperatorWrapper requires vector space "
                     "blocks to be SPMD vector spaces");


  const SerialComm<int>* serialComm 
    = dynamic_cast<const SerialComm<int>*>(spmd->getComm().get());

#ifdef HAVE_MPI
  const MpiComm<int>* mpiComm 
    = dynamic_cast<const MpiComm<int>*>(spmd->getComm().get());

  TEST_FOR_EXCEPTION(mpiComm==0 && serialComm==0, runtime_error, 
                     "SPMD vector space has a communicator that is "
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
  TEST_FOR_EXCEPTION(serialComm==0, runtime_error, 
                     "SPMD vector space has a communicator that is "
                     "neither a serial comm nor an MPI comm");
  rtn = rcp(new Epetra_SerialComm());
  
#endif

  TEST_FOR_EXCEPTION(rtn.get()==0, runtime_error, "null communicator created");
  return rtn;
}


} // namespace Thyra


Teuchos::RefCountPtr<const Thyra::LinearOpBase<double> > 
Thyra::makeEpetraWrapper(const ConstLinearOperator<double>& thyraOp)
{
  RefCountPtr<const LinearOpBase<double> > rtn 
    = rcp(new EpetraLinearOp(rcp(new EpetraOperatorWrapper(thyraOp))));
  return rtn;
}
