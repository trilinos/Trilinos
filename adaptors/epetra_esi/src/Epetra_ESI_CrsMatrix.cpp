#ifndef _Epetra_ESI_CrsMatrix_cpp_
#define _Epetra_ESI_CrsMatrix_cpp_

#include "Epetra_ESI_CrsMatrix.h"

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
epetra_esi::CrsMatrix<Scalar,Ordinal>::
CrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph& graph)
   : epetra_esi::Object(),
     setupHasBeenCalled_(false),
     ispace_(NULL),
     comm_(NULL),
     rowMap_(NULL),
     whichConstructor_(0)
{
   epetra_crsmatrix_ = new Epetra_CrsMatrix(CV, graph);
   int err = addInterface("esi::Object", (void*)((esi::Object*)this) );
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::Object) error");

#ifdef EPETRA_MPI
   err = -1;
   const Epetra_MpiComm* empicomm =
          dynamic_cast<const Epetra_MpiComm*>(&(graph.Comm()));
   if (empicomm != NULL) {
     mpicomm_ = empicomm->Comm();
     mpicommPtr_ = &mpicomm_;
     err = setRunTimeModel("MPI", (void*)(&mpicomm_));
   }
   if (err) msgAbort("epetra_esi::CrsMatrix ctor setRunTimeModel error");
#endif

   err = addInterface("esi::Operator",
                      (void*)((esi::Operator<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::Operator) error");

   err = addInterface("esi::OperatorTranspose",
                      (void*)((esi::OperatorTranspose<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::OperatorTranspose) error");

   err = addInterface("esi::MatrixData",
                      (void*)((esi::MatrixData<Ordinal>*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::MatrixData) error");

   err = addInterface("esi::MatrixRowWriteAccess",
                      (void*)((esi::MatrixRowWriteAccess<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::MatrixRowWriteAccess) error");

   err = addInterface("esi::MatrixRowReadAccess",
                      (void*)((esi::MatrixRowReadAccess<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::MatrixRowReadAccess) error");

   err = addInterface("epetra_esi::CrsMatrix", (void*)this);
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(epetra_esi::CrsMatrix) error");

   err = addInterface("Epetra_CrsMatrix",
                      (void*)((Epetra_CrsMatrix*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(Epetra_CrsMatrix) error");
}

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
epetra_esi::CrsMatrix<Scalar,Ordinal>::
CrsMatrix(Epetra_DataAccess CV,
          const epetra_esi::IndexSpace<Ordinal>& indexspace,
          Ordinal estimatedNumEntriesPerRow)
   : epetra_esi::Object(),
     setupHasBeenCalled_(false),
     ispace_(const_cast<epetra_esi::IndexSpace<Ordinal>*>(&indexspace)),
     comm_(NULL),
     rowMap_(NULL),
     whichConstructor_(1)
{
   epetra_crsmatrix_ = new Epetra_CrsMatrix(CV, indexspace, estimatedNumEntriesPerRow);
   int err = addInterface("esi::Object", (void*)((esi::Object*)this) );
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::Object) error");

  ispace_->addReference();
#ifdef EPETRA_MPI
   err = ispace_->getRunTimeModel("MPI", (void*&)mpicommPtr_);
   if (err) msgAbort("epetra_esi::CrsMatrix ctor setRunTimeModel error");
   err = setRunTimeModel("MPI", (void*)(mpicommPtr_));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor setRunTimeModel error");
#endif

   err = addInterface("esi::Operator",
                      (void*)((esi::Operator<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::Operator) error");

   err = addInterface("esi::OperatorTranspose",
                      (void*)((esi::OperatorTranspose<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::OperatorTranspose) error");

   err = addInterface("esi::MatrixData",
                      (void*)((esi::MatrixData<Ordinal>*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::MatrixData) error");

   err = addInterface("esi::MatrixRowWriteAccess",
                      (void*)((esi::MatrixRowWriteAccess<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::MatrixRowWriteAccess) error");

   err = addInterface("esi::MatrixRowReadAccess",
                      (void*)((esi::MatrixRowReadAccess<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(esi::MatrixRowReadAccess) error");

   err = addInterface("epetra_esi::CrsMatrix", (void*)this);
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(epetra_esi::CrsMatrix) error");

   err = addInterface("Epetra_CrsMatrix",
                      (void*)((Epetra_CrsMatrix*)this));
   if (err) msgAbort("epetra_esi::CrsMatrix ctor addInterface(Epetra_CrsMatrix) error");
}

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
epetra_esi::CrsMatrix<Scalar,Ordinal>::
~CrsMatrix()
{
  if (ispace_ != NULL) ispace_->deleteReference();
  delete comm_;
  delete rowMap_;
  delete epetra_crsmatrix_;
}

#endif

