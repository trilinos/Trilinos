#ifndef _Epetra_ESI_IndexSpace_cpp_
#define _Epetra_ESI_IndexSpace_cpp_

//------------------------------------------------------------------------------
template<class Ordinal>
epetra_esi::IndexSpace<Ordinal>::
IndexSpace(Ordinal globalSize, Ordinal localSize, Ordinal indexBase,
                      const Epetra_Comm& comm)
    : epetra_esi::Object(),
      Epetra_Map(globalSize, localSize, indexBase, comm)
{
   int err = addInterface("esi::Object", (void*)((esi::Object*)this) );
   if (err) msgAbort("epetra_esi::IndexSpace ctor addInterface(esi::Object) error");

   err = addInterface("esi::IndexSpace", (void*)((esi::IndexSpace<Ordinal>*)this) );
   if (err) msgAbort("epetra_esi::IndexSpace ctor addInterface(esi::IndexSpace) error");

   err = addInterface("Epetra_Map", (void*)((Epetra_Map*)this) );
   if (err) msgAbort("epetra_esi::IndexSpace ctor addInterface(Epetra_Map) error");

   err = addInterface("epetra_esi::IndexSpace", (void*)((epetra_esi::IndexSpace<Ordinal>*)this) );
   if (err) msgAbort("epetra_esi::IndexSpace ctor addInterface(Epetra_Map) error");

#ifdef EPETRA_MPI
   err = -1;
   const Epetra_MpiComm* mpicomm = dynamic_cast<const Epetra_MpiComm*>(&comm);
   if (mpicomm != NULL) {
     comm_ = mpicomm->Comm();
     err = setRunTimeModel("MPI", (void*)(&comm_));
   }
   if (err) msgAbort("epetra_esi::IndexSpace ctor setRunTimeModel(MPI) error");
#endif
}

//------------------------------------------------------------------------------
template<class Ordinal>
epetra_esi::IndexSpace<Ordinal>::
~IndexSpace()
{
}

//------------------------------------------------------------------------------
template<class Ordinal>
esi::ErrorCode epetra_esi::IndexSpace<Ordinal>::
getGlobalSize(Ordinal& globalSize)
{
  globalSize = NumGlobalPoints();
  return(0);
}

#endif

