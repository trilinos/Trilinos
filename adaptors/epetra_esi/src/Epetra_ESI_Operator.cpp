#ifndef _Epetra_ESI_Operator_cpp_
#define _Epetra_ESI_Operator_cpp_

#include "Epetra_ESI_Operator.h"
#include "Epetra_Import.h"

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
epetra_esi::Operator<Scalar,Ordinal>::
Operator(esi::Operator<Scalar, Ordinal>& esi_op)
   : epetra_esi::Object(),
     esi_op_(esi_op), esi_op_trans_(NULL),
     setupHasBeenCalled_(false), haveTranspose_(false),
     haveMatrixData_(false), haveRowReadAccess_(false), localOffset_(0)
{
   int err = addInterface("esi::Object", (void*)((esi::Object*)this) );
   if (err) msgAbort("epetra_esi::Operator ctor addInterface(esi::Object) error");

   err = addInterface("esi::Operator",
                      (void*)((esi::Operator<Scalar,Ordinal>*)this));
   if (err) msgAbort("epetra_esi::Operator ctor addInterface(esi::Operator) error");

   esi_op_trans_ = dynamic_cast<esi::OperatorTranspose<Scalar, Ordinal>*>(&esi_op_);
   if (esi_op_trans_ != NULL) {
     haveTranspose_ = true;
     err = addInterface("esi::OperatorTranspose",
                      (void*)((esi::OperatorTranspose<Scalar,Ordinal>*)this));
     if (err) msgAbort("epetra_esi::Operator ctor addInterface(esi::OperatorTranspose) error");
   }

   err = esi_op.getInterface("esi::MatrixData", (void*&)esi_matrix_data_);
   if (err == 0) {
     haveMatrixData_ = true;
     err = esi_matrix_data_->getGlobalSizes(numGlobalRows_, numGlobalCols_);
     err += esi_matrix_data_->getLocalSizes(numLocalRows_, numLocalCols_);
     if (err) msgAbort("epetra_esi::Operator ctor esi::MatrixData available, but can't answer size queries.");

     esi::IndexSpace<Ordinal> *ispace = NULL, *dummy = NULL;
     err = esi_matrix_data_->getIndexSpaces(ispace, dummy);
     if (err) msgAbort("epetra_esi::Operator ctor couldn't get IndexSpaces");

     err = ispace->getLocalPartitionOffset(localOffset_);
     if (err) msgAbort("epetra_esi::Operator ctor couldn't get local offset");
   }

   err = esi_op.getInterface("esi::MatrixRowReadAccess", (void*&)esi_row_read_);
   if (err == 0) {
     haveRowReadAccess_ = true;
   }

   err = addInterface("epetra_esi::Operator", (void*)this);
   if (err) msgAbort("epetra_esi::Operator ctor addInterface(epetra_esi::Operator)error");

   err = addInterface("Epetra_RowMatrix",
                      (void*)((Epetra_RowMatrix*)this));
   if (err) msgAbort("epetra_esi::Operator ctor addInterface(Epetra_RowMatrix) error");

#ifdef EPETRA_MPI
   MPI_Comm* mpicomm;
   err = esi_op.getRunTimeModel("MPI", (void*&)mpicomm);
   if (err == 0) {
     petra_comm_ = new Epetra_MpiComm(*mpicomm);
   }
   else {
     cout << __FILE__ <<", "<<__LINE__
      <<" epetra_esi::Operator constructor ERROR, unable to obtain an"
      << " MPI communicator from the input esi::Operator." << endl;
     petra_comm_ = new Epetra_SerialComm;
   }
#else
   petra_comm_ = new Epetra_SerialComm;
#endif

   petra_row_map_ =
      new Epetra_Map(numGlobalRows_, numLocalRows_, 0, *petra_comm_);
   petra_import_map_ =
      new Epetra_Map(numGlobalRows_, numGlobalRows_, 0, *petra_comm_);
   petra_import_ = new Epetra_Import(*petra_import_map_, *petra_row_map_);
   petra_domain_map_ =
	   petra_import_map_; // By default the map for the domain is the same as the import map?
   petra_range_map_ =
	   petra_row_map_;    // By default the map for the range is the same as the row map?

   if (petra_row_map_ == NULL || petra_import_map_ == NULL ||
       petra_import_ == NULL) {
     msgAbort("epetra_esi::Operator ctor failed to allocate map or import");
   }
}

#endif

