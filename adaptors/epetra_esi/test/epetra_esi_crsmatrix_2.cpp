#include <iostream.h>

//==============================================================================
//This program demonstrates using the epetra_esi::CrsMatrix, which implements
//the esi::Operator, esi::MatrixData, esi::MatrixRowReadAccess and
//esi::MatrixRowWriteAccess interfaces.
//
//In particular, this demonstration uses the epetra_esi::CrsMatrix constructor
//which does not require a fully initialized Epetra_CrsGraph object containing
//the structure of the matrix. Instead, it requires only a esi::IndexSpace
//object describing the distribution of rows across processors.
//
//Note that this program uses the 'copyIntoRow' method for putting data into
//the matrix. If the 'sumIntoRow' method is used, the specified column-indices
//must already have been specified via a call to 'copyIntoRow'. i.e., the
//'sumIntoRow' method may not be used to provide coefficients for "new"
//column positions.
//==============================================================================
//
//Epetra_ESI.h includes all of the Epetra_ESI_* headers, which in turn
//include the Epetra headers that they depend on.
//
#include "Epetra_ESI.h"

//========= utility function prototypes ========

int create_ESI_IndexSpace(int numGlobal, int numLocal,
                   Epetra_Comm& comm, esi::IndexSpace<int>*& esi_ispc);

int create_ESI_Vector(esi::IndexSpace<int>& ispc,
                      esi::Vector<double,int>*& vec);

int create_ESI_Operator(esi::IndexSpace<int>& indexspace,
                        esi::Operator<double,int>*& esiop);

int fill_ESI_RowMatrix(esi::MatrixRowWriteAccess<double,int>* esimat);

int out_ESIMatrix(esi::MatrixRowReadAccess<double,int>* esimat,
                       ostream& out);

int out_ESI_Vector(esi::Vector<double,int>* esivec,
                   ostream& out);

//========= main program =======================

int main(int argc, char** argv) {
   Epetra_Comm* comm = NULL;

#ifdef EPETRA_MPI
   if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      cerr << "error in MPI_Init." << endl; return(-1);
   }
   comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
   comm = new Epetra_SerialComm();
#endif

   if (comm == NULL) {
      cerr << "error allocating Epetra_Comm." << endl; return(-1);
   }

   int numProcs = comm->NumProc();
   int localProc = comm->MyPID();
   int numLocal = 5;
   int numGlobal = numProcs*numLocal;
   int localOffset = numLocal*localProc;

   esi::IndexSpace<int>* esiIndexSpace = NULL;
   CHK_ERR( create_ESI_IndexSpace(numGlobal, numLocal, *comm, esiIndexSpace) );

   //now let's create an esi::Operator.
   esi::Operator<double,int>* esiOp = NULL;
   CHK_ERR( create_ESI_Operator(*esiIndexSpace, esiOp) );

   //now let's get a write-access interface from the operator.
   esi::MatrixRowWriteAccess<double,int>* esiwrt = NULL;
   CHK_ERR( esiOp->getInterface("esi::MatrixRowWriteAccess", (void*&)esiwrt));

   //pass it to our fill function where we'll exercise the copyIntoRow functions.
   CHK_ERR( fill_ESI_RowMatrix(esiwrt) );

   //now let's get an esi MatrixRowReadAccess interface from the operator.
   esi::MatrixRowReadAccess<double,int>* esiread = NULL;
   CHK_ERR( esiOp->getInterface("esi::MatrixRowReadAccess", (void*&)esiread));

   //now let's output the matrix to an ostream.
   CHK_ERR( out_ESIMatrix(esiread, cout) );

   CHK_ERR( esiOp->setup() );

   //now let's create a vector x
   esi::Vector<double,int>* x = NULL;
   CHK_ERR( create_ESI_Vector(*esiIndexSpace, x) );

   //and a vector y which is a clone of x...
   esi::Vector<double,int>* y = NULL;
   CHK_ERR( x->clone(y) );

   //fill x with 1's
   CHK_ERR( x->put(1.0) );

   //for the matrix-vector product y = esiOp * x
   CHK_ERR( esiOp->apply(*x, *y) );

   //now let's output the vector y to an ostream.
   CHK_ERR( out_ESI_Vector(y, cout) );

   delete x;
   delete y;
   delete esiOp;
   delete esiIndexSpace;
   delete comm;

#ifdef EPETRA_MPI
   MPI_Finalize();
#endif

   return(0);
}

//----------------------------------------------
int create_ESI_IndexSpace(int numGlobal, int numLocal,
                   Epetra_Comm& comm, esi::IndexSpace<int>*& esi_ispc)
{
  esi_ispc = new epetra_esi::IndexSpace<int>(numGlobal, numLocal, 0, comm);
  if (esi_ispc == NULL) return(-1);

  return(0);
}

//----------------------------------------------
int create_ESI_Vector(esi::IndexSpace<int>& ispc,
                      esi::Vector<double,int>*& vec)
{
//This function insists that the incoming IndexSpace have run-time type
//epetra_esi::IndexSpace.
  epetra_esi::IndexSpace<int>* pispc =
       dynamic_cast<epetra_esi::IndexSpace<int>*>(&ispc);
  if (pispc == NULL) {
    cerr << "create_ESI_Vector ERROR dynamic_cast failed" << endl;
    EPETRA_ESI_ERR_BEHAVIOR(-1);
  }

  vec = new epetra_esi::Vector<double,int>(*pispc);
  if (vec == NULL) return(-1);
  return(0);
}

//----------------------------------------------
int create_ESI_Operator(esi::IndexSpace<int>& indexspace,
                        esi::Operator<double,int>*& esiop)
{
   epetra_esi::IndexSpace<int>* petra_ispace=
      dynamic_cast<epetra_esi::IndexSpace<int>*>(&indexspace);
   if (petra_ispace == NULL) {
    cerr << "create_ESI_Operator ERROR dynamic_cast failed" << endl;
    EPETRA_ESI_ERR_BEHAVIOR(-1);
   }

   esiop = new epetra_esi::CrsMatrix<double,int>(Copy, *petra_ispace, 0);
   if (esiop == NULL) return(-1);
   return(0);
}

//----------------------------------------------
int fill_ESI_RowMatrix(esi::MatrixRowWriteAccess<double,int>* esimat)
{
   esi::IndexSpace<int>* rowIS = NULL, *colIS = NULL;
   CHK_ERR(esimat->getIndexSpaces(rowIS, colIS) );

   //we need to know how many local equations there are, and what the first
   //local equation is.
   int globalSize, localSize, localOffset;
   CHK_ERR( rowIS->getGlobalSize(globalSize) );
   CHK_ERR( rowIS->getLocalSize(localSize) );
   CHK_ERR( rowIS->getLocalPartitionOffset(localOffset) );

   int colIndices[3];
   double coefs[3];

   for(int i=0; i<localSize; i++) {
      int numCols = 3;

      if (i+localOffset == 0) {
        numCols = 2;
        colIndices[0] = 0; colIndices[1] = 1;
        coefs[0] = 2.0; coefs[1] = -1.0;
      }
      else if (i+localOffset == globalSize-1) {
        numCols = 2;
        colIndices[0] = i+localOffset-1; colIndices[1] = i+localOffset;
        coefs[0] = -1.0; coefs[1] = 2.0;
      }
      else {
        for(int j=0; j<numCols; j++) colIndices[j] = i+localOffset-1 + j;
        coefs[0] = -1.0; coefs[1] = 2.0; coefs[2] = -1.0;
      }

      CHK_ERR( esimat->copyIntoRow(localOffset+i, coefs,
                                 colIndices, numCols) );
   }

   return(0);
}

//----------------------------------------------
int out_ESIMatrix(esi::MatrixRowReadAccess<double,int>* esimat,
  ostream& out)
{
   esi::IndexSpace<int>* rowIS = NULL, *colIS = NULL;
   CHK_ERR(esimat->getIndexSpaces(rowIS, colIS) );

   //we need to know how many local equations there are, and what the first
   //local equation is. (Since we happen to know that the index-spaces were
   //built with an indexBase of 0,
   //then firstLocalEqn == 'getLocalPartitionOffset'.)
   int localSize, firstLocalEqn, localRank;
   CHK_ERR( rowIS->getLocalSize(localSize) );
   CHK_ERR( rowIS->getLocalPartitionOffset(firstLocalEqn) );
   CHK_ERR( rowIS->getLocalPartitionRank(localRank) );

   int numRows, numCols;
   CHK_ERR( esimat->getGlobalSizes(numRows, numCols) );
   out << localRank << ": global rows " << numRows << ", global cols "
     << numCols << endl;

   Epetra_Array<int> colIndices;
   Epetra_Array<double> coefs;

   for(int i=0; i<localSize; i++) {
      //first, make sure our colIndices and coefs arrays are the right length.
      //(This operation is a no-op if the array is already the right size or
      //bigger.)
      int rowLen;
      int row = firstLocalEqn+i;
      CHK_ERR( esimat->getRowNonzeros(row, rowLen) );

      CHK_ERR( colIndices.resize(rowLen) );
      CHK_ERR( coefs.resize(rowLen) );

      CHK_ERR( esimat->copyOutRow(row, coefs.dataPtr(), colIndices.dataPtr(),
                                  coefs.length(), rowLen) );
      out << localRank << ": row " << row << ": ";
      for(int j=0; j<rowLen; j++) {
         out << "("<<colIndices[j]<<","<<coefs[j]<<") ";
      }
      out << endl;
   }

   return(0);
}

//----------------------------------------------
int out_ESI_Vector(esi::Vector<double,int>* esivec, ostream& out)
{
   esi::IndexSpace<int>* indexSpace = NULL;
   CHK_ERR(esivec->getIndexSpace(indexSpace) );

   //we need to know how many local equations there are, and what the first
   //local equation is. (Since we happen to know that the index-space was built
   //with an indexBase of 0, then firstLocalEqn == 'getLocalPartitionOffset'.)
   int localSize, firstLocalEqn, localRank;
   CHK_ERR( indexSpace->getLocalSize(localSize) );
   CHK_ERR( indexSpace->getLocalPartitionOffset(firstLocalEqn) );
   CHK_ERR( indexSpace->getLocalPartitionRank(localRank) );

   double* coefs;
   CHK_ERR( esivec->getCoefPtrReadLock(coefs) );

   for(int i=0; i<localSize; i++) {
      out << localRank << ": "<<firstLocalEqn+i<<", "<<coefs[i]<<endl;
   }

   CHK_ERR( esivec->releaseCoefPtrLock(coefs) );

   return(0);
}

