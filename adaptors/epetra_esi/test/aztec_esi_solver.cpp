#include <iostream.h>

//
//Epetra_ESI.h includes all of the Epetra_ESI_* headers, which in turn
//include the Epetra headers that they depend on.
//
#include "Epetra_ESI.h"

//========= utility function prototypes ========

int create_ESI_IndexSpace(int numGlobal, int numLocal, Epetra_Comm& comm, 
		   esi::IndexSpace<int>*& esi_indexspace);

int create_ESI_Vector(esi::IndexSpace<int>* indexspace, esi::Vector<double,int>*& vec);

int create_ESI_Operator(int numGlobal, int numLocal,
			esi::IndexSpace<int>* indexspace,
			Epetra_CrsGraph*& graph,
                        esi::Operator<double,int>*& esiop);

int create_ESI_Solver(esi::Operator<double,int>* A,
		      esi::Solver<double,int>*& solver);

int create_ESI_Solver_Epetra(esi::Operator<double,int>* A,
		      esi::Solver<double,int>*& solver);

int cout_ESI_RowMatrix(esi::MatrixRowReadAccess<double,int>* esimat);

int cout_ESI_Vector(esi::Vector<double,int>* esivec);

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

   //This is a simple program to test the aztecoo_esi::Solver, which is an
   //implementation of the abstract esi::Solver interface.
   //
   //First we'll set up a trivially simple decomposition of an equation space
   //with 50 equations per processor.
   //
   int numProcs = comm->NumProc();
   int localProc = comm->MyPID();
   int numLocal = 50;
   int numGlobal = numProcs*numLocal;


   //Now create a indexspace to hold the description of our decomposition. The run-time
   //type of the esi::IndexSpace will be epetra_esi::IndexSpace, but our main program
   //scope doesn't need to be aware of that.
   //
   esi::IndexSpace<int>* esiIndexSpace = NULL;
   CHK_ERR( create_ESI_IndexSpace(numGlobal, numLocal, *comm, esiIndexSpace) );

   MPI_Comm* mpicomm;
   CHK_ERR( esiIndexSpace->getRunTimeModel("MPI", (void*&)mpicomm));


   //Now create an esi::Operator. The run-time type of the esi::Operator will
   //be epetra_esi::CrsMatrix. We need a Epetra_CrsGraph to construct our
   //matrix with. The 'create_ESI_Operator' function will not only define the
   //structure of our sparse matrix, it will also fill it with some coefficient
   //values.
   //
   Epetra_CrsGraph* graph = NULL;
   esi::Operator<double,int>* A = NULL;
   CHK_ERR( create_ESI_Operator(numGlobal, numLocal, esiIndexSpace, graph, A) );

   CHK_ERR( A->getRunTimeModel("MPI", (void*&)mpicomm));

   //now let's create a vector x and then clone it to create a vector b
   esi::Vector<double,int> *x = NULL, *b = NULL;

   CHK_ERR( create_ESI_Vector(esiIndexSpace, x) );

   CHK_ERR( x->clone(b) );

   //fill x with 1's
   CHK_ERR( x->put(1.0) );

   //Lets fill our right-hand-side b with A*x.
   CHK_ERR( A->apply(*x, *b) );

   //now re-zero x.
   CHK_ERR( x->put(0.0) );

   //now let's create an esi::Solver.
   //
   esi::Solver<double,int>* solver = NULL;

   CHK_ERR( create_ESI_Solver_Epetra(A, solver) );


   //now we're going to try out the parameters method. We happen to know that
   //the esi::Solver is actually an aztecoo_esi::Solver, so we'll give it some
   //Aztec-specific parameters.
   //
   int numParams = 8;
   char** paramStrings = new char*[numParams];
   for(int i=0; i<numParams; i++) paramStrings[i] = new char[128];
   sprintf(paramStrings[0], "AZ_solver AZ_gmres");
   sprintf(paramStrings[1], "AZ_tol 1.e-7");
   sprintf(paramStrings[2], "AZ_max_iter 200");
   sprintf(paramStrings[3], "AZ_output AZ_all");
//   sprintf(paramStrings[4], "AZ_precond AZ_Neumann");
   sprintf(paramStrings[4], "AZ_precond AZ_dom_decomp");
   sprintf(paramStrings[5], "AZ_subdomain_solver AZ_ilut");
   sprintf(paramStrings[6], "AZ_ilut_fill 0.8");
   sprintf(paramStrings[7], "AZ_drop 0.01");

   CHK_ERR( solver->parameters(numParams, paramStrings) );


   //Now we're finally ready to solve Ax=b.
   CHK_ERR( solver->solve((esi::Vector<double,int>&)(*b),
                          (esi::Vector<double,int>&)(*x)) );

   //Now let's check the solution. We happen to know that it should be a
   //vector of ones, so let's just see whether the 1-norm of the solution is
   //equal to the global size of the problem.
   double nrm1 = 0.0;
   CHK_ERR( x->norm1(nrm1) );

   int returnValue = 0;
   if (fabs(nrm1 - 1.0*numGlobal) > 1.e-5) {
     //obviously this is not a general check. If AZ_tol above is set to a
     //large value, then the difference between nrm1 and numGlobal could
     //easily be greater than 1.e-5
     if (localProc == 0) {
       cerr << "Solution vector doesn't seem to be correct." << endl;
       cerr << "1-norm of soln: " << nrm1 << ". It should be "<<numGlobal<<endl;
     }
     returnValue = -1;
   }

   for(int j=0; j<numParams; j++) delete [] paramStrings[j];
   delete [] paramStrings;

   delete solver;
   delete x;
   delete b;
   delete A;
   delete graph;
   delete esiIndexSpace;
   delete comm;

#ifdef EPETRA_MPI
   MPI_Finalize();
#endif

   return(returnValue);
}

//----------------------------------------------
int create_ESI_IndexSpace(int numGlobal, int numLocal,
                   Epetra_Comm& comm, esi::IndexSpace<int>*& esi_indexspace)
{
   //we're using indexBase = 0...
   esi_indexspace = new epetra_esi::IndexSpace<int>(numGlobal, numLocal, 0, comm);

   //return error-code -1 if the allocation failed
   if (esi_indexspace == NULL) return(-1);

   return(0);
}

//----------------------------------------------
int create_ESI_Vector(esi::IndexSpace<int>* indexspace,
                      esi::Vector<double,int>*& vec)
{

  //In order to construct a epetra_esi::Vector, we need an actual
  //epetra_esi::IndexSpace...
  //
  epetra_esi::IndexSpace<int>* petraindexspace = NULL;
  CHK_ERR( indexspace->getInterface("epetra_esi::IndexSpace", (void*&)petraindexspace) );

  vec = new epetra_esi::Vector<double,int>(*petraindexspace);
  if (vec == NULL) return(-1);
  return(0);
}

//----------------------------------------------
int create_ESI_Operator(int numGlobal, int numLocal,
			esi::IndexSpace<int>* indexspace, Epetra_CrsGraph*& graph,
                        esi::Operator<double,int>*& esiop)
{
  //The goal of this function is to create an esi::Operator. But since the
  //run-time type of the esi::Operator will be epetra_esi::CrsMatrix, we first
  //have to create and fill a Epetra_CrsGraph (the container that defines the
  //structure of the matrix.
  //
  int i;

  int localOffset = 0;
  CHK_ERR( indexspace->getLocalPartitionOffset(localOffset) );

  //This will be a very simple structure.
  //We want our matrix to have a diagonal, and the first off-diagonals both
  //above and below the main diagonal.
  //
  int* rowLengths = new int[numLocal];
  if (rowLengths == NULL) return(-1);

  for(i=0; i<numLocal; i++) rowLengths[i] = 3;
  if (localOffset == 0) rowLengths[0] = 2;
  if (localOffset+numLocal == numGlobal) rowLengths[numLocal-1] = 2;

  //We need a native Epetra_Map to construct the Epetra_CrsGraph with.
  Epetra_Map* petramap = NULL;
  CHK_ERR( indexspace->getInterface("Epetra_Map", (void*&)petramap) );

  graph = new Epetra_CrsGraph(Copy, *petramap, rowLengths);

  delete[] rowLengths;

  //Now we've allocated the 'shell' of the graph, so let's proceed with
  //setting the column-indices...

  int colIndices[3];

  for(i=0; i<numLocal; i++) {
    int numCols = 3;

    if (i+localOffset == 0) {
      numCols = 2;
      colIndices[0] = 0; colIndices[1] = 1;
    }
    else if (i+localOffset == numGlobal-1) {
      numCols = 2;
      colIndices[0] = i+localOffset-1; colIndices[1] = i+localOffset;
    }
    else {
      for(int j=0; j<numCols; j++) colIndices[j] = i+localOffset-1 + j;
    }

    int err = graph->InsertGlobalIndices(i+localOffset, numCols, colIndices);
    if (err != 0) return(err);
  }

  int err = graph->TransformToLocal();
  if (err != 0) return(-1);

  //Now we've got a fully initialized graph holding the structure of our
  //matrix, so we're ready to go ahead and construct the matrix.
  //
  epetra_esi::CrsMatrix<double,int>* petramat = 
    new epetra_esi::CrsMatrix<double,int>(Copy, *graph);
  if (petramat == NULL) return(-1);

  //Now let's run through and supply coefficients to go with the structure
  //we built above.
  //
  double coefs[3];

  for(i=0; i<numLocal; i++) {
      int numCols = 0;
      CHK_ERR( petramat->getRowAllocatedLength(i+localOffset, numCols) );

      if (i == 0 && numCols == 2) {
	//If this is global row 0, then we use coefs -1 and 2
        colIndices[0] = 0; colIndices[1] = 1;
        coefs[0] = 2.0; coefs[1] = -1.0;
      }
      else if (numCols == 2) {
        //If this is the last global row, then we use coefs 2 and -1
        colIndices[0] = i+localOffset-1; colIndices[1] = i+localOffset;
        coefs[0] = -1.0; coefs[1] = 2.0;
      }
      else {
	//for most rows, we use three coefs: -1, 2, -1
        for(int j=0; j<numCols; j++) colIndices[j] = i+localOffset-1 + j;
        coefs[0] = -1.0; coefs[1] = 2.0; coefs[2] = -1.0;
      }

      CHK_ERR( petramat->copyIntoRow(localOffset+i, coefs,
				   colIndices, numCols) );
   }

  CHK_ERR( petramat->getInterface("esi::Operator", (void*&)esiop) );

  //Finally, we finish off by calling esi::Operator::setup, which is where
  //Epetra_CrsMatrix::TransformToLocal gets called.
  //
  CHK_ERR( esiop->setup() );

  return(0);
}


//----------------------------------------------
int create_ESI_Solver(esi::Operator<double,int>* A,
		      esi::Solver<double,int>*& solver)
{
  solver = new aztecoo_esi::Solver<double,int>(A);
  if (solver == NULL) return(-1);

  return(0);
}

//----------------------------------------------
int create_ESI_Solver_Epetra(esi::Operator<double,int>* A,
                      esi::Solver<double,int>*& solver)
{
  //The run-time type of the esi::Solver will be aztecoo_esi::Solver. We
  //Need a epetra_esi::CrsMatrix to construct an aztecoo_esi::Solver.
  //
  epetra_esi::CrsMatrix<double,int>* petraA = NULL;
  CHK_ERR( A->getInterface("epetra_esi::CrsMatrix", (void*&)petraA) );

  solver = new aztecoo_esi::Solver<double,int>(petraA);
  if (solver == NULL) return(-1);

  return(0);
}

//----------------------------------------------
int cout_ESI_RowMatrix(esi::MatrixRowReadAccess<double,int>* esimat)
{
   esi::IndexSpace<int>* rowIS = NULL, *colIS = NULL;
   CHK_ERR(esimat->getIndexSpaces(rowIS, colIS) );

   //we need to know how many local equations there are, and what the first
   //local equation is. (Since we happen to know that the index-spaces were
   //built with an indexBase of 0,
   // then firstLocalEqn == 'getLocalPartitionOffset'.)
   int localSize, firstLocalEqn, localRank;
   CHK_ERR( rowIS->getLocalSize(localSize) );
   CHK_ERR( rowIS->getLocalPartitionOffset(firstLocalEqn) );
   CHK_ERR( rowIS->getLocalPartitionRank(localRank) );

   int numRows, numCols;
   CHK_ERR( esimat->getGlobalSizes(numRows, numCols) );
   cout << localRank << ": global rows " << numRows << ", global cols "
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
      cout << localRank << ": row " << row << ": ";
      for(int j=0; j<rowLen; j++) {
         cout << "("<<colIndices[j]<<","<<coefs[j]<<") ";
      }
      cout << endl;
   }

   return(0);
}

//----------------------------------------------
int cout_ESI_Vector(esi::Vector<double,int>* esivec)
{
   esi::IndexSpace<int>* indexspace = NULL;
   CHK_ERR(esivec->getIndexSpace(indexspace) );

   //we need to know how many local equations there are, and what the first
   //local equation is. (Since we happen to know that the index-space was built
   //with an indexBase of 0, then firstLocalEqn == 'getLocalPartitionOffset'.)
   int localSize, firstLocalEqn, localRank;
   CHK_ERR( indexspace->getLocalSize(localSize) );
   CHK_ERR( indexspace->getLocalPartitionOffset(firstLocalEqn) );
   CHK_ERR( indexspace->getLocalPartitionRank(localRank) );

   double* coefs;
   CHK_ERR( esivec->getCoefPtrReadLock(coefs) );

   for(int i=0; i<localSize; i++) {
      cout << localRank << ": "<<firstLocalEqn+i<<", "<<coefs[i]<<endl;
   }

   CHK_ERR( esivec->releaseCoefPtrLock(coefs) );

   return(0);
}

