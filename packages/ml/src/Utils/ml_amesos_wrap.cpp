/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_utils.h"
#include "ml_xyt.h"

#include "Epetra_Map.h" 
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h" 
#include "Epetra_LinearProblem.h" 
#ifdef HAVE_ML_AMESOS
#include "ml_amesos.h"
#include "ml_amesos_wrap.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_Factory.h" 
#include "AmesosClassType.h"

//  Jonathan - I need to convert an ml to an ML_Operator 
//  Did I pick off the right ML_Operator?

#ifdef EPETRA_MPI
#ifndef ML_MPI
   Garbage - ML_MPI and EPETRA_MPI must be the same 
#endif
#include "Epetra_MpiComm.h"
#else
#ifdef ML_MPI
   Garbage - ML_MPI and EPETRA_MPI must be the same 
#endif
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
#include "Amesos_Parameter_List.h"

   extern "C" { 
   double *SendThisToEMV ; 
   }

int ML_Amesos_Gen(ML *ml, int curr_level, int choice,
		  int MaxProcs, void **Amesos_Handle)
{
  double *global_nodes, *global_rows;
  int  allocated = 15, * colInd = NULL;
  double * colVal = NULL;
  int *global_rows_as_int, *global_nodes_as_int;
  int    N_nodes, node_offset, row_offset;
  int    ncnt;
  char str[80];
  ML_Comm *comm;
  int Nnodes_global, Nrows_global;
  int Nghost_nodes;
  int Nrows;
  ML_Operator *Ke = &(ml->Amat[curr_level]);
#ifdef ML_MPI
  MPI_Comm mpi_comm ;
#endif

  comm = Ke->comm;
#ifdef ML_MPI
  mpi_comm = comm->USR_comm; 
  Epetra_MpiComm EpetraComm( mpi_comm ) ; 
#else
  Epetra_SerialComm EpetraComm ; 
#endif  

  if (Ke->getrow->pre_comm == NULL) 
    Nghost_nodes = 0;
  else {
    if (Ke->getrow->pre_comm->total_rcv_length <= 0)
      ML_CommInfoOP_Compute_TotalRcvLength(Ke->getrow->pre_comm);
    Nghost_nodes = Ke->getrow->pre_comm->total_rcv_length;
  }

  int dummy;

  N_nodes = Ke->invec_leng;
  node_offset = ML_gpartialsum_int(N_nodes, comm);
  Nnodes_global = N_nodes;
  ML_gsum_scalar_int(&Nnodes_global, &dummy, comm);

  Nrows = Ke->outvec_leng;
  row_offset = ML_gpartialsum_int(Nrows, comm);
  Nrows_global = Nrows;
  ML_gsum_scalar_int(&Nrows_global, &dummy, comm);

  assert( N_nodes==Nrows );
  assert( Nnodes_global == Nrows_global ) ; 
  int num_global_rows;
  EpetraComm.SumAll( &N_nodes, &num_global_rows, 1 ) ;

  global_nodes  =(double *) ML_allocate(sizeof(double)*(N_nodes+Nghost_nodes));
  global_nodes_as_int  =(int *) ML_allocate(sizeof(int)*(N_nodes+Nghost_nodes));
  global_rows   =(double *) ML_allocate(sizeof(double)*(Nrows));
  global_rows_as_int   =(int *) ML_allocate(sizeof(int)*(Nrows));

  for (int i = 0 ; i < N_nodes; i++) global_nodes[i] = (double) (node_offset + i);
  for (int i = 0 ; i < Nrows; i++) {
    global_rows[i] = (double) (row_offset + i);
    global_rows_as_int[i] = row_offset + i;
  }
  for (int i = 0 ; i < Nghost_nodes; i++) global_nodes[i+N_nodes] = -1;

  Epetra_Map  EpetraMap( num_global_rows, Nrows, global_rows_as_int, 0, EpetraComm ) ; 

  Epetra_CrsMatrix *Amesos_CrsMatrix = new Epetra_CrsMatrix( Copy, EpetraMap, 0 ); 

  ML_exchange_bdry(global_nodes,Ke->getrow->pre_comm, 
 		 Ke->invec_leng,comm,ML_OVERWRITE,NULL);


  for ( int j = 0; j < N_nodes+Nghost_nodes; j++ ) { 
    global_nodes_as_int[j] = (int) global_nodes[j];
  }

  // MS // introduced variable allocation for colInd and colVal
  // MS // improved efficiency in InsertGlobalValues
  
  allocated = 1;
  colInd = new int[allocated];
  colVal = new double[allocated];
  int NumNonzeros;
  int ierr;
  for (int i = 0; i < Nrows; i++) {
    ierr = ML_Operator_Getrow(Ke,1,&i,allocated,colInd,colVal,&ncnt);
    if( ierr == 0 ) {
      while( ierr == 0 ) {
	delete [] colInd;
	delete [] colVal;
	allocated *= 2;
	colInd = new int[allocated];
	colVal = new double[allocated];
	ierr = ML_Operator_Getrow(Ke,1,&i,allocated,colInd,colVal,&ncnt);
      }
    }
    // MS // check out how many nonzeros we have
    // MS // NOTE: this may result in a non-symmetric patter for Amesos_CrsMatrix
    NumNonzeros = 0;
    for (int j = 0; j < ncnt; j++) {
      if (colVal[j] != 0.0) {
	colInd[NumNonzeros] = global_nodes_as_int[colInd[j]];
	colVal[NumNonzeros] = colVal[j];
	NumNonzeros++;
      }
    }
    Amesos_CrsMatrix->InsertGlobalValues( global_rows_as_int[i], NumNonzeros, 
					  colVal, colInd);
  }

  delete [] colInd;
  delete [] colVal;

  // MS // introduce support for Amesos_BaseFactory to
  // MS // allow different Amesos_Solvers
  
  assert(Amesos_CrsMatrix->TransformToLocal()==0);

  Epetra_LinearProblem *Amesos_LinearProblem = new Epetra_LinearProblem;
  Amesos_LinearProblem->SetOperator( Amesos_CrsMatrix ) ; 

  AMESOS::Parameter::List ParamList ;

  AMESOS::Parameter::List SluParamList=ParamList.sublist("Superludist");

  // this is specific to Superludist-2.0
  if( choice == ML_AMESOS_SUPERLUDIST ) {
    
    if( MaxProcs == -2 ) {
      if( Amesos_CrsMatrix->RowMatrixRowMap().LinearMap() == true ) {
	ParamList.setParameter("Redistribute",false);
      } else {
	if( EpetraComm.MyPID() == 0 ) {
	  cout << "*ML*WRN* in Amesos_Smoother, you can set MaxProcs = -1\n"
	       << "*ML*WRN* (that is, matrix will not be redistributed)\n"
	       << "*ML*WRN* ONLY if the matrix map is linear. Now proceeding\n"
	       << "*ML*WRN* with redistribution of the matrix\n"
	       << "*ML*WRN* (file " << __FILE__ << ", line "
	       << __LINE__ << ")\n";
	}
      }
    } else {
      ParamList.setParameter("Redistribute",true);
      SluParamList.setParameter("MaxProcesses",MaxProcs);
    }
  }
  
  Amesos_BaseSolver* A_Base;
  
  Amesos_Factory A_Factory;

  switch( choice ) {

  case ML_AMESOS_UMFPACK:
    if( EpetraComm.MyPID() == 0 )
      cout << "ML_Gen_Smoother_Amesos : building UMFPACK\n";
    A_Base = A_Factory.Create( AMESOS_UMFPACK, *Amesos_LinearProblem, ParamList );
    break;

  case ML_AMESOS_SUPERLUDIST:
    if( EpetraComm.MyPID() == 0 )
      cout << "ML_Gen_Smoother_Amesos : building SUPERLUDIST\n";
    A_Base = A_Factory.Create( AMESOS_SUPERLUDIST, *Amesos_LinearProblem, ParamList );
    break;

  case ML_AMESOS_KLU:
  default:
    if( EpetraComm.MyPID() == 0 )
      cout << "ML_Gen_Smoother_Amesos : building KLU\n";
    A_Base = A_Factory.Create( AMESOS_KLU, *Amesos_LinearProblem, ParamList );
    break;
  }
  
  A_Base->SymbolicFactorization();
  A_Base->NumericFactorization();

  ML_free(global_nodes_as_int);
  ML_free(global_rows_as_int);
  ML_free(global_rows);
  ML_free(global_nodes);

  *Amesos_Handle = (void *) A_Base ;

  return 0;
}

int ML_Amesos_Solve( void *Amesos_Handle, double x[], double rhs[] ) {

  Amesos_BaseSolver *A_Base = (Amesos_BaseSolver *) Amesos_Handle ;
  Epetra_LinearProblem *Amesos_LinearProblem = (Epetra_LinearProblem *) A_Base->GetProblem() ; 

  Epetra_BlockMap map = Amesos_LinearProblem->GetOperator()->OperatorDomainMap() ; 


  Epetra_Vector EV_rhs( Copy, map, rhs ) ;
  Epetra_Vector EV_lhs( View, map, x ) ;

  Amesos_LinearProblem->SetRHS( &EV_rhs ) ; 
  Amesos_LinearProblem->SetLHS( &EV_lhs ) ;

  A_Base->Solve() ; 

}

void ML_Amesos_Destroy(void *Amesos_Handle){

  Amesos_BaseSolver *A_Base = (Amesos_BaseSolver *) Amesos_Handle ;
  const Epetra_LinearProblem *Amesos_LinearProblem;
  Amesos_LinearProblem = A_Base->GetProblem() ; 
  const Epetra_Operator *EO = Amesos_LinearProblem->GetOperator() ; 

  delete Amesos_LinearProblem->GetOperator() ; 

  delete Amesos_LinearProblem ;
  delete A_Base ;
}

#else

int ciao=0;

#endif
