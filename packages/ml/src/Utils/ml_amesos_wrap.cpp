/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_utils.h"
#include "ml_xyt.h"

#include "ml_amesos_wrap.h"
#include "Epetra_Map.h" 
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h" 
#include "Epetra_LinearProblem.h" 
#include "Amesos_Factory.h" 


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

int ML_Amesos_Gen(ML *ml, int curr_level, void **Amesos_Handle)
{
  double *global_nodes, *global_rows, colVal[15];
  int *global_rows_as_int, *global_nodes_as_int;
  int    N_nodes, node_offset, row_offset;
  int colInd[15], ncnt;
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

  //
  //  Ken - Need to improve the efficiency here - Create arrays and make
  //  a single call to InsertGlobalValues 
  //
  for (int i = 0; i < Nrows; i++) {
    int j = ML_Operator_Getrow(Ke,1,&i,15,colInd,colVal,&ncnt);
    //    Amesos_CrsMatrix.InsertGlobalValues( global_rows_as_int[i], ncnt, global_nodes_as_int, 
    for (int j = 0; j < ncnt; j++) {
      if (colVal[j] != 0.0) {
	Amesos_CrsMatrix->InsertGlobalValues( global_rows_as_int[i], 1, 
					      &(colVal[j]),
					      &(global_nodes_as_int[colInd[j]]) );
					    
      }
    }
  }

  
  assert(Amesos_CrsMatrix->TransformToLocal()==0);

  Epetra_LinearProblem *Amesos_LinearProblem = new Epetra_LinearProblem;
  Amesos_LinearProblem->SetOperator( Amesos_CrsMatrix ) ; 

  AMESOS::Parameter::List ParamList ;

  Amesos_BaseSolver* A_Base;
  
  Amesos_Factory A_Factory;
  A_Base = A_Factory.Create( AMESOS_KLU, *Amesos_LinearProblem, ParamList );
  
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

  cout << "In the desctructor" << endl;
   
  
  delete Amesos_LinearProblem->GetOperator() ; 

  delete Amesos_LinearProblem ;
  delete A_Base ;
}

