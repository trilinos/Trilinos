#include"Isorropia_EpetraMatcher.hpp"

int main(int argc, char** argv) {

	if(argc>1)
	{	
		int rc=0;
		int localProc = 0;
		
		#ifdef HAVE_EPETRAEXT
		/*#ifdef HAVE_MPI
			int numProcs;
		 	MPI_Init(&argc, &argv);
		  	MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
		  	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
		  	const Epetra_MpiComm Comm(MPI_COMM_WORLD);
		  	const Epetra_MpiComm Comm;
		#else*/
		  const Epetra_SerialComm Comm;
		//#endif

		  Epetra_CrsMatrix *matrixPtr;
		  rc = EpetraExt::MatrixMarketFileToCrsMatrix(argv[1], Comm, matrixPtr);
		  
		  if (rc < 0){
			 if (localProc==0){
				cout << "error reading input file" << std::endl << "FAIL" << std::endl;
			 }
			 exit(1);
		  }
		  else
		  		cout<<"Crs Matrix Created!!!...."<<endl;
		#else
		  fail = 0;
		  if (localProc == 0){
			 cout << "Test not run because it requires EPETRA_EXT" << std::endl;
		  }
		#endif
		
		Teuchos::ParameterList paramlist;
		paramlist.set("Matching Algorithm",argv[2]);
        Isorropia::Epetra::Isorropia_EpetraMatcher pm(matrixPtr,paramlist);
		
        //Teuchos::RCP<const Epetra_CrsMatrix> r(Teuchos::RCP<const
        //Epetra_CrsMatrix>(matrixPtr,true));
        //Isorropia::Epetra::Isorropia_EpetraMatcher pm(r,paramlist);
        
        pm.match();
	}
	else
		cout<<"Specify input file.."<<endl;
	
	return 0;
}
