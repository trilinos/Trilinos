//
//  PartialFactorization checks to make sure that we clean up after our mess
//  before everything is done.
//  

int PartialFactorization( const char* AmesosClass,
			  const Epetra_Comm &Comm, 
			  bool transpose, 
			  bool verbose, 
			  Teuchos::ParameterList ParamList, 
			  Epetra_CrsMatrix *& Amat, 
			  double Rcond ) ;
