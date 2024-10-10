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

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif
