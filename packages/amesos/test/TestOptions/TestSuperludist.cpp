//
//  Returns the number of failures.
//  Note:  If AMESOS_SUPERLUDIST is not supported, SubTest() will 
//  always return 0
// 
int SubTest( Epetra_Comm &Comm, bool transpose, 
	     AMESOS::Parameter::List ParamList, 
	     Epetra_CrsMatrix *& Amat, 
	     double maxrelresidual )
{
	
  int iam = Comm.MyPID() ; 
  int errors = 0 ; 
  bool verbose = false; 

  Epetra_Map *Map = Amat->RowMap() ; 
  
  Epetra_Vector x2(*Map), x1(*Map), x(*Map), b(*Map), residual(*Map), temp(*Map);

  //
  //  Solve Ax = b using Amesos_SUPERLUDIST via the Amesos_Factory interface
  //
  Epetra_LinearProblem Problem;
  Amesos_BaseSolver* Abase ; 
  Amesos_Factory Afactory;
  Abase = Afactory.Create( AMESOS_SUPERLUDIST, Problem, ParamList ) ; 
  if ( Abase != 0 ) {

    //
    //  Factor A
    //
    Problem.SetOperator( &Amat );
    EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
    EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

    b.Random();
    b.PutScalar(1.0);
    //
    //  Solve A x = b 
    //
    Problem.SetLHS( &x );
    Problem.SetRHS( &b );
    EPETRA_CHK_ERR( Abase->Solve(  ) ); 


    //  if (verbose) cout << " x = " << x << endl ; 
    //
    int ind[1];
    double val[1];
    ind[0] = 0;
    val[0] = 1 ; 
    if ( A.MyGRID( 0 ) )
      A.SumIntoMyValues( 0, 1, val, ind ) ; 

    //  if (verbose) cout << " A' = " << A << endl ; 
    EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

    //
    //  Solve A' x1 = x 
    //
    Problem.SetLHS( &x1 );
    Problem.SetRHS( &x );
    EPETRA_CHK_ERR( Abase->Solve(  ) ); 

    //  if (verbose) cout << " x1 = " << x1 << endl ; 

    if ( A.MyGRID( 0 ) )
      A.SumIntoMyValues( 0, 1, val, ind ) ; 

    //  if (verbose) cout << " A'' = " << A << endl ; 
    EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

    //
    //  Solve A" x2 = x1
    //
    Problem.SetLHS( &x2 );
    Problem.SetRHS( &x1 );
    EPETRA_CHK_ERR( Abase->Solve(  ) ); 

    //  if (verbose) cout << " x2 = " << x2 << endl ; 

    //
    //  Compute the residual: A A' A" x2 - b
    //

    A.Multiply( false, x2, temp ) ; //  temp = A x2

    //  if (verbose) cout << " temp = " << temp << endl ; 

    val[0] = -val[0] ; 
    if ( A.MyGRID( 0 ) )
      A.SumIntoMyValues( 0, 1, val, ind ) ; 
    A.Multiply( false, temp, x2 ) ; //  x2 = A' A" x2



    //  if (verbose) cout << " x2 = " << x2 << endl ; 


    if ( A.MyGRID( 0 ) )
      A.SumIntoMyValues( 0, 1, val, ind ) ; 
    A.Multiply( false, x2, temp ) ; //  temp = A A' A'' x2


    //  if (verbose) cout << " temp = " << temp << endl ; 
    //  if (verbose) cout << " b = " << b << endl ; 




    residual.Update( 1.0, temp, -1.0, b, 0.0 ) ;
    //  if (verbose) cout << " residual = " << residual << endl ; 

    double norm_residual ;
    residual.Norm2( &norm_residual ) ; 

    if (iam == 0 ) {
      if (verbose) cout << " norm2(A A' A'' x-b) = " << norm_residual << endl ; 
      //
      //  This is an ill-conditioned problem
      //
      if ( norm_residual < (1e-15)*(1.0*NumPoints*NumPoints*NumPoints*
				    NumPoints*NumPoints*NumPoints) ) {
	if (verbose) cout << " Test Passed " << endl ;
      } else {
	if (verbose) cout << " TEST FAILED " << endl ;
	errors += 1 ; 
      }
    }




    // #define FACTOR_B
#ifdef FACTOR_B
    //
    //  Now we check to make sure that we can change the problem and 
    //  re-factorize.  
    //






    const int BNumPoints = NumPoints;  // Must be between 2 and 100 (on large matrices,
    // the problem is quite ill-conditioned) 

    // Construct a Map that puts approximately the same number of 
    // equations on each processor.
    Epetra_Map BMap(BNumPoints, 0, Comm);

    //  Create an empty EpetraCrsMatrix 
    Epetra_CrsMatrix B(Copy, BMap, 0);

    //
    //  Populate A with a [-1,2,-1] tridiagonal matrix WITH -1 in the
    //  off diagonal corners.
    //  See CreateTridi.cpp in this directory 
    CreateTridiPlus( B ) ; 

    Epetra_Vector Bx2(BMap), Bx1(BMap), Bx(BMap), Bb(BMap), Bresidual(BMap), Btemp(BMap);

    //


    //
    //  Factor B
    //
    Problem.SetOperator( &B );
    EPETRA_CHK_ERR( Abase->SymbolicFactorization(  ) ); 
    EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

    Bb.Random();
    Bb.PutScalar(1.0);
    //
    //  Solve B x = b 
    //
    Problem.SetLHS( &Bx );
    Problem.SetRHS( &Bb );
    EPETRA_CHK_ERR( Abase->Solve(  ) ); 


    //  if (verbose) cout << " b = " << b << endl ; 
    //  if (verbose) cout << " x = " << x << endl ; 
    //
    ind[0] = 0;
    val[0] = 1 ; 
    if ( B.MyGRID( 0 ) )
      B.SumIntoMyValues( 0, 1, val, ind ) ; 

    //  if (verbose) cout << " B' = " << B << endl ; 
    EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

    //
    //  Solve B' x1 = x 
    //
    Problem.SetLHS( &Bx1 );
    Problem.SetRHS( &Bx );
    EPETRA_CHK_ERR( Abase->Solve(  ) ); 

    //  if (verbose) cout << " x1 = " << x1 << endl ; 

    if ( B.MyGRID( 0 ) )
      B.SumIntoMyValues( 0, 1, val, ind ) ; 

    //  if (verbose) cout << " B'' = " << B << endl ; 
    EPETRA_CHK_ERR( Abase->NumericFactorization(  ) ); 

    //
    //  Solve B" x2 = x1
    //
    Problem.SetLHS( &Bx2 );
    Problem.SetRHS( &Bx1 );
    EPETRA_CHK_ERR( Abase->Solve(  ) ); 

    //  if (verbose) cout << " x2 = " << x2 << endl ; 

    //
    //  Compute the residual: B B' B" x2 - b
    //

    B.Multiply( false, Bx2, Btemp ) ; //  temp = B x2

    //  if (verbose) cout << " temp = " << temp << endl ; 

    val[0] = -val[0] ; 
    if ( B.MyGRID( 0 ) )
      B.SumIntoMyValues( 0, 1, val, ind ) ; 
    B.Multiply( false, Btemp, Bx2 ) ; //  x2 = B' B" x2



    //  if (verbose) cout << " x2 = " << x2 << endl ; 


    if ( B.MyGRID( 0 ) )
      B.SumIntoMyValues( 0, 1, val, ind ) ; 
    B.Multiply( false, Bx2, Btemp ) ; //  temp = B B' B'' x2


    //  if (verbose) cout << " temp = " << temp << endl ; 
    //  if (verbose) cout << " b = " << b << endl ; 




    Bresidual.Update( 1.0, Btemp, -1.0, Bb, 0.0 ) ;
    //  if (verbose) cout << " residual = " << residual << endl ; 

    Bresidual.Norm2( &norm_residual ) ; 

    if (iam == 0 ) {
      if (verbose) cout << " norm2(B B' B'' x-b) = " << norm_residual << endl ; 
      //
      //  This is an ill-conditioned problem
      //
      if ( norm_residual < (1e-15)*(1.0*NumPoints*NumPoints*NumPoints*
				    NumPoints*NumPoints*NumPoints) ) {
	if (verbose) cout << " Test Passed " << endl ;
      } else {
	if (verbose) cout << " TEST FAILED " << endl ;
	errors += 1 ; 
      }
    }
#endif

    delete Abase;
  }

  return errors;
  
}


int TestSuperludist( EpetraCrsMatrix *& Amat, bool transpose, 
		     double &error, double &residual ) {

{
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", true );
   ParamList.setParameter( "AddZeroToDiag", true );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", true );
   SuperludistParams.setParameter( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, transpose, ParamList, Amat ) ; 
 }
  
#if 0
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", true );
   ParamList.setParameter( "AddZeroToDiag", true );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", false );
   SuperludistParams.setParameter( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }

 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", true );
   ParamList.setParameter( "AddZeroToDiag", true );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", true );
   SuperludistParams.setParameter( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", true );
   ParamList.setParameter( "AddZeroToDiag", true );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", false );
   SuperludistParams.setParameter( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", true );
   ParamList.setParameter( "AddZeroToDiag", true );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", false );
   SuperludistParams.setParameter( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", true );
   ParamList.setParameter( "AddZeroToDiag", false );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", true );
   SuperludistParams.setParameter( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", true );
   ParamList.setParameter( "AddZeroToDiag", false );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", true );
   SuperludistParams.setParameter( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", true );
   ParamList.setParameter( "AddZeroToDiag", false );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", false );
   SuperludistParams.setParameter( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", true );
   ParamList.setParameter( "AddZeroToDiag", false );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", false );
   SuperludistParams.setParameter( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  

 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", false );
   ParamList.setParameter( "AddZeroToDiag", true );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", true );
   SuperludistParams.setParameter( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", false );
   ParamList.setParameter( "AddZeroToDiag", true );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", true );
   SuperludistParams.setParameter( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", false );
   ParamList.setParameter( "AddZeroToDiag", true );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", false );
   SuperludistParams.setParameter( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", false );
   ParamList.setParameter( "AddZeroToDiag", true );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", false );
   SuperludistParams.setParameter( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", false );
   ParamList.setParameter( "AddZeroToDiag", false );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", true );
   SuperludistParams.setParameter( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", false );
   ParamList.setParameter( "AddZeroToDiag", false );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", true );
   SuperludistParams.setParameter( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", false );
   ParamList.setParameter( "AddZeroToDiag", false );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", false );
   SuperludistParams.setParameter( "MaxProcesses", 1 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
  
  
 {
   AMESOS::Parameter::List ParamList ;
   ParamList.setParameter( "Redistribute", false );
   ParamList.setParameter( "AddZeroToDiag", false );
   AMESOS::Parameter::List& SuperludistParams = ParamList.sublist("Superludist") ;
   SuperludistParams.setParameter( "ReuseSymbolic", false );
   SuperludistParams.setParameter( "MaxProcesses", 2 );
   //  ParamList.print( cerr, 10 ) ; 
   
   errors += SubTest( Comm, ParamList ) ; 
 }
 #endif
  
}
