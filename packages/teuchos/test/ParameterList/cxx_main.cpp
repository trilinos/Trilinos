#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

using namespace std;
using namespace Teuchos;

void print_break() { cout << "---------------------------------------------------" << endl; }

int main(int argc, char *argv[])
{
  bool verbose = false;
  bool debug = false;
  bool InvalidCmdLineArgs = false;
  int FailedTests = 0;

  // Read options from the command line. 
  CommandLineProcessor  clp(false); // Don't throw exceptions
  clp.setOption( "v", "q", &verbose, "Set if output is printed or not." );
  CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
  if( parse_return != CommandLineProcessor::PARSE_SUCCESSFULL ) return parse_return;

  // Create Main Parameter List / Sublist Structure
  ParameterList PL_Main;
  ParameterList& PL_Direction = PL_Main.sublist("Direction");
  ParameterList& PL_Newton = PL_Direction.sublist("Newton");
  ParameterList& PL_LinSol = PL_Newton.sublist("Linear Solver");
  ParameterList& PL_LineSearch = PL_Main.sublist("Line Search");

  // Check Parameter List Structure
  if (verbose) {
	print_break();
	cout << "Empty Parameter List Structure" << endl;
	print_break();
  	cout<<PL_Main<< endl;
  }
  if (verbose) cout << "Is 'Direction' recognized as a sublist of 'Main' ... ";
  if ( PL_Main.isParameterSublist( "Direction" ) ) {  
	if (verbose) cout << "yes"<< endl;
  } else {
	if (verbose) cout << "no"<< endl;
	FailedTests++;		
  }
  if (verbose) cout << "Is 'Newton' recognized as a sublist of 'Direction' ... ";
  if ( PL_Direction.isParameterSublist( "Newton" ) ) {  
	if (verbose) cout << "yes"<< endl;
  } else {
	if (verbose) cout << "no"<< endl;
	FailedTests++;		
  }
  if (verbose) cout << "Is 'Linear Solver' recognized as a sublist of 'Newton' ... ";
  if ( PL_Newton.isParameterSublist( "Linear Solver" ) ) {  
	if (verbose) cout << "yes"<< endl;
  } else {
	if (verbose) cout << "no"<< endl;
	FailedTests++;		
  }
  if (verbose) cout << "Is 'Line Search' recognized as a sublist of 'Main' ... ";
  if ( PL_Main.isParameterSublist( "Line Search" ) ) {  
	if (verbose) cout << "yes"<< endl;
  } else {
	if (verbose) cout << "no"<< endl;
	FailedTests++;		
  }

  // Fill in Direction Sublist
  PL_Direction.getParameter("Method", "Newton");
  PL_LinSol.setParameter("Tol",1e-5);
  double tol = PL_LinSol.getParameter("Tolerance",1e-10);
  bool RBNS = PL_Newton.getParameter("Rescue Bad Newton Solve", true );

  // Print out Direction Sublist
  if (verbose) {
	print_break();
	cout << "Direction Parameter List" << endl;
	print_break();
  	PL_Direction.print(cout);
  }
  if (verbose) cout << "Is 'Newton' recognized as a parameter of 'Direction' ... ";
  if ( PL_Direction.isParameter( "Newton" ) ) {  
	if (verbose) cout << "yes"<< endl;
  } else {
	if (verbose) cout << "no"<< endl;
	FailedTests++;		
  }
  if (verbose) cout << "Is 'Tolerance' recognized as a parameter of 'Newton' ... ";
  if ( PL_Newton.isParameter( "Tolerance" ) ) {  
	if (verbose) cout << "yes (should be no)"<< endl;
	FailedTests++;
  } else {
	if (verbose) cout << "no (as expected)"<< endl;
  }
  if (verbose) cout << "Is 'Tolerance' recognized as a parameter of 'Linear Solver' ... ";
  if ( PL_LinSol.isParameter( "Tolerance" ) ) {  
	if (verbose) cout << "yes"<< endl;
  } else {
	if (verbose) cout << "no"<< endl;
	FailedTests++;		
  }
  if (verbose) cout << "Is 'Rescue Bad Newton Solve' recognized as a parameter of 'Newton' ... ";
  if ( PL_Newton.isParameter( "Rescue Bad Newton Solve" ) ) {  
	if (verbose) cout << "yes"<< endl;
  } else {
	if (verbose) cout << "no"<< endl;
	FailedTests++;		
  }

  // Line Search Sublist (if there are no failures, this will be constructed and added)
  if (!FailedTests) {
    ParameterList PL_My_LineSearch;
    string ls_method = PL_My_LineSearch.getParameter("Method", "Polynomial");
    ParameterList& PL_Polynomial = PL_My_LineSearch.sublist("Polynomial");
    int ARI = PL_Polynomial.getParameter("Allowed Relative Increase", 100 );
    double alpha_factor = PL_Polynomial.getParameter("Alpha Factor", 0.0001 );
    int default_step = PL_Polynomial.getParameter("Default Step", 1 );
    bool force_interp = PL_Polynomial.getParameter("Force Interpolation", false );
    string interp_type = PL_Polynomial.getParameter("Interpolation Type", "Cubic" );
    double max_bnds_factor = PL_Polynomial.getParameter("Max Bounds Factor", 0.5 );
    PL_Polynomial.setParameter("Max Iters", 3 );
    int max_iter_inc = PL_Polynomial.getParameter("Maximum Iteration for Increase", 0 );
    double min_bnds_factor = PL_Polynomial.getParameter("Min Bounds Factor", 0.1 );
    int rec_step = PL_Polynomial.getParameter("Recovery Step", 1 );
    string rec_step_type = PL_Polynomial.getParameter("Recovery Step Type", "Constant");
    string suff_dec_cond = PL_Polynomial.getParameter("Sufficient Decrease Condition", "Armijo-Goldstein" );
    bool use_cntrs = PL_Polynomial.getParameter("Use Counters", true );

    PL_Main.setParameter("Nonlinear Solver", "Line Search Based"); 

    // Set the Line Search Parameter List equal to the one just constructed
    PL_LineSearch = PL_My_LineSearch;
    ParameterList& PL_My_Polynomial = PL_LineSearch.sublist("Polynomial");
    ParameterList Copied_PL_Main( PL_Main );
    if (verbose) cout<< "Is 'operator=' functional ... ";
    if ( PL_My_Polynomial.isParameter("Recovery Step Type") ) {
	if (verbose) cout<< "yes" << endl;
    } else {
	if (verbose) cout<< "no" << endl;
	FailedTests++;
    }  
    if (verbose) cout<< "Is the copy constructor functional ... ";
    if ( Copied_PL_Main.isParameter("Nonlinear Solver") ) {
	if (verbose) cout<< "yes" << endl;
    } else {
	if (verbose) cout<< "no" << endl;
	FailedTests++;
    }  
    
    // Retrieve some information from the parameter list using templated "get" method.
    int max_iters;
    string nonlin_solver;
    bool tempMeth = true;
    try {
      max_iters = PL_My_Polynomial.template getParameter<int>("Max Iters");
      nonlin_solver = PL_Main.template getParameter<string>("Nonlinear Solver");
    }
    catch( exception& e ) { tempMeth = false; }  
    if (verbose) {
	cout<< "Is the templated 'getParameter' method functional ... "<<endl;
	cout<< "  Can we retrieve information using the CORRECT variable type ... ";
	if (tempMeth && max_iters==3) { cout << "yes" << endl; }
	else { cout << "no" << endl; FailedTests++; }
    }

    // Retrieve some information from the parameter list that we know is a bad "get".
    float mbf;
    tempMeth = false;
    FailedTests++;  // Increment it prematurely, it will get decremented if the test passes.
    try {
	mbf = PL_LinSol.template getParameter<float>( "Tol" );
    }
    catch( exception& e ) {
	tempMeth = true;
	FailedTests--;		
    }
    if (verbose) {
	cout<< "  Can we retrieve information using the WRONG variable type ... ";
	if (tempMeth) { cout << "no" << endl; }
	else { cout << "yes" << endl; }
    }

    // Check templated isParameterType functionality
    bool PT1, PT2, PT3, PT4, PT5;
    PT1 = PL_Polynomial.template isParameterType<int>("Default Step");
    PT2 = PL_Polynomial.template isParameterType<long int>("Default Step");
    PT3 = PL_Polynomial.template isParameterType<string>("Interpolation Type");
    if (verbose) {
	cout<< "Is the templated 'isParameterType' method functional ... "<<endl;
	cout<< "  Is the 'Default Step' of type 'int' ... ";
	if (PT1) { cout<< "yes" << endl; }
	else { cout<< "no" << endl; FailedTests++; }
	cout<< "  Is the 'Default Step' of type 'long int' ... ";
	if (PT2) { cout<< "yes" << endl; FailedTests++; }
	else { cout<< "no (as expected)" << endl; }
	cout<< "  Is the 'Interpolation Type' of type 'string' ... ";
	if (PT3) { cout<< "yes" << endl; }
	else { cout<< "no" << endl; FailedTests++; }
    }
    PT4 = isParameterType<double>(PL_Polynomial, "Max Bounds Factor");
    PT5 = isParameterType<float>(PL_Polynomial, "Max Bounds Factor");    
    if (verbose) {
	cout<< "Is the helper function 'isParameterType' functional ... "<<endl;
	cout<< "  Is the 'Max Bounds Factor' of type 'double' ... ";
	if (PT4) { cout<< "yes" <<endl; }
	else { cout<< "no" << endl; FailedTests++; }
	cout<< "  Is the 'Max Bounds Factor' of type 'float' ... ";
	if (PT5) { cout<< "yes" <<endl; FailedTests++; }
	else { cout<< "no (as expected)" << endl; }
    }	
  }

    // Print out main list
  if (verbose) {
	print_break();
	cout << "The Final Parameter List" << endl;
	print_break();
	PL_Main.print(cout);
	print_break();
	cout << "Number of Failed Tests : " << FailedTests << endl;
	print_break();
  }

  // Return -1 if there are any failed tests, 
  // else 0 will be returned indicating a clean finish!  
  if ( FailedTests > 0 ) { return(-1); }

  return 0;
}

