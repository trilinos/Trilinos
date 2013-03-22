// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#endif

#ifdef HAVE_MPI
#include "mpi.h"
#endif

using Teuchos::CommandLineProcessor;
using Teuchos::ParameterList;
using Teuchos::getParameter;
typedef ParameterList::PrintOptions PLPrintOptions;
using Teuchos::ParameterEntry;
using Teuchos::OSTab;
using Teuchos::rcp;
using Teuchos::inoutArg;

void print_break() { std::cout << "---------------------------------------------------" << std::endl; }
double Plus ( double a, double b ) { return a+b; }

int main( int argc, char *argv[] )
{

  using std::cout;
  using std::endl;

  bool verbose = true;
  int FailedTests = 0;
  bool result;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const int procRank = Teuchos::GlobalMPISession::getRank();

  bool success = true;

  try {

  // Read options from the command line. 
  CommandLineProcessor  clp(false); // Don't throw exceptions
  clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
  CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
  if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
    cout << "Processor "<< procRank <<", parse_return "<< parse_return << std::endl;
    cout << "End Result: TEST FAILED" << std::endl;
    return parse_return;
  }

  // Only print on 0 processor
  if (procRank != 0 && verbose)
    verbose = false;

  if (verbose)
    cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  //-----------------------------------------------------------
  // Create Main Parameter List / Sublist Structure
  //-----------------------------------------------------------

  ParameterList PL_Main("PL_Main");
  const std::string Direction_Doc = "This sublist controls how direction is computed.";
  ParameterList& PL_Direction = PL_Main.sublist("Direction",false,Direction_Doc);
  ParameterList& PL_Newton = PL_Direction.sublist("Newton");
  ParameterList& PL_LinSol = PL_Newton.sublist("Linear Solver");
  ParameterList& PL_LineSearch = PL_Main.sublist("Line Search");

  //-----------------------------------------------------------
  // Check Parameter List Structure
  //-----------------------------------------------------------
  if (verbose) {
    print_break();
    cout << "Empty Parameter List Structure" << std::endl;
    print_break();
    cout<<PL_Main<< std::endl;
  }
  if (verbose) cout << "Is 'Direction' recognized as a sublist of 'Main' ... ";
  if ( PL_Main.isSublist( "Direction" ) ) {  
    if (verbose) cout << "yes"<< std::endl;
  } else {
    if (verbose) cout << "no"<< std::endl;
    FailedTests++;        
  }
  if (verbose) cout << "Is 'Newton' recognized as a sublist of 'Direction' ... ";
  if ( PL_Direction.isSublist( "Newton" ) ) {  
    if (verbose) cout << "yes"<< std::endl;
  } else {
    if (verbose) cout << "no"<< std::endl;
    FailedTests++;        
  }
  if (verbose) cout << "Is 'Linear Solver' recognized as a sublist of 'Newton' ... ";
  if ( PL_Newton.isSublist( "Linear Solver" ) ) {  
    if (verbose) cout << "yes"<< std::endl;
  } else {
    if (verbose) cout << "no"<< std::endl;
    FailedTests++;        
  }
  if (verbose) cout << "Is 'Line Search' recognized as a sublist of 'Main' ... ";
  if ( PL_Main.isSublist( "Line Search" ) ) {  
    if (verbose) cout << "yes"<< std::endl;
  } else {
    if (verbose) cout << "no"<< std::endl;
    FailedTests++;        
  }

  if (verbose) cout << "Is subist documentation std::string maintained ...\n";
  {
    Teuchos::OSTab tab(cout);
    const Teuchos::ParameterEntry
      *paramEntry = PL_Main.getEntryPtr("Direction");
    TEUCHOS_TEST_FOR_EXCEPT(0==paramEntry);
    const std::string extracted_Direction_Doc = paramEntry->docString();
    if (verbose) tab.o() << "Expected doc std::string = \"" << Direction_Doc << "\"\n";
    if (verbose) tab.o() << "Extracted doc std::string = \"" << extracted_Direction_Doc << "\"\n";
    if (extracted_Direction_Doc == Direction_Doc) {
      if (verbose) tab.o() << "passed!  They match :-)\n";
    }
    else {
      if (verbose) tab.o() << "failed!  They do not match :-("<< std::endl;
      FailedTests++;        
    }
  }
  

  //-----------------------------------------------------------
  // Fill in Direction Sublist
  //-----------------------------------------------------------

  double tol = 0.0;
  bool RBNS = false;
  PL_Direction.get("Method", "Newton");
  PL_LinSol.set("Tol",1e-5);
  tol = PL_LinSol.get("Tolerance",1e-10);
  (void)tol; // Not used, bad test!
  RBNS = PL_Newton.get("Rescue Bad Newton Solve", true );
  (void)RBNS; // Not used, bad test!

  //-----------------------------------------------------------
  // Print out Direction Sublist
  //-----------------------------------------------------------
  if (verbose) {
    print_break();
    cout << "Direction Parameter List" << std::endl;
    print_break();
    PL_Direction.print(cout);
  }
  if (verbose) cout << "Is 'Newton' recognized as a parameter of 'Direction' ... ";
  if ( PL_Direction.isParameter( "Newton" ) ) {  
    if (verbose) cout << "yes"<< std::endl;
  } else {
    if (verbose) cout << "no"<< std::endl;
    FailedTests++;        
  }
  if (verbose) cout << "Is 'Tolerance' recognized as a parameter of 'Newton' ... ";
  if ( PL_Newton.isParameter( "Tolerance" ) ) {  
    if (verbose) cout << "yes (should be no)"<< std::endl;
    FailedTests++;
  } else {
    if (verbose) cout << "no (as expected)"<< std::endl;
  }
  if (verbose) cout << "Is 'Tolerance' recognized as a parameter of 'Linear Solver' ... ";
  if ( PL_LinSol.isParameter( "Tolerance" ) ) {  
    if (verbose) cout << "yes"<< std::endl;
  } else {
    if (verbose) cout << "no"<< std::endl;
    FailedTests++;        
  }
  if (verbose) cout << "Is 'Rescue Bad Newton Solve' recognized as a parameter of 'Newton' ... ";
  if ( PL_Newton.isParameter( "Rescue Bad Newton Solve" ) ) {  
    if (verbose) cout << "yes"<< std::endl;
  } else {
    if (verbose) cout << "no"<< std::endl;
    FailedTests++;
  }

  //-----------------------------------------------------------
  // Line Search Sublist 
  // (if there are no failures, this will be constructed and added)
  //-----------------------------------------------------------
  if (!FailedTests) {
    int ARI = 0, default_step = 0, max_iter_inc = 0, rec_step = 0;
    double alpha_factor = 0.0, min_bnds_factor = 0.0, max_bnds_factor = 0.0;
    bool force_interp = true, use_cntrs = false;
    std::string ls_method = "Polynomial";
    // This is to force a char* to be passed into the get method to see if the template
    // specialization for char* is working.
    char* ls_method_char = const_cast<char *>(ls_method.c_str());
    ParameterList PL_My_LineSearch("PL_My_LineSearch");
    ls_method = PL_My_LineSearch.get("Method", ls_method_char);
    ParameterList& PL_Polynomial = PL_My_LineSearch.sublist("Polynomial");
    ARI = PL_Polynomial.get("Allowed Relative Increase", 100 );
    (void)ARI; // Not used, bad test!
    alpha_factor = PL_Polynomial.get("Alpha Factor", 0.0001 );
    (void)alpha_factor; // Not used, bad test!
    default_step = PL_Polynomial.get("Default Step", 1 );
    (void)default_step; // Not used, bad test!
    force_interp = PL_Polynomial.get("Force Interpolation", false );
    (void)force_interp; // Not used, bad test!
    std::string interp_type = PL_Polynomial.get("Interpolation Type", "Cubic" );
    max_bnds_factor = PL_Polynomial.get("Max Bounds Factor", 0.5 );
    (void)max_bnds_factor; // Not used, bad test!
    PL_Polynomial.set("Max Iters", 3 );
    max_iter_inc = PL_Polynomial.get("Maximum Iteration for Increase", 0 );
    (void)max_iter_inc; // Not used, bad test!
    min_bnds_factor = PL_Polynomial.get("Min Bounds Factor", 0.1 );
    (void)min_bnds_factor; // Not used, bad test!
    rec_step = PL_Polynomial.get("Recovery Step", 1 );
    (void)rec_step; // Not used, bad test!
    std::string rec_step_type = PL_Polynomial.get("Recovery Step Type", "Constant");
    std::string suff_dec_cond = PL_Polynomial.get("Sufficient Decrease Condition", "Armijo-Goldstein" );
    use_cntrs = PL_Polynomial.get("Use Counters", true );
    (void)use_cntrs; // Not used, bad test!
    PL_Main.set("Nonlinear Solver", "Line Search Based");

    //-----------------------------------------------------------
    // Set the Line Search Parameter List equal to the one just constructed
    //-----------------------------------------------------------
    PL_LineSearch.setParameters(PL_My_LineSearch);
    ParameterList& PL_My_Polynomial = PL_LineSearch.sublist("Polynomial");
    if (verbose) cout<< "Is 'operator=' functional ... ";
    if ( PL_My_Polynomial.isParameter("Recovery Step Type") ) {
      if (verbose) cout<< "yes" << std::endl;
    } else {
      if (verbose) cout<< "no" << std::endl;
      FailedTests++;
    }

    //-----------------------------------------------------------
    // Set Copying of parameter sublists and names
    //-----------------------------------------------------------

    if (verbose) {
      print_break();
      if (verbose) cout << "Test copying of sublist\n";
      print_break();
      PL_Direction.print(cout);
    }
    {
      const ParameterList
        &linearSolverPL = PL_Main.sublist("Direction").sublist("Newton").sublist("Line Search");
      const ParameterList
        linearSolverPL_copy(linearSolverPL);
      if (verbose) cout << "linearSolverPL.name() = " << linearSolverPL.name() << endl;
      if (verbose) cout << "linearSolverPL_copy.name() = " << linearSolverPL_copy.name() << endl;
      if (verbose) cout << "linearSolverPL_copy == linearSolverPL.name() : ";
      if (linearSolverPL_copy == linearSolverPL.name()) {
        if (verbose) cout << "passed" << endl;
      }
      else {
        if (verbose) cout << "failed" << endl;
        FailedTests++;
      }
    }

    if (verbose) {
      print_break();
      if (verbose) cout << "General tests\n";
      print_break();
      PL_Direction.print(cout);
    }

    ParameterList Copied_PL_Main(PL_Main);

    if (verbose) cout << "Copied_PL_Main.name() == PL_Main.name() : ";
    if (Copied_PL_Main.name() == PL_Main.name()) {
      if (verbose) cout << "passed" << endl;
    }
    else {
      if (verbose) cout << "failed" << endl;
      FailedTests++;
      if (verbose) cout << "Copyed_PL_Main.name() = " << Copied_PL_Main.name() << endl;
    }

    if (verbose) cout<< "Is the copy constructor functional ... ";
    if ( Copied_PL_Main.isParameter("Nonlinear Solver") ) {
      if (verbose) cout<< "yes" << std::endl;
    } else {
      if (verbose) cout<< "no" << std::endl;
      FailedTests++;
    }  

    bool tempMeth = true;

    //-----------------------------------------------------------
    // Check the templated 'get' method.
    //-----------------------------------------------------------
    //
    //-----------------------------------------------------------
    // Retrieve some information from the parameter list using templated "get" method.
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------
    int max_iters = 0, max_iters_again = 0;
    std::string nonlin_solver;
    tempMeth = true;
    try {
      max_iters = PL_My_Polynomial.INVALID_TEMPLATE_QUALIFIER get<int>("Max Iters");
      max_iters_again = Teuchos::getConst(PL_My_Polynomial).INVALID_TEMPLATE_QUALIFIER get<int>("Max Iters");
      nonlin_solver = PL_Main.INVALID_TEMPLATE_QUALIFIER get<std::string>("Nonlinear Solver");
    }
    catch( const Teuchos::Exceptions::InvalidParameter&) { tempMeth = false; }  
    if (verbose) {
      cout<< "Is the templated 'get' method functional ... "<<std::endl;
      cout<< "  Can we retrieve information using the CORRECT variable type ... ";
    }
    if (tempMeth && max_iters==3) { if (verbose) cout << "yes" << std::endl; }
    else { if (verbose) cout << "no" << std::endl; FailedTests++; }
    if (verbose) {
      cout<< "  Can we retrieve const information using the CORRECT variable type ... ";
    }
    if (tempMeth && max_iters_again==3) { if (verbose) cout << "yes" << std::endl; }
    else { if (verbose) cout << "no" << std::endl; FailedTests++; }

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list that we know is a bad "get".
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------
    float mbf = 0.0;
    tempMeth = false;
    try {
      mbf = PL_LinSol.INVALID_TEMPLATE_QUALIFIER get<float>( "Tol" );
      (void)mbf; // Not used, bad test!
      FailedTests++;
    }
    catch( const Teuchos::Exceptions::InvalidParameter&) {
      tempMeth = true;
    }
    if (verbose) {
      cout<< "  Can we retrieve information using the WRONG variable type ... ";
    }
    if (tempMeth) { if (verbose) cout << "no" << std::endl; }
    else { if (verbose) cout << "yes" << std::endl; }

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list using templated "get" method.
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------
    tempMeth = true;
    try {
      max_iters = PL_My_Polynomial.INVALID_TEMPLATE_QUALIFIER get<int>("Max Iters");
      nonlin_solver = PL_Main.INVALID_TEMPLATE_QUALIFIER get<std::string>("Nonlinear Solver");
    }
    catch( const Teuchos::Exceptions::InvalidParameter&) { tempMeth = false; }  
    if (verbose) {
      cout<< "Is the templated 'get' method functional ... "<<std::endl;
      cout<< "  Can we retrieve information using the CORRECT variable type ... ";
    }
    if (tempMeth && max_iters==3) { if (verbose) cout << "yes" << std::endl; }
    else { if (verbose) cout << "no" << std::endl; FailedTests++; }

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list that we know is a bad "get".
    //-----------------------------------------------------------
    tempMeth = false;
    try {
      mbf = PL_LinSol.INVALID_TEMPLATE_QUALIFIER get<float>( "Tol" );
      (void)mbf; // Not used, bad test!
      FailedTests++;
    }
    catch( const Teuchos::Exceptions::InvalidParameter&) {
      tempMeth = true;
    }
    if (verbose) {
      cout<< "  Can we retrieve information using the WRONG variable type ... ";
    }
    if (tempMeth) { if (verbose) cout << "no" << std::endl; }
    else { if (verbose) cout << "yes" << std::endl; }

    //-----------------------------------------------------------
    // Check the templated 'getPtr' method.
    //-----------------------------------------------------------

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list using templated "get" method.
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------
    int *max_iters_ptr = 0;
    const int *max_iters_ptr_again = 0;
    std::string* nonlin_solver_ptr;

    max_iters_ptr = PL_My_Polynomial.INVALID_TEMPLATE_QUALIFIER getPtr<int>("Max Iters");
    max_iters_ptr_again = Teuchos::getConst(PL_My_Polynomial).INVALID_TEMPLATE_QUALIFIER getPtr<int>("Max Iters");
    nonlin_solver_ptr = PL_Main.INVALID_TEMPLATE_QUALIFIER getPtr<std::string>("Nonlinear Solver");

    if (verbose) {
      cout<< "Is the templated 'getPtr' method functional ... "<<std::endl;
      cout<< "  Can we retrieve information using the CORRECT variable type ... ";
    }
    if (max_iters_ptr) {
      if ((*max_iters_ptr)==3) {
        if (verbose) cout << "yes" << std::endl; 
      }
      else { if (verbose) cout << "no" << std::endl; FailedTests++; }
    }
    if (verbose) {
      cout<< "  Can we retrieve const information using the CORRECT variable type ... ";
    }
    if (max_iters_ptr_again) {
      if ((*max_iters_ptr_again)==3) {
        if (verbose) cout << "yes" << std::endl; 
      }
      else { if (verbose) cout << "no" << std::endl; FailedTests++; }
    }

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list that we know is a bad "get".
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------
    float* mbf_ptr = 0;

    mbf_ptr = PL_LinSol.INVALID_TEMPLATE_QUALIFIER getPtr<float>( "Tol" );

    if (mbf_ptr)
      ++FailedTests;        

    if (verbose) {
      cout<< "  Can we retrieve information using the WRONG variable type ... ";
    }
    if (!mbf_ptr) { if (verbose) cout << "no" << std::endl; }
    else { if (verbose) cout << "yes" << std::endl; }

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list using templated "get" method.
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------

    max_iters_ptr = PL_My_Polynomial.INVALID_TEMPLATE_QUALIFIER getPtr<int>("Max Iters");
    nonlin_solver_ptr = PL_Main.INVALID_TEMPLATE_QUALIFIER getPtr<std::string>("Nonlinear Solver");
    (void)nonlin_solver_ptr; // Not used, bad test!

    if (verbose) {
      cout<< "Is the templated 'getPtr' method functional ... "<<std::endl;
      cout<< "  Can we retrieve information using the CORRECT variable type ... ";
    }
    if (max_iters_ptr) {
      if ((*max_iters_ptr)==3) {
        if (verbose) cout << "yes" << std::endl; 
      }
      else { if (verbose) cout << "no" << std::endl; FailedTests++; }
    }

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list that we know is a bad "get".
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------

    mbf_ptr = PL_LinSol.INVALID_TEMPLATE_QUALIFIER getPtr<float>( "Tol" );

    if (mbf_ptr)
      ++FailedTests;        

    if (verbose) {
      cout<< "  Can we retrieve information using the WRONG variable type ... ";
    }
    if (!mbf_ptr) { if (verbose) cout << "no" << std::endl; }
    else { if (verbose) cout << "yes" << std::endl; }

    //-----------------------------------------------------------
    // Check the 'getParameter' helper function.
    //-----------------------------------------------------------
    int def_step = 0;
    double alpha_fact = 0.0;
    tempMeth = true;
    try {
      def_step = Teuchos::getParameter<int>(PL_Polynomial, "Default Step");
      alpha_fact = Teuchos::getParameter<double>(PL_Polynomial, "Alpha Factor");
      (void)alpha_fact; // Not used, bad test!
    }
    catch( const Teuchos::Exceptions::InvalidParameter&) { tempMeth = false; }
    if (verbose && def_step==1) {
      cout<< "Is the helper function 'getParameter' functional ... ";
    }
    if (tempMeth) { if (verbose) cout << "yes" << std::endl; }
    else { if (verbose) cout << "no" << std::endl; FailedTests++; }

    //-----------------------------------------------------------
    // Check templated isType functionality
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------
    bool PT1, PT2, PT3;
    PT1 = PL_Polynomial.INVALID_TEMPLATE_QUALIFIER isType<int>("Default Step");
    PT2 = PL_Polynomial.INVALID_TEMPLATE_QUALIFIER isType<long int>("Default Step");
    PT3 = PL_Polynomial.INVALID_TEMPLATE_QUALIFIER isType<std::string>("Interpolation Type");
    if (verbose) {
      cout<< "Is the templated 'isType' method functional ... "<<std::endl;
      cout<< "  Is the 'Default Step' of type 'int' ... ";
    }
    if (PT1) { if (verbose) cout<< "yes" << std::endl; }
    else { if (verbose) cout<< "no" << std::endl; FailedTests++; }
    if (verbose) {
      cout<< "  Is the 'Default Step' of type 'long int' ... ";
    }
    if (PT2) { if (verbose) cout<< "yes" << std::endl; FailedTests++; }
    else { if (verbose) cout<< "no (as expected)" << std::endl; }
    if (verbose) {
    	cout<< "  Is the 'Interpolation Type' of type 'std::string' ... ";
    }
    if (PT3) { if (verbose) cout<< "yes" << std::endl; }
    else { if (verbose) cout<< "no" << std::endl; FailedTests++; }

    //-----------------------------------------------------------
    // Check the 'isParameterType' helper function.
    //-----------------------------------------------------------
    bool PT4, PT5;
    PT4 = Teuchos::isParameterType<double>(PL_Polynomial, "Max Bounds Factor");
    PT5 = Teuchos::isParameterType<float>(PL_Polynomial, "Max Bounds Factor");    
    if (verbose) {
      cout<< "Is the helper function 'isParameterType' functional ... "<<std::endl;
      cout<< "  Is the 'Max Bounds Factor' of type 'double' ... ";
    }
    if (PT4) { if (verbose) cout<< "yes" <<std::endl; }
    else { if (verbose) cout<< "no" << std::endl; FailedTests++; }
    if (verbose) {
      cout<< "  Is the 'Max Bounds Factor' of type 'float' ... ";
    }
    if (PT5) { if (verbose) cout<< "yes" <<std::endl; FailedTests++; }
    else { if (verbose) cout<< "no (as expected)" << std::endl; }

    //-----------------------------------------------------------
    // Can we pass a pointer to a std::vector to the parameter list.
    //-----------------------------------------------------------
    Teuchos::ArrayRCP<double> tempvec1_arcp = Teuchos::arcp<double>(10);
    {
      double * tempvec1 = tempvec1_arcp.get();
      for (int i=0; i<10; i++) { tempvec1[i] = i; }
      PL_Main.set( "Address of Norm Vector", tempvec1 );
      double* tempvec2 = Teuchos::getParameter<double*>( PL_Main, "Address of Norm Vector" );
      tempvec1[4] = 2.0; tempvec1[6] = 1.0;
      if (verbose) {
        cout<< "Can we pass a pointer to a std::vector to a parameter list ... ";
      }
      if ((tempvec2[4]-tempvec1[4])!=0.0 || (tempvec2[6]-tempvec1[6])!=0.0) {
        if (verbose) { cout<<"no"<<std::endl; }
        FailedTests++;
      } else {
        if (verbose) { cout<<"yes"<<std::endl; }
      }
    }

    //-----------------------------------------------------------
    // We can add Array<T> objects of various types T as std::string parameters!
    //-----------------------------------------------------------
    if(verbose) {
      print_break();
      cout << "Setting int and double array objects as std::string parameters ...\n";
      print_break();
    }
    const Teuchos::Array<int>
      intArray = Teuchos::tuple<int>(0,1,2,3,4,5,6);
    const Teuchos::Array<double>
      doubleArray = Teuchos::tuple<double>(0,1.0,2.0,3.0,4.0,5.0,6.0);
    //const Teuchos::Array<bool>
    //boolArray = Teuchos::tuple<bool>(true,true,false,false);
    Teuchos::setStringParameterFromArray("Int Array",intArray,&PL_Main);
    Teuchos::setStringParameterFromArray("Double Array",doubleArray,&PL_Main);
    //Teuchos::setStringParameterFromArray("Bool Array",boolArray,&PL_Main);
    if(verbose) {
      print_break();
      cout << "Testing retrieval of set array objects ...\n";
      print_break();
    }
    {

      const Teuchos::Array<int>
        readIntArray = Teuchos::getArrayFromStringParameter<int>(PL_Main,"Int Array");
      result = readIntArray == intArray;
      if(!result) ++FailedTests;
      if(verbose)
        cout
          << "readIntArray = " << readIntArray << " == intArray = " << intArray << " ? "
          << (result ? "passed" : "failed")
          << "\n";

      const Teuchos::Array<int>
        readDoubleAsIntArray = Teuchos::getArrayFromStringParameter<int>(PL_Main,"Double Array");
      result = readDoubleAsIntArray == intArray;
      if(!result) ++FailedTests;
      if(verbose)
        cout
          << "readDoubleAsIntArray = " << readDoubleAsIntArray << " == intArray = " << intArray << " ? "
          << (result ? "passed" : "failed")
          << "\n";

/*
      const Teuchos::Array<bool>
        readBoolArray = Teuchos::getArrayFromStringParameter<bool>(PL_Main,"Bool Array");
      result = readBoolArray == boolArray;
      if(!result) ++FailedTests;
      if(verbose)
        cout
          << "readBoolArray = " << readBoolArray << " == boolArray = " << boolArray << " ? "
          << (result ? "passed" : "failed")
          << "\n";
*/

      if(verbose) {
        print_break();
      }
      
    }

    //-----------------------------------------------------------
    // We can directly set the array strings!
    //-----------------------------------------------------------

    if(verbose) {
      print_break();
      cout << "Setting a array of doubles as a std::string parameter directly ...\n";
      print_break();
    }

    PL_Main.set(
      "Double Array"
      ,"  {\n"
      "      0.00000\n"
      "      ,1.0e0\n"
      "      ,0.2e1\n"
      "      ,30.0e-1\n"
      "      ,4\n"
      "      ,5.0000\n"
      "      ,6\n"
      "  }\n  "
      );

    {

      const Teuchos::Array<double>
        readDoubleArray = Teuchos::getArrayFromStringParameter<double>(PL_Main,"Double Array");
      Teuchos::Array<int>
        doubleAsIntArray(readDoubleArray.size());
      for(int i=0;i<static_cast<int>(readDoubleArray.size());++i)
        doubleAsIntArray[i] = static_cast<int>(readDoubleArray[i]);
      result = doubleAsIntArray == intArray;
      if(!result) ++FailedTests;
      if(verbose)
        cout
          << "doubleAsIntArray = " << doubleAsIntArray << " == intArray = " << intArray << " ? "
          << (result ? "passed" : "failed")
          << "\n";

      if(verbose) {
        print_break();
      }
      
    }

    //-----------------------------------------------------------
    // Can we pass a pointer to a function to the parameter list.
    // Use a simple function, pass it in and get it back out ...
    // ( HKT 03/23/2004 This test is not supported on Janus )
    //-----------------------------------------------------------
#ifndef JANUS_STLPORT 
    double (*pt2Function) (double, double);
    PL_Main.set( "Address to Simple Function", &Plus );
    pt2Function = Teuchos::getParameter<double(*)(double,double)>( PL_Main, "Address to Simple Function" ); 
    if (verbose) {
      cout<< "Can we pass a pointer to a function to a parameter list ... ";
    }
    if ( pt2Function( 1.0, 2.0 ) != 3.0 ) {
      if (verbose) cout<<"no"<<std::endl;
      FailedTests++;
    } else {
      if (verbose) cout<<"yes"<<std::endl;
    }    
#endif
  }

  //-----------------------------------------------------------
  // We can store and retrieve void* pointers!
  //-----------------------------------------------------------
  
  {
    ParameterList pl;
    int someInt = 1;
    void *someIntPtr = &someInt;
    pl.set("Some Pointer", someIntPtr);
    void *someIntPtrRtn = getParameter<void*>(pl, "Some Pointer");
    TEUCHOS_TEST_FOR_EXCEPT(someIntPtrRtn != someIntPtr);
    if (verbose)
      cout << "someIntPtrRtn = " << someIntPtrRtn << " == " << someIntPtr << " : ";
    if (someIntPtrRtn == someIntPtr) {
      if (verbose) cout << "passed\n";
    }
    else {
      if (verbose) cout << "failed\n";
      FailedTests++;
    }
  }

  //-----------------------------------------------------------
  // Print using the public iterators
  // KL - 7 August 2004
  //-----------------------------------------------------------
  ParameterList::ConstIterator iter;
  
  if (verbose) 
  {
    print_break();
    cout << " printing using public iterators " 
         << std::endl;
    print_break();
  }
  for (iter = PL_Main.begin(); iter != PL_Main.end(); ++iter)
  {
    const ParameterEntry& val = PL_Main.entry(iter);
    const std::string& name = PL_Main.name(iter);
    if (val.isList())
    {
      if (verbose) cout << name << std::endl;
      const ParameterList& sublist = Teuchos::getValue<ParameterList>(val);
      ParameterList::ConstIterator i;
      for (i=sublist.begin(); i != sublist.end(); ++i)
      {
        const std::string& nm = sublist.name(i);              
        const ParameterEntry& v = sublist.entry(i);
        if (v.isList())
        {
          if (verbose) cout << "  " << nm << std::endl;
          if (verbose) Teuchos::getValue<ParameterList>(v).print(cout, 6);
        }
        else
        {
          if (verbose) cout << "  " << nm << " " << v << std::endl;
        }
      }
    }
    else
    {
      if (verbose) cout << name << " " << val << std::endl;
    }
  }


#if defined(HAVE_TEUCHOS_EXTENDED)

  try {

    if (verbose) {

      print_break();
      cout << "writing to XML std::ostream" << std::endl;
      print_break();
      writeParameterListToXmlOStream(PL_Main,cout);

      print_break();
      cout << "writing to XML file" << std::endl;
      print_break();
      writeParameterListToXmlFile(PL_Main,"PL_Main.xml");

      print_break();
      cout << "reading from XML file" << std::endl;
      print_break();
      ParameterList readBack;
      updateParametersFromXmlFile("PL_Main.xml", inoutArg(readBack));
      if (verbose) readBack.print(cout);

      print_break();
      cout << "reading from XML std::string" << std::endl;
      print_break();
      std::ifstream xmlInFile("PL_Main.xml");
      std::string xmlStr;
      while(!xmlInFile.eof()) {
        std::string line;
        std::getline(xmlInFile,line);
        xmlStr += line + "\n";
      }
      readBack = ParameterList();
      updateParametersFromXmlString(xmlStr, inoutArg(readBack));
      if (verbose) readBack.print(cout);

    }

  }
  catch(const std::exception& e)
  {
    if(verbose) {
      std::cerr << "caught std::exception:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
    FailedTests++;
  }

#endif // defined(HAVE_TEUCHOS_EXTENDED)

  //-----------------------------------------------------------
  // Print out main list
  //-----------------------------------------------------------

  if (verbose) {
    print_break();
    cout << "The Final Parameter List" << std::endl;
    print_break();
    PL_Main.print(cout);
    print_break();
    cout << "The unused parameters" << std::endl;
    PL_Main.unused(cout);
  }

  //-----------------------------------------------------------
  // Show error outputs
  //-----------------------------------------------------------

  if (verbose) {
    print_break();
    cout << "Accessing a sublist using the wrong name (should throw a Teuchos::Exceptions::InvalidParameterName std::exception)...\n";
    print_break();
  }
  try {
    PL_Main.sublist("Direction").sublist("Newton").sublist("Linear Solvers",true);
    if (verbose) cout << "Did not throw std::exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterName &e) {
    std::cerr << "caught expected Teuchos::Exceptions::InvalidParameterName:\n\n";
    OSTab tab(std::cerr);
    std::cerr << e.what() << std::endl;
  }
  if (verbose) {
    print_break();
    cout << "Accessing a parameter using the wrong name (should throw a Teuchos::Exceptions::InvalidParameterName std::exception)...\n";
    print_break();
  }
  try {
    Teuchos::getParameter<int>(PL_Main.sublist("Direction").sublist("Newton").sublist("Linear Solver"),"Tolerances");
    if (verbose) cout << "Did not throw std::exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterName &e) {
    if(verbose) {
      std::cerr << "caught expected Teuchos::Exceptions::InvalidParameterName:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }

  if (verbose) {
    print_break();
    cout << "Accessing a parameter using the wrong parameter type (should throw a Teuchos::Exceptions::InvalidParameterType std::exception)...\n";
    print_break();
  }
  try {
    Teuchos::getParameter<int>(PL_Main.sublist("Direction").sublist("Newton").sublist("Linear Solver"),"Tolerance");
    if (verbose) cout << "Did not throw std::exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterType &e) {
    if(verbose) {
      std::cerr << "caught expected Teuchos::Exceptions::InvalidParameterType:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }

  //-----------------------------------------------------------
  // Validate the parameter list
  //-----------------------------------------------------------

  // Create a validator version of PL_Main that we will validate against!
  Teuchos::ParameterList PL_Main_valid("PL_Main_copy");
  PL_Main_valid.setParameters(PL_Main);

  // Create a validator for the "Nonlinear Solver" parameter
  Teuchos::setStringToIntegralParameter<int>(
    "Nonlinear Solver", "Line Search Based",
    "Selects the type of nonlinear solver to use",
    Teuchos::tuple<std::string>("Line Search Based","Trust Region Based"),
    &PL_Main
    );

/*    
  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    nonlinearSolverValidator = rcp(
    new Teuchos::StringToIntegralParameterEntryValidator<int>(
      Teuchos::tuple<std::string>("Line Search Based","Trust Region Based")
      ,"Nonlinear Solver"
      )
    );
  PL_Main_valid.set(
    "Nonlinear Solver", "Line Search Based"
    ,"Selects the type of nonlinear solver to use"
    ,nonlinearSolverValidator
    );
*/

  // Create a validator for the parameter "Line Search"->"Polynomial"->"Max Iters"
  // that accepts an 'int', a 'double' or a 'std::string' value!
  typedef Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes AcceptedTypes;
  Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
    linesearchMaxItersValiator = rcp(
      new Teuchos::AnyNumberParameterEntryValidator(
        Teuchos::AnyNumberParameterEntryValidator::PREFER_INT, // Not used here!
        AcceptedTypes(false).allowInt(true).allowDouble(true).allowString(true)
        )
      );
  PL_Main_valid.sublist("Line Search").sublist("Polynomial").set(
    "Max Iters",3
    ,"The maximum number of inner linear search iterations allowed."
    ,linesearchMaxItersValiator
    );

  // Create a validator for the parameter "Direction"->"Newton"->"Linear Solver"->"Tol"
  // that accepts a 'double' or a 'std::string' value!
  typedef Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes AcceptedTypes;
  Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
    linSolveTolValidator = rcp(
      new Teuchos::AnyNumberParameterEntryValidator(
        Teuchos::AnyNumberParameterEntryValidator::PREFER_INT, // Not used here!
        AcceptedTypes(false).allowDouble(true).allowString(true)
        )
      );
  PL_Main_valid.sublist("Direction",true).sublist("Newton",true)
    .sublist("Linear Solver",true).set(
      "Tol", double(1e-5)
      ,"Select the linear solve tolerance"
    ,linSolveTolValidator
    );

  if (verbose) {
    print_break();
    cout << "Validating the parameter list against itself (should not throw std::exception)...\n";
    print_break();
  }
  try {
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw std::exception, success!\n\n";
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }

  if (verbose) {
    print_break();
    cout << "Adding an invalid parameter type then validating (should throw a Teuchos::Exceptions::InvalidParameterType std::exception)...\n";
    print_break();
  }
  try {
    PL_Main.sublist("Line Search").sublist("Polynomial").set("Max Iters",(short int)(3)); // Should be an int!
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw std::exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterType &e) {
    if(verbose) {
      std::cerr << "caught expected Teuchos::Exceptions::InvalidParameterType:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }
  PL_Main.sublist("Line Search").sublist("Polynomial").set("Max Iters",3); // Put back the valid int!

  if (verbose) {
    print_break();
    cout << "Adding an invalid parameter name then validating (should throw a Teuchos::Exceptions::InvalidParameterName std::exception)...\n";
    print_break();
  }
  try {
    PL_Main.sublist("Line Search").sublist("Polynomial").set("Max Iter",10);
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw std::exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterName &e) {
    if(verbose) {
      std::cerr << "caught expected Teuchos::Exceptions::InvalidParameterName:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }
  PL_Main.sublist("Line Search").sublist("Polynomial").remove("Max Iter");

  if (verbose) {
    print_break();
    cout << "Adding an invalid parameter type then validating using validator (should throw a Teuchos::Exceptions::InvalidParameterType std::exception)...\n";
    print_break();
  }
  try {
    PL_Main.set("Nonlinear Solver",int(0)); // Should be a std::string!
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw std::exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterType &e) {
    if(verbose) {
      std::cerr << "caught expected Teuchos::Exceptions::InvalidParameterType:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }
  PL_Main.set("Nonlinear Solver","Line Search Based"); // Put back the valid value!

  if (verbose) {
    print_break();
    cout << "Adding an invalid parameter value then validating using validator (should throw a Teuchos::Exceptions::InvalidParameterValue std::exception)...\n";
    print_break();
  }
  try {
    PL_Main.set("Nonlinear Solver","LineSearch Based"); // Should be "Line Search Based"!
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw std::exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterValue &e) {
    if(verbose) {
      std::cerr << "caught expected Teuchos::Exceptions::InvalidParameterValue:\n\n";
      OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }
  PL_Main.set("Nonlinear Solver","Line Search Based"); // Put back the valid value!

  if (verbose) {
    print_break();
    cout << "Use the validator to access integral value (should *not* throw std::exception)...\n";
    print_break();
  }
  try {
    const int
      nonlinearSolverValue = Teuchos::getIntegralValue<int>(PL_Main,"Nonlinear Solver");
    const bool
      l_result = (nonlinearSolverValue == 0);
    cout
      << "Read value = " << nonlinearSolverValue << " == 0 : "
      << ( l_result ? "passed" : "failed") << "\n";
    if(!l_result) ++FailedTests;
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }

  if (verbose) {
    print_break();
    cout << "Use the validator to access std::string value (should *not* throw std::exception)...\n";
    print_break();
  }
  try {
    const std::string
      nonlinearSolverValue = Teuchos::getStringValue<int>(PL_Main,"Nonlinear Solver");
    const bool
      l_result = (nonlinearSolverValue == "Line Search Based");
    cout
      << "Read value = \"" << nonlinearSolverValue << " == \"Line Search Based\" : "
      << ( l_result ? "passed" : "failed") << "\n";
    if(!l_result) ++FailedTests;
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }

  //-----------------------------------------------------------
  // Validate and set defaults
  //-----------------------------------------------------------

  // Set the default parameters for an emplty list
  ParameterList validatedPL;

  if (verbose) {
    print_break();
    cout << "Validating and setting defaults for an empty parameter list (should not throw) ...\n";
    print_break();
  }
  try {
    validatedPL.validateParametersAndSetDefaults(PL_Main_valid);
    if (verbose) cout << "Did not throw std::exception, success!\n\n";
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }

  if (verbose) {
    print_break();
    cout << "Parameter list with defaults set:" << std::endl;
    print_break();
    validatedPL.print(cout,PLPrintOptions().showTypes(true).showDoc(true));
    print_break();
  }

  if (verbose) {
    print_break();
    cout << "Checking that validatedPL and PL_Main_valid have the same values : ";
  }
  result = haveSameValues(validatedPL,PL_Main_valid);
  if(!result)
    ++FailedTests;
  if (verbose) {
    cout << ( result ? "passed" : "failed" ) << "\n";
    print_break();
  }

  //
  // Testing access of numbers using validator where validator is not buried
  // in the parameter
  //

  for( int type_i = 0; type_i < 3; ++type_i ) {

    ParameterList &Polynomial_sublist
      = PL_Main.sublist("Line Search",true).sublist("Polynomial",true);
    
    std::string typeName;

    // Set the input type

    switch(type_i) {
      case 0:
        typeName = "int";
        Polynomial_sublist.set("Max Iters",(int)(3));
        break;
      case 1:
        typeName = "double";
        Polynomial_sublist.set("Max Iters",(double)(3.0));
        break;
      case 2:
        typeName = "std::string";
        Polynomial_sublist.set("Max Iters",(std::string)("3"));
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }

    // Extract using external validator

    if (verbose) {
      print_break();
      cout << "Use the external number validator to access a "<<typeName<<" as an int ...\n";
      print_break();
    }
    try {
      const int
        lineserchMaxIters
        = linesearchMaxItersValiator->getInt(
          PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters",0
          );
      const bool
        l_result = (lineserchMaxIters == int(3));
      cout
        << "Read value = " << lineserchMaxIters << " == 3 : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

    if (verbose) {
      print_break();
      cout << "Use the external number validator to access a "<<typeName<<" as a double ...\n";
      print_break();
    }
    try {
      const double
        lineserchMaxIters
        = linesearchMaxItersValiator->getDouble(
          PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters",0.0
          );
      const bool
        l_result = (lineserchMaxIters == double(3.0));
      cout
        << "Read value = " << lineserchMaxIters << " == 3 : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

    if (verbose) {
      print_break();
      cout << "Use the external number validator to access a "<<typeName<<" as a std::string ...\n";
      print_break();
    }
    try {
      const std::string
        lineserchMaxIters
        = linesearchMaxItersValiator->getString(
          PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters","0"
          );
      const bool
        l_result = (lineserchMaxIters == "3");
      cout
        << "Read value = \"" << lineserchMaxIters << "\" == \"3\" : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

    // Extract using nonmember functions

    if (verbose) {
      print_break();
      cout << "Use the nomember help function to access a "<<typeName<<" as an int ...\n";
      print_break();
    }
    try {
      const int
        lineserchMaxIters
        = Teuchos::getIntParameter(
          PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters"
          );
      const bool
        l_result = (lineserchMaxIters == int(3));
      cout
        << "Read value = " << lineserchMaxIters << " == 3 : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

    if (verbose) {
      print_break();
      cout << "Use the nomember help function to access a "<<typeName<<" as a double ...\n";
      print_break();
    }
    try {
      const double
        lineserchMaxIters
        = Teuchos::getDoubleParameter(
          PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters"
          );
      const bool
        l_result = (lineserchMaxIters == double(3.0));
      cout
        << "Read value = " << lineserchMaxIters << " == 3 : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

    if (verbose) {
      print_break();
      cout << "Use the nomember help function to access a "<<typeName<<" as a std::string ...\n";
      print_break();
    }
    try {
      const std::string
        lineserchMaxIters
        = Teuchos::getNumericStringParameter(
          PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters"
          );
      const bool
        l_result = (lineserchMaxIters == "3");
      cout
        << "Read value = \"" << lineserchMaxIters << "\" == \"3\" : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

  }

  //
  // Testing access of numbers using validator where validator is buried in
  // the parameter
  //

  for( int type_i = 0; type_i < 3; ++type_i ) {

    ParameterList &Polynomial_sublist
      = PL_Main.sublist("Line Search",true).sublist("Polynomial",true);
    
    std::string typeName;

    // Set the input type

    switch(type_i) {
      case 0:
        typeName = "int";
        Teuchos::setIntParameter("Max Iters",3,"",&Polynomial_sublist);
        break;
      case 1:
        typeName = "double";
        Teuchos::setDoubleParameter("Max Iters",3.0,"",&Polynomial_sublist);
        break;
      case 2:
        typeName = "std::string";
        Teuchos::setNumericStringParameter("Max Iters","3","",&Polynomial_sublist);
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }

    // Extract using nonmember functions (which should use the internal validator)

    if (verbose) {
      print_break();
      cout << "Use the nomember help function to access a "<<typeName<<" as an int ...\n";
      print_break();
    }
    try {
      const int
        lineserchMaxIters
        = Teuchos::getIntParameter(
          PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters"
          );
      const bool
        l_result = (lineserchMaxIters == int(3));
      cout
        << "Read value = " << lineserchMaxIters << " == 3 : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

    if (verbose) {
      print_break();
      cout << "Use the nomember help function to access a "<<typeName<<" as a double ...\n";
      print_break();
    }
    try {
      const double
        lineserchMaxIters
        = Teuchos::getDoubleParameter(
          PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters"
          );
      const bool
        l_result = (lineserchMaxIters == double(3.0));
      cout
        << "Read value = " << lineserchMaxIters << " == 3 : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

    if (verbose) {
      print_break();
      cout << "Use the nomember help function to access a "<<typeName<<" as a std::string ...\n";
      print_break();
    }
    try {
      const std::string
        lineserchMaxIters
        = Teuchos::getNumericStringParameter(
          PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters"
          );
      const bool
        l_result = (lineserchMaxIters == "3");
      cout
        << "Read value = \"" << lineserchMaxIters << "\" == \"3\" : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

  }

  //
  // Testing access of numbers where correct number is set using
  // validateParametersAndSetDefaults(...) with no special access.
  //

  for( int type_i = 0; type_i < 3; ++type_i ) {

    ParameterList valid_PL_Main(PL_Main);

    ParameterList &Polynomial_sublist
      = PL_Main.sublist("Line Search",true).sublist("Polynomial",true);
    
    std::string typeName;

    // Set the input type

    switch(type_i) {
      case 0:
        typeName = "int";
        Teuchos::setIntParameter("Max Iters",3,"",&Polynomial_sublist);
        break;
      case 1:
        typeName = "double";
        Teuchos::setDoubleParameter("Max Iters",3.0,"",&Polynomial_sublist);
        break;
      case 2:
        typeName = "std::string";
        Teuchos::setNumericStringParameter("Max Iters","3","",&Polynomial_sublist);
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }

    // Extract using nonmember functions (which should use the internal validator)

    if (verbose) {
      print_break();
      cout << "Use validateParemetersAndSetDefaults(...) to access a "<<typeName<<" as an int ...\n";
      print_break();
    }
    try {
      Teuchos::setIntParameter(
        "Max Iters", 0, "",
        &valid_PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
        );
      ParameterList copied_PL_Main(PL_Main);
      copied_PL_Main.validateParametersAndSetDefaults(valid_PL_Main);
      const int
        lineserchMaxIters
        = Teuchos::getParameter<int>(
          copied_PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters"
          );
      const bool
        l_result = (lineserchMaxIters == int(3));
      cout
        << "Read value = " << lineserchMaxIters << " == 3 : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

    if (verbose) {
      print_break();
      cout << "Use validateParemetersAndSetDefaults(...) to access a "<<typeName<<" as a double ...\n";
      print_break();
    }
    try {
      Teuchos::setDoubleParameter(
        "Max Iters", 0.0, "",
        &valid_PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
        );
      ParameterList copied_PL_Main(PL_Main);
      copied_PL_Main.validateParametersAndSetDefaults(valid_PL_Main);
      const double
        lineserchMaxIters
        = Teuchos::getParameter<double>(
          copied_PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters"
          );
      const bool
        l_result = (lineserchMaxIters == double(3.0));
      cout
        << "Read value = " << lineserchMaxIters << " == 3 : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

    if (verbose) {
      print_break();
      cout << "Use validateParemetersAndSetDefaults(...) to access a "<<typeName<<" as a std::string ...\n";
      print_break();
    }
    try {
      Teuchos::setNumericStringParameter(
        "Max Iters", "0", "",
        &valid_PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
        );
      ParameterList copied_PL_Main(PL_Main);
      copied_PL_Main.validateParametersAndSetDefaults(valid_PL_Main);
      const std::string
        lineserchMaxIters
        = Teuchos::getParameter<std::string>(
          copied_PL_Main.sublist("Line Search",true).sublist("Polynomial",true)
          ,"Max Iters"
          );
      const bool
        l_result = (lineserchMaxIters == "3");
      cout
        << "Read value = \"" << lineserchMaxIters << "\" == \"3\" : "
        << ( l_result ? "passed" : "failed") << "\n";
      if(!l_result) ++FailedTests;
    }
    catch(const std::exception &e) {
      if(verbose) {
        std::cerr << "caught unexpected std::exception:\n\n";
        OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
      }
      ++FailedTests;
    }

  }

  if (verbose) {
    print_break();
    cout << "Adding an invalid sublist then validating (should throw a Teuchos::Exceptions::InvalidParameterName std::exception)...\n";
    print_break();
  }
  try {
    PL_Main.sublist("Line Search").sublist("Polynomials").set("Max Iters",3); // param correct, sublist wrong
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw std::exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterName &e) {
    if(verbose) {
      std::cerr << "caught expected Teuchos::Exceptions::InvalidParameterName:\n\n";
      OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }
  PL_Main.sublist("Line Search").remove("Polynomials");

  if (verbose) {
    print_break();
    cout << "Validating only the top level list (should not throw std::exception)...\n";
    print_break();
  }
  try {
    PL_Main.validateParameters(PL_Main_valid,0);
    if (verbose) cout << "Did not throw std::exception, success!\n\n";
  }
  catch(const std::exception &e) {
    if(verbose) {
      std::cerr << "caught unexpected std::exception:\n\n";
      OSTab tab(std::cerr); std::cerr << e.what() << std::endl;
    }
    ++FailedTests;
  }

  //-----------------------------------------------------------
  // Compare lists
  //-----------------------------------------------------------

  if (verbose) {
    print_break();
    cout << "Checking that PL_Main == PL_Main == true : ";
  }
  result = (PL_Main == PL_Main);
  if(!result)
    ++FailedTests;
  if (verbose) {
    cout << ( result ? "passed" : "failed" ) << "\n";
    print_break();
  }

  if (verbose) {
    print_break();
    cout << "Checking that PL_Main != PL_Main == false : ";
  }
  result = !(PL_Main != PL_Main);
  if(!result)
    ++FailedTests;
  if (verbose) {
    cout << ( result ? "passed" : "failed" ) << "\n";
    print_break();
  }

  if (verbose) {
    print_break();
    cout << "Checking that PL_Main and PL_Main have the same values : ";
  }
  result = haveSameValues(PL_Main,PL_Main);
  if(!result)
    ++FailedTests;
  if (verbose) {
    cout << ( result ? "passed" : "failed" ) << "\n";
    print_break();
  }

  if (verbose) {
    print_break();
    cout << "Create copy PL_Main_copy, change PL_Main_copy, and check that PL_Main != PL_Main == true : ";
  }
  ParameterList PL_Main_copy(PL_Main);
  PL_Main_copy.sublist("Line Search",true).sublist("Polynomial",true).set("Max Iters",100); // Not the default!
  result = (PL_Main_copy != PL_Main);
  if(!result)
    ++FailedTests;
  if (verbose) {
    cout << ( result ? "passed" : "failed" ) << "\n";
    print_break();
  }

  //-----------------------------------------------------------
  // Print out main list showing the types
  //-----------------------------------------------------------

  if (verbose) {
    print_break();
    cout << "The Final Parameter List with Types and Documentation" << std::endl;
    print_break();
    PL_Main.print(cout,PLPrintOptions().showTypes(true).showDoc(true));
    print_break();
    cout << "The unused parameters" << std::endl;
    PL_Main.unused(cout);
    print_break();
    cout << "Number of Failed Tests : " << FailedTests << std::endl;
    print_break();
  }

  //-----------------------------------------------------------
  // Return -1 if there are any failed tests, 
  // else 0 will be returned indicating a clean finish!  
  //-----------------------------------------------------------

  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);
  if(!success) ++FailedTests;

  if ( FailedTests > 0 ) { 
    cout << "End Result: TEST FAILED" << std::endl;
    return 1; // Can't return negative numbers from main()!
  }

  if ( FailedTests == 0 )
    cout << "End Result: TEST PASSED" << std::endl;

  return 0;

}

