// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

#ifdef HAVE_MPI
#include "mpi.h"
#endif

using namespace Teuchos;

void print_break() { cout << "---------------------------------------------------" << endl; }
double Plus ( double a, double b ) { return a+b; }

int main( int argc, char *argv[] )
{
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
    cout << "Processor "<< procRank <<", parse_return "<< parse_return << endl;
    cout << "End Result: TEST FAILED" << endl;
    return parse_return;
  }

  // Only print on 0 processor
  if (procRank != 0 && verbose)
    verbose = false;

  if (verbose)
    cout << Teuchos::Teuchos_Version() << endl << endl;

  //-----------------------------------------------------------
  // Create Main Parameter List / Sublist Structure
  //-----------------------------------------------------------

  ParameterList PL_Main("PL_Main");
  ParameterList& PL_Direction = PL_Main.sublist("Direction");
  ParameterList& PL_Newton = PL_Direction.sublist("Newton");
  ParameterList& PL_LinSol = PL_Newton.sublist("Linear Solver");
  ParameterList& PL_LineSearch = PL_Main.sublist("Line Search");

  //-----------------------------------------------------------
  // Check Parameter List Structure
  //-----------------------------------------------------------
  if (verbose) {
    print_break();
    cout << "Empty Parameter List Structure" << endl;
    print_break();
    cout<<PL_Main<< endl;
  }
  if (verbose) cout << "Is 'Direction' recognized as a sublist of 'Main' ... ";
  if ( PL_Main.isSublist( "Direction" ) ) {  
    if (verbose) cout << "yes"<< endl;
  } else {
    if (verbose) cout << "no"<< endl;
    FailedTests++;        
  }
  if (verbose) cout << "Is 'Newton' recognized as a sublist of 'Direction' ... ";
  if ( PL_Direction.isSublist( "Newton" ) ) {  
    if (verbose) cout << "yes"<< endl;
  } else {
    if (verbose) cout << "no"<< endl;
    FailedTests++;        
  }
  if (verbose) cout << "Is 'Linear Solver' recognized as a sublist of 'Newton' ... ";
  if ( PL_Newton.isSublist( "Linear Solver" ) ) {  
    if (verbose) cout << "yes"<< endl;
  } else {
    if (verbose) cout << "no"<< endl;
    FailedTests++;        
  }
  if (verbose) cout << "Is 'Line Search' recognized as a sublist of 'Main' ... ";
  if ( PL_Main.isSublist( "Line Search" ) ) {  
    if (verbose) cout << "yes"<< endl;
  } else {
    if (verbose) cout << "no"<< endl;
    FailedTests++;        
  }

  //-----------------------------------------------------------
  // Fill in Direction Sublist
  //-----------------------------------------------------------

  double tol = 0.0;
  bool RBNS = false;
  PL_Direction.get("Method", "Newton");
  PL_LinSol.set("Tol",1e-5);
  tol = PL_LinSol.get("Tolerance",1e-10);
  RBNS = PL_Newton.get("Rescue Bad Newton Solve", true );

  //-----------------------------------------------------------
  // Print out Direction Sublist
  //-----------------------------------------------------------
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

  //-----------------------------------------------------------
  // Line Search Sublist 
  // (if there are no failures, this will be constructed and added)
  //-----------------------------------------------------------
  if (!FailedTests) {
    int ARI = 0, default_step = 0, max_iter_inc = 0, rec_step = 0;
    double alpha_factor = 0.0, min_bnds_factor = 0.0, max_bnds_factor = 0.0;
    bool force_interp = true, use_cntrs = false;
    string ls_method = "Polynomial";
    // This is to force a char* to be passed into the get method to see if the template
    // specialization for char* is working.
    char* ls_method_char = const_cast<char *>(ls_method.c_str());
    ParameterList PL_My_LineSearch("PL_My_LineSearch");
    ls_method = PL_My_LineSearch.get("Method", ls_method_char);
    ParameterList& PL_Polynomial = PL_My_LineSearch.sublist("Polynomial");
    ARI = PL_Polynomial.get("Allowed Relative Increase", 100 );
    alpha_factor = PL_Polynomial.get("Alpha Factor", 0.0001 );
    default_step = PL_Polynomial.get("Default Step", 1 );
    force_interp = PL_Polynomial.get("Force Interpolation", false );
    string interp_type = PL_Polynomial.get("Interpolation Type", "Cubic" );
    max_bnds_factor = PL_Polynomial.get("Max Bounds Factor", 0.5 );
    PL_Polynomial.set("Max Iters", 3 );
    max_iter_inc = PL_Polynomial.get("Maximum Iteration for Increase", 0 );
    min_bnds_factor = PL_Polynomial.get("Min Bounds Factor", 0.1 );
    rec_step = PL_Polynomial.get("Recovery Step", 1 );
    string rec_step_type = PL_Polynomial.get("Recovery Step Type", "Constant");
    string suff_dec_cond = PL_Polynomial.get("Sufficient Decrease Condition", "Armijo-Goldstein" );
    use_cntrs = PL_Polynomial.get("Use Counters", true );
    PL_Main.set("Nonlinear Solver", "Line Search Based");
    //-----------------------------------------------------------
    // Set the Line Search Parameter List equal to the one just constructed
    //-----------------------------------------------------------
    PL_LineSearch.setParameters(PL_My_LineSearch);
    ParameterList& PL_My_Polynomial = PL_LineSearch.sublist("Polynomial");
    if (verbose) cout<< "Is 'operator=' functional ... ";
    if ( PL_My_Polynomial.isParameter("Recovery Step Type") ) {
      if (verbose) cout<< "yes" << endl;
    } else {
      if (verbose) cout<< "no" << endl;
      FailedTests++;
    }  
    ParameterList Copied_PL_Main(PL_Main);
    if (verbose) cout<< "Is the copy constructor functional ... ";
    if ( Copied_PL_Main.isParameter("Nonlinear Solver") ) {
      if (verbose) cout<< "yes" << endl;
    } else {
      if (verbose) cout<< "no" << endl;
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
    string nonlin_solver;
    tempMeth = true;
    try {
      max_iters = PL_My_Polynomial.INVALID_TEMPLATE_QUALIFIER get<int>("Max Iters");
      max_iters_again = Teuchos::getConst(PL_My_Polynomial).INVALID_TEMPLATE_QUALIFIER get<int>("Max Iters");
      nonlin_solver = PL_Main.INVALID_TEMPLATE_QUALIFIER get<string>("Nonlinear Solver");
    }
    catch( const Teuchos::Exceptions::InvalidParameter& e ) { tempMeth = false; }  
    if (verbose) {
      cout<< "Is the templated 'get' method functional ... "<<endl;
      cout<< "  Can we retrieve information using the CORRECT variable type ... ";
    }
    if (tempMeth && max_iters==3) { if (verbose) cout << "yes" << endl; }
    else { if (verbose) cout << "no" << endl; FailedTests++; }
    if (verbose) {
      cout<< "  Can we retrieve const information using the CORRECT variable type ... ";
    }
    if (tempMeth && max_iters_again==3) { if (verbose) cout << "yes" << endl; }
    else { if (verbose) cout << "no" << endl; FailedTests++; }

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list that we know is a bad "get".
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------
    float mbf = 0.0;
    tempMeth = false;
    try {
      mbf = PL_LinSol.INVALID_TEMPLATE_QUALIFIER get<float>( "Tol" );
      FailedTests++;
    }
    catch( const Teuchos::Exceptions::InvalidParameter& e ) {
      tempMeth = true;
    }
    if (verbose) {
      cout<< "  Can we retrieve information using the WRONG variable type ... ";
    }
    if (tempMeth) { if (verbose) cout << "no" << endl; }
    else { if (verbose) cout << "yes" << endl; }

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list using templated "get" method.
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------
    tempMeth = true;
    try {
      max_iters = PL_My_Polynomial.INVALID_TEMPLATE_QUALIFIER get<int>("Max Iters");
      nonlin_solver = PL_Main.INVALID_TEMPLATE_QUALIFIER get<string>("Nonlinear Solver");
    }
    catch( const Teuchos::Exceptions::InvalidParameter& e ) { tempMeth = false; }  
    if (verbose) {
      cout<< "Is the templated 'get' method functional ... "<<endl;
      cout<< "  Can we retrieve information using the CORRECT variable type ... ";
    }
    if (tempMeth && max_iters==3) { if (verbose) cout << "yes" << endl; }
    else { if (verbose) cout << "no" << endl; FailedTests++; }

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list that we know is a bad "get".
    //-----------------------------------------------------------
    tempMeth = false;
    try {
      mbf = PL_LinSol.INVALID_TEMPLATE_QUALIFIER get<float>( "Tol" );
      FailedTests++;
    }
    catch( const Teuchos::Exceptions::InvalidParameter& e ) {
      tempMeth = true;
    }
    if (verbose) {
      cout<< "  Can we retrieve information using the WRONG variable type ... ";
    }
    if (tempMeth) { if (verbose) cout << "no" << endl; }
    else { if (verbose) cout << "yes" << endl; }

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
    string* nonlin_solver_ptr;

    max_iters_ptr = PL_My_Polynomial.INVALID_TEMPLATE_QUALIFIER getPtr<int>("Max Iters");
    max_iters_ptr_again = Teuchos::getConst(PL_My_Polynomial).INVALID_TEMPLATE_QUALIFIER getPtr<int>("Max Iters");
    nonlin_solver_ptr = PL_Main.INVALID_TEMPLATE_QUALIFIER getPtr<string>("Nonlinear Solver");

    if (verbose) {
      cout<< "Is the templated 'getPtr' method functional ... "<<endl;
      cout<< "  Can we retrieve information using the CORRECT variable type ... ";
    }
    if (max_iters_ptr) {
      if ((*max_iters_ptr)==3) {
        if (verbose) cout << "yes" << endl; 
      }
      else { if (verbose) cout << "no" << endl; FailedTests++; }
    }
    if (verbose) {
      cout<< "  Can we retrieve const information using the CORRECT variable type ... ";
    }
    if (max_iters_ptr_again) {
      if ((*max_iters_ptr_again)==3) {
        if (verbose) cout << "yes" << endl; 
      }
      else { if (verbose) cout << "no" << endl; FailedTests++; }
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
    if (!mbf_ptr) { if (verbose) cout << "no" << endl; }
    else { if (verbose) cout << "yes" << endl; }

    //-----------------------------------------------------------
    // Retrieve some information from the parameter list using templated "get" method.
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------

    max_iters_ptr = PL_My_Polynomial.INVALID_TEMPLATE_QUALIFIER getPtr<int>("Max Iters");
    nonlin_solver_ptr = PL_Main.INVALID_TEMPLATE_QUALIFIER getPtr<string>("Nonlinear Solver");

    if (verbose) {
      cout<< "Is the templated 'getPtr' method functional ... "<<endl;
      cout<< "  Can we retrieve information using the CORRECT variable type ... ";
    }
    if (max_iters_ptr) {
      if ((*max_iters_ptr)==3) {
        if (verbose) cout << "yes" << endl; 
      }
      else { if (verbose) cout << "no" << endl; FailedTests++; }
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
    if (!mbf_ptr) { if (verbose) cout << "no" << endl; }
    else { if (verbose) cout << "yes" << endl; }

    //-----------------------------------------------------------
    // Check the 'getParameter' helper function.
    //-----------------------------------------------------------
    int def_step = 0;
    double alpha_fact = 0.0;
    tempMeth = true;
    try {
      def_step = getParameter<int>(PL_Polynomial, "Default Step");
      alpha_fact = getParameter<double>(PL_Polynomial, "Alpha Factor");
    }
    catch( const Teuchos::Exceptions::InvalidParameter& e ) { tempMeth = false; }
    if (verbose && def_step==1) {
      cout<< "Is the helper function 'getParameter' functional ... ";
    }
    if (tempMeth) { if (verbose) cout << "yes" << endl; }
    else { if (verbose) cout << "no" << endl; FailedTests++; }

    //-----------------------------------------------------------
    // Check templated isType functionality
    // (This will be tested using the INVALID_TEMPLATE_QUALIFIER which indicates whether a
    //  non-templated code needs ".template" before the method name )
    //-----------------------------------------------------------
    bool PT1, PT2, PT3;
    PT1 = PL_Polynomial.INVALID_TEMPLATE_QUALIFIER isType<int>("Default Step");
    PT2 = PL_Polynomial.INVALID_TEMPLATE_QUALIFIER isType<long int>("Default Step");
    PT3 = PL_Polynomial.INVALID_TEMPLATE_QUALIFIER isType<string>("Interpolation Type");
    if (verbose) {
      cout<< "Is the templated 'isType' method functional ... "<<endl;
      cout<< "  Is the 'Default Step' of type 'int' ... ";
    }
    if (PT1) { if (verbose) cout<< "yes" << endl; }
    else { if (verbose) cout<< "no" << endl; FailedTests++; }
    if (verbose) {
      cout<< "  Is the 'Default Step' of type 'long int' ... ";
    }
    if (PT2) { if (verbose) cout<< "yes" << endl; FailedTests++; }
    else { if (verbose) cout<< "no (as expected)" << endl; }
    if (verbose) {
    	cout<< "  Is the 'Interpolation Type' of type 'string' ... ";
    }
    if (PT3) { if (verbose) cout<< "yes" << endl; }
    else { if (verbose) cout<< "no" << endl; FailedTests++; }

    //-----------------------------------------------------------
    // Check the 'isParameterType' helper function.
    //-----------------------------------------------------------
    bool PT4, PT5;
    PT4 = isParameterType<double>(PL_Polynomial, "Max Bounds Factor");
    PT5 = isParameterType<float>(PL_Polynomial, "Max Bounds Factor");    
    if (verbose) {
      cout<< "Is the helper function 'isParameterType' functional ... "<<endl;
      cout<< "  Is the 'Max Bounds Factor' of type 'double' ... ";
    }
    if (PT4) { if (verbose) cout<< "yes" <<endl; }
    else { if (verbose) cout<< "no" << endl; FailedTests++; }
    if (verbose) {
      cout<< "  Is the 'Max Bounds Factor' of type 'float' ... ";
    }
    if (PT5) { if (verbose) cout<< "yes" <<endl; FailedTests++; }
    else { if (verbose) cout<< "no (as expected)" << endl; }

    //-----------------------------------------------------------
    // Can we pass a pointer to a vector to the parameter list.
    //-----------------------------------------------------------
    double * tempvec1 = new double[10];
    for (int i=0; i<10; i++) { tempvec1[i] = i; }
    PL_Main.set( "Address of Norm Vector", tempvec1 );
    double* tempvec2 = getParameter<double*>( PL_Main, "Address of Norm Vector" );
    tempvec1[4] = 2.0; tempvec1[6] = 1.0;
    if (verbose) {
      cout<< "Can we pass a pointer to a vector to a parameter list ... ";
    }
    if ((tempvec2[4]-tempvec1[4])!=0.0 || (tempvec2[6]-tempvec1[6])!=0.0) {
      if (verbose) { cout<<"no"<<endl; }
      FailedTests++;
    } else {
      if (verbose) { cout<<"yes"<<endl; }
    }	    

    //-----------------------------------------------------------
    // Can we pass a pointer to a function to the parameter list.
    // Use a simple function, pass it in and get it back out ...
    // ( HKT 03/23/2004 This test is not supported on Janus )
    //-----------------------------------------------------------
#ifndef JANUS_STLPORT 
    double (*pt2Function) (double, double);
    PL_Main.set( "Address to Simple Function", &Plus );
    pt2Function = getParameter<double(*)(double,double)>( PL_Main, "Address to Simple Function" ); 
    if (verbose) {
      cout<< "Can we pass a pointer to a function to a parameter list ... ";
    }
    if ( pt2Function( 1.0, 2.0 ) != 3.0 ) {
      if (verbose) cout<<"no"<<endl;
      FailedTests++;
    } else {
      if (verbose) cout<<"yes"<<endl;
    }    
#endif
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
         << endl;
    print_break();
  }
  for (iter = PL_Main.begin(); iter != PL_Main.end(); ++iter)
  {
    const ParameterEntry& val = PL_Main.entry(iter);
    const string& name = PL_Main.name(iter);
    if (val.isList())
    {
      if (verbose) cout << name << endl;
      const ParameterList& sublist = getValue<ParameterList>(val);
      ParameterList::ConstIterator i;
      for (i=sublist.begin(); i != sublist.end(); ++i)
      {
        const string& nm = sublist.name(i);              
        const ParameterEntry& v = sublist.entry(i);
        if (v.isList())
        {
          if (verbose) cout << "  " << nm << endl;
          if (verbose) getValue<ParameterList>(v).print(cout, 6);
        }
        else
        {
          if (verbose) cout << "  " << nm << " " << v << endl;
        }
      }
    }
    else
    {
      if (verbose) cout << name << " " << val << endl;
    }
  }


#if defined(HAVE_TEUCHOS_EXTENDED)

  try {

    if (verbose) {

      print_break();
      cout << "writing to XML ostream" << endl;
      print_break();
      writeParameterListToXmlOStream(PL_Main,cout);

      print_break();
      cout << "writing to XML file" << endl;
      print_break();
      writeParameterListToXmlFile(PL_Main,"PL_Main.xml");

      print_break();
      cout << "reading from XML file" << endl;
      print_break();
      ParameterList readBack;
      updateParametersFromXmlFile("PL_Main.xml",&readBack);
      if (verbose) readBack.print(cout);

      print_break();
      cout << "reading from XML string" << endl;
      print_break();
      std::ifstream xmlInFile("PL_Main.xml");
      std::string xmlStr;
      while(!xmlInFile.eof()) {
        std::string line;
        std::getline(xmlInFile,line);
        xmlStr += line + "\n";
      }
      readBack = ParameterList();
      updateParametersFromXmlString(xmlStr,&readBack);
      if (verbose) readBack.print(cout);

    }

  }
  catch(const std::exception& e)
  {
    if(verbose) {
      cerr << "caught exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    FailedTests++;
  }

#endif // defined(HAVE_TEUCHOS_EXTENDED)

  //-----------------------------------------------------------
  // Print out main list
  //-----------------------------------------------------------

  if (verbose) {
    print_break();
    cout << "The Final Parameter List" << endl;
    print_break();
    PL_Main.print(cout);
    print_break();
    cout << "The unused parameters" << endl;
    PL_Main.unused(cout);
  }

  //-----------------------------------------------------------
  // Show error outputs
  //-----------------------------------------------------------

  if (verbose) {
    print_break();
    cout << "Accessing a sublist using the wrong name (should throw a Teuchos::Exceptions::InvalidParameterName exception)...\n";
    print_break();
  }
  try {
    PL_Main.sublist("Direction").sublist("Newton").sublist("Linear Solvers",true);
    if (verbose) cout << "Did not throw exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterName &e) {
    cerr << "caught expected Teuchos::Exceptions::InvalidParameterName:\n\n";
    OSTab(cerr)() << e.what() << endl;
  }
  if (verbose) {
    print_break();
    cout << "Accessing a parameter using the wrong name (should throw a Teuchos::Exceptions::InvalidParameterName exception)...\n";
    print_break();
  }
  try {
    getParameter<int>(PL_Main.sublist("Direction").sublist("Newton").sublist("Linear Solver"),"Tolerances");
    if (verbose) cout << "Did not throw exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterName &e) {
    if(verbose) {
      cerr << "caught expected Teuchos::Exceptions::InvalidParameterName:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    ++FailedTests;
  }

  if (verbose) {
    print_break();
    cout << "Accessing a parameter using the wrong parameter type (should throw a Teuchos::Exceptions::InvalidParameterType exception)...\n";
    print_break();
  }
  try {
    getParameter<int>(PL_Main.sublist("Direction").sublist("Newton").sublist("Linear Solver"),"Tolerance");
    if (verbose) cout << "Did not throw exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterType &e) {
    if(verbose) {
      cerr << "caught expected Teuchos::Exceptions::InvalidParameterType:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    ++FailedTests;
  }

  //-----------------------------------------------------------
  // Validate the parameter list
  //-----------------------------------------------------------

  // Create a validator version of PL_Main that we will validate against!
  Teuchos::ParameterList PL_Main_valid("PL_Main_copy");
  PL_Main_valid.setParameters(PL_Main);

  Teuchos::RefCountPtr<Teuchos::StringToIntegralParameterEntryValidator<int> >
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

  if (verbose) {
    print_break();
    cout << "Validating the parameter list against itself (should not throw exception)...\n";
    print_break();
  }
  try {
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw exception, success!\n\n";
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    ++FailedTests;
  }

  if (verbose) {
    print_break();
    cout << "Adding an invalid parameter type then validating (should throw a Teuchos::Exceptions::InvalidParameterType exception)...\n";
    print_break();
  }
  try {
    PL_Main.sublist("Line Search").sublist("Polynomial").set("Max Iters",10.0); // Should be an int!
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterType &e) {
    if(verbose) {
      cerr << "caught expected Teuchos::Exceptions::InvalidParameterType:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    ++FailedTests;
  }
  PL_Main.sublist("Line Search").sublist("Polynomial").set("Max Iters",10); // Put back the valid int!

  if (verbose) {
    print_break();
    cout << "Adding an invalid parameter name then validating (should throw a Teuchos::Exceptions::InvalidParameterName exception)...\n";
    print_break();
  }
  try {
    PL_Main.sublist("Line Search").sublist("Polynomial").set("Max Iter",10);
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterName &e) {
    if(verbose) {
      cerr << "caught expected Teuchos::Exceptions::InvalidParameterName:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    ++FailedTests;
  }

  if (verbose) {
    print_break();
    cout << "Adding an invalid parameter type then validating using validator (should throw a Teuchos::Exceptions::InvalidParameterType exception)...\n";
    print_break();
  }
  try {
    PL_Main.set("Nonlinear Solver",int(0)); // Should be a string!
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterType &e) {
    if(verbose) {
      cerr << "caught expected Teuchos::Exceptions::InvalidParameterType:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    ++FailedTests;
  }
  PL_Main.set("Nonlinear Solver","Line Search Based"); // Put back the valid value!

  if (verbose) {
    print_break();
    cout << "Adding an invalid parameter value then validating using validator (should throw a Teuchos::Exceptions::InvalidParameterValue exception)...\n";
    print_break();
  }
  try {
    PL_Main.set("Nonlinear Solver","LineSearch Based"); // Should be "Line Search Based"!
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterValue &e) {
    if(verbose) {
      cerr << "caught expected Teuchos::Exceptions::InvalidParameterValue:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    ++FailedTests;
  }
  PL_Main.set("Nonlinear Solver","Line Search Based"); // Put back the valid value!

  if (verbose) {
    print_break();
    cout << "Use the validator to access integral value (should *not* throw exception)...\n";
    print_break();
  }
  try {
    const int
      nonlinearSolverValue
      = nonlinearSolverValidator->getIntegralValue(PL_Main,"Nonlinear Solver");
    const bool
      result = (nonlinearSolverValue == 0);
    cout
      << "Read value = " << nonlinearSolverValue << " == 0 : "
      << ( result ? "passed" : "failed") << "\n";
    if(!result) ++FailedTests;
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    ++FailedTests;
  }

  if (verbose) {
    print_break();
    cout << "Use the validator to access string value (should *not* throw exception)...\n";
    print_break();
  }
  try {
    const std::string
      nonlinearSolverValue
      = nonlinearSolverValidator->getStringValue(PL_Main,"Nonlinear Solver");
    const bool
      result = (nonlinearSolverValue == "Line Search Based");
    cout
      << "Read value = \"" << nonlinearSolverValue << " == \"Line Search Based\" : "
      << ( result ? "passed" : "failed") << "\n";
    if(!result) ++FailedTests;
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    ++FailedTests;
  }

  if (verbose) {
    print_break();
    cout << "Adding an invalid sublist then validating (should throw a Teuchos::Exceptions::InvalidParameterName exception)...\n";
    print_break();
  }
  try {
    PL_Main.sublist("Line Search").sublist("Polynomials").set("Max Iters",10); // param correct, sublist wrong
    PL_Main.validateParameters(PL_Main_valid);
    if (verbose) cout << "Did not throw exception, error!\n";
    ++FailedTests;
  }
  catch(const Teuchos::Exceptions::InvalidParameterName &e) {
    if(verbose) {
      cerr << "caught expected Teuchos::Exceptions::InvalidParameterName:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
    }
    ++FailedTests;
  }

  if (verbose) {
    print_break();
    cout << "Validating only the top level list (should not throw exception)...\n";
    print_break();
  }
  try {
    PL_Main.validateParameters(PL_Main_valid,0);
    if (verbose) cout << "Did not throw exception, success!\n\n";
  }
  catch(const std::exception &e) {
    if(verbose) {
      cerr << "caught unexpected exception:\n\n";
      OSTab(cerr)() << e.what() << endl;
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
  result = PL_Main == PL_Main;
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
    cout << "Create copy PL_Main_copy, change PL_Main_copy, and check that PL_Main != PL_Main == true : ";
  }
  ParameterList PL_Main_copy(PL_Main);
  PL_Main_copy.sublist("Line Search").sublist("Polynomial").set("Max Iters",100); // Not the default!
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
    cout << "The Final Parameter List with Types" << endl;
    print_break();
    PL_Main.print(cout,0,true);
    print_break();
    cout << "The unused parameters" << endl;
    PL_Main.unused(cout);
    print_break();
    cout << "Number of Failed Tests : " << FailedTests << endl;
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
    cout << "End Result: TEST FAILED" << endl;
    return 1; // Can't return negative numbers from main()!
  }

  if ( FailedTests == 0 )
    cout << "End Result: TEST PASSED" << endl;

  return 0;

}

