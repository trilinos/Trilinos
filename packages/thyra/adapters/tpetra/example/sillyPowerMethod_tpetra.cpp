
#include "sillyPowerMethod.hpp"
#include "createExampleTridiagTpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"

/** This is commented out for now since at the time when this was coded, you
 * can not change a Tpetra::CisMatrix after it is initally constructed.
 */
/*
//
// Increase the first diagonal element of your tridiagonal matrix.
//
template<class Ordinal, class Scalar>
void scaleFirstDiagElement( const Scalar diagScale, Thyra::LinearOpBase<Scalar> *A )
{
  using Teuchos::RCP;
  TEST_FOR_EXCEPT(A==NULL);
  // (A) Get at the underlying Tpetra::Operator object that the TpetraLinearOp
  // object directly maintains.
  const RCP<Tpetra::Operator<Ordinal,Scalar> >
    tpetra_op = Thyra::get_Tpetra_Operator<Ordinal>(*A);
  // (B) Perform a dynamic cast to Tpetra::CisMatrix.
  // Note, the dyn_cast<>() template function will throw std::bad_cast
  // with a nice error message if the cast fails! 
  Tpetra::CisMatrix<Ordinal,Scalar>
    &cisMatrix = Teuchos::dyn_cast<Tpetra::CisMatrix<Ordinal,Scalar> >(*tpetra_op);
  TEST_FOR_EXCEPT(true); // ToDo: Implement below!
} // end scaleFirstDiagElement()
*/

//
// Example showing how to use Tpetra/Thyra adapters.
//
template<class Ordinal, class Scalar>
bool runPowerMethodExample(
  const int    globalDim
  ,const int   maxNumIters
  ,const bool  verbose
  ,const bool  dumpAll
  )
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType    ScalarMag;
  
  bool success = true;
  bool result;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = (verbose ? Teuchos::VerboseObjectBase::getDefaultOStream() : Teuchos::null);

  if(verbose)
    *out << "\n***\n*** Running power method example using scalar type = \'" << ST::name() << "\' ...\n***\n" << std::scientific;

#ifdef HAVE_MPI
  MPI_Comm mpiComm = MPI_COMM_WORLD;
#endif


  //
  // (1) Setup the initial tridiagonal operator
  //
  //       [  2  -1             ]
  //       [ -1   2  -1         ]
  //  A =  [      .   .   .     ]
  //       [          -1  2  -1 ]
  //       [             -1   2 ]
  //
  if(verbose) *out << "\n(1) Constructing tridagonal Tpetra::CisMatrix A of global dimension = " << globalDim << " ...\n";
  RCP<Thyra::LinearOpBase<Scalar> >
    A = createExampleTridiagTpetraLinearOp<Ordinal,Scalar>(
      globalDim
#ifdef HAVE_MPI
      ,mpiComm
#endif
      ,1.0,verbose,*out
      );
  if( verbose && dumpAll ) *out << "\nA =\n" << *A; // This works even in parallel!
  
  //
  // (2) Test the linear operator
  //
  if(verbose) *out << "\n(2) Testing the linear operator A ...\n";
  ScalarMag errTol = ScalarMag(1e2)*Teuchos::ScalarTraits<ScalarMag>::eps();
  Thyra::LinearOpTester<Scalar> linearOpTester;
  linearOpTester.enable_all_tests(false);
  linearOpTester.set_all_warning_tol(ScalarMag(ScalarMag(1e2)*errTol));
  linearOpTester.set_all_error_tol(errTol);
  linearOpTester.check_linear_properties(true);
  result = linearOpTester.check(*A,OSTab(out).get());
  if(!result) success = false;
  
  //
  // (3) Run the power method ANA
  //
  if(verbose) *out << "\n(3) Running the power method on matrix A ...\n";
  Scalar     lambda      = ST::zero();
  ScalarMag  tolerance   = ScalarMag(1e-3)*Teuchos::ScalarTraits<ScalarMag>::one();
  {
    OSTab tab(out);
    result = sillyPowerMethod(*A,maxNumIters,tolerance,&lambda,out.get());
    if(!result) success = false;
    if(verbose) *out << "\nEstimate of dominate eigenvalue lambda = " << lambda << std::endl;
  }

  /** The below code is commented out because when this was coded up, you
   * could not change the entries in a Tpetra::CisMatrix after it is fully
   * constructed.
   */
/*
  //
  // (3) Increase dominance of first eigenvalue
  //
  if(verbose) *out << "\n(3) Scale the diagonal of A by a factor of 10 ...\n";
  scaleFirstDiagElement<Ordinal,Scalar>( 10.0, &*A );
  
  //
  // (4) Run the power method ANA
  //
  if(verbose) *out << "\n(4) Running the power method again on matrix A ...\n";
  {
    OSTab tab(out);
    result = sillyPowerMethod(*A,maxNumIters,tolerance,&lambda,out.get());
    if(!result) success = false;
    if(verbose) *out << "\nEstimate of dominate eigenvalue lambda = " << lambda << std::endl;
  }
*/
  
  return success;

} // end runPowerMethodExample()

//
// Main driver program for serial implementation of the power method.
//
// Note that several different scalar types are used here (float,
// double, std::complex<float>, std::complex<double> etc.).
//
int main(int argc, char *argv[])
{

  using Teuchos::CommandLineProcessor;
 
  bool success = true;
  bool verbose = true;
  bool result;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read in command-line options
    //
    
    int    globalDim    = 500;
    bool   dumpAll      = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
    clp.setOption( "global-dim", &globalDim, "Global dimension of the linear system." );
    clp.setOption( "dump-all", "no-dump", &dumpAll, "Determines if quantities are dumped or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;
    
    TEST_FOR_EXCEPTION( globalDim < 2, std::logic_error, "Error, globalDim=" << globalDim << " < 2 is not allowed!" );
    
    int    maxNumIters  = 10*globalDim;
    
#if defined(HAVE_THYRA_FLOAT)
    // Run using float
    result = runPowerMethodExample<int,float>(globalDim,maxNumIters,verbose,dumpAll);
    if(!result) success = false;
#endif

    // Run using double
    result = runPowerMethodExample<int,double>(globalDim,maxNumIters,verbose,dumpAll);
    if(!result) success = false;

#ifdef HAVE_THYRA_COMPLEX
    
#if defined(HAVE_THYRA_FLOAT)
    // Run using std::complex<float>
    result = runPowerMethodExample<int,std::complex<float> >(globalDim,maxNumIters,verbose,dumpAll);
    if(!result) success = false;
#endif
    
    // Run using std::complex<double>
    result = runPowerMethodExample<int,std::complex<double> >(globalDim,maxNumIters,verbose,dumpAll);
    if(!result) success = false;

#endif

#ifdef HAVE_TEUCHOS_GNU_MP
    
    // Run using mpf_class
    //result = runPowerMethodExample<mpf_class >(globalDim,maxNumIters,verbose,dumpAll);
    //if(!result) success = false;

#ifdef HAVE_THYRA_COMPLEX
    
    // Run using std::complex<mpf_class>
    //result = runPowerMethodExample<std::complex<mpf_class> >(globalDim,maxNumIters,verbose,dumpAll);
    //if(!result) success = false;
    //The above commented-out code throws a floating-point exception?
 
#endif

#endif
    
  }
  catch( const std::exception &excpt ) {
    std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  catch( ... ) {
    std::cerr << "*** Caught an unknown exception\n";
    success = false;
  }
  
  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }
  
  return success ? 0 : 1;

} // end main()
