#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Tpetra_DefaultPlatform.hpp"

#include <sstream>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static Teuchos::ParameterList
parseOptions (int argc, char* argv[], const bool allowedToPrint, bool& printedHelp)
{
  using std::cerr;
  using std::endl;

  printedHelp = false;

  bool verbose = false;
  bool debug = false;
  int nrows = 100000; 
  int ncols = 10;
  int ncores = 1;
  int cache_block_size = 0;

  Teuchos::ParameterList plist;
  plist.set ("verbose", verbose, "Print messages and results");
  plist.set ("debug", debug, "Print debugging information");
  plist.set ("nrows", nrows, "Number of rows (globally) in the test matrix to factor");
  plist.set ("ncols", ncols, "Number of columns in the test matrix to factor");
  plist.set ("ncores", ncores, "Number of cores to use (per MPI process)");
  plist.set ("cache-block-size", cache_block_size, "Cache block size (0 means we set a reasonable default)");

  try {
    Teuchos::CommandLineProcessor cmdLineProc (/* throwExceptions= */ true, /* recognizeAllOptions=*/ true);

    // mfh 08 Jul 2010
    //
    // C++ lacks sufficient introspection mechanisms for cleanly
    // expressing this iteration over options of various types as a
    // loop.  Teuchos::ParameterList has a ConstIterator for iterating
    // over all the parameters, but it stores each parameter as a
    // Teuchos::any (wrapped in a ParameterEntry).  Teuchos::any only
    // gives the type back as an std::type_info, and C++ can't treat
    // that as a template parameter.  Thus, there are only two ways to
    // express this as a loop over command-line arguments:
    //
    // 1. Map the std::type_info of the contents of the Teuchos::any
    // object to a function that invokes the correct overload of
    // setOption() 
    //
    // 2. Use a Smalltalk - style chain of "if type is this do that,
    // else if type is this do that, else if...".  (Stroustrup will
    // hate you for this.)
    //
    // #1 is made more difficult because std::type_info doesn't allow
    // copy construction or assignment; it has a "const char* name()"
    // method, though, so you could map from the type name to the
    // appropriate handler.  #2 is bad in general, but in this case
    // there are only three command-line argument types of interest
    // (bool, int, std::string), so a Smalltalk-style loop would be
    // acceptable.
    cmdLineProc.setOption ("verbose", "quiet", &verbose, plist.getEntry("verbose").docString().c_str());
    cmdLineProc.setOption ("debug", "nodebug", &debug, plist.getEntry("debug").docString().c_str());
    cmdLineProc.setOption ("nrows", &nrows, plist.getEntry("nrows").docString().c_str());
    cmdLineProc.setOption ("ncols", &ncols, plist.getEntry("ncols").docString().c_str());
    cmdLineProc.setOption ("ncores", &ncores, plist.getEntry("ncores").docString().c_str());
    cmdLineProc.setOption ("cache-block-size", &cache_block_size, plist.getEntry("cache-block-size").docString().c_str());
    cmdLineProc.parse (argc, argv);
  } 
  catch (Teuchos::CommandLineProcessor::UnrecognizedOption& e) { 
    if (allowedToPrint)
      cerr << "Unrecognized command-line option: " << e.what() << endl;
    throw e;
  }
  catch (Teuchos::CommandLineProcessor::HelpPrinted& e) { 
    printedHelp = true;
  } 

  // Validate.  TODO (mfh 08 Jul 2010) Figure out how to do this with
  // ParameterList validators.
  if (nrows <= 0)
    throw std::domain_error ("Number of rows must be positive");
  else if (ncols <= 0)
    throw std::domain_error ("Number of columns must be positive");
  else if (nrows < ncols)
    throw std::domain_error ("Number of rows must be >= number of columns");
  else if (ncores < 1)
    throw std::domain_error ("Number of cores must be positive");
  else if (cache_block_size < 0)
    throw std::domain_error ("Cache block size must be nonnegative");

  if (! printedHelp)
    {
      // Fetch the (possibly altered) values of the command-line options.
      plist.set ("verbose", verbose, plist.getEntry("verbose").docString().c_str());
      plist.set ("debug", debug, plist.getEntry("debug").docString().c_str());
      plist.set ("nrows", nrows, plist.getEntry("nrows").docString().c_str());
      plist.set ("ncols", ncols, plist.getEntry("ncols").docString().c_str());
    }
  return plist;
}


int 
main (int argc, char *argv[]) 
{
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  Teuchos::RCP< const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  size_t myRank = comm->getRank();
  //size_t numProc = comm->getSize();
  bool allowedToPrint = (myRank==0);
  bool printedHelp = false;
  
  Teuchos::ParameterList plist = parseOptions (argc, argv, allowedToPrint, printedHelp);
  if (printedHelp)
    return 0;

  if (allowedToPrint) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}


