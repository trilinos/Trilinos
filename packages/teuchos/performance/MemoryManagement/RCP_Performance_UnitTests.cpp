#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Time.hpp"

#ifdef HAVE_TEUCHOS_BOOST
#  include "Teuchos_RCPBoostSharedPtrConversions.hpp"
#endif


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;


double relCpuSpeed = 1e-2;
int maxArraySize = 10000;
double maxRcpCreateDestroyRatio = 10.0;


inline
double adjustTime( const double &time_in )
{
  return ( time_in > 0.0 ? time_in : -1.0 );
}


TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp =
    Teuchos::UnitTestRepository::getCLP();
  clp.setOption(
    "rel-cpu-speed", &relCpuSpeed,
    "The relative speed of the CPU (higher means the machine runs faster)"
    );
  clp.setOption(
    "max-array-size", &maxArraySize,
    "The maximum size of the arrays created"
    );
  clp.setOption(
    "max-rcp-create-destroy-ratio", &maxRcpCreateDestroyRatio,
    "The ratio of the final CPU time ratio of creating and destroying"
    "std::vector<char>(size) objects wrapped in an RCP object verses"
    "using just raw new and delete."
    );

}


TEUCHOS_UNIT_TEST( RCP, createDestroyOverhead )
{
  using std::setw;
  using std::left;
  using std::right;

  const int maxLoopIters = 1000;

  const double relTestCost = 1e-4;

  const double numInnerLoops = relCpuSpeed / relTestCost;

  Teuchos::Time timer("");

  const int dbl_w = 15;
  const std::string dbl_line = "---------------";

  const int int_w = 10;
  const std::string int_line = "----------";

  out << "\n"
      << "Messuring the overhead of createing and destorying objects of different sizes\n"
      << "using raw C++ pointers and using RCP.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n"
      << "  " << setw(int_w) << left << "obj size"
      << "  " << setw(int_w) << left << "num loops"
      << "  " << setw(dbl_w) << left << "raw"
#ifdef HAVE_TEUCHOS_BOOST
      << "  " << setw(dbl_w) << left << "shared_ptr"
#endif
      << "  " << setw(dbl_w) << left << "RCP"
#ifdef HAVE_TEUCHOS_BOOST
      << "  " << setw(dbl_w) << left << "shared_ptr/raw"
#endif
      << "  " << setw(dbl_w) << left << "RCP/raw"
      << "\n"
      << "  " << setw(int_w) << right << int_line // obj size
      << "  " << setw(int_w) << right << int_line // num loops
      << "  " << setw(dbl_w) << right << dbl_line // raw
#ifdef HAVE_TEUCHOS_BOOST
      << "  " << setw(dbl_w) << right << dbl_line // shared_ptr
#endif
      << "  " << setw(dbl_w) << right << dbl_line // RCP
#ifdef HAVE_TEUCHOS_BOOST
      << "  " << setw(dbl_w) << right << dbl_line // shared_ptr/raw
#endif
      << "  " << setw(dbl_w) << right << dbl_line // RCP/raw
      << "\n";

  double finalRcpRatio = -1.0;

  int arraySize = 1;
  for (int test_case_k = 0;
    test_case_k < maxLoopIters && arraySize <= maxArraySize;
    ++test_case_k
    )
  {

    // obj size
    out << "  " << setw(int_w) << right << arraySize;

    // num loops
    const int numActualLoops =
      TEUCHOS_MAX(
        static_cast<int>(
          (numInnerLoops / arraySize)
          * std::log(static_cast<double>(arraySize+1))
          ),
        1);
    out << "  " << setw(int_w) << right << numActualLoops;

    // raw
    timer.reset();
    timer.start();
    for ( int k = 0; k < numActualLoops; ++k ) {
      std::vector<char> *p = new std::vector<char>(arraySize, 1);
      delete p;
    }
    timer.stop();
    const double rawPtrTime = adjustTime(timer.totalElapsedTime()) / numActualLoops;
    out << "  " << setw(dbl_w) << right << rawPtrTime;

#ifdef HAVE_TEUCHOS_BOOST
    // shared_ptr
    timer.reset();
    timer.start();
    {
      typedef boost::shared_ptr<std::vector<char> > shared_ptr_t;
      shared_ptr_t p;
      for ( int k = 0; k < numActualLoops; ++k ) {
        p = shared_ptr_t(new std::vector<char>(arraySize, 1));
      }
    }
    timer.stop();
    const double spTime = adjustTime(timer.totalElapsedTime()) / numActualLoops;
    out << "  " << setw(dbl_w) << right << spTime;
#endif

    // RCP
    timer.reset();
    timer.start();
    {
      RCP<std::vector<char> > p;
      for ( int k = 0; k < numActualLoops; ++k ) {
        p = rcp(new std::vector<char>(arraySize, 1));
      }
    }
    timer.stop();
    const double rcpTime = adjustTime(timer.totalElapsedTime()) / numActualLoops;
    out << "  " << setw(dbl_w) << right << rcpTime;

#ifdef HAVE_TEUCHOS_BOOST
    // shared_ptr/rawPtr
    const double spRatio = spTime / rawPtrTime;
    out << "  " << setw(dbl_w) << right << spRatio;
#endif

    // RCP/rawPtr
    const double rcpRatio = rcpTime / rawPtrTime;
    out << "  " << setw(dbl_w) << right << rcpRatio;

    out << "\n";

    arraySize *= 4;
    finalRcpRatio = rcpRatio;
  }

  out << "\n";
  TEST_COMPARE( finalRcpRatio, <=, maxRcpCreateDestroyRatio );
  out << "\n";

}


} // namespace
