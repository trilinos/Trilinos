#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TabularOutputter.hpp"

#ifdef HAVE_TEUCHOS_BOOST
#  include "Teuchos_RCPBoostSharedPtrConversions.hpp"
#endif


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TabularOutputter;


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

  typedef Teuchos::TabularOutputter TO;

  const int maxLoopIters = 1000;

  const double relTestCost = 1e-4;

  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Messuring the overhead of creating and destorying objects of different sizes\n"
      << "using raw C++ pointers, shared_ptr, and using RCP.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, 8);
  outputter.setFieldTypePrecision(TO::INT, 8);

  outputter.pushField("obj size", TO::INT);
  outputter.pushField("num loops", TO::INT);
  outputter.pushField("raw", TO::DOUBLE);
#ifdef HAVE_TEUCHOS_BOOST
  outputter.pushField("shared_ptr", TO::DOUBLE);
#endif
  outputter.pushField("RCP", TO::DOUBLE);
#ifdef HAVE_TEUCHOS_BOOST
  outputter.pushField("shared_ptr/raw", TO::DOUBLE);
#endif
  outputter.pushField("RCP/raw", TO::DOUBLE);

  outputter.outputHeader();

  double finalRcpRatio = -1.0;

  int arraySize = 1;
  for (int test_case_k = 0;
    test_case_k < maxLoopIters && arraySize <= maxArraySize;
    ++test_case_k
    )
  {

    // obj size
    outputter.outputField(arraySize);

    // num loops
    const int numActualLoops =
      TEUCHOS_MAX(
        static_cast<int>(
          (numInnerLoops / arraySize)
          * std::log(static_cast<double>(arraySize+1))
          ),
        1
        );
    outputter.outputField(numActualLoops);

    // raw
    TEUCHOS_START_PERF_OUTPUT_TIMER(outputter, numActualLoops)
    {
      std::vector<char> *p = new std::vector<char>(arraySize, 1);
      delete p;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);
    
#ifdef HAVE_TEUCHOS_BOOST
    // shared_ptr
    {
      typedef boost::shared_ptr<std::vector<char> > shared_ptr_t;
      shared_ptr_t p;
      TEUCHOS_START_PERF_OUTPUT_TIMER(outputter, numActualLoops)
      {
        p = shared_ptr_t(new std::vector<char>(arraySize, 1));
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, spTime);
#endif

    // RCP
    {
      RCP<std::vector<char> > p;
      TEUCHOS_START_PERF_OUTPUT_TIMER(outputter, numActualLoops)
      {
        p = rcp(new std::vector<char>(arraySize, 1));
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rcpTime);

#ifdef HAVE_TEUCHOS_BOOST
    // shared_ptr/rawPtr
    const double spRatio = spTime / rawPtrTime;
    outputter.outputField(spRatio);
#endif

    // RCP/rawPtr
    const double rcpRatio = rcpTime / rawPtrTime;
    outputter.outputField(rcpRatio);

    outputter.nextRow();
    
    arraySize *= 4;
    finalRcpRatio = rcpRatio;

  }

  out << "\n";
  TEST_COMPARE( finalRcpRatio, <=, maxRcpCreateDestroyRatio );
  out << "\n";

}


} // namespace
