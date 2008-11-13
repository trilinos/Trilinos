#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TabularOutputter.hpp"

#ifdef HAVE_TEUCHOS_BOOST
#  include "boost/shared_ptr.hpp"
#endif


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TabularOutputter;


double relCpuSpeed = 1e-2;
int maxArraySize = 10000;
double maxRcpCreateDestroyRatio = 10.0;
double maxRcpAjustRefCountRatio = 10.0;


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
  clp.setOption(
    "max-rcp-adjust-ref-count-ratio", &maxRcpAjustRefCountRatio,
    "The ratio of the final CPU time ratio for adjusting the reference"
    "count of RCP objects verses boost::shared_ptr objects."
    );

}


TEUCHOS_UNIT_TEST( RCP, createDestroyOverhead )
{

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

  outputter.pushFieldSpec("obj size", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw", TO::DOUBLE);
#ifdef HAVE_TEUCHOS_BOOST
  outputter.pushFieldSpec("shared_ptr", TO::DOUBLE);
#endif
  outputter.pushFieldSpec("RCP", TO::DOUBLE);
#ifdef HAVE_TEUCHOS_BOOST
  outputter.pushFieldSpec("shared_ptr/raw", TO::DOUBLE);
#endif
  outputter.pushFieldSpec("RCP/raw", TO::DOUBLE);

  outputter.outputHeader();

  double finalRcpRatio = 100000.0;

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
    finalRcpRatio = TEUCHOS_MIN(rcpRatio, finalRcpRatio);

  }

  out << "\n";
  TEST_COMPARE( finalRcpRatio, <=, maxRcpCreateDestroyRatio );
  out << "\n";

}


#ifdef HAVE_TEUCHOS_BOOST


TEUCHOS_UNIT_TEST( RCP, referenceCountManipulationOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-5;

  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Messuring the overhead of manipuliating the reference count by\n"
      << "comparing shared_ptr and using RCP.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, 8);
  outputter.setFieldTypePrecision(TO::INT, 8);

  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("shared_ptr", TO::DOUBLE);
  outputter.pushFieldSpec("RCP", TO::DOUBLE);
  outputter.pushFieldSpec("RCP/shared_ptr", TO::DOUBLE);

  outputter.outputHeader();

  double finalRcpRatio = -1.0;

  const int arraySize = 1;
  
  // num loops
  const int numActualLoops = static_cast<int>(numInnerLoops);
  outputter.outputField(numActualLoops);
  
  // shared_ptr
  {
    typedef boost::shared_ptr<std::vector<char> > shared_ptr_t;
    shared_ptr_t p
      = shared_ptr_t(new std::vector<char>(arraySize, 1));
    shared_ptr_t p2;
    TEUCHOS_START_PERF_OUTPUT_TIMER(outputter, numActualLoops)
    {
      p2 = p;
      p2 = shared_ptr_t();
    }
  }
  TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, spTime);
  
  // RCP
  {
    RCP<std::vector<char> > p
      = rcp(new std::vector<char>(arraySize, 1));
    RCP<std::vector<char> > p2;
    TEUCHOS_START_PERF_OUTPUT_TIMER(outputter, numActualLoops)
    {
      p2 = p;
      p2 = null;
    }
  }
  TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rcpTime);
  
  // RCP/shared_ptr
  const double rcpRatio = rcpTime / spTime;
  outputter.outputField(rcpRatio);
  
  outputter.nextRow();

  finalRcpRatio = rcpRatio;

  out << "\n";
  TEST_COMPARE( finalRcpRatio, <=, maxRcpAjustRefCountRatio );
  out << "\n";
  
}


#endif // HAVE_TEUCHOS_BOOST


} // namespace
