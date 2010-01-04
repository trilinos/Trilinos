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
double maxRcpRawCreateDestroyRatio = 10.0;
double maxRcpRawAdjustRefCountRatio = 100.0;
double maxRcpSpAdjustRefCountRatio = 5.0;
double maxRcpRawObjAccessRatio = 10.0;


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
    "max-rcp-create-destroy-ratio", &maxRcpRawCreateDestroyRatio,
    "The ratio of the final CPU time ratio of creating and destroying"
    "std::vector<char>(size) objects wrapped in an RCP object versus"
    "using just raw new and delete."
    );
  clp.setOption(
    "max-rcp-raw-adjust-ref-count-ratio", &maxRcpRawAdjustRefCountRatio,
    "The ratio of the final CPU time ratio for adjusting the reference"
    "count of RCP objects versus a raw pointer."
    );
  clp.setOption(
    "max-rcp-sp-adjust-ref-count-ratio", &maxRcpSpAdjustRefCountRatio,
    "The ratio of the final CPU time ratio for adjusting the reference"
    "count of RCP objects versus boost::shared_ptr objects."
    );
  clp.setOption(
    "max-rcp-raw-obj-access-ratio", &maxRcpRawObjAccessRatio,
    "The ratio of the final CPU time ratio for accessing the object for RCP"
    "versus a raw pointer."
    );

}


TEUCHOS_UNIT_TEST( RCP, createDestroyOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const int maxLoopIters = 1000;
  const double relTestCost = 1e-3;
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

  double finalRcpRawRatio = 100000.0;

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
    {
      std::vector<std::vector<char>*> p_raw_vec(numActualLoops);
      int i = 0;
      TEUCHOS_START_PERF_OUTPUT_TIMER(outputter, numActualLoops)
      {
        p_raw_vec[i] = new std::vector<char>(arraySize, 1);
        delete p_raw_vec[i];
        ++i;
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);
    
#ifdef HAVE_TEUCHOS_BOOST
    // shared_ptr
    {
      typedef boost::shared_ptr<std::vector<char> > shared_ptr_t;
      std::vector<shared_ptr_t > sp_vec(numActualLoops);
      int i = 0;
      TEUCHOS_START_PERF_OUTPUT_TIMER(outputter, numActualLoops)
      {
        sp_vec[i] = shared_ptr_t(new std::vector<char>(arraySize, 1));
        sp_vec[i].reset();
        ++i;
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, spTime);
#endif

    // RCP
    {
      std::vector<RCP<std::vector<char> > > p_vec(numActualLoops);
      int i = 0;
      TEUCHOS_START_PERF_OUTPUT_TIMER(outputter, numActualLoops)
      {
        p_vec[i] = rcp(new std::vector<char>(arraySize, 1));
        p_vec[i] = null;
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
    finalRcpRawRatio = TEUCHOS_MIN(rcpRatio, finalRcpRawRatio);

  }

  out << "\n";
  TEST_COMPARE( finalRcpRawRatio, <=, maxRcpRawCreateDestroyRatio );
  out << "\n";

}


TEUCHOS_UNIT_TEST( RCP, referenceCountManipulationOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 5e-3;
  const int maxLoopIters = 1000;
  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Messuring the overhead of incrementing and deincrementing the reference count\n"
      << "comparing RCP to raw pointer and boost::shared_ptr.\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, 8);
  outputter.setFieldTypePrecision(TO::INT, 8);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw", TO::DOUBLE);
  outputter.pushFieldSpec("shared_ptr", TO::DOUBLE);
  outputter.pushFieldSpec("RCP", TO::DOUBLE);
  outputter.pushFieldSpec("RCP/raw", TO::DOUBLE);
  outputter.pushFieldSpec("RCP/shared_ptr", TO::DOUBLE);

  outputter.outputHeader();

  double finalRcpRawRatio = 100000.0;
  double finalRcpSpRatio = 100000.0;

  int arraySize = 64;

  for (
    int test_case_k = 0;
    test_case_k < maxLoopIters && arraySize <= maxArraySize;
    ++test_case_k
    )
  {

    // array dim
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
    {
      char dummy_char = 'n';
      std::vector<char*> p_raw_vec(arraySize);
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          p_raw_vec[i] = &dummy_char;
          p_raw_vec[i] = 0;
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);
    
#ifdef HAVE_TEUCHOS_BOOST
    // shared_ptr
    {
      typedef boost::shared_ptr<char> shared_ptr_t;
      shared_ptr_t sp(new char('n'));
      std::vector<shared_ptr_t> sp_vec(arraySize);
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          sp_vec[i] = sp;
          sp_vec[i].reset();
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, spTime);
#else
    outputter.outputField("-");
#endif

    // RCP
    {
      RCP<char> p(new char('n'));
      std::vector<RCP<char> > p_vec(arraySize);
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          p_vec[i] = p;
          p_vec[i].reset();
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rcpTime);

    // RCP/raw
    const double rcpRawRatio = rcpTime / rawPtrTime;
    finalRcpRawRatio = TEUCHOS_MIN(rcpRawRatio, finalRcpRawRatio);
    outputter.outputField(rcpRawRatio);

#ifdef HAVE_TEUCHOS_BOOST
    // RCP/shared_ptr
    const double rcpSpRatio = rcpTime / spTime;
    finalRcpSpRatio = TEUCHOS_MIN(rcpSpRatio, finalRcpSpRatio);
    outputter.outputField(rcpSpRatio);
#else
    outputter.outputField("-");
#endif

    outputter.nextRow();
    
    arraySize *= 4;

  }

  out << "\n";
  TEST_COMPARE( finalRcpRawRatio, <=, maxRcpRawAdjustRefCountRatio );
#ifdef HAVE_TEUCHOS_BOOST
  out << "\n";
  TEST_COMPARE( finalRcpSpRatio, <=, maxRcpSpAdjustRefCountRatio );
  out << "\n";
#else
  (void)finalRcpSpRatio;
#endif
  
}


TEUCHOS_UNIT_TEST( RCP, dereferenceOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;
  const int maxLoopIters = 1000;
  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Messuring the overhead of dereferencing RCP, shared_ptr and a raw pointer.\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, 8);
  outputter.setFieldTypePrecision(TO::INT, 8);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw", TO::DOUBLE);
  outputter.pushFieldSpec("shared_ptr", TO::DOUBLE);
  outputter.pushFieldSpec("RCP", TO::DOUBLE);
  outputter.pushFieldSpec("RCP/raw", TO::DOUBLE);
  outputter.pushFieldSpec("RCP/shared_ptr", TO::DOUBLE);

  outputter.outputHeader();

  double finalRcpRawRatio = 100000.0;

  int arraySize = 64;

  for (
    int test_case_k = 0;
    test_case_k < maxLoopIters && arraySize <= maxArraySize;
    ++test_case_k
    )
  {

    // array dim
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
    {
      char dummy_char = 'n';
      std::vector<char*> p_raw_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        p_raw_vec[i] = &dummy_char;
      }
      char dummy_char_out = '\0';
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_char_out = *p_raw_vec[i];
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);
    
#ifdef HAVE_TEUCHOS_BOOST
    // shared_ptr
    {
      typedef boost::shared_ptr<char> shared_ptr_t;
      shared_ptr_t sp(new char('n'));
      std::vector<shared_ptr_t> sp_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        sp_vec[i] = sp;
      }
      char dummy_char_out = '\0';
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_char_out = *sp_vec[i];
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, spTime);
#else
    outputter.outputField("-");
#endif

    // RCP
    {
      RCP<char> p(new char('n'));
      std::vector<RCP<char> > p_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        p_vec[i] = p;
      }
      char dummy_char_out = '\0';
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_char_out = *p_vec[i];
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rcpTime);

    // RCP/raw
    const double rcpRawRatio = rcpTime / rawPtrTime;
    finalRcpRawRatio = TEUCHOS_MIN(rcpRawRatio, finalRcpRawRatio);
    outputter.outputField(rcpRawRatio);

#ifdef HAVE_TEUCHOS_BOOST
    // RCP/shared_ptr
    const double rcpSpRatio = rcpTime / spTime;
    outputter.outputField(rcpSpRatio);
#else
    outputter.outputField("-");
#endif

    outputter.nextRow();
    
    arraySize *= 4;

  }

  out << "\n";
  TEST_COMPARE( finalRcpRawRatio, <=, maxRcpRawObjAccessRatio );
  
}


struct SomeStruct {
  SomeStruct(char member_in) : member(member_in) {}
  char member;
};


TEUCHOS_UNIT_TEST( RCP, memberAccessOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;
  const int maxLoopIters = 1000;
  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Messuring the overhead of dereferencing RCP, shared_ptr and a raw pointer.\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, 8);
  outputter.setFieldTypePrecision(TO::INT, 8);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw", TO::DOUBLE);
  outputter.pushFieldSpec("shared_ptr", TO::DOUBLE);
  outputter.pushFieldSpec("RCP", TO::DOUBLE);
  outputter.pushFieldSpec("RCP/raw", TO::DOUBLE);
  outputter.pushFieldSpec("RCP/shared_ptr", TO::DOUBLE);

  outputter.outputHeader();

  double finalRcpRawRatio = 100000.0;

  int arraySize = 64;

  for (
    int test_case_k = 0;
    test_case_k < maxLoopIters && arraySize <= maxArraySize;
    ++test_case_k
    )
  {

    // array dim
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
    {
      SomeStruct dummy_SomeStruct('n');
      std::vector<SomeStruct*> p_raw_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        p_raw_vec[i] = &dummy_SomeStruct;
      }
      char dummy_SomeStruct_out = '\0';
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_SomeStruct_out = p_raw_vec[i]->member;
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);
    
#ifdef HAVE_TEUCHOS_BOOST
    // shared_ptr
    {
      typedef boost::shared_ptr<SomeStruct> shared_ptr_t;
      shared_ptr_t sp(new SomeStruct('n'));
      std::vector<shared_ptr_t> sp_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        sp_vec[i] = sp;
      }
      char dummy_SomeStruct_out = '\0';
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_SomeStruct_out = sp_vec[i]->member;
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, spTime);
#else
    outputter.outputField("-");
#endif

    // RCP
    {
      RCP<SomeStruct> p(new SomeStruct('n'));
      std::vector<RCP<SomeStruct> > p_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        p_vec[i] = p;
      }
      char dummy_SomeStruct_out = '\0';
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_SomeStruct_out = p_vec[i]->member;
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rcpTime);

    // RCP/raw
    const double rcpRawRatio = rcpTime / rawPtrTime;
    finalRcpRawRatio = TEUCHOS_MIN(rcpRawRatio, finalRcpRawRatio);
    outputter.outputField(rcpRawRatio);

#ifdef HAVE_TEUCHOS_BOOST
    // RCP/shared_ptr
    const double rcpSpRatio = rcpTime / spTime;
    outputter.outputField(rcpSpRatio);
#else
    outputter.outputField("-");
#endif

    outputter.nextRow();
    
    arraySize *= 4;

  }

  out << "\n";
  TEST_COMPARE( finalRcpRawRatio, <=, maxRcpRawObjAccessRatio );
  
}







} // namespace
