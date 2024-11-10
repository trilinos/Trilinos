// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TabularOutputter.hpp"
#ifdef HAVE_TEUCHOSCORE_CXX11
#  include <memory>
#endif // HAVE_TEUCHOSCORE_CXX11

namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TabularOutputter;


double relCpuSpeed = 1e-2;
int maxArraySize = 10000;
double maxRcpRawCreateDestroyRatio = 10.0;
double maxRcpRawAdjustRefCountRatio = 100.0;
#ifdef HAVE_TEUCHOSCORE_CXX11
double maxRcpSpAdjustRefCountRatio = 5.0;
#endif
double maxRcpRawObjAccessRatio = 13.5;

const int intPrec = 8;
const int dblPrec = 6;


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
#ifdef HAVE_TEUCHOSCORE_CXX11
  clp.setOption(
    "max-rcp-sp-adjust-ref-count-ratio", &maxRcpSpAdjustRefCountRatio,
    "The ratio of the final CPU time ratio for adjusting the reference"
    "count of RCP objects versus std::shared_ptr objects."
    );
#endif
  clp.setOption(
    "max-rcp-raw-obj-access-ratio", &maxRcpRawObjAccessRatio,
    "The ratio of the final CPU time ratio for accessing the object for RCP"
    "versus a raw pointer."
    );

}


template<typename T>
struct DeleteDeleter {};


TEUCHOS_UNIT_TEST( RCP, _sizeofObjects )
{
  out << "\nPrinting the size the RCP and RCPNodeImpl objects ...\n";
  TEST_INEQUALITY_CONST(sizeof(bool), 0);
  TEST_INEQUALITY_CONST(sizeof(double), 0);
  TEST_INEQUALITY_CONST(sizeof(double*), 0);
  TEST_INEQUALITY_CONST(sizeof(std::vector<double>), 0);
  TEST_INEQUALITY_CONST(sizeof(Teuchos::RCPNode*), 0);
  TEST_INEQUALITY_CONST(sizeof(Teuchos::ERCPStrength), 0);
  TEST_INEQUALITY_CONST(sizeof(Teuchos::RCPNodeHandle), 0);
  TEST_INEQUALITY_CONST(sizeof(Teuchos::RCP<std::vector<double> >), 0);
  TEST_INEQUALITY_CONST(
    sizeof(Teuchos::RCPNodeTmpl<std::vector<double>,
      Teuchos::DeallocDelete<std::vector<double> > >),
    0);
}


TEUCHOS_UNIT_TEST( RCP, createDestroyOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const int maxLoopIters = 1000;
  const double relTestCost = 1e-3;
  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Messuring the overhead of creating and destorying objects of different sizes\n"
      << "using raw C++ pointers,"
#ifdef HAVE_TEUCHOSCORE_CXX11
      << " shared_ptr,"
#endif
      << " and using RCP.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("obj size", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw", TO::DOUBLE);
#ifdef HAVE_TEUCHOSCORE_CXX11
  outputter.pushFieldSpec("shared_ptr", TO::DOUBLE);
#endif
  outputter.pushFieldSpec("RCP", TO::DOUBLE);
#ifdef HAVE_TEUCHOSCORE_CXX11
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

#ifdef HAVE_TEUCHOSCORE_CXX11
    // shared_ptr
    {
      typedef std::shared_ptr<std::vector<char> > shared_ptr_t;
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

#ifdef HAVE_TEUCHOSCORE_CXX11
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
      << "comparing RCP to raw pointer"
#ifdef HAVE_TEUCHOSCORE_CXX11
      << " and std::shared_ptr"
#endif
      << ".\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw", TO::DOUBLE);
#ifdef HAVE_TEUCHOSCORE_CXX11
  outputter.pushFieldSpec("shared_ptr", TO::DOUBLE);
#endif
  outputter.pushFieldSpec("RCP", TO::DOUBLE);
  outputter.pushFieldSpec("RCP/raw", TO::DOUBLE);
#ifdef HAVE_TEUCHOSCORE_CXX11
  outputter.pushFieldSpec("RCP/shared_ptr", TO::DOUBLE);
#endif

  outputter.outputHeader();

  double finalRcpRawRatio = 100000.0;
#ifdef HAVE_TEUCHOSCORE_CXX11
  double finalRcpSpRatio = 100000.0;
#endif
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

    // Note on std::shared_ptr and modification to the test
    // Originally this test copied a single ptr
    // Added 1 and 2 types ('n' and 'o') so that each copy would be unique
    // std::shared_ptr for gcc (but not clang) will handle the case of setting
    // a = b with b already equal to a in an optimized way and the original
    // test format spent most of it's time in this case.

    // raw
    {
      char dummy_char1 = 'n';
      char dummy_char2 = 'o'; // See above note for std::shared_ptr
      std::vector<char*> p_raw_vec(arraySize);
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          p_raw_vec[i] = &dummy_char1;
          p_raw_vec[i] = &dummy_char2; // See above note for std::shared_ptr
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);

#ifdef HAVE_TEUCHOSCORE_CXX11
    // shared_ptr
    {
      typedef std::shared_ptr<char> shared_ptr_t;
      shared_ptr_t sp1(new char('n'));
      shared_ptr_t sp2(new char('o')); // See above note for std::shared_ptr
      std::vector<shared_ptr_t> sp_vec(arraySize);
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          sp_vec[i] = sp1;
          sp_vec[i] = sp2; // See above note for std::shared_ptr
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, spTime);
#endif

    // RCP
    {
      RCP<char> p1(new char('n'));
      RCP<char> p2(new char('o')); // See above note for std::shared_ptr
      std::vector<RCP<char> > p_vec(arraySize);
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          p_vec[i] = p1;
          p_vec[i] = p2; // See above note for std::shared_ptr
          // NOTE: This assignment operation tests the copy constructor and
          // the swap function.  This calls both bind() and unbind()
          // underneath.
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rcpTime);

    // RCP/raw
    const double rcpRawRatio = rcpTime / rawPtrTime;
    finalRcpRawRatio = TEUCHOS_MIN(rcpRawRatio, finalRcpRawRatio);
    outputter.outputField(rcpRawRatio);

#ifdef HAVE_TEUCHOSCORE_CXX11
    // RCP/shared_ptr
    const double rcpSpRatio = rcpTime / spTime;
    finalRcpSpRatio = TEUCHOS_MIN(rcpSpRatio, finalRcpSpRatio);
    outputter.outputField(rcpSpRatio);
#endif

    outputter.nextRow();

    arraySize *= 4;

  }

  out << "\n";
  TEST_COMPARE( finalRcpRawRatio, <=, maxRcpRawAdjustRefCountRatio );
  out << "\n";
#ifdef HAVE_TEUCHOSCORE_CXX11
  TEST_COMPARE( finalRcpSpRatio, <=, maxRcpSpAdjustRefCountRatio );
  out << "\n";
#endif

}


TEUCHOS_UNIT_TEST( RCP, dereferenceOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;
  const int maxLoopIters = 1000;
  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Measuring the overhead of dereferencing RCP"
#ifdef HAVE_TEUCHOSCORE_CXX11
      << ", shared_ptr"
#endif
      << " and a raw pointer.\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw", TO::DOUBLE);
#ifdef HAVE_TEUCHOSCORE_CXX11
  outputter.pushFieldSpec("shared_ptr", TO::DOUBLE);
#endif
  outputter.pushFieldSpec("RCP", TO::DOUBLE);
  outputter.pushFieldSpec("RCP/raw", TO::DOUBLE);
#ifdef HAVE_TEUCHOSCORE_CXX11
  outputter.pushFieldSpec("RCP/shared_ptr", TO::DOUBLE);
#endif

  outputter.outputHeader();

  double finalRcpRawRatio = 100000.0;
  int arraySize = 64;
  const int dummy_int_val = 1;
  int overall_dummy_int_out = 0;


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

    int dummy_int_out = 0;

    // raw
    {
      int dummy_int = dummy_int_val;
      std::vector<int*> p_raw_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        p_raw_vec[i] = &dummy_int;
      }
      dummy_int_out = 0;
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_int_out += *p_raw_vec[i];
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);
    overall_dummy_int_out += dummy_int_out;

#ifdef HAVE_TEUCHOSCORE_CXX11
    // shared_ptr
    {
      typedef std::shared_ptr<int> shared_ptr_t;
      shared_ptr_t sp(new int(dummy_int_val));
      std::vector<shared_ptr_t> sp_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        sp_vec[i] = sp;
      }
      dummy_int_out = 0;
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_int_out += *sp_vec[i];
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, spTime);
    overall_dummy_int_out += dummy_int_out;
#endif

    // RCP
    {
      RCP<int> p(new int(dummy_int_val));
      std::vector<RCP<int> > p_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        p_vec[i] = p;
      }
      dummy_int_out = 0;
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_int_out += *p_vec[i];
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rcpTime);
    overall_dummy_int_out += dummy_int_out;

    // RCP/raw
    const double rcpRawRatio = rcpTime / rawPtrTime;
    finalRcpRawRatio = TEUCHOS_MIN(rcpRawRatio, finalRcpRawRatio);
    outputter.outputField(rcpRawRatio);

#ifdef HAVE_TEUCHOSCORE_CXX11
    // RCP/shared_ptr
    const double rcpSpRatio = rcpTime / spTime;
    outputter.outputField(rcpSpRatio);
#endif

    outputter.nextRow();

    arraySize *= 4;

  }

  out << "\n";
  TEST_COMPARE( finalRcpRawRatio, <=, maxRcpRawObjAccessRatio );
  out << "\n";

  // This silly variable must be accumulated or compilers like MSVC++ will
  // optimize away the loops!
  if (overall_dummy_int_out == 0)
    success = false;

}


struct SomeStruct {
  SomeStruct(int member_in) : member(member_in) {}
  int member;
};


TEUCHOS_UNIT_TEST( RCP, memberAccessOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;
  const int maxLoopIters = 1000;
  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Measuring the overhead of dereferencing RCP"
#ifdef HAVE_TEUCHOSCORE_CXX11
      << ", shared_ptr"
#endif
      << " and a raw pointer.\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw", TO::DOUBLE);
#ifdef HAVE_TEUCHOSCORE_CXX11
  outputter.pushFieldSpec("shared_ptr", TO::DOUBLE);
#endif
  outputter.pushFieldSpec("RCP", TO::DOUBLE);
  outputter.pushFieldSpec("RCP/raw", TO::DOUBLE);
#ifdef HAVE_TEUCHOSCORE_CXX11
  outputter.pushFieldSpec("RCP/shared_ptr", TO::DOUBLE);
#endif

  outputter.outputHeader();

  double finalRcpRawRatio = 100000.0;
  int arraySize = 64;
  const int dummy_int_val = 1;
  int overall_dummy_int_out = 0;

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

    int dummy_int_out = 0;

    // raw
    {
      SomeStruct dummy_SomeStruct(dummy_int_val);
      std::vector<SomeStruct*> p_raw_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        p_raw_vec[i] = &dummy_SomeStruct;
      }
      dummy_int_out = 0;
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_int_out += p_raw_vec[i]->member;
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);
    overall_dummy_int_out += dummy_int_out;

#ifdef HAVE_TEUCHOSCORE_CXX11
    // shared_ptr
    {
      typedef std::shared_ptr<SomeStruct> shared_ptr_t;
      shared_ptr_t sp(new SomeStruct(dummy_int_val));
      std::vector<shared_ptr_t> sp_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        sp_vec[i] = sp;
      }
      dummy_int_out = 0;
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_int_out += sp_vec[i]->member;
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, spTime);
    overall_dummy_int_out += dummy_int_out;
#endif

    // RCP
    {
      RCP<SomeStruct> p(new SomeStruct(dummy_int_val));
      std::vector<RCP<SomeStruct> > p_vec(arraySize);
      for (int i=0; i < arraySize; ++i) {
        p_vec[i] = p;
      }
      dummy_int_out = 0;
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (int i=0; i < arraySize; ++i) {
          dummy_int_out += p_vec[i]->member;
        }
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rcpTime);
    overall_dummy_int_out += dummy_int_out;

    // RCP/raw
    const double rcpRawRatio = rcpTime / rawPtrTime;
    finalRcpRawRatio = TEUCHOS_MIN(rcpRawRatio, finalRcpRawRatio);
    outputter.outputField(rcpRawRatio);

#ifdef HAVE_TEUCHOSCORE_CXX11
    // RCP/shared_ptr
    const double rcpSpRatio = rcpTime / spTime;
    outputter.outputField(rcpSpRatio);
#endif

    outputter.nextRow();

    arraySize *= 4;

  }

  out << "\n";
  TEST_COMPARE( finalRcpRawRatio, <=, maxRcpRawObjAccessRatio );
  out << "\n";

  // This silly variable must be accumulated or compilers like MSVC++ will
  // optimize away the loops!
  if (overall_dummy_int_out == 0)
    success = false;

}







} // namespace
