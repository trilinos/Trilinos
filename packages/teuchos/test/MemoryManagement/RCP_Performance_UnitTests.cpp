/*
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
*/

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
#ifdef HAVE_TEUCHOS_BOOST
  TEST_INEQUALITY_CONST(sizeof(boost::detail::shared_count), 0);
  TEST_INEQUALITY_CONST(sizeof(boost::shared_ptr<std::vector<double> >), 0);
  TEST_INEQUALITY_CONST(sizeof(boost::detail::sp_counted_impl_p<std::vector<double> >), 0);
  TEST_INEQUALITY_CONST(
    sizeof(boost::detail::sp_counted_impl_pd<std::vector<double>,
      DeleteDeleter<std::vector<double> > >),
    0);
#endif
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
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

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
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

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
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

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
    
    // shared_ptr
#ifdef HAVE_TEUCHOS_BOOST
    {
      typedef boost::shared_ptr<int> shared_ptr_t;
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
#else
    outputter.outputField("-");
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
  out << "\n";

  // This silly varible must be accumulated or compilers like MSVC++ will
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
      << "Messuring the overhead of dereferencing RCP, shared_ptr and a raw pointer.\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

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
    
    // shared_ptr
#ifdef HAVE_TEUCHOS_BOOST
    {
      typedef boost::shared_ptr<SomeStruct> shared_ptr_t;
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
#else
    outputter.outputField("-");
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
  out << "\n";

  // This silly varible must be accumulated or compilers like MSVC++ will
  // optimize away the loops!
  if (overall_dummy_int_out == 0)
    success = false;
  
}







} // namespace
