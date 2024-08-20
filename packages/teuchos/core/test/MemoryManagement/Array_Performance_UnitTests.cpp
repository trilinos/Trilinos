// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TabularOutputter.hpp"

#include "Teuchos_Array.hpp"


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TabularOutputter;
using Teuchos::Ordinal;


double relCpuSpeed = 1e-2;
int maxArraySize = 10000;
double maxArrayBracketRatio =100.0;
double maxArrayIterRatio = 200.0;
double maxArrayRCPSelfIterRatio =200.0;

const int minArraySize = 100;
const int maxLoopIters = 1000;
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
    "max-array-bracket-ratio", &maxArrayBracketRatio,
    "The max allowed CPU timing ratio of the Array[RCP,View] braket operator relative"
    " to the std::vector braket operator."
    );
  clp.setOption(
    "max-array-iter-ratio", &maxArrayIterRatio,
    "The max allowed CPU timing ratio of the Array[RCP,View] iterators relative"
    " to using raw pointers as iterators."
    );
  clp.setOption(
    "max-arrayrcp-self-iter-ratio", &maxArrayRCPSelfIterRatio,
    "The max allowed CPU timing ratio of the ArrayrCP as a self iterator relative"
    " to raw pointer arithmetic."
    );
}


TEUCHOS_UNIT_TEST( Array, braketOperatorOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;

  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Measuring the overhead of the Array braket operator relative to raw pointers.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw ptr", TO::DOUBLE);
  outputter.pushFieldSpec("vector", TO::DOUBLE);
  outputter.pushFieldSpec("Array", TO::DOUBLE);
  outputter.pushFieldSpec("vector/raw", TO::DOUBLE);
  outputter.pushFieldSpec("Array/raw", TO::DOUBLE);

  outputter.outputHeader();

  // Start out really big to make sure it fails if not set correctly!
  double finalArrayBraketRatio = 100000.0;

  Ordinal arraySize = minArraySize;
  for (int test_case_k = 0;
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

    std::vector<double> vec(arraySize);

    // raw ptr
    {
      double *p_raw = &vec[0];
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (Ordinal i=0; i < arraySize; ++i)
          p_raw[i] = 0.0;
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);

    // vector
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      for (Ordinal i=0; i < arraySize; ++i)
        vec[i] = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, vectorTime);

    // Array
    {
      Teuchos::Array<double> a(arraySize);
      TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
      {
        for (Ordinal i=0; i < arraySize; ++i)
          a[i] = 0.0;
      }
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, arrayTime);

    // vector/raw
    const double vectorRatio = vectorTime / rawPtrTime;
    outputter.outputField(vectorRatio);

    // Array/raw
    const double arrayRatio = arrayTime / rawPtrTime;
    outputter.outputField(arrayRatio);

    outputter.nextRow();

    arraySize *= 4;
    finalArrayBraketRatio = TEUCHOS_MIN(arrayRatio, finalArrayBraketRatio);

  }

  out << "\n";
  TEST_COMPARE( finalArrayBraketRatio, <=, maxArrayBracketRatio );
  out << "\n";

}


TEUCHOS_UNIT_TEST( ArrayView, braketOperatorOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;

  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Measuring the overhead of the ArrayView braket operator relative to raw pointers.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw ptr", TO::DOUBLE);
  outputter.pushFieldSpec("ArrayView", TO::DOUBLE);
  outputter.pushFieldSpec("ArrayView/raw", TO::DOUBLE);

  outputter.outputHeader();

  // Start out really big to make sure it fails if not set correctly!
  double finalArrayViewBraketRatio = 100000.0;

  Ordinal arraySize = minArraySize;
  for (int test_case_k = 0;
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

    std::vector<double> vec(arraySize);

    // raw ptr
    double *p_raw = &vec[0];
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      for (Ordinal i=0; i < arraySize; ++i)
        p_raw[i] = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);

    // ArrayView
    Teuchos::Array<double> a(arraySize);
    Teuchos::ArrayView<double> av = a;
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      for (Ordinal i=0; i < arraySize; ++i)
        av[i] = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, arrayviewTime);

    // Array/raw
    const double arrayviewRatio = arrayviewTime / rawPtrTime;
    outputter.outputField(arrayviewRatio);

    outputter.nextRow();

    arraySize *= 4;
    finalArrayViewBraketRatio = TEUCHOS_MIN(arrayviewRatio, finalArrayViewBraketRatio);

  }

  out << "\n";
  TEST_COMPARE( finalArrayViewBraketRatio, <=, maxArrayBracketRatio );
  out << "\n";

}


TEUCHOS_UNIT_TEST( ArrayRCP, braketOperatorOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;

  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Measuring the overhead of the ArrayRCP braket operator relative to raw pointers.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw ptr", TO::DOUBLE);
  outputter.pushFieldSpec("ArrayRCP", TO::DOUBLE);
  outputter.pushFieldSpec("ArrayRCP/raw", TO::DOUBLE);

  outputter.outputHeader();

  // Start out really big to make sure it fails if not set correctly!
  double finalArrayRCPBraketRatio = 100000.0;

  Ordinal arraySize = minArraySize;
  for (int test_case_k = 0;
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

    std::vector<double> vec(arraySize);

    // raw ptr
    double *p_raw = &vec[0];
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      for (Ordinal i=0; i < arraySize; ++i)
        p_raw[i] = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);

    // ArrayRCP
    Teuchos::ArrayRCP<double> arcp = Teuchos::arcp<double>(arraySize);
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      for (Ordinal i=0; i < arraySize; ++i)
        arcp[i] = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, arrayrcpTime);

    // Array/raw
    const double arrayrcpRatio = arrayrcpTime / rawPtrTime;
    outputter.outputField(arrayrcpRatio);

    outputter.nextRow();

    arraySize *= 4;
    finalArrayRCPBraketRatio = TEUCHOS_MIN(arrayrcpRatio, finalArrayRCPBraketRatio);

  }

  out << "\n";
  TEST_COMPARE( finalArrayRCPBraketRatio, <=, maxArrayBracketRatio );
  out << "\n";

}


TEUCHOS_UNIT_TEST( Array, iteratorOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;

  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Measuring the overhead of the Array iterators relative to raw pointers.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw ptr", TO::DOUBLE);
  outputter.pushFieldSpec("vector", TO::DOUBLE);
  outputter.pushFieldSpec("Array", TO::DOUBLE);
  outputter.pushFieldSpec("vector/raw", TO::DOUBLE);
  outputter.pushFieldSpec("Array/raw", TO::DOUBLE);

  outputter.outputHeader();

  // Start out really big to make sure it fails if not set correctly!
  double finalArrayIterRatio = 100000.0;

  Ordinal arraySize = minArraySize;
  for (int test_case_k = 0;
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

    std::vector<double> vec(arraySize);

    // raw ptr
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      double
        *p_raw_itr = &vec[0],
        *p_raw_end = &vec[0] + arraySize;
      for ( ; p_raw_itr < p_raw_end; ++p_raw_itr)
        *p_raw_itr = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);

    // vector
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      std::vector<double>::iterator
        vec_itr = vec.begin(),
        vec_end = vec.end();
      for ( ; vec_itr < vec_end; ++vec_itr)
        *vec_itr = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, vectorTime);

    // Array
    Teuchos::Array<double> a(arraySize);
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      Teuchos::Array<double>::iterator
        a_itr = a.begin(),
        a_end = a.end();
      for ( ; a_itr < a_end; ++a_itr)
        *a_itr = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, arrayTime);

    // vector/raw
    const double vectorRatio = vectorTime / rawPtrTime;
    outputter.outputField(vectorRatio);

    // Array/raw
    const double arrayRatio = arrayTime / rawPtrTime;
    outputter.outputField(arrayRatio);

    outputter.nextRow();

    arraySize *= 4;
    finalArrayIterRatio = TEUCHOS_MIN(arrayRatio, finalArrayIterRatio);

  }

  out << "\n";
  TEST_COMPARE( finalArrayIterRatio, <=, maxArrayIterRatio );
  out << "\n";

}


TEUCHOS_UNIT_TEST( ArrayView, iteratorOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;

  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Measuring the overhead of the ArrayView iterators relative to raw pointers.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw ptr", TO::DOUBLE);
  outputter.pushFieldSpec("ArrayView", TO::DOUBLE);
  outputter.pushFieldSpec("ArrayView/raw", TO::DOUBLE);

  outputter.outputHeader();

  // Start out really big to make sure it fails if not set correctly!
  double finalArrayViewIterRatio = 100000.0;

  Ordinal arraySize = minArraySize;
  for (int test_case_k = 0;
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

    std::vector<double> vec(arraySize);

    // raw ptr
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      double
        *p_raw_itr = &vec[0],
        *p_raw_end = &vec[0] + arraySize;
      for ( ; p_raw_itr < p_raw_end; ++p_raw_itr)
        *p_raw_itr = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);

    // ArrayView
    Teuchos::Array<double> a(arraySize);
    Teuchos::ArrayView<double> av = a;
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      Teuchos::ArrayView<double>::iterator
        av_itr = av.begin(),
        av_end = av.end();
      for ( ; av_itr < av_end ; ++av_itr)
        *av_itr = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, arrayviewTime);

    // ArrayView/raw
    const double arrayviewRatio = arrayviewTime / rawPtrTime;
    outputter.outputField(arrayviewRatio);

    outputter.nextRow();

    arraySize *= 4;
    finalArrayViewIterRatio = TEUCHOS_MIN(arrayviewRatio, finalArrayViewIterRatio);

  }

  out << "\n";
  TEST_COMPARE( finalArrayViewIterRatio, <=, maxArrayIterRatio );
  out << "\n";

}


TEUCHOS_UNIT_TEST( ArrayRCP, iteratorOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;

  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Measuring the overhead of the ArrayRCP iterators relative to raw pointers.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw ptr", TO::DOUBLE);
  outputter.pushFieldSpec("ArrayRCP", TO::DOUBLE);
  outputter.pushFieldSpec("ArrayRCP/raw", TO::DOUBLE);

  outputter.outputHeader();

  // Start out really big to make sure it fails if not set correctly!
  double finalArrayRCPIterRatio = 100000.0;

  Ordinal arraySize = minArraySize;
  for (int test_case_k = 0;
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

    std::vector<double> vec(arraySize);

    // raw ptr
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      double
        *p_raw_itr = &vec[0],
        *p_raw_end = &vec[0] + arraySize;
      for ( ; p_raw_itr < p_raw_end; ++p_raw_itr)
        *p_raw_itr = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);

    // ArrayRCP
    Teuchos::ArrayRCP<double> ap = Teuchos::arcp<double>(arraySize);
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      Teuchos::ArrayRCP<double>::iterator
        ap_itr = ap.begin(),
        ap_end = ap.end();
      for ( ; ap_itr < ap_end; ++ap_itr)
        *ap_itr = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, arrayviewTime);

    // ArrayRCP/raw
    const double arrayviewRatio = arrayviewTime / rawPtrTime;
    outputter.outputField(arrayviewRatio);

    outputter.nextRow();

    arraySize *= 4;
    finalArrayRCPIterRatio = TEUCHOS_MIN(arrayviewRatio, finalArrayRCPIterRatio);

  }

  out << "\n";
  TEST_COMPARE( finalArrayRCPIterRatio, <=, maxArrayIterRatio );
  out << "\n";

}


TEUCHOS_UNIT_TEST( ArrayRCP, selfIteratorOverhead )
{

  typedef Teuchos::TabularOutputter TO;

  const double relTestCost = 1e-4;

  const double numInnerLoops = relCpuSpeed / relTestCost;

  out << "\n"
      << "Measuring the overhead of the ArrayRCP as a self iterataor relative to raw pointers.\n"
      << "\n"
      << "Number of loops = relCpuSpeed/relTestCost = "
      << relCpuSpeed << "/" << relTestCost << " = " << numInnerLoops << "\n"
      << "\n";

  TabularOutputter outputter(out);
  outputter.setFieldTypePrecision(TO::DOUBLE, dblPrec);
  outputter.setFieldTypePrecision(TO::INT, intPrec);

  outputter.pushFieldSpec("array dim", TO::INT);
  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("raw ptr", TO::DOUBLE);
  outputter.pushFieldSpec("ArrayRCP", TO::DOUBLE);
  outputter.pushFieldSpec("ArrayRCP/raw", TO::DOUBLE);

  outputter.outputHeader();

  // Start out really big to make sure it fails if not set correctly!
  double finalArrayRCPIterRatio = 100000.0;

  Ordinal arraySize = minArraySize;
  for (int test_case_k = 0;
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

    std::vector<double> vec(arraySize);

    // raw ptr
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      double
        *p_raw_itr = &vec[0],
        *p_raw_end = &vec[0] + arraySize;
      for ( ; p_raw_itr < p_raw_end; ++p_raw_itr)
        *p_raw_itr = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, rawPtrTime);

    // ArrayRCP
    Teuchos::ArrayRCP<double> ap = Teuchos::arcp<double>(arraySize);
    TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(outputter, numActualLoops, arraySize)
    {
      Teuchos::ArrayRCP<double>
        ap_itr = ap,
        ap_end = ap + arraySize;
      for ( ; ap_itr < ap_end; ++ap_itr)
        *ap_itr = 0.0;
    }
    TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, arrayviewTime);

    // ArrayRCP/raw
    const double arrayviewRatio = arrayviewTime / rawPtrTime;
    outputter.outputField(arrayviewRatio);

    outputter.nextRow();

    arraySize *= 4;
    finalArrayRCPIterRatio = TEUCHOS_MIN(arrayviewRatio, finalArrayRCPIterRatio);

  }

  out << "\n";
  TEST_COMPARE( finalArrayRCPIterRatio, <=, maxArrayRCPSelfIterRatio );
  out << "\n";

}


} // namespace
