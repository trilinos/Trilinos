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
#include "Teuchos_TabularOutputter.hpp"


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TabularOutputter;


TEUCHOS_UNIT_TEST( TabularOutputter, basic1 )
{

  typedef Teuchos::TabularOutputter TO;

  std::stringstream sout;
  sout << "\n";

  TabularOutputter outputter(sout);

  outputter.pushFieldSpec("very long col name", TO::INT);
  outputter.pushFieldSpec("col b", TO::DOUBLE);
  outputter.pushFieldSpec("col cc", TO::STRING, TO::LEFT, TO::GENERAL, 6);
  outputter.pushFieldSpec("col d", TO::DOUBLE);
  outputter.pushFieldSpec("col e", TO::STRING);

  outputter.outputHeader();

  outputter.outputField(1);
  outputter.outputField(1.2);
  outputter.outputField("s13");
  outputter.outputField(1.4);
  outputter.outputField("s15");
  outputter.nextRow();

  outputter.outputField(2);
  outputter.outputField(2.2);
  outputter.outputField("s23");
  outputter.outputField(2.4);
  outputter.outputField("s25");
  outputter.nextRow();

  outputter.outputField(3);
  outputter.outputField(3.2);
  outputter.outputField("s33");
  outputter.outputField(3.4);
  outputter.outputField("s35");
  outputter.nextRow();

  std::stringstream expectedOutput;
  expectedOutput
    << "\n"
    << "  very long col name  col b         col cc  col d         col e\n"
    << "  ------------------  ------------  ------  ------------  -----\n"
    << "                   1    1.2000e+00  s13       1.4000e+00    s15\n"
    << "                   2    2.2000e+00  s23       2.4000e+00    s25\n"
    << "                   3    3.2000e+00  s33       3.4000e+00    s35\n"
    ;

  TEST_EQUALITY_CONST( sout.str(), expectedOutput.str() );

  // 2008/11/12: rabartl: Note: The above test may not be portable because it
  // requires the numeric formatting of the doubles to be the same.  To make
  // this more portable, I may have to do some work.

}


TEUCHOS_UNIT_TEST( TabularOutputter, basic2 )
{

  typedef Teuchos::TabularOutputter TO;

  std::stringstream sout;
  sout << "\n";

  TabularOutputter outputter(Teuchos::rcpFromRef(sout));

  outputter.setFieldTypePrecision(TO::DOUBLE, 8);
  outputter.setFieldTypePrecision(TO::INT, 4);
  outputter.setFieldTypePrecision(TO::STRING, 5);

  outputter.pushFieldSpec("col a", TO::INT);
  outputter.pushFieldSpec("col b", TO::DOUBLE);
  outputter.pushFieldSpec("col cc", TO::STRING, TO::LEFT, TO::GENERAL, 6);
  outputter.pushFieldSpec("col d", TO::DOUBLE);
  outputter.pushFieldSpec("col e", TO::STRING);

  outputter.outputHeader();

  outputter.outputField(1);
  outputter.outputField(1.2);
  outputter.outputField("s13");
  outputter.outputField(1.4);
  outputter.outputField("s15");
  outputter.nextRow();

  outputter.outputField(2);
  outputter.outputField(2.2);
  outputter.outputField("s23");
  outputter.outputField(2.4);
  outputter.outputField("s25");
  outputter.nextRow();

  outputter.outputField(3);
  outputter.outputField(3.2);
  outputter.outputField("s33");
  outputter.outputField(3.4);
  outputter.outputField("s35");
  outputter.nextRow();

  std::stringstream expectedOutput;
  expectedOutput
    << "\n"
    << "  col a  col b             col cc  col d             col e\n"
    << "  -----  ----------------  ------  ----------------  -----\n"
    << "      1    1.20000000e+00  s13       1.40000000e+00    s15\n"
    << "      2    2.20000000e+00  s23       2.40000000e+00    s25\n"
    << "      3    3.20000000e+00  s33       3.40000000e+00    s35\n"
    ;

  TEST_EQUALITY_CONST( sout.str(), expectedOutput.str() );

  // 2008/11/12: rabartl: Note: See the comment in the basic1 test above!

}


TEUCHOS_UNIT_TEST( TabularOutputter, perfTiming )
{

  typedef Teuchos::TabularOutputter TO;

  std::stringstream sout;
  sout << "\n";

  TabularOutputter outputter(sout);

  outputter.pushFieldSpec("num loops", TO::INT);
  outputter.pushFieldSpec("vecTime", TO::DOUBLE);
  outputter.pushFieldSpec("dequeTime", TO::DOUBLE);

  outputter.outputHeader();

  const int numLoops = 15;

  // num loops
  outputter.outputField(numLoops);

  // vecTime
  TEUCHOS_START_PERF_OUTPUT_TIMER(outputter, numLoops)
  {
    std::vector<int> a(numLoops);
    std::vector<int> b(numLoops);
    a = b;
  }
  TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, vecTime);
  TEST_INEQUALITY_CONST(vecTime, 0);

  // dequeTime
  TEUCHOS_START_PERF_OUTPUT_TIMER(outputter, numLoops)
  {
    std::deque<int> a(numLoops);
    std::deque<int> b(numLoops);
    a = b;
  }
  TEUCHOS_END_PERF_OUTPUT_TIMER(outputter, dequeTime);
  TEST_INEQUALITY_CONST(dequeTime, 0);

  outputter.nextRow();

  std::stringstream expectedOutput;
  expectedOutput
    << "\n"
    << "Nothing\n"
    ;

  TEST_INEQUALITY_CONST( sout.str(), expectedOutput.str() );

  // 2008/11/12: rabartl: Above, this is not the greatest test but it would be
  // hard to produce the exact same formatted output since it involves timing
  // results.

}


#ifdef TEUCHOS_DEBUG


TEUCHOS_UNIT_TEST( TabularOutputter, nullOStream )
{

  typedef Teuchos::TabularOutputter TO;

  TabularOutputter outputter(out);

  TEST_THROW(
    outputter.setOStream(Teuchos::null),
    Teuchos::NullReferenceError
    );

}


TEUCHOS_UNIT_TEST( TabularOutputter, invalidFieldSpecError )
{

  typedef Teuchos::TabularOutputter TO;

  TabularOutputter outputter(out);

  outputter.setFieldTypePrecision(TO::DOUBLE, 8);
  outputter.setFieldTypePrecision(TO::INT, 4);
  outputter.setFieldTypePrecision(TO::STRING, 3);

  outputter.pushFieldSpec("col d", TO::DOUBLE);

  TEST_THROW(
    outputter.pushFieldSpec(
      "very long field name", TO::INT, TO::LEFT, TO::GENERAL, 4),
    TO::InvalidFieldSpecError
    );

}


TEUCHOS_UNIT_TEST( TabularOutputter, missingHeaderError )
{

  typedef Teuchos::TabularOutputter TO;

  TabularOutputter outputter(out);

  outputter.pushFieldSpec("col a", TO::INT);
  outputter.pushFieldSpec("col b", TO::DOUBLE);
  outputter.pushFieldSpec("col c", TO::STRING);
  outputter.pushFieldSpec("col d", TO::DOUBLE);

  TEST_THROW(outputter.outputField(1), TO::MissingHeaderError);

}


TEUCHOS_UNIT_TEST( TabularOutputter, missingNextRowError )
{

  typedef Teuchos::TabularOutputter TO;

  TabularOutputter outputter(out);

  outputter.pushFieldSpec("col a", TO::INT);
  outputter.pushFieldSpec("col b", TO::DOUBLE);
  outputter.pushFieldSpec("col c", TO::STRING);
  outputter.pushFieldSpec("col d", TO::DOUBLE);

  outputter.outputHeader();

  outputter.outputField(1);
  outputter.outputField(1.2);
  outputter.outputField("s13");
  outputter.outputField(1.4);

  // Missing nextRow()!

  TEST_THROW(outputter.outputField(2), TO::InvalidFieldOutputError);

}


TEUCHOS_UNIT_TEST( TabularOutputter, missingFieldOutputError )
{

  typedef Teuchos::TabularOutputter TO;

  TabularOutputter outputter(out);

  outputter.pushFieldSpec("col a", TO::INT);
  outputter.pushFieldSpec("col b", TO::DOUBLE);
  outputter.pushFieldSpec("col c", TO::STRING);
  outputter.pushFieldSpec("col d", TO::DOUBLE);

  outputter.outputHeader();

  outputter.outputField(1);
  outputter.outputField(1.2);
  outputter.outputField("s13");

  // Missing a call to outputField(...);
  
  out << "\n\n";

  TEST_THROW(outputter.nextRow(), TO::InvalidFieldOutputError);

}


TEUCHOS_UNIT_TEST( TabularOutputter, missingFieldOutputOkay )
{

  typedef Teuchos::TabularOutputter TO;

  TabularOutputter outputter(out);

  outputter.pushFieldSpec("col a", TO::INT);
  outputter.pushFieldSpec("col b", TO::DOUBLE);
  outputter.pushFieldSpec("col c", TO::STRING);
  outputter.pushFieldSpec("col d", TO::DOUBLE);

  outputter.outputHeader();

  outputter.outputField(1);
  outputter.outputField(1.2);
  outputter.outputField("s13");

  // Missing a call to outputField(...);

  outputter.nextRow(true); // Just fine!

}


TEUCHOS_UNIT_TEST( TabularOutputter, missingFields )
{

  typedef Teuchos::TabularOutputter TO;

  std::ostringstream sout;
  TabularOutputter outputter(sout);

  TEST_THROW(outputter.outputHeader(), TO::MissingFieldsError);

}


#endif // TEUCHOS_DEBUG


} // namespace
