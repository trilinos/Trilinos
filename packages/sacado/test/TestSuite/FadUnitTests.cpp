// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "FadUnitTests.hpp"
#include "FadUnitTests2.hpp"

#include "Sacado_Fad_SimpleFad.hpp"

void FAD::error(const char *msg) {
  std::cout << msg << endl;
}

typedef FadOpsUnitTest<Sacado::Fad::DFad<double>,double> DFadDoubleTest;
typedef FadOpsUnitTest<Sacado::Fad::SFad<double,5>,double> SFadDoubleTest;
typedef FadOpsUnitTest<Sacado::Fad::SLFad<double,10>,double> SLFadDoubleTest;
typedef FadOpsUnitTest<Sacado::Fad::SimpleFad<double>,double> SimpleFadDoubleTest;

typedef RealFadOpsUnitTest2<Sacado::Fad::DFad<double>,double> DFadDoubleTest2;
typedef RealFadOpsUnitTest2<Sacado::Fad::SFad<double,5>,double> SFadDoubleTest2;
typedef RealFadOpsUnitTest2<Sacado::Fad::SLFad<double,10>,double> SLFadDoubleTest2;
typedef RealFadOpsUnitTest2<Sacado::Fad::SimpleFad<double>,double> SimpleFadDoubleTest2;

CPPUNIT_TEST_SUITE_REGISTRATION(DFadDoubleTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SFadDoubleTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SLFadDoubleTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SimpleFadDoubleTest);

CPPUNIT_TEST_SUITE_REGISTRATION(DFadDoubleTest2);
CPPUNIT_TEST_SUITE_REGISTRATION(SFadDoubleTest2);
CPPUNIT_TEST_SUITE_REGISTRATION(SLFadDoubleTest2);
CPPUNIT_TEST_SUITE_REGISTRATION(SimpleFadDoubleTest2);

#ifdef HAVE_SACADO_COMPLEX
typedef FadOpsUnitTest2<Sacado::Fad::DFad<std::complex<double> >,
			std::complex<double> > DFadComplexDoubleTest;
typedef FadOpsUnitTest2<Sacado::Fad::SFad<std::complex<double>,5>,
			std::complex<double> > SFadComplexDoubleTest;
typedef FadOpsUnitTest2<Sacado::Fad::SLFad<std::complex<double>,10>,
			std::complex<double> > SLFadComplexDoubleTest;
typedef FadOpsUnitTest2<Sacado::Fad::SimpleFad<std::complex<double> >,
			std::complex<double> > SimpleFadComplexDoubleTest;
CPPUNIT_TEST_SUITE_REGISTRATION(DFadComplexDoubleTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SFadComplexDoubleTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SLFadComplexDoubleTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SimpleFadComplexDoubleTest);
#endif
