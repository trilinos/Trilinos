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

#include "FadUnitTests2.hpp"

#include "Sacado_Fad_SimpleFad.hpp"

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
