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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "FadUnitTests2.hpp"
#include "Sacado_CacheFad_DFad.hpp"
#include "Sacado_CacheFad_SFad.hpp"
#include "Sacado_CacheFad_SLFad.hpp"

#ifdef HAVE_SACADO_COMPLEX
typedef FadOpsUnitTest2<Sacado::CacheFad::DFad<std::complex<double> >,
			std::complex<double> > DFadComplexDoubleTest;
typedef FadOpsUnitTest2<Sacado::CacheFad::SFad<std::complex<double>,5>,
			std::complex<double> > SFadComplexDoubleTest;
typedef FadOpsUnitTest2<Sacado::CacheFad::SLFad<std::complex<double>,10>,
			std::complex<double> > SLFadComplexDoubleTest;
CPPUNIT_TEST_SUITE_REGISTRATION(DFadComplexDoubleTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SFadComplexDoubleTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SLFadComplexDoubleTest);
#endif
