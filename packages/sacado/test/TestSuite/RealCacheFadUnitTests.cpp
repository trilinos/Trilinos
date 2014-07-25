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

typedef RealFadOpsUnitTest2<Sacado::CacheFad::DFad<double>,double> DFadDoubleTest2;
typedef RealFadOpsUnitTest2<Sacado::CacheFad::SFad<double,5>,double> SFadDoubleTest2;
typedef RealFadOpsUnitTest2<Sacado::CacheFad::SLFad<double,10>,double> SLFadDoubleTest2;

CPPUNIT_TEST_SUITE_REGISTRATION(DFadDoubleTest2);
CPPUNIT_TEST_SUITE_REGISTRATION(SFadDoubleTest2);
CPPUNIT_TEST_SUITE_REGISTRATION(SLFadDoubleTest2);

