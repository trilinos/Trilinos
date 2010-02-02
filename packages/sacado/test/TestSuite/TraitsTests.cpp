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

#include "TraitsTests.hpp"

#include "Sacado_Fad_SimpleFad.hpp"
#include "Sacado_CacheFad_DFad.hpp"
#include "Sacado_CacheFad_SFad.hpp"
#include "Sacado_CacheFad_SLFad.hpp"
#include "Sacado_Tay_CacheTaylor.hpp"
#ifdef HAVE_SACADO_STOKHOS
#include "Sacado_PCE_OrthogPoly.hpp"
#endif

template <>
Sacado::Fad::MemPool* Sacado::Fad::MemPoolStorage<double>::defaultPool_ = NULL;
template <>
Sacado::Fad::MemPool* Sacado::Fad::MemPoolStorage< Sacado::Fad::DMFad<double> >::defaultPool_ = NULL;

typedef TraitsTests< Sacado::Fad::DFad<double> > DFadTest;
typedef TraitsTests< Sacado::Fad::SFad<double,5> > SFadTest;
typedef TraitsTests< Sacado::Fad::SLFad<double,10> > SLFadTest;
typedef TraitsTests< Sacado::Fad::SimpleFad<double> > SimpleFadTest;
typedef TraitsTests< Sacado::Fad::DMFad<double> > DMFadTest;
typedef TraitsTests< Sacado::Fad::DVFad<double> > DVFadTest;

CPPUNIT_TEST_SUITE_REGISTRATION(DFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SLFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SimpleFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(DMFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(DVFadTest);

typedef TraitsTests< Sacado::ELRFad::DFad<double> > ELRDFadTest;
typedef TraitsTests< Sacado::ELRFad::SFad<double,5> > ELRSFadTest;
typedef TraitsTests< Sacado::ELRFad::SLFad<double,10> > ELRSLFadTest;

CPPUNIT_TEST_SUITE_REGISTRATION(ELRDFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(ELRSFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(ELRSLFadTest);

typedef TraitsTests< Sacado::CacheFad::DFad<double> > CacheDFadTest;
typedef TraitsTests< Sacado::CacheFad::SFad<double,5> > CacheSFadTest;
typedef TraitsTests< Sacado::CacheFad::SLFad<double,10> > CacheSLFadTest;

CPPUNIT_TEST_SUITE_REGISTRATION(CacheDFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(CacheSFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(CacheSLFadTest);

typedef TraitsTests< Sacado::ELRCacheFad::DFad<double> > ELRCacheDFadTest;
typedef TraitsTests< Sacado::ELRCacheFad::SFad<double,5> > ELRCacheSFadTest;
typedef TraitsTests< Sacado::ELRCacheFad::SLFad<double,10> > ELRCacheSLFadTest;

CPPUNIT_TEST_SUITE_REGISTRATION(ELRCacheDFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(ELRCacheSFadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(ELRCacheSLFadTest);

typedef TraitsTests< Sacado::LFad::LogicalSparse<double,bool> > LSFadTest;
CPPUNIT_TEST_SUITE_REGISTRATION(LSFadTest);

typedef TraitsTests< Sacado::FlopCounterPack::ScalarFlopCounter<double> > SFCTest;
CPPUNIT_TEST_SUITE_REGISTRATION(SFCTest);

typedef TraitsTests< Sacado::Tay::Taylor<double> > TaylorTest;
typedef TraitsTests< Sacado::Tay::CacheTaylor<double> > CacheTaylorTest;
CPPUNIT_TEST_SUITE_REGISTRATION(TaylorTest);
CPPUNIT_TEST_SUITE_REGISTRATION(CacheTaylorTest);

typedef TraitsTests< Sacado::Rad::ADvar<double> > RadTest;
typedef TraitsTests< Sacado::Rad2::ADvar<double> > Rad2Test;
typedef TraitsTests< Sacado::RadVec::ADvar<double> > RadVecTest;
CPPUNIT_TEST_SUITE_REGISTRATION(RadTest);
CPPUNIT_TEST_SUITE_REGISTRATION(Rad2Test);
CPPUNIT_TEST_SUITE_REGISTRATION(RadVecTest);

#ifdef HAVE_SACADO_STOKHOS
typedef TraitsTests< Sacado::PCE::OrthogPoly<double> > OrthogPolyTest;
CPPUNIT_TEST_SUITE_REGISTRATION(OrthogPolyTest);
#endif
