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

#include "Teuchos_SimpleObjectTable.hpp"
#include "Teuchos_RCP.hpp"

#include "TestClasses.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::RangeError;
using Teuchos::NullReferenceError;
using Teuchos::m_bad_cast;
using Teuchos::SimpleObjectTable;


/* SimpleObjectTable::storeNew() */

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleObjectTable, storeNew, T )
{
  ECHO(SimpleObjectTable<T> sot);
  ECHO(int id = sot.storeNew(new T));
  TEST_EQUALITY_CONST(id, 0);
  TEST_EQUALITY_CONST(nonnull(sot.getRCP(id)), true);
  TEST_EQUALITY_CONST(is_null(sot.getRCP(id)), false);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, storeNewNull )
{
  ECHO(SimpleObjectTable<A> sot);
  TEST_THROW(sot.storeNew(NULL), NullReferenceError); 
}


/* SimpleObjectTable::storeRCP() */

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleObjectTable, storeRCP, T )
{
  ECHO(SimpleObjectTable<T> sot);
  ECHO(RCP<T> rcpT = rcp(new T));
  TEST_EQUALITY_CONST(nonnull(rcpT), true);
  TEST_EQUALITY_CONST(is_null(rcpT), false);
  ECHO(int id = sot.storeRCP(rcpT));
  TEST_EQUALITY_CONST(id, 0);
  ECHO(RCP<T> rcpT2 = sot.getRCP(id));
  TEST_EQUALITY_CONST(nonnull(rcpT2), true);
  TEST_EQUALITY_CONST(is_null(rcpT2), false);
  TEST_EQUALITY(rcpT.get(), rcpT2.get());
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, storeRCPNull1 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(RCP<A> rcpA);
  TEST_THROW(sot.storeRCP(rcpA), NullReferenceError);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, storeRCPNull2 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(A *a=NULL);
  TEST_THROW(sot.storeRCP(rcp(a)), NullReferenceError);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, storeRCPNull3 )
{
  ECHO(SimpleObjectTable<A> sot);
  TEST_THROW(sot.storeRCP(Teuchos::null), NullReferenceError);
}


/* SimpleObjectTable::removeRCP() */

TEUCHOS_UNIT_TEST( SimpleObjectTable, removeRCPFromNew )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id, 0);
  ECHO(sot.removeRCP(id));
  TEST_EQUALITY_CONST(id, -1);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, removeRCPFromRCP )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(RCP<A> rcpA = rcp(new A));
  ECHO(int id = sot.storeRCP(rcpA));
  TEST_EQUALITY_CONST(id, 0);
  ECHO(sot.removeRCP(id));
  TEST_EQUALITY_CONST(id, -1);
}

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( SimpleObjectTable, removeRCPInvalid1 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = -1);
  TEST_THROW(sot.removeRCP(id), RangeError);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, removeRCPInvalid2 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = -2);
  TEST_THROW(sot.removeRCP(id), RangeError);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, removeRCPInvalid3 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = 0);
  TEST_THROW(sot.removeRCP(id), RangeError);
}

#endif /* TEUCHOS_DEBUG */


/* SimpleObjectTable::getRCP() */

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( SimpleObjectTable, getRCPInvalid1 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = -1);
  TEST_THROW(sot.getRCP(id), RangeError);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, getRCPInvalid2 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = -2);
  TEST_THROW(sot.getRCP(id), RangeError);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, getRCPInvalid3 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = 0);
  TEST_THROW(sot.getRCP(id), RangeError);
}

#endif /* TEUCHOS_DEBUG */

TEUCHOS_UNIT_TEST( SimpleObjectTable, getRCPInvalid4 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id, 0);
  ECHO(int id2 = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id2, 1);
  ECHO(int id3 = id);
  ECHO(sot.removeRCP(id));
  TEST_THROW(sot.getRCP(id3), RangeError);
}


/* SimpleObjectTable::storeCastedRCP() */

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( SimpleObjectTable, storeCastedRCP, T1, T2 )
{
  ECHO(SimpleObjectTable<T2> sot);
  ECHO(RCP<T1> rcpT1 = rcp(new T1));
  ECHO(T2 *pT2 = dynamic_cast<T2*>(rcpT1.get()));
  if (pT2 == NULL) {
    TEST_THROW(sot.storeCastedRCP(rcpT1), m_bad_cast);
  } else {
    ECHO(int id = sot.storeCastedRCP(rcpT1));
    TEST_EQUALITY_CONST(id, 0);
    TEST_EQUALITY_CONST(nonnull(sot.getRCP(id)), true);
    TEST_EQUALITY_CONST(rcpT1.shares_resource(sot.getRCP(id)), true);
  }
}


/* SimpleObjectTable::purge() */

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( SimpleObjectTable, purge )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = sot.storeNew(new A));
  TEST_EQUALITY_CONST(nonnull(sot.getRCP(id)), true);
  ECHO(sot.purge());
  TEST_THROW(sot.getRCP(id), RangeError);
}

#endif /* TEUCHOS_DEBUG */


/* SimpleObjectTable's freedIndices table */

TEUCHOS_UNIT_TEST( SimpleObjectTable, recycleIndex1 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id, 0);
  ECHO(sot.removeRCP(id));
  ECHO(int id2 = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id2, 0);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, recycleIndex2 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id, 0);
  ECHO(int id2 = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id2, 1);
  ECHO(sot.removeRCP(id));
  ECHO(int id3 = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id3, 0);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, recycleIndex3 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id, 0);
  ECHO(int id2 = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id2, 1);
  ECHO(sot.removeRCP(id2));
  ECHO(int id3 = sot.storeNew(new A));
  TEST_EQUALITY_CONST(id3, 1);
}


/* SimpleObjectTable's RCP counts */

TEUCHOS_UNIT_TEST( SimpleObjectTable, rcpNewShared1 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(A a);
  ECHO(int id = sot.storeNew(&a, false));

  ECHO(RCP<A> rcpA = sot.getRCP(id));
  TEST_EQUALITY(rcpA.get(), &a);
  TEST_EQUALITY_CONST(rcpA.has_ownership(), false);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 2);

  ECHO(rcpA = null);
  ECHO(int cnt = sot.removeRCP(id));
  TEST_EQUALITY_CONST(cnt, 0);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, rcpNewShared2 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(A a);
  ECHO(int id = sot.storeNew(&a, false));

  ECHO(RCP<A> rcpA = sot.getRCP(id));
  TEST_EQUALITY(rcpA.get(), &a);
  TEST_EQUALITY_CONST(rcpA.has_ownership(), false);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 2);

  ECHO(int cnt = sot.removeRCP(id));
  TEST_EQUALITY_CONST(cnt, 1);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 1);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, rcpNewOwned1 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = sot.storeNew(new A));

  ECHO(RCP<A> rcpA = sot.getRCP(id));
  TEST_EQUALITY_CONST(rcpA.has_ownership(), true);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 2);

  ECHO(rcpA = null);
  ECHO(int cnt = sot.removeRCP(id));
  TEST_EQUALITY_CONST(cnt, 0);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, rcpNewOwned2 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(int id = sot.storeNew(new A));

  ECHO(RCP<A> rcpA = sot.getRCP(id));
  TEST_EQUALITY_CONST(rcpA.has_ownership(), true);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 2);

  ECHO(int cnt = sot.removeRCP(id));
  TEST_EQUALITY_CONST(cnt, 1);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 1);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, rcpRCPOwned1 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(RCP<A> rcpA = rcp(new A));
  TEST_EQUALITY_CONST(rcpA.has_ownership(), true);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 1);

  ECHO(int id = sot.storeRCP(rcpA));
  TEST_EQUALITY_CONST(rcpA.has_ownership(), true);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 2);

  ECHO(RCP<A> rcpA2 = sot.getRCP(id));
  TEST_EQUALITY(rcpA2.get(), rcpA.get());
  TEST_EQUALITY_CONST(rcpA2.has_ownership(), true);
  TEST_EQUALITY_CONST(rcpA2.strong_count(), 3);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, rcpRCPOwned2 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(RCP<A> rcpA = rcp(new A));
  TEST_EQUALITY_CONST(rcpA.has_ownership(), true);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 1);

  ECHO(int id = sot.storeRCP(rcpA));
  TEST_EQUALITY_CONST(rcpA.strong_count(), 2);

  ECHO(rcpA = null);
  ECHO(int cnt = sot.removeRCP(id));
  TEST_EQUALITY_CONST(cnt, 0);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, rcpRCPOwned3 )
{
  ECHO(SimpleObjectTable<A> sot);
  ECHO(RCP<A> rcpA = rcp(new A));
  TEST_EQUALITY_CONST(rcpA.has_ownership(), true);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 1);

  ECHO(int id = sot.storeRCP(rcpA));
  TEST_EQUALITY_CONST(rcpA.strong_count(), 2);

  ECHO(int cnt = sot.removeRCP(id));
  TEST_EQUALITY_CONST(cnt, 1);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 1);
}

TEUCHOS_UNIT_TEST( SimpleObjectTable, rcpDestructTable )
{
  ECHO(SimpleObjectTable<A> *psot = new SimpleObjectTable<A>);
  ECHO(A *pA = new A);
  ECHO(int id = psot->storeNew(pA));

  ECHO(RCP<A> rcpA = psot->getRCP(id));
  TEST_EQUALITY_CONST(rcpA.strong_count(), 2);
  TEST_EQUALITY(rcpA.get(), pA);

  ECHO(delete psot);
  TEST_EQUALITY_CONST(rcpA.strong_count(), 1);
  TEST_EQUALITY(rcpA.get(), pA);
}


//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SimpleObjectTable, storeNew, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SimpleObjectTable, storeRCP, T ) \
  DEBUG_UNIT_TEST_GROUP( T )

#define UNIT_TEST_GROUP_PAIR( T1, T2 ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( SimpleObjectTable, storeCastedRCP, T1, T2 )

#define UNIT_TEST_GROUP_PAIR_SYM( T1, T2 ) \
  UNIT_TEST_GROUP_PAIR( T1, T2 ) \
  UNIT_TEST_GROUP_PAIR( T2, T1 )


UNIT_TEST_GROUP(A)
UNIT_TEST_GROUP(B1)
UNIT_TEST_GROUP(B2)
UNIT_TEST_GROUP(C)

UNIT_TEST_GROUP_PAIR(A, A)
UNIT_TEST_GROUP_PAIR(B1, B1)
UNIT_TEST_GROUP_PAIR(B2, B2)
UNIT_TEST_GROUP_PAIR(C, C)

UNIT_TEST_GROUP_PAIR_SYM(A, B1)
UNIT_TEST_GROUP_PAIR_SYM(A, B2)
UNIT_TEST_GROUP_PAIR_SYM(A, C)
UNIT_TEST_GROUP_PAIR_SYM(B1, B2)
UNIT_TEST_GROUP_PAIR_SYM(B1, C)
UNIT_TEST_GROUP_PAIR_SYM(B2, C)


} // namespace
