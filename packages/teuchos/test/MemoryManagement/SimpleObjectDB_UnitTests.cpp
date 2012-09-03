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

#include "Teuchos_SimpleObjectDB.hpp"
#include "Teuchos_getConst.hpp"

#include "TestClasses.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::RangeError;
using Teuchos::NullReferenceError;
using Teuchos::m_bad_cast;
using Teuchos::SimpleObjectDB;
using Teuchos::getConst;


//
// SimpleObjectDB::SimpleObjectDB()
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleObjectDB, defaultConstruct, T )
{
  ECHO(SimpleObjectDB<T> sot);
  TEST_EQUALITY_CONST(sot.tableSize(), 0);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(sot.numObjects(), 0);
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW(sot.getNonconstObjRCP(0), RangeError);
  TEST_THROW(sot.getConstObjRCP(0), RangeError);
  TEST_THROW(sot.getNonconstObjPtr(0), RangeError);
  TEST_THROW(sot.getConstObjPtr(0), RangeError);
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
}


//
// createSimpleObjectDB()
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleObjectDB, createSimpleObjectDB, T )
{
  ECHO(RCP<SimpleObjectDB<T> > sot = Teuchos::createSimpleObjectDB<T>());
  TEST_EQUALITY_CONST(sot->numObjects(), 0);
}


//
// SimpleObjectDB::storeNonconstObj()
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleObjectDB, storeNonconstObj, T )
{
  ECHO(SimpleObjectDB<T> sot);
  ECHO(const int id = sot.storeNonconstObj(T::create()));
  TEST_EQUALITY_CONST(sot.tableSize(), 1);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(id, 0);
  TEST_EQUALITY_CONST(nonnull(sot.getNonconstObjRCP(id)), true);
}


//
// SimpleObjectDB::get[Nonconst,Const]Obj[RCP,Ptr]()
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleObjectDB, getNonconstObjRCP, T )
{
  ECHO(SimpleObjectDB<T> sot);
  ECHO(const RCP<T> obj = T::create());
  ECHO(const int id = sot.storeNonconstObj(obj));
  TEST_EQUALITY(obj.get(), sot.getNonconstObjRCP(id).get());
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleObjectDB, getConstObjRCP, T )
{
  ECHO(SimpleObjectDB<T> sot);
  ECHO(const RCP<T> obj = T::create());
  ECHO(const int id = sot.storeNonconstObj(obj));
  TEST_EQUALITY(obj.get(), getConst(sot).getConstObjRCP(id).get());
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleObjectDB, getNonconstObjPtr, T )
{
  ECHO(SimpleObjectDB<T> sot);
  ECHO(const RCP<T> obj = T::create());
  ECHO(const int id = sot.storeNonconstObj(obj));
  TEST_EQUALITY(obj.get(), sot.getNonconstObjPtr(id).get());
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleObjectDB, getConstObjPtr, T )
{
  ECHO(SimpleObjectDB<T> sot);
  ECHO(const RCP<T> obj = T::create());
  ECHO(const int id = sot.storeNonconstObj(obj));
  TEST_EQUALITY(obj.get(), getConst(sot).getConstObjPtr(id).get());
}


//
// SimpleObjectDB::storeConstObj(), getNonconstObjRCP()
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleObjectDB, storeConstObj, T )
{
  ECHO(SimpleObjectDB<T> sot);
  ECHO(RCP<const T> obj = T::create());
  ECHO(const int id = sot.storeConstObj(obj));
  TEST_EQUALITY_CONST(id, 0);
  TEST_EQUALITY_CONST(sot.tableSize(), 1);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY(sot.getConstObjRCP(id).get(), obj.get());
  TEST_THROW(sot.getNonconstObjRCP(id), Teuchos::NonconstAccessError);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, storeNonconstObjNull1 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(RCP<A> rcpA);
  TEST_THROW(sot.storeNonconstObj(rcpA), NullReferenceError);
  TEST_EQUALITY_CONST(sot.tableSize(), 0);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, storeNonconstObjNull2 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(A *a=NULL);
  TEST_THROW(sot.storeNonconstObj(rcp(a)), NullReferenceError);
  TEST_EQUALITY_CONST(sot.tableSize(), 0);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, storeNonconstObjNull3 )
{
  ECHO(SimpleObjectDB<A> sot);
  TEST_THROW(sot.storeNonconstObj(Teuchos::null), NullReferenceError);
  TEST_EQUALITY_CONST(sot.tableSize(), 0);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
}


//
// SimpleObjectDB::remove[Nonconst,Const]Obj()
//


TEUCHOS_UNIT_TEST( SimpleObjectDB, removeObj_storeNonconstObj_1_0 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(const RCP<A> obj = A::create());
  ECHO(const int id1 = sot.storeNonconstObj(obj));
  TEST_EQUALITY(sot.getNonconstObjRCP(id1).ptr(), obj.ptr());
  ECHO(sot.removeObj(id1));
  TEST_EQUALITY_CONST(sot.numObjects(), 0);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, removeObj_storeConstObj_1_0 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(const RCP<A> obj = A::create());
  ECHO(const int id1 = sot.storeConstObj(obj));
  TEST_EQUALITY(sot.getConstObjRCP(id1).ptr(), obj.ptr());
  ECHO(sot.removeObj(id1));
  TEST_EQUALITY_CONST(sot.numObjects(), 0);
  ECHO(const int id2 = sot.storeConstObj(obj));
  TEST_EQUALITY_CONST(id2, 0);
  TEST_EQUALITY(sot.getConstObjRCP(id2).ptr(), obj.ptr());
  TEST_EQUALITY_CONST(sot.numObjects(), 1);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, removeNonconstObj_storeNonconstObj_1_0 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(const RCP<A> obj = A::create());
  ECHO(const int id1 = sot.storeNonconstObj(obj));
  TEST_EQUALITY(sot.getNonconstObjRCP(id1).ptr(), obj.ptr());
  ECHO(const RCP<A> obj2 = sot.removeNonconstObj(id1));
  TEST_EQUALITY(obj2.ptr(), obj.ptr());
  TEST_EQUALITY_CONST(sot.numObjects(), 0);
  ECHO(const int id2 = sot.storeNonconstObj(obj));
  TEST_EQUALITY_CONST(id2, 0);
  TEST_EQUALITY(sot.getNonconstObjRCP(id2).ptr(), obj.ptr());
  TEST_EQUALITY_CONST(sot.numObjects(), 1);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, removeConstObj_storeConstObj_1_0 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(const RCP<const A> obj = A::create());
  ECHO(const int id1 = sot.storeConstObj(obj));
  TEST_EQUALITY(sot.getConstObjRCP(id1).ptr(), obj.ptr());
  ECHO(const RCP<const A> obj2 = sot.removeConstObj(id1));
  TEST_EQUALITY(obj2.ptr(), obj.ptr());
  TEST_EQUALITY_CONST(sot.numObjects(), 0);
  ECHO(const int id2 = sot.storeConstObj(obj));
  TEST_EQUALITY_CONST(id2, 0);
  TEST_EQUALITY(sot.getConstObjRCP(id2).ptr(), obj.ptr());
  TEST_EQUALITY_CONST(sot.numObjects(), 1);
}


#ifdef TEUCHOS_DEBUG


TEUCHOS_UNIT_TEST( SimpleObjectDB, removeNonconstObjInvalid1 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(int id = -1);
  TEST_THROW(sot.removeNonconstObj(id), RangeError);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, removeNonconstObjInvalid2 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(int id = -2);
  TEST_THROW(sot.removeNonconstObj(id), RangeError);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, removeNonconstObjInvalid3 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(int id = 0);
  TEST_THROW(sot.removeNonconstObj(id), RangeError);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, removeObjTwice )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(const int id = sot.storeNonconstObj(A::create()));
  ECHO(sot.removeObj(id));
  TEST_THROW(sot.removeObj(id), NullReferenceError);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, removeNonconstObjTwice )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(const int id = sot.storeNonconstObj(A::create()));
  ECHO(sot.removeNonconstObj(id));
  TEST_THROW(sot.removeNonconstObj(id), NullReferenceError);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, removeConstObjTwice )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(const int id = sot.storeNonconstObj(A::create()));
  ECHO(sot.removeConstObj(id));
  TEST_THROW(sot.removeConstObj(id), NullReferenceError);
}


#endif // TEUCHOS_DEBUG


//
// SimpleObjectDB::getNonconstObjRCP()
//


#ifdef TEUCHOS_DEBUG


TEUCHOS_UNIT_TEST( SimpleObjectDB, getNonconstObjRCPInvalid1 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(int id = -1);
  TEST_THROW(sot.getNonconstObjRCP(id), RangeError);
}

TEUCHOS_UNIT_TEST( SimpleObjectDB, getNonconstObjRCPInvalid2 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(int id = -2);
  TEST_THROW(sot.getNonconstObjRCP(id), RangeError);
}

TEUCHOS_UNIT_TEST( SimpleObjectDB, getNonconstObjRCPInvalid3 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(int id = 0);
  TEST_THROW(sot.getNonconstObjRCP(id), RangeError);
}


#endif // TEUCHOS_DEBUG


TEUCHOS_UNIT_TEST( SimpleObjectDB, getNonconstObjRCPInvalid4 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(int id = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id, 0);
  ECHO(int id2 = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id2, 1);
  ECHO(sot.removeNonconstObj(id));
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW(sot.getNonconstObjRCP(id), NullReferenceError);
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
}


//
// SimpleObjectDB::storeCastedNonconstObj()
//


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( SimpleObjectDB, storeCastedNonconstObj, T1, T2 )
{
  ECHO(SimpleObjectDB<T2> sot);
  ECHO(RCP<T1> rcpT1 = rcp(new T1));
  ECHO(T2 *pT2 = dynamic_cast<T2*>(rcpT1.get()));
  if (pT2 == NULL) {
    TEST_THROW(sot.storeCastedNonconstObj(rcpT1), m_bad_cast);
  } else {
    ECHO(int id = sot.storeCastedNonconstObj(rcpT1));
    TEST_EQUALITY_CONST(id, 0);
    TEST_EQUALITY_CONST(nonnull(sot.getNonconstObjRCP(id)), true);
    TEST_EQUALITY_CONST(rcpT1.shares_resource(sot.getNonconstObjRCP(id)), true);
  }
}


//
// SimpleObjectDB::purge()
//


#ifdef TEUCHOS_DEBUG


TEUCHOS_UNIT_TEST( SimpleObjectDB, purge )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(const RCP<A> a(new A));
  ECHO(int id = sot.storeNonconstObj(a));
  ECHO(int id2 = sot.storeNonconstObj(a));
  TEST_EQUALITY_CONST(nonnull(sot.getNonconstObjRCP(id)), true);
  TEST_EQUALITY_CONST(nonnull(sot.getNonconstObjRCP(id2)), true);
  TEST_EQUALITY_CONST(sot.tableSize(), 2);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(sot.numObjects(), 2);
  ECHO(sot.removeNonconstObj(id));
  TEST_EQUALITY_CONST(sot.tableSize(), 2);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 1);
  TEST_EQUALITY_CONST(sot.numObjects(), 1);
  ECHO(sot.purge());
  TEST_EQUALITY_CONST(sot.tableSize(), 0);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(sot.numObjects(), 0);
  TEST_THROW(sot.getNonconstObjRCP(id), RangeError);
  TEST_EQUALITY_CONST(a.strong_count(), 1); // sot gave up its RCP?
}


#endif // TEUCHOS_DEBUG


//
// SimpleObjectDB's freedIndices table
//


TEUCHOS_UNIT_TEST( SimpleObjectDB, recycleIndex1 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(int id = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id, 0);
  TEST_EQUALITY_CONST(sot.tableSize(), 1);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(sot.numObjects(), 1);
  ECHO(sot.removeNonconstObj(id));
  TEST_EQUALITY_CONST(sot.tableSize(), 1);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 1);
  TEST_EQUALITY_CONST(sot.numObjects(), 0);
  ECHO(int id2 = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id2, 0);
  TEST_EQUALITY_CONST(sot.tableSize(), 1);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(sot.numObjects(), 1);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, recycleIndex2 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(int id = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id, 0);
  ECHO(int id2 = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id2, 1);
  TEST_EQUALITY_CONST(sot.tableSize(), 2);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(sot.numObjects(), 2);
  ECHO(sot.removeNonconstObj(id));
  TEST_EQUALITY_CONST(sot.tableSize(), 2);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 1);
  TEST_EQUALITY_CONST(sot.numObjects(), 1);
  ECHO(int id3 = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id3, 0);
  TEST_EQUALITY_CONST(sot.tableSize(), 2);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(sot.numObjects(), 2);
}


TEUCHOS_UNIT_TEST( SimpleObjectDB, recycleIndex3 )
{
  ECHO(SimpleObjectDB<A> sot);
  ECHO(int id = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id, 0);
  ECHO(int id2 = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id2, 1);
  ECHO(int id3 = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id3, 2);
  TEST_EQUALITY_CONST(sot.tableSize(), 3);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(sot.numObjects(), 3);
  ECHO(sot.removeNonconstObj(id2));
  TEST_EQUALITY_CONST(sot.tableSize(), 3);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 1);
  TEST_EQUALITY_CONST(sot.numObjects(), 2);
  ECHO(int id4 = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id4, 1);
  TEST_EQUALITY_CONST(sot.tableSize(), 3);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(sot.numObjects(), 3);
  ECHO(int id5 = sot.storeNonconstObj(A::create()));
  TEST_EQUALITY_CONST(id5, 3);
  TEST_EQUALITY_CONST(sot.tableSize(), 4);
  TEST_EQUALITY_CONST(sot.numFreeIndexes(), 0);
  TEST_EQUALITY_CONST(sot.numObjects(), 4);
}


//
// Template Instantiations
//


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SimpleObjectDB, defaultConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SimpleObjectDB, createSimpleObjectDB, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SimpleObjectDB, storeNonconstObj, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SimpleObjectDB, getNonconstObjRCP, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SimpleObjectDB, getConstObjRCP, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SimpleObjectDB, getNonconstObjPtr, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SimpleObjectDB, getConstObjPtr, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SimpleObjectDB, storeConstObj, T )

#define UNIT_TEST_GROUP_PAIR( T1, T2 ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( SimpleObjectDB, storeCastedNonconstObj, T1, T2 )

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
