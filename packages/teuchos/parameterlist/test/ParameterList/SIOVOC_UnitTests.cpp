// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_StringIndexedOrderedValueObjectContainer.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace Teuchos {


TEUCHOS_UNIT_TEST( OrdinalIndex, defaultConstruct )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  ECHO(SIOVOCB::OrdinalIndex idx);
  TEST_EQUALITY_CONST(idx.idx, -1); // Depends on implementation choice!
}


TEUCHOS_UNIT_TEST( OrdinalIndex, construct )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  ECHO(SIOVOCB::OrdinalIndex idx(5));
  TEST_EQUALITY_CONST(idx.idx, 5);
}


TEUCHOS_UNIT_TEST( KeyObjectPair, defaultConstruct )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  ECHO(SIOVOCB::KeyObjectPair<SIOVOCB::OrdinalIndex> kop);
  TEST_EQUALITY_CONST(kop.first, "");
  TEST_EQUALITY_CONST(kop.second.idx, -1);
  TEST_EQUALITY_CONST(kop.isActive(), true);
}


TEUCHOS_UNIT_TEST( KeyObjectPair, construct_default_isActive )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  ECHO(SIOVOCB::KeyObjectPair<SIOVOCB::OrdinalIndex> kop("key_name1", 7));
  TEST_EQUALITY_CONST(kop.first, "key_name1");
  TEST_EQUALITY_CONST(kop.second.idx, 7);
  TEST_EQUALITY_CONST(kop.isActive(), true);
}


TEUCHOS_UNIT_TEST( KeyObjectPair, makeInvalid )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  ECHO(SIOVOCB::KeyObjectPair<SIOVOCB::OrdinalIndex> kop =
    SIOVOCB::KeyObjectPair<SIOVOCB::OrdinalIndex>::makeInvalid());
  TEST_EQUALITY_CONST(kop.first, "");
  TEST_EQUALITY_CONST(kop.second.idx, -1);
  TEST_EQUALITY_CONST(kop.isActive(), false);
}


TEUCHOS_UNIT_TEST( KeyObjectPair, construct_set_isActive )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  ECHO(SIOVOCB::KeyObjectPair<SIOVOCB::OrdinalIndex> kop("key_name2", 8, false));
  TEST_EQUALITY_CONST(kop.first, "key_name2");
  TEST_EQUALITY_CONST(kop.second.idx, 8);
  TEST_EQUALITY_CONST(kop.isActive(), false);
}


TEUCHOS_UNIT_TEST( KeyObjectPair, copyConstruct )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  ECHO(SIOVOCB::KeyObjectPair<SIOVOCB::OrdinalIndex> kop1("key_name", 3, false));
  ECHO(SIOVOCB::KeyObjectPair<SIOVOCB::OrdinalIndex> kop2(kop1));
  TEST_EQUALITY_CONST(kop2.first, "key_name");
  TEST_EQUALITY_CONST(kop2.second.idx, 3);
  TEST_EQUALITY_CONST(kop2.isActive(), false);
}


TEUCHOS_UNIT_TEST( KeyObjectPair, assign )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  ECHO(SIOVOCB::KeyObjectPair<SIOVOCB::OrdinalIndex> kop1("key_name", 3, false));
  ECHO(SIOVOCB::KeyObjectPair<SIOVOCB::OrdinalIndex> kop2);
  TEST_EQUALITY_CONST(kop2.isActive(), true);
  ECHO(kop2 = kop1);
  TEST_EQUALITY_CONST(kop2.first, "key_name");
  TEST_EQUALITY_CONST(kop2.second.idx, 3);
  TEST_EQUALITY_CONST(kop2.isActive(), false);
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, defaultConstruct )
{
  StringIndexedOrderedValueObjectContainer<int> oc;
  TEST_EQUALITY_CONST(oc.numObjects(), 0);
  TEST_EQUALITY_CONST(oc.numStorage(), 0);
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, basic_set_get )
{

  //typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB; // unused
  //typedef SIOVOCB::Ordinal Ordinal;
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);

  ECHO(const Ordinal my_int_1_idx1 = oc.setObj("my_int_1", 3));
  TEST_EQUALITY_CONST(my_int_1_idx1, 0);
  TEST_EQUALITY_CONST(oc.numObjects(), 1);
  TEST_EQUALITY_CONST(oc.numStorage(), 1);
  ECHO(const Ordinal my_int_1_idx2 = oc.getObjOrdinalIndex(("my_int_1")));
  TEST_EQUALITY(my_int_1_idx2, my_int_1_idx2);
  TEST_EQUALITY_CONST(*oc.getNonconstObjPtr(my_int_1_idx1), 3);
  TEST_EQUALITY_CONST(*oc.getObjPtr(my_int_1_idx1), 3);
  TEST_EQUALITY_CONST(*oc.getNonconstObjPtr("my_int_1"), 3);
  TEST_EQUALITY_CONST(*oc.getObjPtr("my_int_1"), 3);

  ECHO(const Ordinal my_int_2_idx1 = oc.setObj("my_int_2", 4));
  TEST_EQUALITY_CONST(my_int_2_idx1, 1);
  TEST_EQUALITY_CONST(oc.numObjects(), 2);
  TEST_EQUALITY_CONST(oc.numStorage(), 2);
  ECHO(const Ordinal my_int_2_idx2 = oc.getObjOrdinalIndex(("my_int_2")));
  TEST_EQUALITY(my_int_2_idx2, my_int_2_idx2);
  TEST_EQUALITY_CONST(*oc.getNonconstObjPtr(my_int_2_idx1), 4);
  TEST_EQUALITY_CONST(*oc.getObjPtr(my_int_2_idx1), 4);
  TEST_EQUALITY_CONST(*oc.getNonconstObjPtr("my_int_2"), 4);
  TEST_EQUALITY_CONST(*oc.getObjPtr("my_int_2"), 4);

}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, set_two_keep_ref )
{
  // Test test makes sure that objects keep the same address when adding new
  // objects.
  //typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB; // unused
  //typedef SIOVOCB::Ordinal Ordinal;
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);
  ECHO(const Ordinal my_int_1_idx = oc.setObj("my_int_1", 3));
  ECHO(const int &my_int_1 = *oc.getObjPtr(my_int_1_idx));
  ECHO(oc.setObj("my_int_2", 4));
  TEST_EQUALITY_CONST(*oc.getObjPtr(my_int_1_idx), 3);
  TEST_EQUALITY(&my_int_1, oc.getObjPtr(my_int_1_idx).get());
  TEST_EQUALITY_CONST(my_int_1, 3);
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, getObjOrdinalIndex )
{

  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  //typedef SIOVOCB::Ordinal Ordinal;
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);
  ECHO(oc.setObj("my_int_1", 3));
  TEST_EQUALITY_CONST(oc.getObjOrdinalIndex("my_int_1"), 0);
  ECHO(oc.setObj("my_int_2", 3));
  TEST_EQUALITY_CONST(oc.getObjOrdinalIndex("my_int_2"), 1);
  TEST_EQUALITY_CONST(oc.getObjOrdinalIndex("does_not_exist"), SIOVOCB::getInvalidOrdinal());
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, set_and_set_again )
{

  //typedef StringIndexedOrderedValueObjectContainerBase::Ordinal Ordinal;
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);

  ECHO(const Ordinal my_int_1_idx1 = oc.setObj("my_int_1", 3));
  TEST_EQUALITY_CONST(my_int_1_idx1, 0);
  ECHO(const Ordinal my_int_1_idx2 = oc.getObjOrdinalIndex(("my_int_1")));
  TEST_EQUALITY_CONST(my_int_1_idx2, my_int_1_idx1);
  TEST_EQUALITY_CONST(*oc.getNonconstObjPtr("my_int_1"), 3);

  ECHO(const Ordinal my_int_1_idx3 = oc.setObj("my_int_1", 4));
  TEST_EQUALITY_CONST(my_int_1_idx3, 0);
  ECHO(const Ordinal my_int_1_idx4 = oc.getObjOrdinalIndex(("my_int_1")));
  TEST_EQUALITY_CONST(my_int_1_idx3, my_int_1_idx4);
  TEST_EQUALITY_CONST(*oc.getNonconstObjPtr("my_int_1"), 4);

}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, basicNonconstIterators )
{
  typedef StringIndexedOrderedValueObjectContainer<int>::Iterator Iterator;
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);
  ECHO(oc.setObj("c", 1));
  ECHO(oc.setObj("a", 2));
  ECHO(oc.setObj("b", 3));
  ECHO(Iterator itr = oc.nonconstBegin());
  TEST_EQUALITY_CONST(itr->first, "c");
  TEST_EQUALITY_CONST(itr->second, 1);
  ECHO(itr->second = 5);
  TEST_EQUALITY_CONST(itr->second, 5);
  ECHO(++itr);
  TEST_EQUALITY_CONST(itr->first, "a");
  TEST_EQUALITY_CONST(itr->second, 2);
  ECHO(itr->second = 6);
  TEST_EQUALITY_CONST(itr->second, 6);
  ECHO(++itr);
  TEST_EQUALITY_CONST(itr->first, "b");
  TEST_EQUALITY_CONST(itr->second, 3);
  ECHO(itr->second = 7);
  TEST_EQUALITY_CONST(itr->second, 7);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, oc.nonconstEnd());
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, basicConstIterators )
{
 typedef StringIndexedOrderedValueObjectContainer<int>::ConstIterator ConstIterator;
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);
  ECHO(oc.setObj("c", 1));
  ECHO(oc.setObj("a", 2));
  ECHO(oc.setObj("b", 3));
  ECHO(ConstIterator itr = oc.begin());
  TEST_EQUALITY_CONST(itr->first, "c");
  TEST_EQUALITY_CONST(itr->second, 1);
  ECHO(++itr);
  TEST_EQUALITY_CONST(itr->first, "a");
  TEST_EQUALITY_CONST(itr->second, 2);
  ECHO(++itr);
  TEST_EQUALITY_CONST(itr->first, "b");
  TEST_EQUALITY_CONST(itr->second, 3);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, oc.end());
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, removeObj_idx_first )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  typedef SIOVOCB::OrdinalIndex OI;
  typedef StringIndexedOrderedValueObjectContainer<OI>::ConstIterator ConstIterator;
  ECHO(StringIndexedOrderedValueObjectContainer<OI> oc);
  ECHO(oc.setObj("c", 1));
  ECHO(oc.setObj("a", 2));
  ECHO(oc.setObj("b", 3));
  TEST_EQUALITY_CONST(oc.numObjects(), 3);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  ECHO(const OI &b = *oc.getObjPtr("b"));
  ECHO(oc.removeObj(0));
  TEST_EQUALITY_CONST(oc.numObjects(), 2);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  TEST_EQUALITY(&b, oc.getObjPtr("b").get());
  ECHO(ConstIterator itr = oc.begin());
  TEST_EQUALITY_CONST(itr->first, "a");
  TEST_EQUALITY_CONST(itr->second.idx, 2);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_EQUALITY_CONST(itr->first, "b");
  TEST_EQUALITY_CONST(itr->second.idx, 3);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, oc.end());
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, removeObj_key_first )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  typedef SIOVOCB::OrdinalIndex OI;
  typedef StringIndexedOrderedValueObjectContainer<OI>::ConstIterator ConstIterator;
  ECHO(StringIndexedOrderedValueObjectContainer<OI> oc);
  ECHO(oc.setObj("c", 1));
  ECHO(oc.setObj("a", 2));
  ECHO(oc.setObj("b", 3));
  TEST_EQUALITY_CONST(oc.numObjects(), 3);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  ECHO(const OI &b = *oc.getObjPtr("b"));
  ECHO(oc.removeObj("c"));
  TEST_EQUALITY_CONST(oc.numObjects(), 2);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  TEST_EQUALITY(&b, oc.getObjPtr("b").get());
  ECHO(ConstIterator itr = oc.begin());
  TEST_EQUALITY_CONST(itr->first, "a");
  TEST_EQUALITY_CONST(itr->second.idx, 2);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_EQUALITY_CONST(itr->first, "b");
  TEST_EQUALITY_CONST(itr->second.idx, 3);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, oc.end());
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, removeObj_idx_middle )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  typedef SIOVOCB::OrdinalIndex OI;
  typedef StringIndexedOrderedValueObjectContainer<OI>::ConstIterator ConstIterator;
  ECHO(StringIndexedOrderedValueObjectContainer<OI> oc);
  ECHO(oc.setObj("c", 1));
  ECHO(oc.setObj("a", 2));
  ECHO(oc.setObj("b", 3));
  TEST_EQUALITY_CONST(oc.numObjects(), 3);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  ECHO(const OI &c = *oc.getObjPtr("c"));
  ECHO(oc.removeObj(1));
  TEST_EQUALITY(&c, oc.getObjPtr("c").get());
  ECHO(ConstIterator itr = oc.begin());
  TEST_EQUALITY_CONST(oc.numObjects(), 2);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  TEST_EQUALITY_CONST(itr->first, "c");
  TEST_EQUALITY_CONST(itr->second.idx, 1);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_EQUALITY_CONST(itr->first, "b");
  TEST_EQUALITY_CONST(itr->second.idx, 3);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, oc.end());
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, removeObj_key_middle )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  typedef SIOVOCB::OrdinalIndex OI;
  typedef StringIndexedOrderedValueObjectContainer<OI>::ConstIterator ConstIterator;
  ECHO(StringIndexedOrderedValueObjectContainer<OI> oc);
  ECHO(oc.setObj("c", 1));
  ECHO(oc.setObj("a", 2));
  ECHO(oc.setObj("b", 3));
  TEST_EQUALITY_CONST(oc.numObjects(), 3);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  ECHO(const OI &c = *oc.getObjPtr("c"));
  ECHO(oc.removeObj("a"));
  TEST_EQUALITY(&c, oc.getObjPtr("c").get());
  TEST_EQUALITY_CONST(oc.numObjects(), 2);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  ECHO(ConstIterator itr = oc.begin());
  TEST_EQUALITY_CONST(itr->first, "c");
  TEST_EQUALITY_CONST(itr->second.idx, 1);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_EQUALITY_CONST(itr->first, "b");
  TEST_EQUALITY_CONST(itr->second.idx, 3);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, oc.end());
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, removeObj_idx_last )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  typedef SIOVOCB::OrdinalIndex OI;
  typedef StringIndexedOrderedValueObjectContainer<OI>::ConstIterator ConstIterator;
  ECHO(StringIndexedOrderedValueObjectContainer<OI> oc);
  ECHO(oc.setObj("c", 1));
  ECHO(oc.setObj("a", 2));
  ECHO(oc.setObj("b", 3));
  TEST_EQUALITY_CONST(oc.numObjects(), 3);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  ECHO(const OI &a = *oc.getObjPtr("a"));
  ECHO(oc.removeObj(2));
  TEST_EQUALITY_CONST(oc.numObjects(), 2);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  TEST_EQUALITY(&a, oc.getObjPtr("a").get());
  ECHO(ConstIterator itr = oc.begin());
  TEST_EQUALITY_CONST(itr->first, "c");
  TEST_EQUALITY_CONST(itr->second.idx, 1);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_EQUALITY_CONST(itr->first, "a");
  TEST_EQUALITY_CONST(itr->second.idx, 2);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, oc.end());
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, removeObj_key_last )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  typedef SIOVOCB::OrdinalIndex OI;
  typedef StringIndexedOrderedValueObjectContainer<OI>::ConstIterator ConstIterator;
  ECHO(StringIndexedOrderedValueObjectContainer<OI> oc);
  ECHO(oc.setObj("c", 1));
  ECHO(oc.setObj("a", 2));
  ECHO(oc.setObj("b", 3));
  TEST_EQUALITY_CONST(oc.numObjects(), 3);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  ECHO(const OI &a = *oc.getObjPtr("a"));
  ECHO(oc.removeObj("b"));
  TEST_EQUALITY_CONST(oc.numObjects(), 2);
  TEST_EQUALITY_CONST(oc.numStorage(), 3);
  TEST_EQUALITY(&a, oc.getObjPtr("a").get());
  ECHO(ConstIterator itr = oc.begin());
  TEST_EQUALITY_CONST(itr->first, "c");
  TEST_EQUALITY_CONST(itr->second.idx, 1);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_EQUALITY_CONST(itr->first, "a");
  TEST_EQUALITY_CONST(itr->second.idx, 2);
  TEST_EQUALITY_CONST(itr->isActive(), true);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, oc.end());
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, getNonconstObjPtr_idx_invalid )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  //typedef StringIndexedOrderedValueObjectContainer<int>::Iterator Iterator; // unused
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);
  ECHO(oc.setObj("a", 4));
  TEST_EQUALITY_CONST(*oc.getNonconstObjPtr(0), 4);
  TEST_THROW(oc.getNonconstObjPtr(-1), SIOVOCB::InvalidOrdinalIndexError);
  TEST_THROW(oc.getNonconstObjPtr(1), SIOVOCB::InvalidOrdinalIndexError);
  ECHO(oc.removeObj(0));
  TEST_THROW(oc.getNonconstObjPtr(0), SIOVOCB::InvalidOrdinalIndexError);
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, getObjPtr_idx_invalid )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  //typedef StringIndexedOrderedValueObjectContainer<int>::Iterator Iterator; // unused
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);
  ECHO(oc.setObj("a", 4));
  TEST_EQUALITY_CONST(*oc.getObjPtr(0), 4);
  TEST_THROW(oc.getObjPtr(-1), SIOVOCB::InvalidOrdinalIndexError);
  TEST_THROW(oc.getObjPtr(1), SIOVOCB::InvalidOrdinalIndexError);
  ECHO(oc.removeObj(0));
  TEST_THROW(oc.getObjPtr(0), SIOVOCB::InvalidOrdinalIndexError);
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, removeObj_idx_invalid )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  //typedef StringIndexedOrderedValueObjectContainer<int>::Iterator Iterator; // unused
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);
  ECHO(oc.setObj("a", 4));
  TEST_EQUALITY_CONST(*oc.getObjPtr(0), 4);
  TEST_THROW(oc.removeObj(-1), SIOVOCB::InvalidOrdinalIndexError);
  TEST_THROW(oc.removeObj(1), SIOVOCB::InvalidOrdinalIndexError);
  TEST_EQUALITY_CONST(oc.numObjects(), 1);
  ECHO(oc.removeObj(0));
  TEST_EQUALITY_CONST(oc.numObjects(), 0);
  TEST_THROW(oc.removeObj(0), SIOVOCB::InvalidOrdinalIndexError);
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, getNonconstObjPtr_key_invalid )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  //typedef StringIndexedOrderedValueObjectContainer<int>::Iterator Iterator; // unused
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);
  ECHO(oc.setObj("a", 4));
  TEST_EQUALITY_CONST(*oc.getNonconstObjPtr("a"), 4);
  TEST_THROW(oc.getNonconstObjPtr("does_not_exist"), SIOVOCB::InvalidKeyError);
  ECHO(oc.removeObj("a"));
  TEST_THROW(oc.getNonconstObjPtr("a"), SIOVOCB::InvalidKeyError);
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, getObjPtr_key_invalid )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  //typedef StringIndexedOrderedValueObjectContainer<int>::Iterator Iterator; // unused
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);
  ECHO(oc.setObj("a", 4));
  TEST_EQUALITY_CONST(*oc.getObjPtr("a"), 4);
  TEST_THROW(oc.getObjPtr("does_not_exist"), SIOVOCB::InvalidKeyError);
  ECHO(oc.removeObj("a"));
  TEST_THROW(oc.getObjPtr("a"), SIOVOCB::InvalidKeyError);
}


TEUCHOS_UNIT_TEST( StringIndexedOrderedValueObjectContainer, removeObj_key_invalid )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  //typedef StringIndexedOrderedValueObjectContainer<int>::Iterator Iterator; // unused
  ECHO(StringIndexedOrderedValueObjectContainer<int> oc);
  ECHO(oc.setObj("a", 4));
  TEST_EQUALITY_CONST(*oc.getObjPtr("a"), 4);
  TEST_THROW(oc.removeObj("does_not_exist"), SIOVOCB::InvalidKeyError);
  TEST_EQUALITY_CONST(oc.numObjects(), 1);
  ECHO(oc.removeObj("a"));
  TEST_EQUALITY_CONST(oc.numObjects(), 0);
  TEST_THROW(oc.removeObj("a"), SIOVOCB::InvalidKeyError);
}

// ToDo: Test dangling object references!


} // namespace Teuchos
