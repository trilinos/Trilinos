#ifndef CTHULHU_UNIT_TEST_HELPERS_HPP
#define CTHULHU_UNIT_TEST_HELPERS_HPP

#include "Teuchos_UnitTestHelpers.hpp"

/** \brief Basic unit test creation macro for templated code on five template parameters. */
#define TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(TEST_GROUP, TEST_NAME, TYPE1, TYPE2, TYPE3, TYPE4, TYPE5) \
  template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5>         \
  class TEST_GROUP##_##TEST_NAME##_UnitTest : public Teuchos::UnitTestBase \
	{ \
  public: \
    TEST_GROUP##_##TEST_NAME##_UnitTest( \
      const std::string& type1Name, \
      const std::string& type2Name, \
      const std::string& type3Name, \
      const std::string& type4Name, \
      const std::string& type5Name \
       ) \
      :Teuchos::UnitTestBase( \
         std::string(#TEST_GROUP)+"_"+type1Name+"_"+type2Name+"_"+type3Name+"_"+type4Name+"_"+type5Name, #TEST_NAME ) \
    {} \
    void runUnitTestImpl( Teuchos::FancyOStream &out, bool &success ) const; \
    virtual std::string unitTestFile() const { return __FILE__; } \
    virtual long int unitTestFileLineNumber() const { return __LINE__; } \
  }; \
  \
  template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5>         \
  void TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1,TYPE2,TYPE3,TYPE4, TYPE5>::runUnitTestImpl( \
    Teuchos::FancyOStream &out, bool &success ) const \

/** \brief Template instantiation for five templated types. */
#define TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(TEST_GROUP, TEST_NAME, TYPE1, TYPE2, TYPE3, TYPE4, TYPE5) \
  \
  template class TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3, TYPE4, TYPE5 >; \
  TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3, TYPE4, TYPE5 >     \
  instance_##TEST_GROUP##_##TYPE1##_##TYPE2##_##TYPE3##_##TYPE4##_##TYPE5##_##TEST_NAME##_UnitTest(#TYPE1,#TYPE2,#TYPE3,#TYPE4,#TYPE5);

#endif  // CTHULHU_UNIT_TEST_HELPERS_HPP
