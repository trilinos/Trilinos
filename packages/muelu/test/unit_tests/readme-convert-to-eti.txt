Instructions for converting a unit test to use ETI.

1) At top of file delete lines

   #include "MueLu_UseDefaultTypes.hpp"
   and
   #include <MueLu_UseShortNames.hpp>

   and add

   #include <Teuchos_ScalarTraits.hpp>

2)
   Move

   For each unit test body, add lines

#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

3) Replace
     TEUCHOS_UNIT_TEST
   by
     TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL or
     TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL if no scalar template arg

4) At bottom of file, add

   #define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(group,test 1,Scalar,LocalOrdinal,GlobalOrdinal,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(group,test 2,Scalar,LocalOrdinal,GlobalOrdinal,Node) \
      etc.

#    include <MueLu_ETI_4arg.hpp>
