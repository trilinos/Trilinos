Instructions for converting a unit test to use ETI.

1) At top of file delete

   #include "MueLu_UseDefaultTypes.hpp"

   and add 

   #include <Teuchos_ScalarTraits.hpp>


2)
   Move
     #include <MueLu_UseShortNames.hpp>

   into unit test bodies.  Add line

     MUELU_LIMIT_EPETRA_TESTING_SCOPE(Scalar,GlobalOrdinal,Node);

3) Replace
     TEUCHOS_UNIT_TEST
   by
     TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL or
     TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL if no scalar template arg

4) At bottom of file, add

   #define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(group,test 1,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(group,test 2,Scalar,LO,GO,Node) \
      etc.

   #include <MueLu_ETI_4arg.hpp>
