#include <Tpetra_MultiVectorFiller.hpp>

namespace Tpetra {
  namespace Test {

    void
    testSortAndMergeIn ()
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Tpetra::Details::sortAndMergeIn;

      Array<int> allEntries (5);
      for (Array<int>::size_type k = 0; k < 5; ++k) {
        allEntries[k] = 2 * k;
      }
      Array<int> newEntries (4);
      newEntries[0] = -1;
      newEntries[1] = 3;
      newEntries[2] = 3;
      newEntries[3] = 11;

      ArrayView<int> result =
        sortAndMergeIn<int> (allEntries, allEntries.view (0, 5), newEntries());
      TEUCHOS_TEST_FOR_EXCEPTION(
        result.size() != 8, std::logic_error,
        "Returned ArrayView should have size 8, but instead has size "
        << result.size() << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        allEntries.size() < 8, std::logic_error,
        "Input/output Array argument should have size at least 8, but instead has "
        "size " << allEntries.size() << ".");

      bool success = true;
      ArrayView<int>::size_type firstBadIndex = -1; // size_type is signed
      for (ArrayView<int>::size_type k = 0; k < result.size(); ++k) {
        if (result[k] != allEntries[k]) {
          success = false;
          firstBadIndex = k;
          break;
        }
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        success, std::logic_error, "Returned ArrayView and the input/output "
        "Array argument don't match.  First nonmatching array index is "
        << firstBadIndex << ".");

      Array<int> expectedEntries (8);
      expectedEntries[0] = -1;
      expectedEntries[1] = 0;
      expectedEntries[2] = 2;
      expectedEntries[3] = 3;
      expectedEntries[4] = 4;
      expectedEntries[5] = 6;
      expectedEntries[6] = 8;
      expectedEntries[7] = 11;
      for (ArrayView<int>::size_type k = 0; k < result.size(); ++k) {
        if (expectedEntries[k] != result[k]) {
          success = false;
          firstBadIndex = k;
          break;
        }
      }
      TEUCHOS_TEST_FOR_EXCEPTION(success, std::logic_error, "Returned ArrayView "
        "and the expected results don't match.  First nonmatching array index is "
        << firstBadIndex << ".");
    }

  } // namespace Test
} // namespace Tpetra
