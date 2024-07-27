// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Map.hpp"
#include "Tpetra_Core.hpp"
#include "Teuchos_OrdinalTraits.hpp"

#include <vector>
#include <unordered_map>

// Test program exhibiting github issue #2431
// Submitted by Denis Ridzal, 3/20/18
// Modified and augmentd by Karen Devine, 3/21/18

//////////////////////////////////////////////////////////////////////////////
// Tie-break function that assigns shared IDs to the lowest process that
// has a copy.
namespace {
template <typename LO, typename GO>
class GreedyTieBreak :
      public Tpetra::Details::TieBreak<LO,GO>
{
public:
  GreedyTieBreak() { }

  virtual bool mayHaveSideEffects() const {
    return true;
  }

  virtual std::size_t selectedIndex(
    GO GID,
    const std::vector<std::pair<int,LO> > & pid_and_lid) const
  {
    // always choose index of pair with smallest pid
    auto numLids = pid_and_lid.size();
    decltype(numLids) idx = 0;
    auto minpid = pid_and_lid[0].first;
    decltype(minpid) minidx = 0;
    for (idx = 0; idx < numLids; ++idx) {
      if (pid_and_lid[idx].first < minpid) {
        minpid = pid_and_lid[idx].first;
        minidx = idx;
      }
    }
    return minidx;
  }
};

}

//////////////////////////////////////////////////////////////////////////////
// Given input IDs vecP0, vecP1, vecP2, vecP3, build a (probably overlapping)
// map with these IDs on the respective processors P0-P3.
// Then create one-to-one maps from the overlapping map, with and without
// use of the tie-break function.
// Compare the number of unique IDs in the three maps; the test passes if
// the number of unique IDs matches.

template <typename LO, typename GO, typename NO>
int runTest(
  const char *message,
  std::ostream &outStream,   // allows varying levels of output
  Teuchos::RCP<const Teuchos::Comm<int>> &comm,
  std::vector<GO> &vecP0,
  std::vector<GO> &vecP1,
  std::vector<GO> &vecP2,
  std::vector<GO> &vecP3
)
{
  int errorFlag = 0;

  try {
    Teuchos::Array<GO> arrP0(vecP0), arrP1(vecP1), arrP2(vecP2), arrP3(vecP3);

    typedef Tpetra::Map<LO,GO,NO> map_t;
    Teuchos::RCP<map_t> overlapMap;
    GreedyTieBreak<LO,GO> greedy_tie_break;

    auto pid = comm->getRank();
    auto dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

    if (pid == 0) {
      overlapMap = Teuchos::rcp(new map_t(dummy, arrP0(), 0, comm));
    }
    else if (pid == 1) {
      overlapMap = Teuchos::rcp(new map_t(dummy, arrP1(), 0, comm));
    }
    else if (pid == 2) {
      overlapMap = Teuchos::rcp(new map_t(dummy, arrP2(), 0, comm));
    }
    else if (pid == 3) {
      overlapMap = Teuchos::rcp(new map_t(dummy, arrP3(), 0, comm));
    }

    std::cout << message
              << ": Before Tpetra::createOneToOne on " << "Proc "
              << overlapMap->getComm()->getRank()
              << "; nGids = " << overlapMap->getLocalNumElements() << "\n";

    auto myidx_o =  overlapMap->getMyGlobalIndices();
    outStream << "IDS ON PROC " << comm->getRank() << ": ";
    for (std::size_t idx=0; idx<myidx_o.size(); ++idx) {
      outStream << myidx_o[idx] << ", ";
    }
    outStream << "\n";

    // Create non-overlap Map with TieBreak function
    Teuchos::RCP<const map_t> nonOverlapMapTB =
             Tpetra::createOneToOne<LO,GO,NO>(overlapMap, greedy_tie_break);

    std::cout << message
              << ": After Tpetra::createOneToOne with TieBreak on Proc "
              << overlapMap->getComm()->getRank()
              << "; nGids = " << nonOverlapMapTB->getLocalNumElements() << "\n";

    auto myidx_notb =  nonOverlapMapTB->getMyGlobalIndices();
    outStream << "IDS ON PROC " << comm->getRank() << ": ";
    for (std::size_t idx=0; idx<myidx_notb.size(); ++idx) {
      outStream << myidx_notb[idx] << ", ";
    }
    outStream << "\n";

    // Create non-overlap Map without TieBreak function
    Teuchos::RCP<const map_t> nonOverlapMap =
             Tpetra::createOneToOne<LO,GO,NO>(overlapMap);

    std::cout << message
              << ": After Tpetra::createOneToOne without TieBreak on Proc "
              << overlapMap->getComm()->getRank()
              << "; nGids = " << nonOverlapMap->getLocalNumElements() << "\n";

    auto myidx_no =  nonOverlapMap->getMyGlobalIndices();
    outStream << "IDS ON PROC " << comm->getRank() << ": ";
    for (std::size_t idx=0; idx<myidx_no.size(); ++idx) {
      outStream << myidx_no[idx] << ", ";
    }
    outStream << "\n";

    // Discover how many unique GIDs we had initially
    std::unordered_map<GO,int> uniqueGids;

    for (auto i = arrP0.begin(); i != arrP0.end(); i++) {
      if (uniqueGids.find(*i) != uniqueGids.end())
        uniqueGids[*i]++;
      else
        uniqueGids[*i] = 1;
    }
    for (auto i = arrP1.begin(); i != arrP1.end(); i++) {
      if (uniqueGids.find(*i) != uniqueGids.end())
        uniqueGids[*i]++;
      else
        uniqueGids[*i] = 1;
    }
    for (auto i = arrP2.begin(); i != arrP2.end(); i++) {
      if (uniqueGids.find(*i) != uniqueGids.end())
        uniqueGids[*i]++;
      else
        uniqueGids[*i] = 1;
    }
    for (auto i = arrP3.begin(); i != arrP3.end(); i++) {
      if (uniqueGids.find(*i) != uniqueGids.end())
        uniqueGids[*i]++;
      else
        uniqueGids[*i] = 1;
    }

    GO ncopies = 0;
    for (auto i = uniqueGids.begin(); i != uniqueGids.end(); i++)
      if (i->second > 1) ncopies += (i->second - 1);

    if (pid == 0) {
      std::cout << "\n\n" << message
                << ": Before Tpetra::createOneToOne, there are "
                << uniqueGids.size() << " ids, with "
                << ncopies << " copies.\n";
      std::cout << message
                << ": After Tpetra::createOneToOne with TieBreak, there are "
                << nonOverlapMapTB->getGlobalNumElements() << " ids.\n";
      std::cout << message
                << ": After Tpetra::createOneToOne without TieBreak, there are "
                << nonOverlapMap->getGlobalNumElements() << " ids.\n";
      std::cout << "\n\n";
    }

    // Check the results; number of unique GIDs should be the same in all Maps.
    if (uniqueGids.size() != nonOverlapMapTB->getGlobalNumElements())
      errorFlag = -1;
    if (uniqueGids.size() != nonOverlapMap->getGlobalNumElements())
      errorFlag = -1;
  }
  catch (std::logic_error &err) {
    outStream << err.what() << "\n";
    errorFlag = -1000;
  };

  return errorFlag;
}

//////////////////////////////////////////////////////////////////////////////

int main(int narg, char *arg[]) {

  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type NO;

  Tpetra::ScopeGuard tpetraScope(&narg, &arg);
  auto comm = Tpetra::getDefaultComm();

  if (comm->getSize() != 4) {
    if (comm->getRank() == 0)
      std::cout << "TEST FAILED: This test is written for four processes only. "
                << "You are running on " << comm->getSize() << " processes."
                << std::endl;
    return EXIT_FAILURE;
  }

  int errorFlag  = 0;

  // This little trick lets us print to std::cout only
  // if a (dummy) command-line argument is provided.
  int iprint     = narg - 1;
  Teuchos::oblackholestream bhs; // outputs nothing
  std::ostream &outStream(iprint > 0 ? std::cout : bhs);

  // Sparse test that uses hash tables in directory
  {
    std::vector<GO> vecP0 =
      { 0, 1, 3, 4, 9, 10, 12, 13, 18, 19, 21, 22, 27, 31, 36, 37, 38, 39, 40,
        41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
        59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76,
        77, 78, 79, 80, 81, 82};
    std::vector<GO> vecP1 =
      { 1, 2, 4, 5, 10, 11, 13, 14, 19, 20, 22, 23, 28, 32, 37, 46, 55, 56,
        58, 67, 76, 77, 79, 996, 997, 998, 999, 1000, 1001, 1002, 1004, 1005,
        1006, 1007, 1009, 1010, 1013, 1015, 1016, 1017, 1018, 1019, 1020, 1021,
        1022, 1023, 1024, 1025, 1026, 1027, 1028, 1030, 1031, 1034, 1036, 1037,
        1038, 1039, 1040, 1041, 1042};
    std::vector<GO> vecP2 =
      { 3, 4, 6, 7, 12, 13, 15, 16, 21, 22, 24, 25, 29, 33, 42, 54, 58, 59,
        60, 75, 79, 80, 81, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964,
        1967, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980,
        1981, 1982, 1983, 1984, 1985, 1988, 1991, 1992, 1993, 1994, 1995, 1996,
        1997, 1998, 1999, 2000, 2001, 2002};
    std::vector<GO> vecP3 =
      { 4, 5, 7, 8, 13, 14, 16, 17, 22, 23, 25, 26, 30, 34, 58, 79, 1002,
        1018, 1019, 1020, 1039, 1040, 1041, 1957, 1975, 1976, 1978, 1996, 1997,
        1999, 2917, 2918, 2919, 2920, 2921, 2922, 2924, 2927, 2930, 2933, 2935,
        2936, 2937, 2938, 2939, 2940, 2941, 2942, 2943, 2944, 2945, 2948, 2951,
        2954, 2956, 2957, 2958, 2959, 2960, 2961, 2962};

    errorFlag += runTest<LO,GO,NO>("sparseTest", outStream, comm,
                                   vecP0, vecP1, vecP2, vecP3);

    // Make sure it works if some process has no data
    std::vector<GO> empty = { } ;
    errorFlag += runTest<LO,GO,NO>("sparseTestEmptyP2", outStream, comm,
                                   vecP0, vecP1, empty, vecP3);
  }

  // Dense test that does not use hash tables in directory.
  // Keep same number of IDs and structure of overlap, but
  // narrow the range of global ID values so that processors more than
  // 0.1 * (max ID - min ID).
  {
    std::vector<GO> vecP0 =
      { 0, 1, 3, 4, 9, 10, 12, 13, 18, 19, 21, 22, 27, 31, 36, 37, 38, 39, 40,
        41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
        59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76,
        77, 78, 79, 80, 81, 82};
    std::vector<GO> vecP1 =
      { 1, 2, 4, 5, 10, 11, 13, 14, 19, 20, 22, 23, 28, 32, 37, 46, 55, 56,
        58, 67, 76, 77, 79, 396, 397, 398, 399, 400, 401, 402, 404, 405,
        406, 407, 409, 410, 413, 415, 416, 417, 418, 419, 420, 421,
        422, 423, 424, 425, 426, 427, 428, 430, 431, 434, 436, 437,
        438, 439, 440, 441, 442};
    std::vector<GO> vecP2 =
      { 3, 4, 6, 7, 12, 13, 15, 16, 21, 22, 24, 25, 29, 33, 42, 54, 58, 59,
        60, 75, 79, 80, 81, 557, 558, 559, 560, 561, 562, 563, 564,
        567, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580,
        581, 582, 583, 584, 585, 588, 591, 592, 593, 594, 595, 596,
        597, 598, 599, 600, 601, 602};
    std::vector<GO> vecP3 =
      { 4, 5, 7, 8, 13, 14, 16, 17, 22, 23, 25, 26, 30, 34, 58, 79, 402,
        418, 419, 420, 439, 440, 441, 557, 575, 576, 578, 596, 597,
        599, 617, 618, 619, 620, 621, 622, 624, 627, 630, 633, 635,
        636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 648, 651,
        654, 656, 657, 658, 659, 660, 661, 662};

    errorFlag += runTest<LO,GO,NO>("denseTest", outStream, comm,
                                   vecP0, vecP1, vecP2, vecP3);

    // Make sure it works if some process has no data
    std::vector<GO> empty = { } ;
    errorFlag += runTest<LO,GO,NO>("denseTestEmptyP2", outStream, comm,
                                   vecP0, vecP1, empty, vecP3);
  }

  if (errorFlag != 0) {
    std::cout << "End Result: TEST FAILED" << std::endl;
    return EXIT_FAILURE;
  }
  else {
    std::cout << "End Result: TEST PASSED" << std::endl;
    return EXIT_SUCCESS;
  }
}
