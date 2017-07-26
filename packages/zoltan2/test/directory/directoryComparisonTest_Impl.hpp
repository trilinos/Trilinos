// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <unordered_map>
#include <unordered_set>
#include <Zoltan2_TPLTraits.hpp>

#ifdef HAVE_MPI

#if !defined(TEMP_TRIAL_USER_ARRAY_TYPE) || !defined(CONVERT_DIRECTORY_TPETRA)
#include "Zoltan2_Directory_Impl.hpp"
#endif

#ifdef CONVERT_DIRECTORY_ORIGINAL
#include <zoltan_dd_cpp.h>
#endif

// This is temporary
#if !defined(TEMP_TRIAL_USER_ARRAY_TYPE) || !defined(CONVERT_DIRECTORY_TPETRA)

enum Mode {
  Replace = 0,
  Add,
  Aggregate
};

// Temporary measure to determine find did not get gid - used for validating
// the remove mechanism - how the system behaves when we request an unfound
// id will need some discussion. This -1 works because right now we just set
// positive values.
#define NOT_FOUND_VALUE -1

template <typename gid_t,typename lid_t,typename user_t>
class IDs {
public:
  IDs(gid_t totalIds, gid_t idBase_, int idStride_,
    Teuchos::RCP<const Teuchos::Comm<int> > &comm_, Mode mode_) :
      idBase(idBase_), contiguous(idStride_ == 1), idStride(idStride_),
      comm(comm_), maxPrintSize(10), mode(mode_)
  {
    for(int n = 0; n < totalIds; ++n) {
      bool bAtLeastOne = false;
      gid_t gid = idBase + n * idStride;
      for(int proc = 0; proc < comm->getSize(); ++proc) {
        if(writeID(gid, proc)) {
          bAtLeastOne = true;
        }
      }
      if(!bAtLeastOne) {
        throw std::logic_error("The writeID method must genrate an algorithm"
          " which returns true for at least 1 proce for all IDs. If this is"
          " changed then the useReadID algorithm should not return true for any"
          " of those missing IDs. This is a sanity check.");
      }
    }
    // create ids list - the count here is arbitrary and sets the complexity
    // of the test
    for(int n = 0; n < totalIds; ++n) {
      id_t gid = idBase + n * idStride;

      for(int i = 0; i < writeID(n, comm->getRank()); ++i) {
        writeIds.push_back(gid);
      }

      for(int i = 0; i < readID(n, comm->getRank()); ++i) {
        readIds.push_back(gid);
      }

      for(int i = 0; i < removeID(n, comm->getRank()); ++i) {
        removeIds.push_back(gid);
      }
    }

    // set up all the initial writes
    writeUser.resize(writeIds.size());
    for(size_t n = 0; n < writeUser.size(); ++n) {
      writeUser[n] = getInitialValue(writeIds[n], comm->getRank());
    }
  }

  bool trueForAtLeastOneProc(int index, int rank) const {
    return ((index % comm->getSize()) == rank);
  }

  // some arbitrary subset
  bool subset1(int index, int rank) const {
    return ((index % (comm->getSize()+1)) == rank);
  }

  // some arbitrary subset
  bool subset2(int index, int rank) const {
    return ((index % (comm->getSize()+3)) == rank);
  }

  // some arbitrary subset
  bool subset3(int index, int rank) const {
    return ((index % (comm->getSize()+7)) == rank);
  }

  int writeID(int index, int rank) const {
    // this method should be guaranteed to return true for 1 of the procs
    // it can optionally return true for more than one to make a nice mix
    // the return is the number of repeats as we want the directory to work
    // in such a case
    int some_count = 3;
    // Note Tpetra requires unique IDs so readID must be 1 for it to work
#ifdef CONVERT_DIRECTORY_TPETRA
    some_count = 1;
#endif
    if(trueForAtLeastOneProc(index, rank) || subset1(index, rank)) {
      return subset1(index, rank) ? some_count : 1; // return different to test repeats
    }
    else {
      return 0; // do not write this ID
    }
  }

  int readID(int index, int rank) const {
    // this method can be anything - things should work even if this is empty;
    int some_count = 3;
    // Note Tpetra requires unique IDs so readID must be 1 for it to work
#ifdef CONVERT_DIRECTORY_TPETRA
    some_count = 1;
#endif
    if(index < 20) { // to make sure we always get some good duplicate coverage
      return some_count;
    }
    if(subset1(index, rank)) {
      return subset2(index, rank) ? some_count : 1; // return different to test repeats
    }
    else {
      return 0; // do not read this ID
    }
  }

  int removeID(int index, int rank) const {
    return 0; // TODO - This is not working so disable any removed for now

    // this should include some readID values but not all, but also other
    // values for good testing - we want to remove more than we wrote locally
    // and check for proper removal
    // Tpetra not supporting remove right now
#ifdef CONVERT_DIRECTORY_TPETRA
    return 0;
#endif
    int some_count = 3;
    if(subset2(index, rank)) {
      return subset3(index, rank) ? some_count : 1; // return different to test repeats
    }
    else {
      return 0; // do not remove this ID
    }
  }

  bool removedIDGlobally(gid_t gid) const {
    int index = (gid-idBase)/idStride; // convert back to 0,1,2,...indexing
    for(int proc = 0; proc < comm->getSize(); ++proc) {
      if(removeID(index, proc)) {
        return true;
      }
    }
    return false;
  }

  int sharedCount(gid_t gid) const {
    gid_t index = (gid-idBase)/idStride; // convert back to 0,1,2,...indexing
    int count = 0;
    for(int proc = 0; proc < comm->getSize(); ++proc) {
      count += writeID(index, proc);
    }
    return count;
  }

  void printVector(const std::vector<gid_t> &printIds, std::string name) const {
    std::cout << "  " << printIds.size() << " " << name << ": ";
    for (size_t i = 0; i < printIds.size() && i < maxPrintSize; i++) {
      std::cout << printIds[i] << " ";
    }
    if(printIds.size() > maxPrintSize) {
      std::cout << "... ";
    }
    std::cout << std::endl;
  }

  void print(std::string name)
  {
    for(int proc = 0; proc < comm->getSize(); ++proc) {
      comm->barrier();
      if(proc == comm->getRank()) {
        if(proc == 0) {
          std::cout << std::endl <<
            "############ Output: " << name
            << " ############" << std::endl;
        }
        std::cout << "Rank: " << proc << std::endl;
        printVector(writeIds,  "Write Ids"  );
        printVector(readIds,   "Read Ids"   );
        printVector(removeIds, "Remove Ids" );

#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
        for(size_t n = 0; n < writeIds.size(); ++n) {
          std::cout << "  Write array for GID (" << writeIds[n] << "): ";
          for(size_t a = 0; a < writeUser[n].size(); ++a) {
            std::cout << writeUser[n][a] << " ";
          }
          std::cout << std::endl;
        }
        for(size_t n = 0; n < readIds.size(); ++n) {
          std::cout << "  Read array for GID (" << readIds[n] << "): ";
          for(size_t a = 0; a < readUser[n].size(); ++a) {
            std::cout << readUser[n][a] << " ";
          }
          std::cout << std::endl;
        }
#endif
        if(proc == comm->getSize()-1) {
          std::cout << std::endl;
        }
      }
    }
    comm->barrier();
  }

#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
  // for development and eventually must go away
  // this is used by the tests to generate the variable array sizes
  // but also used by find in one place that needs to be resolved
  size_t test_create_array_length(gid_t gid, int proc) const {
#if defined(TEMP_USED_FIXED_SIZE) || defined(TEMP_CONSTANT_ARRAY_SIZE)
    return CONSTANT_ARRAY_SIZE;
#else
    if(mode == Aggregate) {
      // we want different procs with same gid to generate different lengths
      return (gid%7)+proc; // varies per gid and proc
    }
    else {
      // we want different procs with same gid to generate same lengths
      // at least for ADD ... I need to make this different for Replace to make
      // the test have better coverage of the possible usage cases
      return (gid%7); // 1..7 inclusive
    }
#endif
  }
#endif

  void ZoltanDirectoryTest();
  bool evaluateTests(const std::string& name) const;

private:
  user_t getInitialValue(gid_t gid, int rank) const;
  user_t getExpectedValue(gid_t gid) const;

  gid_t idBase;
  bool contiguous;              // Flag indicating whether IDs are contiguous
  int idStride;
  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  const size_t maxPrintSize;       // print only this number for debugging
  Mode mode;
  std::vector<gid_t> writeIds;   // Ids generated on this proc
  std::vector<user_t> writeUser; // The data we initially updated with
  std::vector<gid_t> readIds;    // Ids to read
  std::vector<user_t> readUser;  // User data we read
  std::vector<gid_t> removeIds;  // Ids to remove
};

template <typename gid_t,typename lid_t, typename user_t>
user_t IDs<gid_t,lid_t,user_t>::getInitialValue(gid_t gid, int rank) const {
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
    // determine length of array
    size_t modLength = test_create_array_length(gid, rank);
    user_t array(modLength);
    for(size_t n = 0; n < array.size(); ++n) {
      switch(mode) {
        case Mode::Replace:
          array[n] = gid + 3; // all elements same for now
          break;
        case Mode::Add:
          array[n] = 1; // all elements equal n - will sum
          break;
        case Mode::Aggregate:
          // Now we want some mix so that, for example, gid 10 will have
          // different but overlapping values for each array element n
          // For example proc 1 could update with gid=10 array={5,7,9}
          //       while proc 2 could update with gid=10 array={3,5,7}
          // Then we expect the result to be {3,5,7,9}
          array[n] = n + rank*2; // this creates overlapping regions
          break;
      }
    }
    return array;
#else
    switch(mode) {
      case Mode::Replace:
        return gid + 3; // arbitrary - just make it something good to test
        break;
      case Mode::Add:
        return 1; // we will be testing if they sum to shared count
        break;
      case Mode::Aggregate:
        throw std::logic_error("Aggregate requested for non array setup.");
        break;
    }
#endif
}

template <typename gid_t,typename lid_t, typename user_t>
user_t IDs<gid_t,lid_t,user_t>::getExpectedValue(gid_t gid) const {
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
  switch(mode) {
    case Mode::Add:
      return user_t(test_create_array_length(gid, comm->getRank()),
        sharedCount(gid)); // should have summed
      break;
    case Mode::Replace:
      // TODO - need to consider if another proc replaces our original
      // value with a new value and test for that - the current testing assigns
      // a fixed result per gid but it could be more complex
      return getInitialValue(gid, comm->getRank());
      break;
    case Mode::Aggregate:
      // for this we need the union of all the procs
      user_t final_aggregated;
      // loop through all possible updaters
      for(int proc = 0; proc < comm->getSize(); ++proc) {
        int index = (gid-idBase)/idStride; // convert back to 0,1,2,...indexing
        if(writeID(index, proc)) { // did this proc update?
          user_t proc_input = getInitialValue(gid, proc); // get original input
          for(size_t i = 0; i < proc_input.size(); ++i) { // scan elements
            auto val = proc_input[i]; // get the array element
            if(final_aggregated.size() == 0 ||
              val > final_aggregated[final_aggregated.size()-1]) {
              final_aggregated.push_back(val); // add first element or tail end
            }
            else { // loop and insert to keep ordering unless we find a match
              for(auto itr = final_aggregated.begin();
                itr != final_aggregated.end(); ++itr) {
                if((*itr) == val) {
                  break; // don't add - already added
                }
                else if((*itr) > val) {
                  final_aggregated.insert(itr, val);
                  break; // insert here to keep ordering
                }
              }
            }
          }
        }
      }
      return final_aggregated; // temp
      break;
  }
#else
  switch(mode) {
    case Mode::Add:
      return sharedCount(gid); // should have summed
      break;
    case Mode::Replace:
      // TODO - need to consider if another proc replaces our original
      // value with a new value and test for that - the current testing assigns
      // a fixed result per gid but it could be more complex
      return getInitialValue(gid, comm->getRank());
      break;
    case Mode::Aggregate:
      throw std::logic_error("Unexpected aggregate mode for non array.");
      break;
  }
#endif
}

template <typename gid_t,typename lid_t, typename user_t>
bool IDs<gid_t,lid_t,user_t>::evaluateTests(const std::string& name) const {

#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
  #define READ_VALUE readUser[i][arrayIndex]
  #define EXPECTED_VALUE expected_value[arrayIndex]
#else
  #define READ_VALUE readUser[i]
  #define EXPECTED_VALUE expected_value
#endif

  // check the results
  bool pass = true;
  for(int proc = 0; proc < comm->getSize(); ++proc) {
    bool passRank = true;
    comm->barrier();
    if(proc == comm->getRank()) {
      for(size_t i = 0; i < readIds.size(); ++i) {
        id_t gid = readIds[i];

        // verify removal - eventually I think we may have the directory
        // just throw an error if we try to find on a gid which was removed.
        // Or some return value will indicate. For now we pass in all values
        // of NOT_FOUND_VALUE (-1) which will still be set to that if not found.
        if(removedIDGlobally(gid)) {
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
          // should be an array of size 1 with element NOT_FOUND_VALUE
          if(readUser[i].size() != 1 || readUser[i][0] != NOT_FOUND_VALUE) {
            passRank = false;
            std::cout << "Removed array for gid: " << gid <<
              " but something set the user data which is incorrect. "
              "Remove FAILED." << std::endl;
          }
          else {
            std::cout << "Successfully removed array for gid: "
              << gid << std::endl;
          }
#else
          if(READ_VALUE != NOT_FOUND_VALUE) {
            passRank = false;
            std::cout << "Removed gid: " << gid <<
              " but got value: " << readUser[i] << " when we expected"
              " to read the unset value of " << NOT_FOUND_VALUE << ". " <<
              " This is incorrect. Remove FAILED." << std::endl;
          }
          else {
            std::cout << "Successfully removed gid: " << gid << std::endl;
          }
#endif
        }

        if(!passRank) {
          break; // we don't need to check further
        }

        user_t expected_value = getExpectedValue(gid);
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
        // first validate the array length is correct
        if(readUser[i].size() != expected_value.size()) {
          std::cout << "Rank : " << proc << " array size is incorrect for gid: "
            << gid << ". Expected size: " << expected_value.size() <<
            " and got size: " << readUser[i].size() << std::endl;
          passRank = false;
          break;
        }
        // now loop the elements and validate each individual element
        for(size_t arrayIndex = 0;
          arrayIndex < readUser[i].size() && passRank; ++arrayIndex) {
          if(READ_VALUE != EXPECTED_VALUE) {
            passRank = false;
            std::cout << "Failed read for global ID: " << gid
              << ". Expected: " << EXPECTED_VALUE << " at array index " <<
              arrayIndex << ". Got: " << READ_VALUE << std::endl;
          }
        }
#else
        if(READ_VALUE != EXPECTED_VALUE) {
          passRank = false;
          std::cout << "Failed read for global ID: " << gid
            << ". Expected: " << EXPECTED_VALUE <<
            " Got: " << READ_VALUE << std::endl;
        }
#endif
      }

      if(!passRank) {
        std::cout << "Checked rank: " << comm->getRank() << " with nIds: " <<
          readIds.size() << " which " << "FAILED" << std::endl;
      }

      if(!passRank) {
        pass = false;
      }
    }
  }
  comm->barrier();

  return pass;
}

template <typename gid_t,typename lid_t, typename user_t>
void IDs<gid_t,lid_t,user_t>::ZoltanDirectoryTest()
{
#ifdef CONVERT_DIRECTORY_ORIGINAL

  if (writeIds.size() > std::numeric_limits<int>::max())
    throw std::runtime_error("Problem too large for zoltan");
  int nIdEnt = Zoltan2::TPL_Traits<ZOLTAN_ID_PTR, id_t>::NUM_ID;
  MPI_Comm mpicomm = Teuchos::getRawMpiComm(*comm); // TODO comm fixes
  Zoltan2_Directory_Clock constructClock("construct", comm);
  Zoltan_DD zz(mpicomm, nIdEnt, 0, sizeof(user_t), writeIds.size(), 0);
  constructClock.complete(); // to match with new directory clocks...

#else

  const bool bUseLocalIDs = false;
#ifdef CONVERT_DIRECTORY_RELIC
  Zoltan2::Zoltan2_Directory<gid_t,lid_t,user_t>
    zz(comm, bUseLocalIDs, writeIds.size(), 0, contiguous);
#else
  Zoltan2::Zoltan2_Directory<gid_t,lid_t,user_t>
    zz(comm, bUseLocalIDs, 0, contiguous);
#endif

#endif // CONVERT_DIRECTORY_ORIGINAL

  // Create data space for reading user

  // TODO - Eventually should just be like following line
  // However using NOT_FOUND_VALUE to track success/failed remove events
  // Need to decide how to handle error messaging for find called on removed id
  //   readUser.resize(readIds.size());
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
  // initially arrays of size 1 with NOT_FOUND_VALUE
  readUser = std::vector<user_t>(readIds.size(), user_t(1, NOT_FOUND_VALUE));
#else
  // initially all set to NOT_FOUND_VALUE
  readUser = std::vector<user_t>(readIds.size(), NOT_FOUND_VALUE);
#endif


#ifdef CONVERT_DIRECTORY_ORIGINAL
  // CONVERT_DIRECTORY_ORIGINAL does not have a working add feature
  Zoltan2_Directory_Clock updateClock("update", comm);
  zz.Update(reinterpret_cast<ZOLTAN_ID_PTR>(&(writeIds[0])), NULL,
    writeIds.size() ? reinterpret_cast<char*>(&(writeUser[0])) : NULL,
    NULL, writeIds.size());
  updateClock.complete(); // to match with new directory clocks...

  Zoltan2_Directory_Clock findClock("find", comm);
  zz.Find(reinterpret_cast<ZOLTAN_ID_PTR>(&(readIds[0])), NULL,
    readIds.size() ? reinterpret_cast<char*>(&(readUser[0])) : NULL,
    NULL, readIds.size(), NULL);
  findClock.complete(); // to match with new directory clocks...

  Zoltan2_Directory_Clock removeClock("remove", comm);
  zz.Remove(reinterpret_cast<ZOLTAN_ID_PTR>(&(removeIds[0])),removeIds.size());
  removeClock.complete(); // to match with new directory clocks...
#else
  // these modes all have working add features
  std::vector<lid_t> ignore_lid; // TODO: decide how to best handle this API
  std::vector<int> ignore_int;   // TODO: decide how to best handle this API

  // convert generic mode to Zoltan2Directory - this awkward step exists because
  // the original mode doesn't have a directory - will improve.
  typedef Zoltan2::Zoltan2_Directory<gid_t,lid_t,user_t> directory_t;
  auto directoryMode = directory_t::Update_Mode::Add;
  switch(mode) {
    case Add:
      directoryMode = directory_t::Update_Mode::Add;
      break;
    case Replace:
      directoryMode = directory_t::Update_Mode::Replace;
      break;
    case Aggregate:
      directoryMode = directory_t::Update_Mode::Aggregate;
      break;
  }

  zz.update(writeIds, ignore_lid, writeUser, ignore_int, directoryMode);
  zz.remove(removeIds);
  zz.find(readIds, ignore_lid, readUser, ignore_int, ignore_int);
#endif
  comm->barrier();
}

class TestManager {
  public:
    TestManager(Teuchos::RCP<const Teuchos::Comm<int> > comm_, int totalIds_) :
      comm(comm_), totalIds(totalIds_), success(true) {}

    template<typename gid_t, typename lid_t, typename user_t>
    void runTest(const std::string& name, Mode mode,
      size_t idBase, int idStride) {
      IDs<gid_t,lid_t,user_t> myIds(totalIds, idBase, idStride, comm, mode);
      std::string base_message = name + " (" + mode_to_string(mode) + ") ";
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
      base_message = "Array User Type: " + base_message;
#else
      base_message = "Simple User Type: " + base_message;
#endif
      myIds.ZoltanDirectoryTest();

      // evaluate it correct
      bool test_passsed = myIds.evaluateTests(name);

      // print output
      // myIds.print(name);

      for(int proc = 0; proc < comm->getSize(); ++proc) {
        comm->barrier();
        if(proc == comm->getRank()) {
          if(!test_passsed) {
            std::cout << "FAILED: " << base_message << " on rank "
              << comm->getRank() << std::endl;
            success = false; // failed
          }
          else if(comm->getRank() == 0) {
            std::cout << "Passed: " << base_message << std::endl;
          }
        }
      }
    }

    bool is_success() const { return success; }

  private:
    std::string mode_to_string(Mode mode) const {
      switch(mode) {
        case Mode::Add: return "Add"; break;
        case Mode::Replace: return "Replace"; break;
        case Mode::Aggregate: return "Aggregate"; break;
        default: throw std::logic_error("Bad mode."); break;
      }
    }
    Teuchos::RCP<const Teuchos::Comm<int> > comm;
    int totalIds;
    bool success;
};

// convert totalIds to a loop count so total test time is of order 30s
// this is only for performance testing
int getTotalLoops(int totalIds) {
  switch(totalIds) {
    case 1:        return 10000; break;
    case 10:       return 10000; break;
    case 100:      return 10000; break;
    case 1000:     return 10000; break;
    case 10000:    return  3000; break;
    case 100000:   return   500; break;
    case 1000000:  return    50; break;
    case 10000000: return     5; break;
    default: throw std::logic_error("Not set up."); break;
  }
}

int runDirectoryComparisonTest(int narg, char **arg) {
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
  typedef std::vector<int> user_t; // user data since we are counting occurrences
#else
  typedef int user_t;
#endif

  Teuchos::GlobalMPISession mpiSession(&narg,&arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  // Note that setting this to 10 will hit an error related to the other issues
  // I will be working on next in execute_wait of Zoltan2_Directory_Comm.
  // When this count is low processes may set no_send_buff 1 in the
  // Zoltan2_Directory_Constructor and then MPI_Waitall fails in execute_wait.
  // I expect the other issues I have ongoing in that method are directly related
  // and all will be resolved together.
  const id_t totalIds = 1045;

  bool bQuickTest = true; // false is for performance testing

  int runLoops = bQuickTest ? 1 : getTotalLoops(totalIds);
  for(int loop = 0; loop < runLoops; ++loop) {
    comm->barrier();

    TestManager manager(comm, totalIds);

    manager.runTest<int, int, user_t>("contiguous int", Mode::Replace, 0, 1);

    // Do Replace tests
    manager.runTest<int, int, user_t>("non-contiguous int", Mode::Replace, 20, 3);
#ifndef CONVERT_DIRECTORY_ORIGINAL // long long not actually set up for original
    manager.runTest<long long, int, user_t>("long long", Mode::Replace, 200, 4);
#endif

    // Do Add tests
#ifndef CONVERT_DIRECTORY_ORIGINAL
    manager.runTest<int, int, user_t>("contiguous int", Mode::Add, 0, 1);
    manager.runTest<int, int, user_t>("non-contiguous int", Mode::Add, 20, 3);
    manager.runTest<long long, int, user_t>("long long", Mode::Add, 200, 4);
#endif

    // Do aggregate tests
#ifndef CONVERT_DIRECTORY_TPETRA
#ifndef CONVERT_DIRECTORY_ORIGINAL
#ifndef CONVERT_DIRECTORY_RELIC // TODO - add this once Kokkos is working
    manager.runTest<int, int, user_t>("contiguous int", Mode::Aggregate, 0, 1);
    manager.runTest<int, int, user_t>("non-contiguous int", Mode::Aggregate, 20, 3);
    manager.runTest<long long, int, user_t>("long long", Mode::Aggregate, 200, 4);
#endif
#endif
#endif

    // Check all processed
    int globalReturnCode = 0;
    int returnCode = manager.is_success() ? 0 : 1;
    Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM, returnCode,
    Teuchos::outArg(globalReturnCode));

    if(globalReturnCode != 0) {
      return globalReturnCode; // one of the ranks failed this test
    }
  }

  return 0;
}

#endif

#else // HAVE_MPI

int runDirectoryComparisonTest(int narg, char **arg) {
  throw std::logic_error( "No imeplementation for serial." );
}

#endif
