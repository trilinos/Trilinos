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
#include <Zoltan2_TPLTraits.hpp>

#ifdef CONVERT_DIRECTORY_ORIGINAL
#include <zoltan_dd_cpp.h>
#include "Zoltan2_Directory_Clock.hpp"
#else
#include "Zoltan2_Directory_Impl.hpp"
#endif

// The directory also has modes but currently working with an Original mode
// as well, which doesn't have it - so created this to have a universal
// reference to use. Eventually will delete all but Kokkos and then eliminate
// this declaration.
enum TestMode {
  Replace = 0,
  Add,
  Aggregate,
  TestMode_Max  // exists to provide a loop through these modes
};

// Do gid with multiple sub gids
// Note that the testing code (things in this file) needs knowledge
// about this structure so it's a hard coded assumption but the directory
// class itself does not - that is why there are sub_gid references in the
// below classes but no special handling for this in the directory.
// as long as sizeof works and it can be returned as a function
// type it should be ok to use as a gid_t.
#define GID_SET_LENGTH 3 // arbitrary for testing
struct gid_set_t {
  int sub_gid[GID_SET_LENGTH];
};

// same as gid above
#define LID_SET_LENGTH 3 // arbitrary for testing
struct lid_set_t {
  int sub_lid[LID_SET_LENGTH];
};

// a utility function to print messages out for each rank or optionally
// just for rank 0. Used to prevent mangled output.
void print_proc_safe(const std::string& message,
  Teuchos::RCP<const Teuchos::Comm<int> > comm, bool only_rank0 = true) {
  for(int proc = 0; proc < comm->getSize(); ++proc) {
    comm->barrier();
    if(proc == comm->getRank() && (!only_rank0 || proc == 0)) {
      std::cout << message << std::endl;
    }
  }
  comm->barrier();
}

// This is checked when remove is called and then find is called on a removed
// gid. The Original mode will update 0 to missing gid but the new code has an
// optional flag to indicate whether to throw (default true) but set false for
// this test. For the new code we use -1 as the 'unset' value to make it easier
// to see it worked but it's arbitrary. TODO: Discuss error behaviors for
// remove when not found.
#ifdef CONVERT_DIRECTORY_ORIGINAL
#define NOT_FOUND_VALUE 0
#else
#define NOT_FOUND_VALUE -1 // easier to tell it's the unset value as -1
#endif

// This class manages the IDs that we will update to the directory and which
// ids to remove and find. The class also handles error checking at the end
// to validate the find got the proper user data we expected.
//
// Derived classes handle vector user_t versus single user_t and
// also multiple gids (gid = [1,2,8,3]) or single gid (gid = 10).
// All this machinery is just for testing/evaluation and nothing really to
// do with the directory class itself.

// Because this testing code started to get complex due covering all possible
// behaviors of the directory, another test was created called
// directoryTest_KokkosSimple. That doesn't reference any of this stuff and is
// a more natural example of the directory being used (but doesn't thoroughly
// test the directory as is done here.
//
// The class structure is this:
//   IDs
//     Single_GID    - code for gid is simple type
//     Multiple_GID  - code for gid such as gid = [1,2,8,3]
//     Single_User   - code for when user is simple single type
//     Vector_User   - code for when user_t is a std::vector
//
// Then to run tests we use multiple inheritance of above classes.
//   Single_User_Single_GID   - inherits Single_GID + Single_User
//   Single_User_Multiple_GID - inherits Multiple_GID + Single_User
//   Vector_User_Single_GID   - inherits Single_GID + Vector_User
//   Vector_User_Multiple_GID - inherits Multiple_GID + Vector_User
template <typename gid_t,typename lid_t,typename user_t>
class IDs {
  public:
    // gcc requires this (not clang) - however it will not be used
    // I think clang is handling this properly based on below public virtual
    // inheritance but maybe will change this if the different compilers
    // will handle it in non-uniform ways
    IDs() {}

    /*! \brief Construct IDs
     *   \param[in] totalIds_ is total gid count across all procs
     *   \param[in] idBase_ is initial gid value
     *   \param[in] idStride_ is gid spacing
     *   \param[in] comm_ is the process communicator
     *   \param[in] mode_ is above enum (Replace, Add, or Aggregate)
     *   \param[in] test_name_ is helper string used for logging info
     *
     *  Note that in the current setup the derived class is responsible
     *  for calling setup() for proper initialization.
     */
    IDs(size_t totalIds_, size_t idBase_, size_t idStride_,
      Teuchos::RCP<const Teuchos::Comm<int> > &comm_, int mode_,
        const std::string& test_name_, bool print_detailed_output_,
        bool performance_test_, bool bUseLocalIDs_) :
        totalIds(totalIds_), idBase(idBase_), contiguous(idStride_ == 1),
        idStride(idStride_), comm(comm_), maxPrintSize(10), mode(mode_),
        did_setup(false), test_name(test_name_),
        print_detailed_output(print_detailed_output_),
        performance_test(performance_test_), bUseLocalIDs(bUseLocalIDs_) {
    }

    // prints all the subsets used as utlity to create interesting overlaps
    // prints all the update/remove/find decisions
    virtual void debug_print_subsets_and_decisions() {
        // map out the behaviors of the subsets
        comm->barrier();
        if(comm->getRank() == 0) {
          std::cout << "Sub1: " << std::endl;
          for(int proc = 0; proc < comm->getSize(); ++proc) {
            std::cout << " rank " << proc << ": ";
            for(size_t n = 0; n < totalIds; ++n) {
              std::cout << subset1(n, proc);
            }
            std::cout << std::endl;
          }

          std::cout << "Sub2: " << std::endl;
          for(int proc = 0; proc < comm->getSize(); ++proc) {
            std::cout << " rank " << proc << ": ";
            for(size_t n = 0; n < totalIds; ++n) {
              std::cout << subset2(n, proc);
            }
            std::cout << std::endl;
          }

          // update
          std::cout << "update: " << std::endl;
          for(int proc = 0; proc < comm->getSize(); ++proc) {
            std::cout << " rank " << proc << ": ";
            for(size_t n = 0; n < totalIds; ++n) {
              std::cout << proc_update_gid(convert_index_to_gid(n), proc);
            }
            std::cout << std::endl;
          }

          // remove
          std::cout << "remove: " << std::endl;
          for(int proc = 0; proc < comm->getSize(); ++proc) {
            std::cout << " rank " << proc << ": ";
            for(size_t n = 0; n < totalIds; ++n) {
              std::cout << proc_remove_gid(convert_index_to_gid(n), proc);
            }
            std::cout << std::endl;
          }

          // find
          std::cout << "find: " << std::endl;
          for(int proc = 0; proc < comm->getSize(); ++proc) {
            std::cout << " rank " << proc << ": ";
            for(size_t n = 0; n < totalIds; ++n) {
              std::cout << proc_find_gid(convert_index_to_gid(n), proc);
            }
            std::cout << std::endl;
          }
        }
        comm->barrier();
    }

    void printResultMessage(bool pass) const {
      std::string base_message = get_test_name() + " mode: " + get_mode_name();
      base_message = get_test_style() + " " + base_message;
      std::string message = pass ?
        "  Passed: " + base_message :
        "  FAILED: " + base_message + " on rank " +
        std::to_string(comm->getRank());
      print_proc_safe(message, comm, pass); // only print passes on rank 0
    }

    const std::string & get_test_name() const { return test_name; }

    /*! \brief evaluateTests - determine if test worked
     */
    virtual bool evaluateTests() const = 0;

    /*! \brief getMode - Replace, Add, or Aggregate
     */
    int getMode() const { return mode; }

    /*! \brief did_test_pass - did test pass
     */
    bool did_test_pass() const { return passed; }

    /*! \brief get_test_style
     *  the test is either vector user_t or single user_t
     *  the test is either multiple gids or single gid
     *  so 4 combinations are possible - this just describes it
     */
    virtual std::string get_test_style() const = 0;

    /*! \brief detailed notes on update IDs, find IDs, etc
     */
    void print() const
    {
      for(int proc = 0; proc < comm->getSize(); ++proc) {
        comm->barrier();
        if(proc == comm->getRank()) {
          if(proc == 0) {
            std::cout << std::endl <<
              "############ Output: " << get_test_name()
              << " ############" << std::endl;
          }
          std::cout << "Rank: " << proc << std::endl;
          print_gids(update_gids,  "Update gids"  );
          if(bUseLocalIDs) {
            print_lids(update_lids,  "Update lids"  );
          }
          print_gids(remove_gids, "Remove gids" );
          print_gids(find_gids,   "Find gids"   );
          if(bUseLocalIDs) {
            print_lids(find_lids,   "Find lids"  );
          }
          print_user_data();  // specific to single user or vector user type
          if(bUseLocalIDs) {
            print_lid_data(); // the lids we found, if bUseLocalIDs = true
          }

          if(proc == comm->getSize()-1) {
            std::cout << std::endl;
          }
        }
      }
      comm->barrier();
    }

  protected:
    std::string get_mode_name() const {
      switch(mode) {
        case TestMode::Add: return "Add"; break;
        case TestMode::Replace: return "Replace"; break;
        case TestMode::Aggregate: return "Aggregate"; break;
        default: throw std::logic_error("Bad mode."); break;
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

    virtual int convert_gid_to_index(const gid_t& gid) const = 0;

    virtual gid_t convert_index_to_gid(int index) const = 0;

    int proc_update_gid(gid_t gid, int rank) const {
      // this method should be guaranteed to return true for 1 of the procs
      // it can optionally return true for more than one to make a nice mix
      // the return is the number of repeats as we want the directory to work
      // in such a case
      int index = convert_gid_to_index(gid); // back to 0,1,2,... indexing

      // note not allowing duplicate gid for update - that will throw an error
      // however we do allow duplicate find gid where each shoudl get the same
      // value.
      if(trueForAtLeastOneProc(index, rank) || subset1(index, rank)) {
        return 1;
      }
      else {
        return 0; // do not update this ID
      }
    }

    int proc_find_gid(gid_t gid, int rank) const {
      // this method can be anything - things should work even if this is empty;
      int index = convert_gid_to_index(gid); // convert back to 0,1,2,...indexing

      // for a few gid's duplicate - I keep this a small percentage because
      // I want performance testing to be based on mostly unique ids.
      int test_duplicate_count = 2;
      // Note Tpetra requires unique IDs so proc_find_gid must be 1 for it to work
  #ifdef CONVERT_DIRECTORY_TPETRA
      test_duplicate_count = 1;
  #endif
      if(subset1(index, rank)) {
        return (index < 100) ? test_duplicate_count : 1;
      }
      else {
        return 0; // do not read this ID
      }
    }

    int proc_remove_gid(gid_t gid, int rank) const {
      // this should include some proc_find_gid values but not all, but also other
      // values for good testing - we want to remove more than we wrote locally
      // and check for proper removal
      int index = convert_gid_to_index(gid); // convert back to 0,1,2,...indexing

  #ifdef CONVERT_DIRECTORY_TPETRA
      return 0; // Tpetra not supporting remove right now
  #endif

      int test_duplicate_count = 2;
      if(subset2(index, rank)) {
        return (index < 100) ? test_duplicate_count : 1;
      }
      else {
        return 0; // do not remove this ID
      }
    }

    virtual user_t get_not_found_user() const = 0;

    virtual lid_t get_not_found_lid() const = 0;

    virtual void initialize_with_not_found_user() {
      find_user = std::vector<user_t>(find_gids.size(), get_not_found_user());
    }

    virtual void initialize_with_not_found_lid() {
      find_lids = std::vector<lid_t>(find_gids.size(), get_not_found_lid());
    }

    bool removedIDGlobally(gid_t gid) const {
      for(int proc = 0; proc < comm->getSize(); ++proc) {
        if(proc_remove_gid(gid, proc)) {
          return true;
        }
      }
      return false;
    }

    int sharedCount(gid_t gid) const {
      int count = 0;
      for(int proc = 0; proc < comm->getSize(); ++proc) {
        count += proc_update_gid(gid, proc);
      }
      return count;
    }

    void print_gids(const std::vector<gid_t> &printIds, std::string name) const {
      std::cout << "  " << printIds.size() << " " << name << ": ";
      for (size_t i = 0; i < printIds.size() && i < maxPrintSize; i++) {
        std::cout << gid_to_string(printIds[i]) << " ";
      }
      if(printIds.size() > maxPrintSize) {
        std::cout << "... ";
      }
      std::cout << std::endl;
    }

    void print_lids(const std::vector<lid_t> &printIds, std::string name) const {
      std::cout << "  " << printIds.size() << " " << name << ": ";
      for (size_t i = 0; i < printIds.size() && i < maxPrintSize; i++) {
        std::cout << lid_to_string(printIds[i]) << " ";
      }
      if(printIds.size() > maxPrintSize) {
        std::cout << "... ";
      }
      std::cout << std::endl;
    }

    virtual void print_user_data() const = 0;

    virtual void print_lid_data() const {
      std::cout << "Find LIDs ";
      for(size_t n = 0; n < this->find_gids.size(); ++n) {
        std::cout << this->gid_to_string(this->find_gids[n]) << ":" <<
          this->lid_to_string(this->find_lids[n]) << " ";
      }
      std::cout << std::endl;
    }

    virtual size_t gid_seed_value(const gid_t& gid) const = 0;

    virtual std::string gid_to_string(gid_t gid) const = 0;

    virtual std::string lid_to_string(lid_t lid) const = 0;

    virtual bool check_lid_equal(const lid_t & a, const lid_t & b) const = 0;

    virtual user_t getExpectedValue(gid_t gid) const = 0;

    virtual lid_t getExpectedlid(gid_t gid) const = 0;

    virtual user_t getInitialValue(gid_t gid, int rank) const = 0;

    virtual void execute() {
      setup();   // get things ready
      test();    // run the actual directory code
      analyze(); // analyze the results
      output();  // print results
      // debug_print_subsets_and_decisions();
    }

    virtual void setup() {
      // this clock is not part of the actual directory work, just test setup
      Zoltan2_Directory_Clock setupClock("setup");

      if(did_setup) {
        throw std::logic_error(
          "setup already called - execute once in highest level class only.");
      }
      did_setup = true;

      for(size_t n = 0; n < totalIds; ++n) {
        gid_t gid = convert_index_to_gid(n);
        bool bAtLeastOne = false;
        for(int proc = 0; proc < comm->getSize(); ++proc) {
          if(proc_update_gid(gid, proc)) {
            bAtLeastOne = true;
          }
        }
        if(!bAtLeastOne) {
          throw std::logic_error("The proc_update_gid method must generate an"
            " algorithm which returns true for at least 1 proce for all IDs. If"
            " this is changed then the findID algorithm should not return true"
            " for any of those missing IDs. This is a sanity check.");
        }
      }

      // create ids list - the count here is arbitrary and sets the complexity
      // of the test
      for(size_t n = 0; n < totalIds; ++n) {
        gid_t gid = convert_index_to_gid(n);

        int count_proc_update_gid = proc_update_gid(gid, comm->getRank());
        for(int i = 0; i < count_proc_update_gid; ++i) {
          update_gids.push_back(gid);
          update_lids.push_back(getExpectedlid(gid));
        }

        for(int i = 0; i < proc_find_gid(gid, comm->getRank()); ++i) {
          find_gids.push_back(gid);
        }

        for(int i = 0; i < proc_remove_gid(gid, comm->getRank()); ++i) {
          remove_gids.push_back(gid);
        }
      }

      // set up all the initial writes
      this->update_user.resize(update_gids.size());
      for(size_t n = 0; n < update_user.size(); ++n) {
        update_user[n] = getInitialValue(update_gids[n], comm->getRank());
      }
    }

    virtual void test() = 0;

#ifndef CONVERT_DIRECTORY_ORIGINAL
    template<typename directory_t>
#endif
    void test_implement() {
      Zoltan2_Directory_Clock allClock("all");

      const int debug_level = 0;

#ifdef CONVERT_DIRECTORY_ORIGINAL
      if (update_gids.size() > std::numeric_limits<int>::max())
        throw std::runtime_error("Problem too large for zoltan");
      int gid_length = Zoltan2::TPL_Traits<ZOLTAN_ID_PTR, gid_t>::NUM_ID;
      int lid_length = Zoltan2::TPL_Traits<ZOLTAN_ID_PTR, lid_t>::NUM_ID;
      MPI_Comm mpicomm = Teuchos::getRawMpiComm(*comm); // TODO comm fixes
      Zoltan2_Directory_Clock constructClock("construct");
      Zoltan_DD zz(mpicomm, gid_length, lid_length, sizeof(user_t),
        update_gids.size(), debug_level);
      constructClock.complete(); // to match with new directory clocks...
#else

    #ifdef CONVERT_DIRECTORY_RELIC
      directory_t zz(comm, bUseLocalIDs, update_gids.size(),
        debug_level, contiguous);
    #else
      directory_t zz(comm, bUseLocalIDs, debug_level, contiguous);
    #endif

#endif // CONVERT_DIRECTORY_ORIGINAL

      // usually this step is not necessary - for this test we initialize all
      // the find_user values to a specific number and then check that it's not
      // changed when we call find on an ID which was previously removed. The
      // formal behavior for this in a normal setup is probably to just have an
      // error but in this special case we override the throw when we call find.
      initialize_with_not_found_user();

      // for lid note that the original mode will leave a not found gid such
      // that the lid is set to the gid. This is how the find_local worked. In
      // the new directory the value will be left untouched so this will be
      // preserved. Also original mode requires this to be set but for new
      // modes with std::vector we can make this fillled automatically.
      if(bUseLocalIDs) { // if false it shouldn't matter if we call this
        initialize_with_not_found_lid();
      }

#ifdef CONVERT_DIRECTORY_ORIGINAL
      // Original zoltan code only supports the replace mode
      if(mode != Replace) {
        throw std::logic_error("Bad setup. Original mode can only run replace.");
      }
      Zoltan2_Directory_Clock updateClock("update");
      zz.Update(reinterpret_cast<ZOLTAN_ID_PTR>(&(update_gids[0])),
        bUseLocalIDs ? reinterpret_cast<ZOLTAN_ID_PTR>(&(update_lids[0])) : NULL,
        update_gids.size() ? reinterpret_cast<char*>(&(update_user[0])) : NULL,
        NULL, update_gids.size());
      updateClock.complete(); // to match with new directory clocks...

      Zoltan2_Directory_Clock removeClock("remove");
      zz.Remove(reinterpret_cast<ZOLTAN_ID_PTR>(&(remove_gids[0])),
        remove_gids.size());
      removeClock.complete(); // to match with new directory clocks...

      Zoltan2_Directory_Clock findClock("find");
      zz.Find(reinterpret_cast<ZOLTAN_ID_PTR>(&(find_gids[0])),
        bUseLocalIDs ? reinterpret_cast<ZOLTAN_ID_PTR>(&(find_lids[0])) : NULL,
        find_gids.size() ? reinterpret_cast<char*>(&(find_user[0])) : NULL,
        NULL, find_gids.size(), NULL);
     findClock.complete(); // to match with new directory clocks...
#else
      // convert generic mode to Zoltan2Directory - this awkward step exists
      // because the original mode doesn't have a directory - will improve.
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

      // perhaps change this back to a ptr so NULL can be considered off
      std::vector<int> ignore_int; // TODO: decide how to best handle this API
      std::vector<lid_t> ignore_lids; // TODO: how to best handle this API

      zz.update(update_gids,
        bUseLocalIDs ? update_lids : ignore_lids,
        update_user, ignore_int, directoryMode);

      #ifndef CONVERT_DIRECTORY_TPETRA // tpetra won't set remove_gids currently
      zz.remove(remove_gids);
      #endif

      zz.find(find_gids,
        bUseLocalIDs ? find_lids : ignore_lids,
        find_user, ignore_int, ignore_int, false);
#endif
      comm->barrier();
    }

    void analyze() {
      passed = evaluateTests();
    }

    void output() {
      if(print_detailed_output) {
        print();
      }
      if(!performance_test) {
        printResultMessage(passed);
      }
    }

    size_t totalIds;
    size_t idBase;
    bool contiguous;              // Flag indicating whether IDs are contiguous
    size_t idStride;
    Teuchos::RCP<const Teuchos::Comm<int> > comm;
    size_t maxPrintSize;             // print only this number for debugging
    int mode;
    std::vector<gid_t> update_gids;  // gids generated on this proc
    std::vector<lid_t> update_lids;  // corresponding lids generated on this proc
    std::vector<user_t> update_user; // user data we initially updated with
    std::vector<gid_t> find_gids;    // gids to read
    std::vector<lid_t> find_lids;    // lids will be filled
    std::vector<user_t> find_user;   // user data we find
    std::vector<gid_t> remove_gids;  // gids to remove
    bool did_setup;                  // error check setup not called twice
    std::string test_name;           // used for logging
    bool passed;                     // did we pass
    bool print_detailed_output;      // log all the details
    bool performance_test;           // is this being run as a performance test
    bool bUseLocalIDs;               // should local ids be tested/applied
};

template <typename gid_t,typename lid_t, typename user_t>
class Single_GID : public virtual IDs<gid_t,lid_t,user_t>
{
  public:
    Single_GID() {}

  protected:
    virtual lid_t get_not_found_lid() const {
      return NOT_FOUND_VALUE;
    }

    virtual std::string gid_to_string(gid_t gid) const {
      return "(" + std::to_string(gid) + ")";
    };

    virtual std::string lid_to_string(lid_t lid) const {
      return "(" + std::to_string(lid) + ")";
    };

    virtual bool check_lid_equal(const lid_t & a, const lid_t & b) const {
      return (a == b);
    }

    virtual size_t gid_seed_value(const gid_t& gid) const {
      return gid;
    }

    virtual int convert_gid_to_index(const gid_t& gid) const {
      return (gid-this->idBase)/this->idStride;
    }

    virtual gid_t convert_index_to_gid(int index) const {
      return this->idBase + index * this->idStride;
    }

    virtual lid_t getExpectedlid(gid_t gid) const {
      return gid + 1; // this is for testing - make the lid equal to gid+1
    }
};

template <typename gid_t,typename lid_t, typename user_t>
class Multiple_GID : public virtual IDs<gid_t,lid_t,user_t>
{
  public:
    Multiple_GID(size_t gid_length_, size_t lid_length_) :
      gid_length(gid_length_), lid_length(lid_length_) {
      // currently we're making the testing for multiple gid cover
      // multiple lid as well - further combinations would be a gid multiple
      // type such as [1,4,5,8] and a simple lid type ( such as int ).
    }

  protected:
    virtual lid_t get_not_found_lid() const {
      lid_t not_found;
      for(size_t n = 0; n < lid_length; ++n) {
        not_found.sub_lid[n] = NOT_FOUND_VALUE;
      }
      return not_found;
    }

    virtual std::string gid_to_string(gid_t gid) const {
      std::string output_string = "(";
      for(size_t n = 0; n < gid_length; ++n) {
        if(n!=0) output_string += " ";
        output_string += std::to_string(gid.sub_gid[n]);
      }
      output_string += ")";
      return output_string;
    };

    virtual std::string lid_to_string(lid_t lid) const {
      std::string output_string = "(";
      for(size_t n = 0; n < lid_length; ++n) {
        if(n!=0) output_string += " ";
        output_string += std::to_string(lid.sub_lid[n]);
      }
      output_string += ")";
      return output_string;
    };

    virtual bool check_lid_equal(const lid_t & a, const lid_t & b) const {
      for(size_t n = 0; n < lid_length; ++n) {
        if(a.sub_lid[n] != b.sub_lid[n]) {
          return false;
        }
      }
      return true;
    }

    virtual size_t gid_seed_value(const gid_t& gid) const {
      return gid.sub_gid[0]; // just uses first to generate values
    }

    virtual int convert_gid_to_index(const gid_t& gid) const {
      // here we are testing gid with multiple elements
      // the first element is set using the same rules as for single gid
      // so just read that and convert back to index
      return (gid.sub_gid[0]-this->idBase)/this->idStride;
    }

    virtual gid_t convert_index_to_gid(int index) const {
      gid_t gid;
      // here we are testing gid with multiple elements
      // set the first element as if it was a single gid - same as other tests
      // set all subsequent gid sub ids in increasing order just to have change
      // the values don't have any impact on the behavior of the directory
      int val = this->idBase + index * this->idStride;
      for(size_t n = 0; n < gid_length; ++n) {
        gid.sub_gid[n] = val + n;
      }
      return gid;
    }

    virtual lid_t getExpectedlid(gid_t gid) const {
      lid_t result;
      for(size_t n = 0; n < lid_length; ++n) {
        result.sub_lid[n] = n + gid.sub_gid[0] + 1; // lid will be gid[0]+1, gid[1]+2, gid[2]+3...
      }
      return result;
    }

  private:
    size_t gid_length;
    size_t lid_length;
};

template <typename gid_t,typename lid_t, typename user_t>
class Single_User : public virtual IDs<gid_t,lid_t,user_t>
{
  public:
    Single_User() {}

  protected:
    virtual void test() {
#ifdef CONVERT_DIRECTORY_ORIGINAL
      this->test_implement(); // no templated type
#else
      typedef typename Zoltan2::Zoltan2_Directory_Simple<gid_t,lid_t,user_t>
        directory_t;
      this->template test_implement<directory_t>();
#endif
    }

    virtual user_t get_not_found_user() const {
      return NOT_FOUND_VALUE;
    }

    virtual user_t getExpectedValue(gid_t gid) const {
      switch(this->mode) {
        case TestMode::Replace:
          return this->getInitialValue(gid, this->comm->getRank());
          break;
        case TestMode::Add:
          return this->sharedCount(gid); // should have summed
          break;
        case TestMode::Aggregate:
          throw std::logic_error("Unexpected aggregate mode for non array.");
          break;
        default:
          throw std::logic_error("Unexpected mode index.");
          break;
      }
    }

    virtual user_t getInitialValue(gid_t gid, int rank) const {
      switch(this->mode) {
        case TestMode::Replace:
          // arbitrary - just make it something good to test
          // note that replace is a tricky case because two procs could both
          // try to update the same gid with different values. Then the result
          // could be arbitrary based on how the messages are passed.
          // For the testing code here the value is fixed to be the same for
          // all procs but the question to resolve is should the directory
          // handle conflicting replace calls with an error or simply assume the
          // user is responsible for making it logically consistent.
          return this->gid_seed_value(gid) + 3;
          break;
        case TestMode::Add:
          return 1; // we will be testing if they sum to shared count
          break;
        case TestMode::Aggregate:
          throw std::logic_error("Aggregate requested for non array setup.");
          break;
        default:
          throw std::logic_error("Unexpected mode index.");
          break;
      }
    }

    virtual bool evaluateTests() const {

      // check the results
      bool pass = true;
      for(int proc = 0; proc < this->comm->getSize(); ++proc) {
        bool passRank = true;
        this->comm->barrier();
        if(proc == this->comm->getRank()) {
          for(size_t i = 0; i < this->find_gids.size(); ++i) {
            gid_t gid = this->find_gids[i];

            // verify removal - eventually I think we may have the directory
            // just throw an error if we try to find on a gid which was removed.
            // Or some return value will indicate. For now we pass in all values
            // of NOT_FOUND_VALUE which will still be set to that if not found.
            if(this->removedIDGlobally(gid)) {
              if(this->find_user[i] != NOT_FOUND_VALUE) {
                passRank = false;
                std::cout << "Removed gid: " << this->gid_to_string(gid) <<
                  " but got value: " << this->find_user[i] <<
                  " when we expected to read the unset value of " <<
                  NOT_FOUND_VALUE << ". " << " This is incorrect. " <<
                  "Remove FAILED." << std::endl;
              }
            }
            else {
              user_t expected_value = this->getExpectedValue(gid);
              if(this->find_user[i] != expected_value) {
                passRank = false;
                std::cout << "Failed read user data for global ID: " <<
                  this->gid_to_string(gid) << ". Expected data: " <<
                  expected_value << " Got data: " <<
                  this->find_user[i] << std::endl;
              }

              if(this->bUseLocalIDs) {
                const lid_t & find_lid = this->find_lids[i];
                lid_t expected_lid = this->getExpectedlid(gid);
                if(!this->check_lid_equal(find_lid, expected_lid)) {
                  passRank = false;
                  std::cout << "Failed read lid for global ID: " <<
                    this->gid_to_string(gid) << ". Expected lid: " <<
                    this->lid_to_string(expected_lid) << " Got lid: " <<
                    this->lid_to_string(find_lid) << std::endl;
                }
              }
            }

            if(!passRank) {
              break; // we don't need to check further
            }
          }

          if(!passRank) {
            std::cout << "Checked rank: " << this->comm->getRank() <<
              " with nIds: " << this->find_gids.size() << " which " <<
              "FAILED" << std::endl;
          }

          if(!passRank) {
            pass = false;
          }
        }
      }
      this->comm->barrier();

      return pass;
    }

    virtual void print_user_data() const {
      std::cout << "Write GID user data ";
      for(size_t n = 0; n < this->update_gids.size(); ++n) {
        std::cout << this->gid_to_string(this->update_gids[n]) << ":" <<
          this->update_user[n] << " ";
      }
      std::cout << std::endl;

      std::cout << "Find GID user data ";
      for(size_t n = 0; n < this->find_gids.size(); ++n) {
        std::cout << this->gid_to_string(this->find_gids[n]) << ":" <<
          this->find_user[n] << " ";
      }
      std::cout << std::endl;
    }
};

template <typename gid_t,typename lid_t, typename user_t>
class Vector_User : public virtual IDs<gid_t,lid_t,user_t>
{
  public:
    Vector_User() {}

  protected:
    virtual void test() {
#ifdef CONVERT_DIRECTORY_ORIGINAL
      this->test_implement(); // no templated type
#else
      typedef typename Zoltan2::Zoltan2_Directory_Vector<gid_t,lid_t,user_t>
        directory_t;
      this->template test_implement<directory_t>();
#endif
    }

    virtual user_t get_not_found_user() const {
      return user_t(1, NOT_FOUND_VALUE);
    }

    virtual user_t getExpectedValue(gid_t gid) const {
      switch(this->mode) {
        case TestMode::Add:
          return user_t(test_create_array_length(gid, this->comm->getRank()),
            this->sharedCount(gid)); // should have summed
          break;
        case TestMode::Replace:
          return getInitialValue(gid, this->comm->getRank());
          break;
        case TestMode::Aggregate:
          {
            // for this we need the union of all the procs
            user_t final_aggregated;
            // loop through all possible updaters
            for(int proc = 0; proc < this->comm->getSize(); ++proc) {
              if(this->proc_update_gid(gid, proc)) { // did this proc update?
                user_t proc_input = getInitialValue(gid, proc); // get original
                for(size_t i = 0; i < proc_input.size(); ++i) { // scan elements
                  auto val = proc_input[i]; // get the array element
                  if(final_aggregated.size() == 0 ||
                    val > final_aggregated[final_aggregated.size()-1]) {
                    // add first element or tail end
                    final_aggregated.push_back(val);
                  }
                  else { // loop and insert to keep ordering until match found
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
            return final_aggregated;
          }
          break;
        default:
          throw std::logic_error("Unexpected mode index.");
          break;
      }
    }

    virtual user_t getInitialValue(gid_t gid, int rank) const {
      // determine length of array
      size_t modLength = test_create_array_length(gid, rank);
      user_t array(modLength);
      for(size_t n = 0; n < array.size(); ++n) {
        switch(this->mode) {
          case TestMode::Replace:
            array[n] = this->gid_seed_value(gid) + 3; // 3 is arbitrary
            break;
          case TestMode::Add:
            array[n] = 1; // all elements sum to comm->getSize()
            break;
          case TestMode::Aggregate:
            // Now we want some mix so that, for example, gid 10 will have
            // different but overlapping values for each array element n
            // For example proc 1 could update with gid=10 array={5,7,9}
            //       while proc 2 could update with gid=10 array={3,5,7}
            // Then we expect the result to be {3,5,7,9}
            array[n] = n + rank*2; // this creates overlapping regions
            break;
          default:
            throw std::logic_error("Unexpected mode index.");
            break;
        }
      }
      return array;
    }

    virtual bool evaluateTests() const {

      // check the results
      bool pass = true;
      for(int proc = 0; proc < this->comm->getSize(); ++proc) {
        bool passRank = true;
        this->comm->barrier();
        if(proc == this->comm->getRank()) {
          for(size_t i = 0; i < this->find_gids.size(); ++i) {
            gid_t gid = this->find_gids[i];

            // verify removal - eventually I think we may have the directory
            // just throw an error if we try to find on a gid which was removed.
            // Or some return value will indicate. For now we pass in all values
            // of NOT_FOUND_VALUE which will still be set to that if not found.
            if(this->removedIDGlobally(gid)) {
              // should be an array of size 1 with element NOT_FOUND_VALUE
              if(this->find_user[i].size() != 1 ||
                this->find_user[i][0] != NOT_FOUND_VALUE) {
                passRank = false;
                std::cout << "Removed array for gid: " <<
                  this->gid_to_string(gid) <<
                  " but something set the user data which is incorrect. "
                  "Remove FAILED." << std::endl;
              }
            }
            else {
              user_t expected_value = getExpectedValue(gid);
              // first validate the array length is correct
              if(this->find_user[i].size() != expected_value.size()) {
                std::cout << "  Rank: " << proc << " array size is incorrect for"
                " gid: " << this->gid_to_string(gid) << ". Expected size: " <<
                expected_value.size() << " and got size: " <<
                this->find_user[i].size() << std::endl;
                passRank = false;
                break;
              }

              // TODO: Fix this code duplicated in the other evaluateTests
              // Generall these two methdos can perhaps be merged better now
              if(this->bUseLocalIDs) {
                const lid_t & find_lid = this->find_lids[i];
                lid_t expected_lid = this->getExpectedlid(gid);
                if(!this->check_lid_equal(find_lid, expected_lid)) {
                  passRank = false;
                  std::cout << "Failed read lid for global ID: " <<
                    this->gid_to_string(gid) << ". Expected lid: " <<
                    this->lid_to_string(expected_lid) << " Got lid: " <<
                    this->lid_to_string(find_lid) << std::endl;
                }
              }

              // now loop the elements and validate each individual element
              for(size_t arrayIndex = 0; arrayIndex < this->find_user[i].size()
                && passRank; ++arrayIndex) {
                if(this->find_user[i][arrayIndex] != expected_value[arrayIndex]) {
                  passRank = false;
                  std::cout << "  Failed vector read for global ID: " <<
                    this->gid_to_string(gid)
                    << ". Expected: " << expected_value[arrayIndex] <<
                    " at array index " << arrayIndex << ". Got: " <<
                    this->find_user[i][arrayIndex] << std::endl;
                }
              }
            }

            if(!passRank) {
              break; // we don't need to check further
            }
          }

          if(!passRank) {
            std::cout << "  Checked rank: " << this->comm->getRank() <<
              " with num find gids: " << this->find_gids.size() << " which " <<
              "FAILED" << std::endl;
          }

          if(!passRank) {
            pass = false;
          }
        }
      }
      this->comm->barrier();

      return pass;
    }

    virtual void print_user_data() const {
      for(size_t n = 0; n < this->update_gids.size(); ++n) {
        std::cout << "  Write array for GID " <<
          this->gid_to_string(this->update_gids[n]) << ": ";
        for(size_t a = 0; a < this->update_user[n].size(); ++a) {
          std::cout << this->update_user[n][a] << " ";
        }
        std::cout << std::endl;
      }
      for(size_t n = 0; n < this->find_gids.size(); ++n) {
        std::cout << "  Read array for GID " <<
          this->gid_to_string(this->find_gids[n]) << ": ";
        for(size_t a = 0; a < this->find_user[n].size(); ++a) {
          std::cout << this->find_user[n][a] << " ";
        }
        std::cout << std::endl;
      }
    }

    virtual size_t test_create_array_length(gid_t gid, int proc) const {
      switch(this->mode) {
        case Replace:
          // replace is in some ways the trickiest because Add and Aggregate
          // are not order dependent. If two different procs both try to update
          // the same gid with replace the final result may be arbitrary ordering
          // and not well defined. Should the directory be responsible for
          // detecting a logic error from the user? Note this issue is not a
          // vector issue, but a general issue with replace. For now we assume
          // the user is consistent and never sends conflicting calls. Replace
          // uses the same value for type or vector and is independent of rank.
          return (this->gid_seed_value(gid)%7); // 1..7 inclusive same as Add
          break;
        case Add:
          // in this case all vector lengths must be identical for each proc
          // or the directory gives an error - cannot add length 2 and length 8
          return (this->gid_seed_value(gid)%7); // 1..7 inclusive
          break;
        case Aggregate:
          // in this case we want different proc but same gid to make different
          // vectors - so for example mix a length 3 with a length 8
          return (this->gid_seed_value(gid)%7)+proc; // varies per gid and proc
          break;
        default:
          throw std::logic_error("test_create_array_length bad mode.");
      }
    }
};

template <typename gid_t,typename lid_t, typename user_t>
class Single_User_Single_GID :
  public Single_User<gid_t,lid_t,user_t>, public Single_GID<gid_t,lid_t,user_t>
{
  public:
    Single_User_Single_GID(size_t totalIds_, size_t idBase_, size_t idStride_,
      Teuchos::RCP<const Teuchos::Comm<int> > &comm_, int mode_,
      const std::string& name_, bool print_detailed_output_,
      bool performance_test_, bool bUseLocalIDs_) :
        IDs<gid_t, lid_t, user_t>(totalIds_, idBase_, idStride_, comm_, mode_,
          name_, print_detailed_output_, performance_test_, bUseLocalIDs_),
        Single_User<gid_t, lid_t, user_t>(),
        Single_GID<gid_t, lid_t, user_t>() {
      this->execute();
    }

    virtual std::string get_test_style() const {
      return "Single_User_Single_GID";
    }
};

template <typename gid_t,typename lid_t, typename user_t>
class Single_User_Multiple_GID :
  public Single_User<gid_t,lid_t,user_t>, public Multiple_GID<gid_t,lid_t,user_t>
{
  public:
    Single_User_Multiple_GID(size_t gid_length_, size_t lid_length_,
      size_t totalIds_, size_t idBase_, size_t idStride_,
      Teuchos::RCP<const Teuchos::Comm<int> > &comm_, int mode_,
      const std::string& name_, bool print_detailed_output_,
      bool performance_test_, bool bUseLocalIDs_) :
        IDs<gid_t, lid_t, user_t>(totalIds_, idBase_, idStride_, comm_, mode_,
          name_, print_detailed_output_, performance_test_, bUseLocalIDs_),
        Single_User<gid_t, lid_t, user_t>(),
        Multiple_GID<gid_t, lid_t, user_t>(gid_length_, lid_length_) {
      this->execute();
    }

    virtual std::string get_test_style() const {
      return "Single_User_Multiple_GID";
    }
};

template <typename gid_t,typename lid_t, typename user_t>
class Vector_User_Single_GID :
  public Vector_User<gid_t,lid_t,user_t>, public Single_GID<gid_t,lid_t,user_t>
{
  public:
    Vector_User_Single_GID(size_t totalIds_, size_t idBase_, size_t idStride_,
      Teuchos::RCP<const Teuchos::Comm<int> > &comm_, int mode_,
      const std::string& name_, bool print_detailed_output_,
      bool performance_test_, bool bUseLocalIDs_) :
        IDs<gid_t, lid_t, user_t>(totalIds_, idBase_, idStride_, comm_, mode_,
          name_, print_detailed_output_, performance_test_, bUseLocalIDs_),
        Vector_User<gid_t, lid_t, user_t>(),
        Single_GID<gid_t, lid_t, user_t>() {
      this->execute();
    }

    virtual std::string get_test_style() const {
      return "Vector_User_Single_GID";
    }
};

template <typename gid_t,typename lid_t, typename user_t>
class Vector_User_Multiple_GID :
  public Vector_User<gid_t,lid_t,user_t>, public Multiple_GID<gid_t,lid_t,user_t>
{
  public:
    Vector_User_Multiple_GID(size_t gid_length_, size_t lid_length_,
      size_t totalIds_, size_t idBase_, size_t idStride_,
      Teuchos::RCP<const Teuchos::Comm<int> > &comm_, int mode_,
      const std::string& name_, bool print_detailed_output_,
      bool performance_test_, bool bUseLocalIDs_) :
        IDs<gid_t, lid_t, user_t>(totalIds_, idBase_, idStride_, comm_, mode_,
          name_, print_detailed_output_, performance_test_, bUseLocalIDs_),
        Vector_User<gid_t, lid_t, user_t>(),
        Multiple_GID<gid_t, lid_t, user_t>(gid_length_, lid_length_) {
      this->execute();
    }

    virtual std::string get_test_style() const {
      return "Vector_User_Multiple_GID";
    }
};

// The TestManager holds some common values and provides 4 API calls for running
// tests using the different options. The Multiple GID option is just for
// setting up tests as the directory class itself works automatically to handle
// gid_t of int or a struct. For Vector User this is currently handled by
// a different named directory class but eventually may merge these back into
// a single class and just let the API distinguish. However I'm not sure about
// how to make that work cleanly with the templating. The 4 testing modes are:
//   Single User + Single GID      (ex. user_t int, gid=8)
//   Single User + Multiple GID    (ex. user_t int, gid=[0,1,8])
//   Vector User + Single GID      (ex. user_t std::vector<int>, gid=8)
//   Vector User + Multiple GID    (ex. user_t std::vector<int>, gid=[0,1,8]
// The point of this is to make sure this unit test actually tries all these
// combinations and fails if one is missing. Only Kokkos mode supports everything.
class TestManager {
  public:
    TestManager(Teuchos::RCP<const Teuchos::Comm<int> > comm_, int totalIds_,
      bool print_detailed_output_, bool performance_test_, bool bUseLocalIDs_) :
        comm(comm_), totalIds(totalIds_), all_pass(true),
        print_detailed_output(print_detailed_output_),
        performance_test(performance_test_), bUseLocalIDs(bUseLocalIDs_) {
    }

    // Single User + Single GID
    template<typename gid_t, typename lid_t, typename user_t>
    void run_single_user_single_gid(const std::string& name_,
      int mode_, size_t idBase_, int idStride_) {
      Single_User_Single_GID<gid_t,lid_t,user_t> ids(totalIds,
        idBase_, idStride_, comm, mode_, name_, print_detailed_output,
        performance_test, bUseLocalIDs);
      if(!ids.did_test_pass()) { all_pass = false; }
    }

    // Single User + Multiple GID
    template<typename gid_t, typename lid_t, typename user_t>
    void run_single_user_multiple_gid(size_t gid_length_, size_t lid_length_,
      const std::string& name_, int mode_, size_t idBase_, int idStride_) {
      Single_User_Multiple_GID<gid_t,lid_t,user_t> ids(gid_length_, lid_length_,
        totalIds, idBase_, idStride_, comm, mode_, name_, print_detailed_output,
        performance_test, bUseLocalIDs);
      if(!ids.did_test_pass()) { all_pass = false; }
    }

    // Vector User + Single GID
    template<typename gid_t, typename lid_t, typename user_t>
    void run_vector_user_single_gid(const std::string& name_,
      int mode_, size_t idBase_, int idStride_) {
      Vector_User_Single_GID<gid_t,lid_t,user_t> ids(totalIds,
        idBase_, idStride_, comm, mode_, name_, print_detailed_output,
        performance_test, bUseLocalIDs);
      if(!ids.did_test_pass()) { all_pass = false; }
    }

    // Vector User + Multiple GID
    template<typename gid_t, typename lid_t, typename user_t >
    void run_vector_user_multiple_gid(size_t gid_length_, size_t lid_length_,
      const std::string& name_, int mode_, size_t idBase_, int idStride_) {
      Vector_User_Multiple_GID<gid_t,lid_t,user_t> ids(gid_length_, lid_length_,
        totalIds, idBase_, idStride_, comm, mode_, name_, print_detailed_output,
        performance_test, bUseLocalIDs);
      if(!ids.did_test_pass()) { all_pass = false; }
    }

    bool did_all_pass() const { return all_pass; }

  private:
    Teuchos::RCP<const Teuchos::Comm<int> > comm;
    int totalIds;
    bool all_pass;
    bool print_detailed_output;
    bool performance_test;
    bool bUseLocalIDs;
};


int runPerformanceTest(Teuchos::RCP<const Teuchos::Comm<int> > comm,
  bool bVectorMode, TestMode mode) {

  int err = 0;

  // just run 10000 and 100000000 for now. Tpetra has a lot of overhead so
  // looping that many times is too slow. This gives a good 1st performanct test.
  std::vector<size_t> test_totalIds = { 10000000 };
  for(auto itr_totalIds = test_totalIds.begin();
    itr_totalIds != test_totalIds.end(); ++itr_totalIds) {
    const size_t & totalIds = (*itr_totalIds);

    // can loop these more for more accurate times - just use 1 for now
    size_t total_loops = 1;

    const bool print_output = false; // noisy output with values for each gid
    const bool performance_test = true;
    const bool bUseLocalIDs = false;
    TestManager manager(comm, totalIds, print_output,
      performance_test, bUseLocalIDs);

    for(size_t n = 0; n < total_loops; ++n) {
      if(bVectorMode) {
#ifdef CONVERT_DIRECTORY_KOKKOS
        manager.run_vector_user_single_gid<int, int, std::vector<int>>
          ("contiguous int", mode, 0, 1);
#else
        throw std::logic_error("This non-kokkos mode doesn't support vectors.");
#endif
      }
      else {
        manager.run_single_user_single_gid<int, int, int>
          ("contiguous int", mode, 0, 1);
      }
      if(!manager.did_all_pass()) {
        err = 1; // if we fail at some point just drop out
      }
    }

    Zoltan2_Directory_Clock::logAndClearTimes(
      std::string("Times for Total IDs: ") + std::to_string(totalIds), comm);
  }

  // Get the global err results so we fail properly if any proc failed
  int errGlobal;
  Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM, err,
    Teuchos::outArg(errGlobal));
  return errGlobal; // only 0 if all tests and all proc return 0
}

// The performance test is used to compare times of different sub sections in
// the zoltan2. This is probably temporary and the 'true' unit test is
// runDirectoryTests below, which grinds through all combinations and validates
// results are returned correctly.
int runPerformanceDirectoryTests(int narg, char **arg) {
  Teuchos::GlobalMPISession mpiSession(&narg,&arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int err = 0;

  // first run a simple replace - all modes can do this
  err += runPerformanceTest(comm, false, Replace);

#ifdef CONVERT_DIRECTORY_KOKKOS
  // then we can compare with a vector Aggregate
  err += runPerformanceTest(comm, true, Aggregate);
#endif

  // Get the global err results so we fail properly if any proc failed
  int errGlobal;
  Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM, err,
    Teuchos::outArg(errGlobal));
  return errGlobal; // only 0 if all tests and all proc return 0
}

// This is a systematic grind through all the possible options the directory
// can do and makes use of above classes to hit on various features like
// vector user type versus single user type, multiple gid, different types,
// etc. This is intended to be a thorough test of the directory but is
// not transparent as a starting point so another test was created called
// directoryTest_KokkosSimple which demonstrates a more natural usage of the
// directory - that provides examples of how to use the directory which is
// independent of all the above stuff.
//
// Setting print_output true below will print a lot of information about the
// gid lists for update, remove, find, as well as the user data and associated
// lids after running find.
int runDirectoryTests(int narg, char **arg) {

  Teuchos::GlobalMPISession mpiSession(&narg,&arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  // run the tests through a range of totalIds
  // we pick 0 and 1 and some low values for edge cases, then just
  // try to hit on some random spots
  std::vector<size_t> run_with_totalIds = {0, 1, 2, 3, 5, 27, 63, 456, 1093};
  int err = 0;

  const bool print_output = false; // noisy output with values for each gid
  const bool performance_test = false;

  for(int run_with_local_ids = 0; run_with_local_ids <= 1; ++run_with_local_ids) {
    #ifdef CONVERT_DIRECTORY_TPETRA
      if(run_with_local_ids == 1) {
        continue; // all modes support local ids except tpetra
      }
    #endif

    for(size_t n = 0; n < run_with_totalIds.size(); ++n) {

      print_proc_safe("Testing totalIds: " + std::to_string(run_with_totalIds[n]),
        comm, true);

      comm->barrier();

      TestManager manager(comm, run_with_totalIds[n], print_output,
        performance_test, run_with_local_ids ? true : false);

      // loop the modes: Replace, Add, Aggregate and then do various tests on
      // each which vary gid_t, single or vector user type, single or multipe gid
      // some of the non-kokkos modes don't support everything and skip tests but
      // eventually they will all be deleted leaving only Kokkos
      for(int test_mode = 0; test_mode < TestMode_Max; ++test_mode) {

  #ifndef CONVERT_DIRECTORY_KOKKOS
        if(test_mode == Aggregate) {
          continue; // only Kokkos supports Aggregate right now
        }
  #endif
  #ifdef CONVERT_DIRECTORY_ORIGINAL
        if(test_mode == Add) {
          continue; // Original mode does not support Add
        }
  #endif

        // Aggregate mode is for vector user type only
        if(test_mode != Aggregate) {
          manager.run_single_user_single_gid<int, int, int>
            ("contiguous int", test_mode, 0, 1);
          manager.run_single_user_single_gid<int, int, int>
            ("non-contiguous int", test_mode, 20, 3);

  #ifndef CONVERT_DIRECTORY_ORIGINAL // long long not set up for original mode
          manager.run_single_user_single_gid<long long, int, int>
            ("long long", test_mode, 200, 4);
  #endif

        }

  #ifdef CONVERT_DIRECTORY_KOKKOS
        manager.run_vector_user_single_gid<int, int, std::vector<int>>
          ("contiguous int", TestMode::Aggregate, 0, 1);
        manager.run_vector_user_single_gid<int, int, std::vector<int>>
          ("non-contiguous int", TestMode::Aggregate, 20, 3);
        manager.run_vector_user_single_gid<long long, int, std::vector<int>>
          ("non-contiguous  long long", TestMode::Aggregate, 200, 4);

        // Aggregate mode is for vector user type only
        if(test_mode != Aggregate) {
          manager.run_single_user_multiple_gid<gid_set_t, lid_set_t, int>
            (GID_SET_LENGTH, LID_SET_LENGTH,  "contiguous int", test_mode, 0, 1);
          manager.run_single_user_multiple_gid<gid_set_t, lid_set_t, int>
            (GID_SET_LENGTH, LID_SET_LENGTH, "non-contiguous int", test_mode, 20, 3);
        }

        manager.run_vector_user_multiple_gid<gid_set_t, lid_set_t, std::vector<int>>
          (GID_SET_LENGTH, LID_SET_LENGTH, "contiguous int", test_mode, 0, 1);
  #endif

        if(!manager.did_all_pass()) {
          err = 1; // if we fail at some point just drop out
        }
      }
    }
  }

  // Get the global err results so we fail properly if any proc failed
  int errGlobal;
  Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM, err,
    Teuchos::outArg(errGlobal));
  return errGlobal; // only 0 if all tests and all proc return 0
}
