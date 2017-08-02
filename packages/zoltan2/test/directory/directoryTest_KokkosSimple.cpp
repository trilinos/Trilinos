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

// Eventually this macro will be removed and only the Kokkos mode will exist
// The other modes (Original, Tpetra, Relic) are for development and
// performance testing purposes
#define CONVERT_DIRECTORY_KOKKOS
#include "Zoltan2_Directory_Impl.hpp"

int main(int narg, char **arg)
{
  Teuchos::GlobalMPISession mpiSession(&narg,&arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  // will reduce err on all ranks
  int err = 0;

  // define a few constants/helpers
  const bool bUseLocalIDs = false;
  const bool contiguous = false;
  typedef int lid_t; // fix this for now - not dealing with the lid testing yet

  // setup some empty vectors (this will be passed to mean not used)
  // however maybe will go back to a ptr API which allows for NULL as an 'off'
  // this seems akward to have
  std::vector<lid_t> ignore_lid; // TODO: decide how to best handle this API
  std::vector<int> ignore_int;   // TODO: decide how to best handle this API

  // for this first test just do a simple update and find
  {
    // first set up some types
    typedef int gid_t;
    typedef long user_t;

    // set up the directory type - this is the single user mode (not vector user)
    typedef Zoltan2::Zoltan2_Directory_Simple<gid_t,lid_t,user_t> directory_t;

    // create the directory
    directory_t directory(comm, bUseLocalIDs, 0, contiguous);

    // now create some gids to write
    // to keep things interesting we'll have just proc 0 write 4 Ids and all the
    // other procs right the same gid (3).
    std::vector<gid_t> writeIds;

    if(comm->getRank() == 0) {
      writeIds =  { 1, 5, 7, 10 };
    }
    else {
      writeIds = { 3 };
    }

    // now create some user values associated with above gids
    // for easy reference just make them gid x 10
    std::vector<user_t> writeUser;
    if(comm->getRank() == 0) {
      writeUser = { 10, 50, 70, 100 };
    }
    else {
      writeUser = { 30 };
    }

    // call update on the directory
    directory.update(writeIds, ignore_lid, writeUser, ignore_int,
      directory_t::Replace);

    // now pick some gids to find
    std::vector<gid_t> findIds = { 3, 5 };

    // now create a user space to accept the values
    // Setting this empty will turn off this option. The original directory uses
    // ptrs and null which has some advantages so this API choice might need some
    // better setup.
    std::vector<user_t> findUser(findIds.size());

    // now call find which will fill findUser
    // for special case of mpi proc 1 the find of gid 3 will give a throw
    // because rank 1 never existed to write it. This test is currently setup
    // so it will pass for any mpi count so for this special case we turn this
    // into a test to verify the throw happened. If mpi count is not 1, then
    // a throw will not happen (or if it does it's a real error).
    // However maybe it's better to change this so the simple test example
    // doesn't have this complication ... to do.
    try {
      directory.find(findIds, ignore_lid, findUser, ignore_int, ignore_int);

      // now check element 0 in the array - make sure it matches writeUser above
      // all procs except 0 sent this one
      if(findUser[0] != 30) {
        std::cout << "We should have gotten 30 for this value." << std::endl;
        ++err;
      }

      // only proc 0 sent the second one - make sure all procs can read it
      if(findUser[1] != 50) {
        std::cout << "We should have gotten 50 for this value." << std::endl;
        ++err;
      }
    }
    catch(std::logic_error) {
      if(comm->getSize() != 1) {
        err = 1;
      }
      else {
        std::cout << "Successfully detected throw since find failed." << std::endl;
      }
    }
  }

  // now do an example using the variable vector length
  {
    // first set up some types
    typedef long gid_t;
    typedef std::vector<long long> user_t; // now user is a vector of long long

    // set up the directory type - this is the vector user mode
    // the way things are setup right now a std::vector user_t must use this
    // class - maybe eventually this could all be handling by templating but
    // I think this might be cleaner and more performant
    typedef Zoltan2::Zoltan2_Directory_Vector<gid_t,lid_t,user_t> directory_t;

    // create the directory
    directory_t directory(comm, bUseLocalIDs, 0, contiguous);

    // now create some gids to write based on the rank
    std::vector<gid_t> writeIds;
    switch(comm->getRank()) {
      case 0:
        writeIds = { 3, 5 }; // rank 0 writes only 3 and 5
        break;
      case 1:
        writeIds = { 3, 7 }; // rank 1 writes only 3 and 7
        break;
      default:
        writeIds = { 5 }; // all other ranks write only 5
        break;
    }

    // now create some user values associated with above gids
    // in this case we can make the values different lengths
    std::vector<user_t> writeUser;
    switch(comm->getRank()) {
      case 0:
        writeUser = { {1,2}, {4,8,10} }; // rank 0 writes only 3 and 5
        break;
      case 1:
        writeUser = { {2,10}, {6,8} }; // rank 1 writes only 3 and 7
        break;
      default:
        writeUser = { {1,2,10} }; // all other ranks write only 5
        break;
    }

    // call update on the directory
    directory.update(writeIds, ignore_lid, writeUser, ignore_int,
      directory_t::Aggregate);

    // now pick some gids to find - could make this rank specific
    std::vector<gid_t> findIds = { 3, 5 };

    // now create a user space to accept the values
    // Setting this empty will turn off this option. The original directory uses
    // ptrs and null which has some advantages so this API choice might need some
    // better setup. Note that findUser will be a set of empty std::vector
    // in this case and each will be individually filled by the directory with
    // the correct aggregated result
    std::vector<user_t> findUser(findIds.size());

    // now call find which will fill findUser
    directory.find(findIds, ignore_lid, findUser, ignore_int, ignore_int);

    // now check element 0 in findUser
    // that was for gid = 3 so from the above we expect the following:
    //   rank 0 provided user value {1,2}
    //   rank 1 provided user value {2,10}
    //   all other ranks provided nothing for gid 3
    // Therefore the aggregated result should be {1,2} + {2,10} = {1,2,10}
    // For a little variation run as MPI 1 proc, then rank 1 will never write
    //   and the expected value is just {1,2}.
    user_t expectedValue0 =
      (comm->getSize() == 1) ? user_t({1,2}) : user_t({1,2,10});
    if(findUser[0] != expectedValue0) {
      std::cout << "findUser[0] did not match expected." << std::endl;
      ++err;
    }

    // now check element 1 in findUser
    // that was for gid = 5 so from the above we expect the following:
    //   rank 0 provided user value {4,8,10}
    //   rank 1 provided nothing
    //   all other ranks provided {1,2,10}
    // Therefore the aggregated result should be {4,8,10} + {1,2,10} = {1,2,4,8,10}
    // Again for variation can run this with proc count < 2, in which case
    //   the expected resut will be just {4,8}.
    user_t expectedValue1 =
      (comm->getSize() <= 2) ? user_t({4,8,10}) : user_t({1,2,4,8,10});
    if(findUser[1] != expectedValue1) {
      std::cout << "findUser[1] did not match expected." << std::endl;
      ++err;
    }
  }

  // now do an Add example and try gid with sub elements gid = {a,b,c,...}
  {
    // first set up some types - make the gid type be a full structure with
    // sub elements - as long as the directory can call sizeof(gid_t) and
    // do things like a = b for gid_t, it's fine. It can't contain any memory
    // allocations (so it could not contain a std::vector for example).
    class gid_struct {
      public:
      // Kokkos will need an empty constructor so either provide one or
      // make a simple struct with no constructors will also work
      gid_struct() {}
      gid_struct(int v0, int v1, int v2, int v3) {
        val[0] = v0;
        val[1] = v1;
        val[2] = v2;
        val[3] = v3;
      }
      int val[4]; // can be any length but this can't have mem alloc (like std::vector)
    };
    typedef gid_struct gid_t;
    typedef int user_t;

    // set up the directory type - this is the single user mode (not vector user)
    typedef Zoltan2::Zoltan2_Directory_Simple<gid_t,lid_t,user_t> directory_t;

    // create the directory
    directory_t directory(comm, bUseLocalIDs, 0, contiguous);

    // set up some test gids
    // in this case the entire set of values is the gid
    std::vector<gid_t> writeIds;
    gid_t gid1(1, 8, 7, 3);
    gid_t gid2(1, 8, 7, 4); // this is a completely different gid from the prior

    // let's have rank 0 write gid1 and gid2 while other ranks write gid2 only
    if(comm->getRank() == 0) {
      writeIds =  { gid1, gid2 };
    }
    else {
      writeIds = { gid2 };
    }

    // now create some user values associated with above gids
    // since this will be an add test, just write 1 for all of them
    // we expect gid1 to end up unchanged (only rank 0 writes it)
    // but we expect gid2 to end up with value equal to comm->getSize()
    // because all ranks will contribute 1 to the sum
    std::vector<user_t> writeUser;
    if(comm->getRank() == 0) {
      writeUser = { 1, 1 }; // two values because rank 0 write both gid1 and gid2
    }
    else {
      writeUser = { 1 }; // one value because other ranks only write gid2
    }

    // call update on the directory using add mode
    directory.update(writeIds, ignore_lid, writeUser, ignore_int,
      directory_t::Add);

    // now check them on all ranks - this could be different for different ranks
    std::vector<gid_t> findIds = { gid1, gid2 };

    // create a user space to accept the values
    // Setting this empty will turn off this option. The original directory uses
    // ptrs and null which has some advantages so this API choice might need some
    // better setup.
    std::vector<user_t> findUser(findIds.size());

    // now call find which will fill findUser
    directory.find(findIds, ignore_lid, findUser, ignore_int, ignore_int);

    // now check element 0 in the array which should have value 1
    if(findUser[0] != 1) {
      std::cout << "We should have gotten 1 for gid1. Only rank 0 wrote to gid1"
        " so the sum should have been 1." << std::endl;
      ++err;
    }

    // now check element 1 in the array should have value comm->getSize()
    // that is because each rank contributed 1 and we are in Add mode
    if(findUser[1] != comm->getSize()) {
      std::cout << "We should have gotten the proc count " << comm->getSize()
        << " since every rank wrote to gid2." << std::endl;
      ++err;
    }
  }

  // Get the global err results so we fail properly if any proc failed
  int errGlobal;
  Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM, err,
    Teuchos::outArg(errGlobal));

  // this proc is ok
  comm->barrier();
  if(comm->getRank() == 0) {
    if(errGlobal == 0) {
      std::cout << "Passed" << std::endl;
    }
    else {
      std::cout << "FAILED!" << std::endl;
    }
  }
  return errGlobal;
}