/*
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan2 Directory for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of theremove_local
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

// TODO: A note on the file structure
// We currently have Zoltan2_Directory.hpp and Zoltan2_Directory_Impl.hpp
// This evolved and the point of the second file is now doubtful so probably
// should be all merged into one, especially since all templated now.
// However there might be utility in keeping the API separate.
#include "Zoltan2_Directory.hpp"

// original points back to zoltan code - block this out completely to make
// sure we don't have any dependencies - eventually this gets deleted
#ifndef CONVERT_DIRECTORY_ORIGINAL

// Relic mode uses the new directory class but points back to the original
// zoltan communicator. It's a halfway point between the new Kokkos mode
// and the original mode which was useful for debugging some stuff.
#ifdef CONVERT_DIRECTORY_RELIC
  #include <comm.h>
#endif

// Only Kokkos mode uses the new Zoltan2_Directory_Comm class which is pulled
// from the original zoltan code and modified to handle variable length
// user data.
#ifdef CONVERT_DIRECTORY_KOKKOS
  #include "Zoltan2_Directory_Comm.hpp"
#endif

// Supplies the current hash - rolled over from zoltan
// Right now the hash_proc method is still being used for the new Kokkos mode
// which determines which proc will own a gid. TODO: Make this use a Kokkos
// setup. The Kokkos mode is using the Kokkos unordered map for the node data.
// I added the declaration here as inline to avoid the warnings generated
// otherwise. All this should eventually go away and either use Kokkos methods
// or implement the hash here so we don't have the links to zoltan
static inline void MurmurHash3_x86_32 (
  const void * key, int len, uint32_t seed, void * out);
#include "murmur3.c"

namespace Zoltan2 {

// These macros were rolled over from the original zoltan code. I've preserved
// them for now until we have further discussion on how we want debug levels
// and logging to work in the new code. This should just be debug related and
// not impact the behavior of the code.
#define ZOLTAN2_PRINT_INFO(proc,yo,str) \
  printf("ZOLTAN2 (Processor %d) %s: %s\n", (proc), (yo), \
         ((str) != NULL ? (char *)(str) : " "));

#define ZOLTAN2_TRACE(proc,where,yo,str) \
  printf("ZOLTAN (Processor %d) %s %s  %s\n", (proc), (where), (yo), \
         ((str) != NULL ? (char *)(str) : " "));

#define ZOLTAN2_TRACE_IN(proc,yo,str) \
  ZOLTAN2_TRACE((proc),"Entering",(yo),(str));

#define ZOLTAN2_TRACE_OUT(proc,yo,str) \
  ZOLTAN2_TRACE((proc),"Exiting",(yo),(str));

// Tags for MPI communications.  These need unique values. Arbitrary
// These were rolled over from the original zoltan code and still used.
#define ZOLTAN2_DD_FIND_MSG_TAG     29137  /* needs 3 consecutive values */
#define ZOLTAN2_DD_UPDATE_MSG_TAG   29140  /* needs 2 consecutive values */
#define ZOLTAN2_DD_REMOVE_MSG_TAG   29142  /* needs 2 consecutive values */
#define ZOLTAN2_DD_RESIZE_MSG_TAG   29150  /*  */

template <typename gid_t,typename lid_t,typename user_t>
void Zoltan2_Directory<gid_t,lid_t,user_t>::allocate()
{
  const char * yo = "Zoltan2_Directory::allocate";

#ifdef CONVERT_DIRECTORY_RELIC
  if (table_length < 0)  {
    throw std::invalid_argument(
     "Zoltan2_Directory() Invalid table_length argument");
  }
#endif

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

  /* insure all processors are using the same GID, LID, USER lengths */
  size_t array[3], max_array[3], min_array[3];
  array[0] = sizeof(lid_t);
  array[1] = sizeof(gid_t);
  array[2] = sizeof(user_t);
  Teuchos::reduceAll<int,size_t>(
    *comm, Teuchos::REDUCE_MAX, 3, array, max_array);
  Teuchos::reduceAll<int,size_t>(
    *comm, Teuchos::REDUCE_MIN, 3, array, min_array);
  if (max_array[0] != min_array[0] || max_array[1] != min_array[1]
    || max_array[2] != min_array[2])  {
    throw std::invalid_argument(
      "Zoltan2_Directory() LID, GID, USER data lengths differ globally");
  }

#ifdef CONVERT_DIRECTORY_RELIC
  if(!table_length) {
    table_length = 100000; // relic value for hash table count
  }
  table_length = recommended_hash_size(table_length);

  // Note that in the original Zoltan implementation memory was allocated for
  // the main struct with additional memory for the table - so this step mimics
  // that part of the code where we allocate memory for the table
  table = Teuchos::arcp(new relice_idx_t[table_length], 0, table_length, true);
  for(Teuchos::ArrayRCP<relice_idx_t>::size_type i = 0; i < table.size(); ++i) {
    table[i] = -1;
  }

  nodelistlen = 0;
  nextfreenode = -1;
#endif

  // get the base size of the user data component
  // for simple mode, this is going to just be the sizeof(user_t)
  // for vector mode, this is sizeof(size_t) for the std::vector length
  // the actual data will be variable
  size_t user_base_size = is_Zoltan2_Directory_Vector() ? sizeof(size_t) :
    size_of_value_type();

  // set up the base size to include the gid and the lid (if used) + user size
  size_t size = sizeof(gid_t) + (use_lid?sizeof(lid_t):0) + user_base_size;

#ifdef CONVERT_DIRECTORY_RELIC
  nodedata_size = size;
#endif

  // now calculate the full update msg size
  update_msg_size = size + sizeof(Zoltan2_DD_Update_Msg<gid_t,lid_t>);

  // for remove message we just need the gid, not any other info
  size = sizeof(gid_t);
  remove_msg_size = size + sizeof(Zoltan2_DD_Remove_Msg<gid_t,lid_t>);

  // Current form of find_local is passed so gid_t is the input and
  // lid_t is the output so this guarantees that ptr is sufficient to
  // cover both (uses the max). This handling is consistent with the original
  // zoltan code but I initially found it very confusing.
  // I suggest searching for all of the places where max_id_size is used and
  // that represents the relevant code. This may be the most optimal way of
  // doing this. I did not get to assessing it in detail.
  //
  // NOTE: To really test this thoroughly we would want to play around with
  // the values GID_SET_LENGTH and LID_SET_LENGTH in the directoryTest_Impl.hpp
  // setup. I made sure this works for gid larger than lid and the opposite.
  // Probably should extend the unit tests to cover all these cases at once as
  // originally I had some bugs where this would work fine except when lid
  // size was greater than gid size. Now this all should be ok.
  max_id_size = std::max(sizeof(gid_t),sizeof(lid_t));

  size = max_id_size + user_base_size;
  find_msg_size = size + sizeof(Zoltan2_DD_Find_Msg<gid_t,lid_t>);

  // TODO: Alignment is currently off for all of the modes for performance
  // comparisons. I was not sure yet the implications of this and how it would
  // best integrate if we have variable length user data
/*
#ifdef CONVERT_DIRECTORY_RELIC
    nodedata_size   = align_size_t(nodedata_size);
#endif
    update_msg_size = align_size_t(update_msg_size);
    remove_msg_size = align_size_t(remove_msg_size);
    find_msg_size   = align_size_t(find_msg_size);
*/

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }
}

// Copy Block
template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::copy(
  const Zoltan2_Directory<gid_t,lid_t,user_t> & src) {

#ifdef CONVERT_DIRECTORY_RELIC
  nodelist = src.nodelist;
  nodedata = src.nodedata;
#endif

#ifdef CONVERT_DIRECTORY_KOKKOS
  node_map = src.node_map;
#endif

#ifdef CONVERT_DIRECTORY_TPETRA
  oto_idMap = src.oto_idMap;
  oto_idVec = src.oto_idVec;
#endif

  return 0;
}

// Not yet implemented - was in progress and will be part of future work.
// Two unit tests are disabled in CMakeLists.txt for tests/directory
//  directoryTest_KokkosOptimized
//  directoryTest_PerformanceKokkosOptimized
// Those tests turn on this flag CONVERT_DIRECTORY_KOKKOS_OPTIMIZED
// Eventually this would become the final version - this was in development.
// Currently this code doesn't compile since I moved update_local back to being
// protected. It needs to be available here as public and I didn't finish
// figuring out how that should be setup. I had started setting up a first test
// case where threads would optimize the application of update_local. Sort of
// got it running but there was some issue with the data and I expect I'm making
// bad assumptions about what is thread safe. The goal here was to get a simple
// application of functors which would speed up the update_local calls.
#ifdef CONVERT_DIRECTORY_KOKKOS_OPTIMIZED
template <typename gid_t,typename lid_t,typename user_t>
class update_local_functor {
public:
  typedef Kokkos::UnorderedMap<gid_t,
    Zoltan2_Directory_Node<gid_t,lid_t,user_t>> node_map_t;
  typedef Zoltan2_DD_Update_Msg<gid_t,lid_t> msg_t;
  // Initialize all members
  update_local_functor(
    Teuchos::ArrayRCP<msg_t*> pMessagesIn,
    Zoltan2_Directory<gid_t,lid_t,user_t> * pDirectoryIn,
    const std::vector<user_t> * pUserIn) :
      pMessages(pMessagesIn),
      pDirectory(pDirectoryIn),
      pUser(pUserIn) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {
    msg_t *ptr = pMessages[i];

    user_t * puser = (ptr->user_flag) ?
      reinterpret_cast<user_t*>(reinterpret_cast<char*>(ptr->adjData) +
        sizeof(gid_t) + (pDirectory->is_use_lid()?sizeof(lid_t):0)) : NULL;

    pDirectory->update_local(ptr->adjData,
      (ptr->lid_flag) ? reinterpret_cast<lid_t*>(ptr->adjData + 1) : NULL,
      puser,
      (ptr->partition_flag) ? (ptr->partition) : -1,  // illegal partition
      ptr->owner);
  }

  Teuchos::ArrayRCP<msg_t*> pMessages;
  Zoltan2_Directory<gid_t,lid_t,user_t> * pDirectory;
  const std::vector<user_t> * pUser;
};
#endif

template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::update(
  const std::vector<gid_t>& gid, const std::vector<lid_t>& lid,
  const std::vector<user_t>& user, const std::vector<int>& partition,
  Update_Mode update_mode)
{
  // for conveniece store but maybe this should be a construct property
  // should a directory allow mixed modes (Replace, then Aggregate... for ex)
  this->mode = update_mode;

  Zoltan2_Directory_Clock clock("update");

  const char * yo = "Zoltan2_Directory::update";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN(comm->getRank(), yo, NULL);
  }

  Zoltan2_Directory_Clock update_setup_clock("update_setup", 1);

  // part of initializing the error checking process
  // for each linked list head, walk its list resetting errcheck

#ifdef CONVERT_DIRECTORY_KOKKOS
  if(debug_level) {
    for(size_t n = 0; n < node_map.size(); ++n) {
      node_map.value_at(n).errcheck = -1; // not possible processor
    }
  }
#endif

#ifdef CONVERT_DIRECTORY_RELIC
  if(debug_level) {
    for (int i = 0; i < table_length; i++) {
      relice_idx_t nodeidx;
      for (nodeidx = table[i]; nodeidx != -1;
        nodeidx = nodelist[nodeidx].next) {
        nodelist[nodeidx].errcheck = -1;  // not possible processor
      }
    }
  }
#endif


  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After reset errcheck");
  }

  update_setup_clock.complete();

  int err = 0;

#ifdef CONVERT_DIRECTORY_TPETRA

  // The purpose of this code was to implement a Tpetra form of the directory
  // so we could do some performance comparisons. The Tpetra mode will support
  // Add and Replace (not Aggregate) and can be run using the same API calls.

  // Compute global indexBase
  gid_t minId = gid.size() ?
    (*std::min_element(std::begin(gid), std::end(gid))) : 0;
  gid_t indexBase;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &minId, &indexBase);

  // Build Tpetra::Map for the given local ids
  size_t dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  rcp_map_t idMap = Teuchos::rcp(new map_t(dummy, gid, indexBase, comm), true);

  // Create Vector using this map.
  rcp_vector_t idVec = Teuchos::rcp(new vector_t(idMap, 0.));

  // Now set the tpetra data to match the input data
  vectordata_t idData = idVec->getDataNonConst();
  for (size_t i = 0; i < idVec->getLocalLength(); i++) {
    idData[i] = user[i];
  }

  // I pulled this for now as I didn't want to have the API requiring the
  // contiguous option just for Tpetra which is just a side show here anyways.
  // If we want to restore this optimization we might scan the ids to auto
  // determine if they are contiguous (or make a different API for Tpetra)
  // I used to pass this in from the unit tests but eliminated that.
  /*
  if (contiguous) {
    // For contigous ids, can use Tpetra default Map
    gid_t gnUnique = idMap->getMaxAllGlobalIndex() - indexBase + 1;
    oto_idMap = Teuchos::rcp(new map_t(gnUnique, indexBase, comm));
  }
  else
  */
  {
    // Since ids may not be contiguous, cannot use default Tpetra Map
    oto_idMap = Tpetra::createOneToOne(idMap);
  }

  oto_idMap = Tpetra::createOneToOne(idMap);

  oto_idVec = Teuchos::rcp(new vector_t(oto_idMap));

  // Create an exporter between the two maps
  rcp_export_t idExporter = Teuchos::rcp(new export_t(idMap, oto_idMap));

  // Now we can use Tpetra for Replace or Add (doesn't support Aggregate)
  switch(mode) {
    case Update_Mode::Replace:
      oto_idVec->doExport(*idVec, *idExporter, Tpetra::REPLACE);
      break;
    case Update_Mode::Add:
      oto_idVec->doExport(*idVec, *idExporter, Tpetra::ADD);
      break;
    case Update_Mode::Aggregate:
      throw std::logic_error("Tpetra mode won't support Aggregate.");
      break;
  }

#else

  Zoltan2_Directory_Clock update_alloc_clock("update_alloc", 1);

  // allocate memory for list of processors to contact
  Teuchos::ArrayRCP<int> procs;
  if(gid.size() > 0) {
    procs = Teuchos::arcp(new int[gid.size()], 0, gid.size(), true);
  }

  // set up the msg sizes for vector mode only
  Teuchos::ArrayRCP<int> msg_sizes;
  int sum_msg_size = 0;
  if(is_Zoltan2_Directory_Vector() && gid.size() > 0) {
    msg_sizes = Teuchos::arcp(new int[gid.size()], 0, gid.size(), true);
    for (size_t i = 0; i < gid.size(); i++) {
      size_t msg_size = get_update_msg_size(user[i]);
      sum_msg_size += msg_size;
      msg_sizes[i] = msg_size;
    }
  }
  else {
    sum_msg_size = update_msg_size * gid.size(); // simple case
  }

  update_alloc_clock.complete();

  Zoltan2_Directory_Clock update_alloc_sbuff_clock("update_alloc_sbuff", 1);

  Teuchos::ArrayRCP<char> sbuff;
  if(sum_msg_size) {
    sbuff = Teuchos::arcp(new char[sum_msg_size], 0, sum_msg_size, true);
  }

  update_alloc_sbuff_clock.complete();

  typedef Zoltan2_DD_Update_Msg<gid_t,lid_t> msg_t;

  Zoltan2_Directory_Clock update_build_raw_clock("update_build_raw", 1);

  // build update messages
  int track_offset = 0;
  char * trackptr = sbuff.getRawPtr();
  for (size_t i = 0; i < gid.size(); i++) {
    // hash the gid to determine which proc will own it
    procs[i] = hash_proc(gid[i]);

    // this ptr is the beginning of the message which can be different
    // lengths for different gids if using Zoltan2_Directory_Vector
    msg_t *ptr = reinterpret_cast<msg_t*>(trackptr);

    // write in all my info
    ptr->lid_flag       = (lid.size())  ? 1 : 0;
    ptr->user_flag      = (user.size()) ? 1 : 0;
    ptr->partition_flag = (partition.size()) ? 1 : 0;
    ptr->partition      = (partition.size()) ? partition[i] :  -1;
    ptr->owner          = comm->getRank();

    // now deal with the gid
    gid_t * pgid = ptr->adjData;
    *pgid = gid[i];

    // optionally write the lid if we are using that
    if(use_lid) {
      if(!lid.size()) {
        throw std::logic_error(
          "Did not pass lid values but directory was created to use them!");
      }
      lid_t * plid = reinterpret_cast<lid_t*>(ptr->adjData + 1);
      *plid = lid[i];
    }
    else {
      if(lid.size()) {
        throw std::logic_error(
          "Passed lid values but directory was created not to use them!");
      }
    }

    // find the spot where the user data begins
    user_t * puser = reinterpret_cast<user_t*>(
      reinterpret_cast<char*>(ptr->adjData) + sizeof(gid_t) +
        (use_lid?sizeof(lid_t):0));

    // write in the user data - for Zoltan2_Directory_Simple this is a trival
    // copy but for Zoltan2_Directory_Vector we write length and then write all
    // of the elements in the vector.
    if (user.size()) {
      user_to_raw(user[i], puser);
    }
    else {
      // The update msg contains space for the vector length (sizeof(size_t))
      // if it's a Zoltan2_Directory_Vector
      *puser = user_t(); // create an empty result
    }

    // vector will have different lengths but the simple mode will have all
    // lengths just update_msg_size
    size_t new_update_msg_size =
      is_Zoltan2_Directory_Vector() ? msg_sizes[i] : update_msg_size;
    track_offset += new_update_msg_size;
    trackptr += new_update_msg_size;
  }

  // this check just makes sure our internal logic above is correct
  // we should have looped to total size sum_msg_size
  if(track_offset != sum_msg_size) {
    throw std::logic_error("Bad summing!");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After fill contact list");
  }

  update_build_raw_clock.complete();

  // now create efficient communication plan

  Zoltan2_Directory_Clock update_build_plan_clock("update_build_plan", 1);

#ifdef CONVERT_DIRECTORY_KOKKOS
  // Kokkos mode uses the new refactored C++ communicator class
  Zoltan2_Directory_Comm directoryComm(gid.size(), procs, comm,
    ZOLTAN2_DD_UPDATE_MSG_TAG);
  int nrec = directoryComm.getNRec();
#else
  // The relic mode will be working with the original zoltan code
  ZOLTAN_COMM_OBJ *plan  = NULL;   // for efficient MPI communication
  int nrec = 0;       // number of receives to expect
  err = Zoltan_Comm_Create (&plan, gid.size(), procs.getRawPtr(),
    Teuchos::getRawMpiComm(*comm), ZOLTAN2_DD_UPDATE_MSG_TAG, &nrec);
#endif

  if (err) {
    throw std::logic_error("Zoltan2_Directory::update() Comm_Create error");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Create");
  }

  update_build_plan_clock.complete();

  Zoltan2_Directory_Clock update_resize_clock("update_resize", 1);

  int sum_recv_sizes = 0;
  if(is_Zoltan2_Directory_Vector()) {
    // Only vector mode has to use the resizing options for getting receive info
#ifdef CONVERT_DIRECTORY_KOKKOS
    err = directoryComm.resize(msg_sizes,
      ZOLTAN2_DD_RESIZE_MSG_TAG, &sum_recv_sizes);
#else
    err = Zoltan_Comm_Resize(plan, msg_sizes.getRawPtr(),
      ZOLTAN2_DD_RESIZE_MSG_TAG, &sum_recv_sizes);
#endif
  }
  else {
    sum_recv_sizes = update_msg_size * nrec;
  }

  if (err) {
    throw std::logic_error("directoryComm.execute_resize error");
  }

  update_resize_clock.complete();

  // If dd has no nodes allocated (e.g., first call to DD_Update;
  // create the nodelist and freelist
  Zoltan2_Directory_Clock build_hash_clock("update_build_hash", 1);

#ifdef CONVERT_DIRECTORY_RELIC
  if (nrec && nodelistlen == 0) {
    allocate_node_list((relice_idx_t) nrec, 0.);
  }
#endif

#ifdef CONVERT_DIRECTORY_KOKKOS
  // TODO upgrade nrec as size_t and all corresponding changes...
  if(nrec && static_cast<int>(node_map.size()) < nrec) {
    // some things to consider here if we will have multiple update calls in
    // series... how will subsequent calls optimally manage this list since the
    // new gid set may be partially overlapped with the original. Currently the
    // update_local has a mechanism to rehash and increase this size if we run
    // out so skipping this call would be logically ok (but not best performance)
    node_map.rehash(nrec);
  }
#endif

  build_hash_clock.complete();

  // allocate receive buffer for nrec DD_Update_Msg structures
  Teuchos::ArrayRCP<char> rbuff;
  if(sum_recv_sizes > 0) {
    rbuff = Teuchos::arcp(
      new char[sum_recv_sizes], 0, sum_recv_sizes, true);   // receive buffer
  }

  // send my update messages & receive updates directed to me
  //if resizing we send size 1 because the sizes will be built individually
  const int nbytes = is_Zoltan2_Directory_Vector() ? 1 : update_msg_size;

  Zoltan2_Directory_Clock update_forward_clock("update_forward", 1);

#ifdef CONVERT_DIRECTORY_KOKKOS
  err = directoryComm.do_forward(ZOLTAN2_DD_UPDATE_MSG_TAG+1,
    sbuff, nbytes, rbuff);
#else
  err = Zoltan_Comm_Do(plan, ZOLTAN2_DD_UPDATE_MSG_TAG+1,
    sbuff.getRawPtr(), nbytes, rbuff.getRawPtr());
#endif

  update_forward_clock.complete();

  if (err) {
    throw std::logic_error("Zoltan2_Directory::update() Comm_Do error");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Do");
  }

  Zoltan2_Directory_Clock update_locals_clock("update_locals", 1);

  int errcount = 0;

#ifdef CONVERT_DIRECTORY_KOKKOS_OPTIMIZED
  // Note: This code was in progress - not live yet
  // See note above where CONVERT_DIRECTORY_KOKKOS_OPTIMIZED is first used.

  // build a map of message to buff location
  // I wanted to pass this to a functor and then have the functor process the
  // update_local calls in parallel using threads - this would be a first
  // example of optimizing the directory. The parallel operations can happen
  // any order so we need to precalculate the buffer locations since each
  // msg can be a different size. This wasn't working yet.
  Teuchos::ArrayRCP<msg_t*> message_locations(nrec);
  trackptr = rbuff.getRawPtr();
  for(int n = 0; n < nrec; ++n) {
    msg_t *ptr = reinterpret_cast<msg_t*>(trackptr);

    user_t * puser = (ptr->user_flag) ?
      (user_t*)(reinterpret_cast<char*>(ptr->adjData) +
        sizeof(gid_t) + (use_lid?sizeof(lid_t):0)) : NULL;
    size_t delta_msg_size = get_update_msg_size(puser);
    message_locations[n] = reinterpret_cast<msg_t*>(trackptr);
    trackptr += delta_msg_size;
  }

  Kokkos::parallel_for (nrec,
    update_local_functor<gid_t,lid_t,user_t>(message_locations, this, &user));
#else

  // for each message rec'd, update local directory information
  track_offset = 0;
  trackptr = rbuff.getRawPtr();
  for (int i = 0; i < nrec; i++) {
    msg_t *ptr = reinterpret_cast<msg_t*>(trackptr);

    user_t * puser = (ptr->user_flag) ?
      (user_t*)(reinterpret_cast<char*>(ptr->adjData) +
        sizeof(gid_t) + (use_lid?sizeof(lid_t):0)) : NULL;

    err = update_local(ptr->adjData,
      (ptr->lid_flag) ? reinterpret_cast<lid_t*>(ptr->adjData + 1) : NULL,
      puser,
      (ptr->partition_flag) ? (ptr->partition) : -1,  // illegal partition
      ptr->owner);

    if (err)
      ++errcount;

    // in this case we are reading the raw data so we calculate the message
    // size from it. For vector mode, this will find the length of the vector
    // and calculate the full message size from that.
    size_t delta_msg_size = get_update_msg_size(puser);
    trackptr += delta_msg_size;
    track_offset += delta_msg_size;
  }

  // safety check - if all this internal logic is correct the ptr offset should
  // not be exactly at the end of the recv buffer.
  if(track_offset != sum_recv_sizes) {
    throw std::logic_error("Did not sum!");
  }
#endif

  update_locals_clock.complete();

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Local update");
  }

  err = 0;

  if (debug_level) {  // overwrite error return if extra checking is on
    err = (errcount) ? 1 : 0;
  }

#ifdef CONVERT_DIRECTORY_RELIC
  Zoltan_Comm_Destroy (&plan);
#endif

  if (debug_level)  {
    char str[100];      // used to build message string
    sprintf (str, "Processed %lu GIDs (%d local), %d GID errors", gid.size(),
#ifdef CONVERT_DIRECTORY_RELIC
      nrec,
#else
      directoryComm.getNRec(),
#endif
      errcount);
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);
  }

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT(comm->getRank(), yo, NULL);
  }

#endif
  return err;
}

#ifndef CONVERT_DIRECTORY_TPETRA
template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::update_local(
 gid_t* gid,                 /* GID to update (in)                        */
 lid_t* lid,                 /* gid's LID (in), NULL if not needed        */
 user_t *user,               /* gid's user data (in), NULL if not needed  */
 int partition,              /* gid's partition (in), -1 if not used      */
 int owner)                  /* gid's current owner (proc number) (in)    */
{
  const char * yo = "Zoltan2_Directory::update_local";

  // input sanity checking
  if (owner < 0) {
    throw std::logic_error(
      "Zoltan2_Directory::update_local() owner < 0");
  }
  if (owner >= comm->getSize()) {
    throw std::logic_error(
      "Zoltan2_Directory::update_local() owner >= comm->getSize()");
  }
  if (gid == NULL) {
    throw std::logic_error(
      "Zoltan2_Directory::update_local() gid == NULL");
  }

  if (debug_level > 5) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

  // compute offset into hash table to find head of linked list
#ifdef CONVERT_DIRECTORY_RELIC
  int index = hash_table(*gid);
  // walk linked list until end looking for matching gid
  for (relice_idx_t nodeidx = table[index]; nodeidx != -1;
    nodeidx = nodelist[nodeidx].next) {
    Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node = nodelist[nodeidx];
    if (*gid == *node.gid) {
#else
  size_t node_index = node_map.find(*gid);
  if(node_map.valid_at(node_index)) {
    Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node =
      node_map.value_at(node_index);
#endif

#ifdef CONVERT_DIRECTORY_RELIC
      // found match, update directory information
      if (lid) {
        lid_t* plid = reinterpret_cast<lid_t*>(node.gid + 1);
        *plid = lid[0];
      }
      if (user) {
        user_t * puser = reinterpret_cast<user_t*>(
          reinterpret_cast<char*>(node.gid) +
          sizeof(gid_t) + (use_lid?sizeof(lid_t):0));
        switch(mode) {
          case Update_Mode::Replace:
            *puser = *user;
            break;
          case Update_Mode::Add:
            *puser += *user;
            break;
          case Update_Mode::Aggregate:
            throw std::logic_error("Aggregate doesn't have meaning for single type.");
            break;
        }
      }
#else
      // found match, update directory information
      if (lid) {
        node.lid = *lid;
      }
      if (user) {
        update_local_user(user, node.userData);
      }
#endif

      node.owner = owner;
      if (partition != -1) {
        node.partition = partition;
      }

      // errcheck is reset to -1 at beginning of update for debug mode
      // then it will get set to the provider of the data when the node is
      // created or on the first call to be updated locally.
      // So if errcheck -1, then we do nothing - it's the first action.
      if(debug_level) {
        // if node.errcheck is -1 we are detecting the first update to a node
        // which already existed at the beginning of the update call so it's ok.
        // if mode is not Replace then duplicate calls are ok (expected) since
        // Add and Aggregate combine from many procs are the outcome is not
        // order dependent.
        if(node.errcheck != -1 && mode == Replace) {
          // here we have detected two actions on a single node in the same
          // update call for replace.
          // The actions could have been:    create node -> change node
          // or if the node already existed: change node -> change node

          // in Replace mode, two actions are potentially problematic.
          // If two procs update the same gid with different data it will
          // be order dependent.

          // For performance testing it's convenient to allow
          // Replace to be called from different procs but expect the data
          // to always be identical. That means we can compare Replace and
          // Aggregate more meaningfully since the gid setup for update is
          // the same.

          // To discuss - should one proc be allowed to submit duplicate
          // gids in an update call using Replace mode. Should two different
          // procs be allowed to submit duplicate gids with the same data
          // for a replace call, in which case the outcome is determined
          // regardless of order but we would be relying on the user to have
          // accurate data submission. Then we might consider a debug check
          // here to validate the data is coming in matches the local data.

          // TODO: Probably add this throw in, but currently the testing
          // framework will feed multiple Replace calls with the same data
          // from different gids - so we would have to change the tests
          // to guarantee Replace was a unique gid lists for each proc.
          // That is easy to do but not ideal at the moment for performance
          // testing reasons. I'd like the pattern of calls to be the same
          // as Add and Aggregate as a reference.

          //   throw std::logic_error( "Two replace calls were detected on."
          //     " the same gid which can be an undefined results.");
        }
      }

      node.errcheck = owner;

      if (debug_level > 5) {
        ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
      }

      // TODO: Change logic flow to avoid this internal return
      return 0;          // ignore all errors
#ifdef CONVERT_DIRECTORY_RELIC
    }
#endif
  }

  // gid not found.
  // Create new Zoltan2_Directory_Node<gid_t,lid_t> and fill it in
#ifdef CONVERT_DIRECTORY_RELIC
  relice_idx_t nodeidx = allocate_node();
  Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node = nodelist[nodeidx];
#else
  Zoltan2_Directory_Node<gid_t,lid_t,user_t> node; // will add to hash at end
#endif

#ifdef CONVERT_DIRECTORY_RELIC
  *node.gid = *gid;

  lid_t *plid = reinterpret_cast<lid_t*>(node.gid + 1);
  if (lid) {
    *plid = *lid;
  }
  else {
    *plid = 0;
  }

  user_t * puser = reinterpret_cast<user_t*>(
    reinterpret_cast<char*>(node.gid) + sizeof(gid_t)
      + (use_lid?sizeof(lid_t):0));
  if (user) {
    *puser = *user;
  }
  else {
    *puser = user_t();
  }
#else
  node.lid = lid ? (*lid) : lid_t();

  // this is more or less doing what above relic commands did except for
  // vector mode there is special handling to account for the variable length
  // of the std::vector
  if(user) {
    raw_to_user(user, node.userData);
  }
  else {
    node.userData = user_t();
  }

#endif

  node.partition = partition;
  node.owner = owner;
  node.errcheck = owner;

#ifdef CONVERT_DIRECTORY_RELIC
  // Add node to the linked list
  node.next = table[index];
  table[index] = nodeidx;
#else

  if(node_map.insert(*gid, node).failed()) {
    // Need more nodes (that's the assumption here if we have failure)
    // A new update has added more to our local list and we need more capacity.
    // TODO: Decide most efficient scheme. Here we bump to at least 10 or if
    // we're already at the level, increase by 10% increments. I think this will
    // be less efficient for small scale problems, when we probably care less,
    // but more efficient as we scale up.
    size_t new_guess_size = (node_map.size() < 10) ? 10 :
      ( node_map.size() + node_map.size()/10); // adds a minimum of 1
    node_map.rehash(new_guess_size);
    if(node_map.insert(*gid, node).failed()) {
      throw std::logic_error("Hash insert failed. Mem sufficient?....");
    }
  }

#endif

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, "Created new directory item");
  }

  if (debug_level > 5) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }

  return 0;
}
#endif

template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::find(
 const std::vector<gid_t>& gid, /* Incoming list of GIDs to get owners proc   */
 std::vector<lid_t>& lid,       /* Outgoing corresponding list of LIDs        */
 std::vector<user_t>& user,     /* Outgoing optional corresponding user data  */
 std::vector<int>& partition,   /* Outgoing optional partition information    */
 std::vector<int>& owner,       /* Outgoing optional list of data owners      */
 bool throw_if_missing)         /* default true - throw if gid is not found.  */
{
  Zoltan2_Directory_Clock clock("find");

  Zoltan2_Directory_Clock find_init_clock("  find_init");

  const char * yo = "Zoltan2_Directory::find";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN(comm->getRank(), yo, NULL);
  }

  int err = 0;

#ifdef CONVERT_DIRECTORY_TPETRA
  // Compute global indexBase
  gid_t minId =
    gid.size() ? (*std::min_element(std::begin(gid), std::end(gid))) : 0;
  gid_t indexBase;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &minId, &indexBase);

  // Build Tpetra::Map for the given ids
  size_t dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  rcp_map_t idMap = Teuchos::rcp(new map_t(dummy, gid, indexBase, comm), true);

  // Create Vector using this map.
  // This vector will store number of occurrences of each id across procs
  rcp_vector_t idVec = Teuchos::rcp(new vector_t(idMap, 0.));

  // Create an exporter between the two maps
  rcp_export_t idExporter = Teuchos::rcp(new export_t(idMap, oto_idMap));

  idVec->doImport(*oto_idVec, *idExporter, Tpetra::REPLACE);

  // Copy the result - TODO - make it pretty - make it fast
  vectordata_t idData = idVec->getDataNonConst();
  for (size_t i = 0; i < idVec->getLocalLength(); i++) {
    user[i] = idData[i];
  }
#else
  /* allocate memory for processors to contact for directory info */
  Teuchos::ArrayRCP<int> procs;
  if(gid.size() > 0) {
    procs = Teuchos::arcp(
      new int[gid.size()], 0, gid.size(), true); // processors to contact
  }

  /* Setup procs list */
  for (size_t i = 0; i < gid.size(); i++) {
    procs[i] = hash_proc(gid[i]); // determines the owner
  }

  // create efficient communication plan
#ifdef CONVERT_DIRECTORY_RELIC
  ZOLTAN_COMM_OBJ *plan  = NULL;     // efficient MPI communication
  int nrec;                          // number of messages to receive
  err = Zoltan_Comm_Create (&plan, gid.size(), procs.getRawPtr(),
    Teuchos::getRawMpiComm(*comm), ZOLTAN2_DD_FIND_MSG_TAG, &nrec);
#else
  Zoltan2_Directory_Comm directoryComm(gid.size(), procs, comm,
    ZOLTAN2_DD_FIND_MSG_TAG);
  int nrec = directoryComm.getNRec();
#endif

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Create");
  }

  if (err) {
    throw std::logic_error("Zoltan2_Directory::find() error");
  }

  find_init_clock.complete();

  Zoltan2_Directory_Clock find_sbuff_alloc_clock("  find_sbuff_alloc");

  Teuchos::ArrayRCP<char> sbuff;
  if(gid.size() > 0) {
    sbuff = Teuchos::arcp(new char[find_msg_size*gid.size()],
      0, find_msg_size*gid.size(), true); // send buffer
  }
  find_sbuff_alloc_clock.complete();

  Zoltan2_Directory_Clock find_prep_clock("  find_prep");

  typedef Zoltan2_DD_Find_Msg<gid_t,lid_t> msg_t;

  /* for each GID, fill DD_Find_Msg buffer and contact list */
  char *trackptr = sbuff.getRawPtr();
  for (size_t i = 0; i < gid.size(); i++)  {
    msg_t *ptr = reinterpret_cast<msg_t*>(trackptr);
    ptr->index = i;
    ptr->proc  = procs[i];
    *(ptr->adjData) = gid[i];
    trackptr += find_msg_size;
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After fill");
  }

  if (err) {
    throw std::logic_error("directoryComm.execute_resize error");
  }

  find_prep_clock.complete();

  Zoltan2_Directory_Clock find_rbuff_alloc_clock("  find_rbuff_alloc");

  // allocate receive buffer
  Teuchos::ArrayRCP<char> rbuff;
  if(nrec > 0) {
    rbuff = Teuchos::arcp(new char[nrec * find_msg_size],
      0, nrec * find_msg_size, true);
  }

  find_rbuff_alloc_clock.complete();

  Zoltan2_Directory_Clock find_comm_clock("  find_comm");

  const int nbytes = find_msg_size; // just getting length not full vector data

#ifdef CONVERT_DIRECTORY_RELIC
  err = Zoltan_Comm_Do (plan, ZOLTAN2_DD_FIND_MSG_TAG+1, sbuff.getRawPtr(),
    nbytes, rbuff.getRawPtr());
#else
  err = directoryComm.do_forward(ZOLTAN2_DD_FIND_MSG_TAG+1, sbuff,
    nbytes, rbuff);
#endif

  if (err) {
    throw std::logic_error("Zoltan2_Directory::find() error");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Do");
  }

  find_comm_clock.complete();

  Zoltan2_Directory_Clock find_build_sizes_clock("  find_build_sizes");

  // get find messages directed to me and determine total recv size
  // and a list of individual sizes (rmsg_sizes_resized)
  Teuchos::ArrayRCP<int> rmsg_sizes_resized;
  if(nrec > 0) {
    rmsg_sizes_resized = Teuchos::arcp(new int[nrec], 0, nrec, true);
  }
  Teuchos::ArrayRCP<int>::size_type sum_rmsg_sizes_resized = 0;

  char *track_ptr = rbuff.getRawPtr();
  for (int i = 0; i < nrec; i++) {
    msg_t *msg = reinterpret_cast<msg_t*>(track_ptr);
    track_ptr += find_msg_size;
    size_t find_rmsg_size_resized = get_local_find_msg_size(msg->adjData,
      throw_if_missing);
    rmsg_sizes_resized[i] = find_rmsg_size_resized;
    sum_rmsg_sizes_resized += find_rmsg_size_resized;
  }

  find_build_sizes_clock.complete();

  Zoltan2_Directory_Clock find_cpy_clock("  find_cpy");

  Teuchos::ArrayRCP<char>::size_type track_offset_resized = 0;

  Teuchos::ArrayRCP<char> rbuff_resized_build;

  if(is_Zoltan2_Directory_Vector()) {

    if(sum_rmsg_sizes_resized > 0) {
      rbuff_resized_build = Teuchos::arcp(new char[sum_rmsg_sizes_resized], 0, sum_rmsg_sizes_resized, true);
    }

    track_ptr = rbuff.getRawPtr();
    char * track_ptr_resized = rbuff_resized_build.getRawPtr();
    for (int i = 0; i < nrec; i++) {
      memcpy(track_ptr_resized, track_ptr, find_msg_size);
      track_ptr += find_msg_size;
      track_ptr_resized += rmsg_sizes_resized[i];
    }
  }

  // for performance, when not using variable sized we can optimize this step
  // and use the original rbuff (there is no resizing to consider)
  Teuchos::ArrayRCP<char> rbuff_resized = is_Zoltan2_Directory_Vector() ?
    rbuff_resized_build : rbuff;

  find_cpy_clock.complete();

  Zoltan2_Directory_Clock find_locals_clock("  find_locals");

  // Fill it with the true array data
  track_offset_resized = 0; // for debugging
  track_ptr = rbuff_resized.getRawPtr();
  for (int i = 0; i < nrec; i++) {
    typedef Zoltan2_DD_Find_Msg<gid_t,lid_t> find_msg_t;
    find_msg_t *ptr = reinterpret_cast<find_msg_t*>(track_ptr);
    user_t * puser = reinterpret_cast<user_t*>(
      reinterpret_cast<char*>(ptr->adjData) + max_id_size);

    // In original DD_Find_Local the first two values (gid and lid) are
    // passed as the same, we send in gid and collect lid if it's used.
    // that is why the find_msg_size is setup originally using max of
    // sizeof(gid_t) and sizeof(lid_t)
    err = find_local(ptr->adjData, (lid_t*)ptr->adjData,
      puser, &ptr->partition, &ptr->proc, throw_if_missing);

    const size_t & size_shift = rmsg_sizes_resized[i];
    track_offset_resized += size_shift;
    track_ptr += size_shift;
  }

  if(track_offset_resized != sum_rmsg_sizes_resized) {
    throw std::logic_error("Bad sum!");
  }

  find_locals_clock.complete();

  Zoltan2_Directory_Clock find_reverse_clock("  find_reverse");

  // This section is handled differently if it's variable sized array
  size_t size_scale = is_Zoltan2_Directory_Vector() ? 1 : find_msg_size;

  Teuchos::ArrayRCP<char> sbuff_resized;
  if(!is_Zoltan2_Directory_Vector()) {
    sbuff_resized = sbuff; // vector mode will set this and fill it
  }

  if(is_Zoltan2_Directory_Vector()) {
#ifdef CONVERT_DIRECTORY_RELIC
    err = Zoltan_Comm_Do_Reverse(plan, ZOLTAN2_DD_FIND_MSG_TAG+2,
      rbuff_resized.getRawPtr(), size_scale, rmsg_sizes_resized.getRawPtr(),
      sbuff_resized.getRawPtr());
#else
    err = directoryComm.do_reverse(ZOLTAN2_DD_FIND_MSG_TAG+2,
      rbuff_resized, size_scale, rmsg_sizes_resized, sbuff_resized);
#endif
  }
  else {
#ifdef CONVERT_DIRECTORY_RELIC
    err = Zoltan_Comm_Do_Reverse(plan, ZOLTAN2_DD_FIND_MSG_TAG+2,
     rbuff_resized.getRawPtr(), size_scale, NULL, sbuff_resized.getRawPtr());
#else
    err = directoryComm.do_reverse(ZOLTAN2_DD_FIND_MSG_TAG+2,
      rbuff_resized, size_scale, Teuchos::null, sbuff_resized);
#endif
  }

  if (err) {
    throw std::logic_error("Zoltan2_Directory::find() do reverse failed");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Reverse");
  }

  find_reverse_clock.complete();

  Zoltan2_Directory_Clock find_collect_reverse_clock("  find_reverse");

  // fill in user supplied lists with returned information
  track_offset_resized = 0;

  // it's not going to be in order... so don't use gid[i]
  char * trackptr_resized = sbuff_resized.getRawPtr();
  for (size_t i = 0; i < gid.size(); i++) {

    if(track_offset_resized >= sbuff_resized.size()) {
      printf( "%d has gid.size() %d track_offset_resized: %d sbuff_resized: %d\n", comm->getRank(),
        (int) gid.size(), (int) track_offset_resized, (int) sbuff_resized.size());
      throw std::logic_error("Bad buffer overflow! Internal error.");
    }

    msg_t *ptr = reinterpret_cast<msg_t*>(trackptr_resized);

    if (owner.size())
      owner[ptr->index] = ptr->proc;
    if (partition.size())
      partition[ptr->index] = ptr->partition ;
    if (lid.size()) {
      lid[ptr->index] = *(reinterpret_cast<lid_t*>(ptr->adjData));
    }

    user_t * pRead = reinterpret_cast<user_t*>(
      reinterpret_cast<char*>(ptr->adjData) + max_id_size);

    // if find_local failed proc is set to -1. Then we can leave the data
    // untouched - the default behavior is to throw but the unit tests are
    // set up to track and verify each remove id was properly taken out. To do
    // this the test overrides with an optional flag on find() and says do not
    // throw - and expects the data to remain untouched - then validates the
    // data is not changed.
    if(ptr->proc != -1  && user.size()) {
      raw_to_user(pRead, user[ptr->index]);
    }

    // don't use smsg_sizes_resized here because we don't get them back
    // in the same order
    size_t incoming_size = get_incoming_find_msg_size(ptr);
    trackptr_resized += incoming_size;
    track_offset_resized += incoming_size;
  }

  if(track_offset_resized != sbuff_resized.size()) {
    throw std::logic_error("Bad buffer sum!");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After fill return lists");
  }

  find_collect_reverse_clock.complete();

  Zoltan2_Directory_Clock find_cleanup_clock("  find_cleanup");

  int errcount = 0;
  Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_SUM, 1, &errcount, &err);
  err = (err) ? 1 : 0;

  // if at least one GID was not found, potentially notify caller of error
  if (debug_level > 0) {
    char str[100];      /* diagnostic message string */
    sprintf(str, "Processed %lu GIDs, GIDs not found: %d", gid.size(), errcount);
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);
  }

  find_cleanup_clock.complete();

#ifdef CONVERT_DIRECTORY_RELIC
  Zoltan_Comm_Destroy (&plan);
#endif

#endif

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT(comm->getRank(), yo, NULL);
  }

  return err;
}

#ifndef CONVERT_DIRECTORY_TPETRA
template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::find_local(
  gid_t* gid,                 /* incoming GID to locate (in)            */
  lid_t* lid,                 /* gid's LID (out)                        */
  user_t *user,               /* gid's user data (out)                  */
  int *partition,             /* gid's partition number (out)           */
  int *owner,                 /* gid's owner (processor number) (out)   */
  bool throw_if_missing) const
{
  const char * yo = "Zoltan2_Directory::find_local";

  /* input sanity check */
  if (owner == NULL || gid == NULL)  {
    throw std::logic_error("Zoltan2_Directory::find_local() Invalid input");
  }

  if (debug_level > 5) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

  /* compute offset into hash table to find head of linked list */
#ifdef CONVERT_DIRECTORY_KOKKOS
  // Note performance is better if we first get index, then check valid_at,
  // and use index to call value_at. Alternative is to call exists(*gid) and
  // then node_map.value_at(node_map.find(*gid))) which is slower.
  // TODO: Can this be optimized further?
  size_t node_index = node_map.find(*gid);
  if(node_map.valid_at(node_index))
  {
    Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node =
      node_map.value_at(node_index);
#else
  int index = hash_table(*gid);
  /* walk link list until end looking for matching global ID */
  for (relice_idx_t nodeidx = table[index]; nodeidx != -1;
    nodeidx = nodelist[nodeidx].next) {
    const Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node = nodelist[nodeidx];
    if (*gid == *node.gid) {
#endif
      /* matching global ID found! Return gid's information */
#ifdef CONVERT_DIRECTORY_RELIC
      if (lid) {
        lid_t *plid = reinterpret_cast<lid_t*>(node.gid + 1);
        *lid = *plid;
      }
      if (user) {
        user_t * puser = reinterpret_cast<user_t*>(
          reinterpret_cast<char*>(node.gid) + sizeof(gid_t) +
          (use_lid?sizeof(lid_t):0));
        *user = *puser;
      }
#else
      if(lid) {
        *lid = node.lid;
      }

      if (user) {
        user_to_raw(node.userData, user);
      }
#endif

      if (owner)     *owner     = node.owner;
      if (partition) *partition = node.partition;

      if (debug_level > 5) {
        ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
      }
      return 0; // success point
#ifdef CONVERT_DIRECTORY_RELIC
    }
#endif
  }

  // failure point
  if (owner != NULL)
    *owner = -1;    /* JDT Added -1 owner not found */

  if (debug_level > 5) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }

  if(throw_if_missing) {
    // TODO: This used to be linked to debug_level but wanted to discuss as it
    // seems we would want an error always if this failed.
    throw std::logic_error("find_local did not succeed");
  }

  return 0;
}
#endif

#ifdef CONVERT_DIRECTORY_RELIC
template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::allocate_node_list(
  relice_idx_t count,         /* Number of GIDs in update list  */
  float overalloc             /* Percentage to extra nodes to
                               allocate (for future dynamic
                               additions).                    */
)
{
  // Allocate node memory and initialize it to be all free.
  // Return error code if memory alloc fails.
  relice_idx_t len = count * (1. + overalloc);
  nodelistlen = len;
  if (len > 0) {
    nodelist = Teuchos::arcp<Zoltan2_Directory_Node<gid_t,lid_t,user_t>>(
      new Zoltan2_Directory_Node<gid_t,lid_t,user_t>[len], 0, len, true);
    nodedata = Teuchos::arcp<char>(
      new char[nodedata_size * len], 0, nodedata_size * len, true);
    nextfreenode = 0;
    for(relice_idx_t nodeidx = 0; nodeidx < len; ++nodeidx) {
      Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node = nodelist[nodeidx];
      node.next = nodeidx + 1;
      node.gid = reinterpret_cast<gid_t*>(&(nodedata[nodeidx*nodedata_size]));
      node.free = 1;
    }
    nodelist[len-1].next = -1;  // NULL value at end of list
  }
  return 0;
}

template <typename gid_t,typename lid_t,typename user_t>
relice_idx_t Zoltan2_Directory<gid_t,lid_t,user_t>::allocate_node()
{
  // "allocate" a node by removing it from the free list and returning
  // its NodeIdx.
  if (nextfreenode == -1) {
    // No more room for new nodes in the node list.
    // Realloc memory and set up a new freelist.
    relice_idx_t newlen = nodelistlen * 2;

    nodelist.resize(newlen);
    nodedata.resize(nodedata_size * newlen);

    // Reinitialize the gid pointers in the realloc'ed nodelist.
    for(relice_idx_t nodeidx = 0; nodeidx < newlen; ++nodeidx) {
      Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node = nodelist[nodeidx];
      node.gid = (gid_t*)(&(nodedata[nodeidx*nodedata_size]));
    }

    // Initialize free list in the newly extended memory
    nextfreenode = nodelistlen;
    for (relice_idx_t i = nodelistlen; i < newlen-1; i++) {
      nodelist[i].next = i+1;
      nodelist[i].free = 1;
    }
    nodelist[newlen-1].next = -1;
    nodelist[newlen-1].free = 1;
    nodelistlen = newlen;
  }

  relice_idx_t returnnode = nextfreenode;
  nextfreenode = nodelist[returnnode].next;
  nodelist[returnnode].next = -1;
  nodelist[returnnode].free = 0;
  return returnnode;
}

template <typename gid_t,typename lid_t,typename user_t>
void Zoltan2_Directory<gid_t,lid_t,user_t>::free_node(
  relice_idx_t freenode)
{
  // "free" a node by returning it to the free list.
  // TODO Error check:  freenode should be < nodelistlen
  nodelist[freenode].next = nextfreenode;
  nodelist[freenode].free = 1;
  nextfreenode = freenode;
}
#endif

// Print block
template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::print() const
{
  const char * yo = "Zoltan2_Directory::print";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

#ifdef CONVERT_DIRECTORY_TPETRA
  printf( "Tpetra mode does not have print() method implemented.\n" );
#else

#ifdef CONVERT_DIRECTORY_RELIC
  /* walk linked list printing each node */
  for (int i = 0; i < table_length; i++) {
    for (relice_idx_t nodeidx = table[i]; nodeidx != -1;
      nodeidx = nodelist[nodeidx].next) {
      const Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node = nodelist[nodeidx];
      printf ("ZOLTAN DD Print(%d): \tList %3d, \tGID ", comm->getRank(), i);
      printf("(");
      printf("%zu ", static_cast<size_t>(node.gid[0]));
#endif

#ifdef CONVERT_DIRECTORY_KOKKOS
  for (size_t i = 0; i < node_map.size(); i++) {
      Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node =
        node_map.value_at(i);
      printf ("ZOLTAN DD Print(%d): \tList, \tGID ", comm->getRank());
      printf("(");
      //  the issue here is that gid could be a simple type (int/long) or it
      //  could be an arbitrary structure such as:
      //  struct some_type {
      //    int x[4];
      //  }
      // Directory works for such arbitrary type but only knows sizeof,
      // not details of the internal implementation. In the testing framework
      // the length of above x array is stored so it can be printed.
      // TODO: Decide how to handle this - if we want directory to be able to
      // print nice output of the gid (and lid) then perhaps we need to require
      // that such structs define a to_string() method.
      printf( "TODO: Decide how to print gid of arbitrary type.");
#endif

      printf(") ");
      if (use_lid) {
        printf("\tLID (");
#ifdef CONVERT_DIRECTORY_RELIC
        lid_t * ptr = reinterpret_cast<lid_t*>(node.gid + 1);
        printf( "%zu ", static_cast<size_t>(*ptr));
#endif

#ifdef CONVERT_DIRECTORY_KOKKOS
        // see above note on the gid and printing
        printf( "TODO: Decide how to print lid of arbitrary type.");
#endif
        printf(") ");
      }
      printf ("\tPart %d\n", node.partition);
      printf ("\tOwner %d\n", node.owner);
#ifdef CONVERT_DIRECTORY_RELIC
    }
#endif

  }
#endif // #ifndef CONVERT_DIRECTORY_TPETRA

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }
  return 0;
}

// Stats block
template <typename gid_t,typename lid_t,typename user_t>
void Zoltan2_Directory<gid_t,lid_t,user_t>::stats() const
{
  const char    *yo = "Zoltan2_Directory::stats";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

#ifdef CONVERT_DIRECTORY_TPETRA
  printf( "Tpetra mode does not have stats() method implemented.\n" );
#else

  char str[100];      // used to build message string

#ifdef CONVERT_DIRECTORY_RELIC
  /* walk down each list in hash table to find every Node */
  int node_count = 0;     /* counts Nodes in local directory      */
  int maxlength  = 0;     /* length of longest linked list        */
  int list_count = 0;     /* number of linked lints in hash table */
  for (int i = 0; i < table_length; i++) {
    int length = 0;                    /* reset length for next count */
    if (table[i] != -1)
      list_count++;               /* count of distict linked lists */
    for (relice_idx_t nodeidx = table[i]; nodeidx != -1;
      nodeidx = nodelist[nodeidx].next) {

      if (debug_level > 6) {
        char str2[100];      // used to build message string
   //     sprintf(str, "GID %zu, Owner %d, Table Index %d.",
   //       (size_t) *nodelist[nodeidx].gid, nodelist[nodeidx].owner, i);
        ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str2);
      }

      length++;                  /* linked list length */
      node_count++;              /* count of Nodes */

      if (length > maxlength)
        maxlength = length;        /* save length of longest linked list */
    }
  }

  sprintf(str, "Hash table size %d, %d nodes on %d lists, max list length %d.",
          table_length, node_count, list_count, maxlength);
#else
  // not much to do here for equivalent to stats
  // TODO: Consider removing stats()
  sprintf(str, "Kokkos unordered map %d nodes.", node_map.size());
#endif

  ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);

#endif // #ifndef CONVERT_DIRECTORY_TPETRA

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }
}

template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::remove(
  const std::vector<gid_t>& gid)   /* Incoming list of GIDs to remove */
{
  Zoltan2_Directory_Clock clock("remove");

  const char * yo = "Zoltan2_Directory::remove";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

  int err = 0;

#ifdef CONVERT_DIRECTORY_TPETRA
  throw std::logic_error("Tpetra does not support remove yet.");
#else

  // allocate memory for processor contact list
  Teuchos::ArrayRCP<int> procs;
  Teuchos::ArrayRCP<char> sbuff;
  if(gid.size() > 0) {
    procs = Teuchos::arcp( // list of processors to contact
      new int[gid.size()], 0, gid.size(), true);
    sbuff = Teuchos::arcp( // send buffer
      new char[gid.size()*remove_msg_size], 0, gid.size()*remove_msg_size, true);
  }

  typedef Zoltan2_DD_Remove_Msg<gid_t,lid_t> msg_t;

  // for each GID, fill in contact list and then message structure
  char * trackptr = sbuff.getRawPtr();
  for (size_t i = 0; i < gid.size(); i++)  {
    procs[i] = hash_proc(gid[i]);
    msg_t *ptr = reinterpret_cast<msg_t*>(trackptr);
    ptr->owner = comm->getRank();
    *(ptr->adjData) = gid[i];
    trackptr += remove_msg_size;
  }

  // now create efficient communication plan
#ifdef CONVERT_DIRECTORY_KOKKOS
  Zoltan2_Directory_Comm directoryComm(gid.size(), procs, comm,
    ZOLTAN2_DD_UPDATE_MSG_TAG);
  int nrec = directoryComm.getNRec();
#else
  ZOLTAN_COMM_OBJ *plan  = NULL;   // for efficient MPI communication
  int nrec;       // number of receives to expect
  err = Zoltan_Comm_Create (&plan, gid.size(), procs.getRawPtr(),
    Teuchos::getRawMpiComm(*comm), ZOLTAN2_DD_REMOVE_MSG_TAG, &nrec);
#endif

  if (err) {
    throw std::logic_error("Zoltan_Comm_Create failed");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, "After Zoltan_Comm_Create");
  }

  // allocate receive buffer for nrec DD_Remove_Msg structures
  Teuchos::ArrayRCP<char> rbuff;
  if(nrec > 0) {
    rbuff = Teuchos::arcp(new char[nrec*remove_msg_size],
      0, nrec*remove_msg_size, true);   // receive buffer
  }

  // send my update messages & receive updates directed to me
#ifdef CONVERT_DIRECTORY_KOKKOS
  err = directoryComm.do_forward(ZOLTAN2_DD_UPDATE_MSG_TAG+1,
    sbuff, remove_msg_size, rbuff);
#else
  err = Zoltan_Comm_Do (plan, ZOLTAN2_DD_UPDATE_MSG_TAG+1, sbuff.getRawPtr(),
    remove_msg_size, rbuff.getRawPtr());
#endif

  if (err) {
    throw std::logic_error("Comm_Do error");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, "After Zoltan_Comm_Do");
  }

  /* for each message rec'd,  remove local directory info */
  int errcount = 0;

  // Note calling begin_erase() and end_erase() without actually erasing
  // something leads to confusing seg faults. Hence nrec>0 check.
  // Pulling begin_erase and end_erase out of the loop is important for
  // performance. Since the actual erase is fairly simple we may consider
  // eliminating remove_local and just calling here but will keep for now to
  // keep symmetry with other modes.
  if(nrec > 0) {

#ifdef CONVERT_DIRECTORY_KOKKOS
  node_map.begin_erase();
#endif

  for (int i = 0; i < nrec; i++)  {
    msg_t *ptr = reinterpret_cast<msg_t*>(&(rbuff[i*remove_msg_size]));
    err = remove_local(ptr->adjData);
    if (err == 1) { // TODO eliminate warns (1) make all errors
      ++errcount;
    }
  }

#ifdef CONVERT_DIRECTORY_KOKKOS
  node_map.end_erase();
#endif

  } // if nrec > 0

  err = 0;

  if (debug_level) {
    char str[100];      // used to build message string
    sprintf (str, "Processed %zu GIDs (%d local), %d GIDs not found",
      gid.size(), nrec, errcount);
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);
    err = (errcount) ? 1 : 0;
  }
#endif

  return err;
}

#ifndef CONVERT_DIRECTORY_TPETRA
template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::remove_local(
  gid_t* gid)                /* GID to be removed (in)  */
{
  const char * yo = "Zoltan2_Directory::remove_local";

  /* input sanity checking */
  if (gid == NULL) {
    throw std::logic_error("Invalid input argument");
  }

  if (debug_level > 5) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

#ifdef CONVERT_DIRECTORY_KOKKOS
  if(node_map.exists(*gid)) {
    node_map.erase(*gid);
    return 0;
  }
#else
  /* compute offset into hash table to find head of linked list */
  int index = hash_table(*gid);

  /* walk linked list until end looking for matching gid (key) */
  relice_idx_t prevnodeidx = -1;
  for (relice_idx_t nodeidx = table[index]; nodeidx != -1;
    nodeidx = nodelist[nodeidx].next) {
    if (*gid == *nodelist[nodeidx].gid) {
      /* found node to remove, need to preserve its next ptr */
        if (prevnodeidx != -1)
          nodelist[prevnodeidx].next = nodelist[nodeidx].next;
        else
          table[index] = nodelist[nodeidx].next;
        free_node(nodeidx);       /* now OK to delete node */

        if (debug_level > 5) {
          ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
        }
        return 0;
    }
    prevnodeidx = nodeidx;
  }
#endif

  /* We get here only if the global ID has not been found */
  if (debug_level > 5) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }

  return 1;
}
#endif // CONVERT_DIRECTORY_TPETRA

#ifdef CONVERT_DIRECTORY_RELIC
template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::equal_id(int n, gid_t* a, gid_t* b) const
{
  /*
   * Returns 1 if a == b; 0 otherwise.
   * a == b if for all i, a[i] == b[i].
   */
  for (int i = 0; i < n; i++) {
    if (a[i] != b[i]) {
      return(0);
    }
  }
  return(1);
}
#endif

template <typename gid_t,typename lid_t,typename user_t>
unsigned int Zoltan2_Directory<gid_t,lid_t,user_t>::hash_proc(const gid_t & gid) const
{
  uint32_t k;
  MurmurHash3_x86_32((void *)(&gid), sizeof(gid_t), 14, (void *)&k);
  return(k % comm->getSize());
}

#ifdef CONVERT_DIRECTORY_RELIC
template <typename gid_t,typename lid_t,typename user_t>
unsigned int Zoltan2_Directory<gid_t,lid_t,user_t>::hash_table(const gid_t& key) const
{
  uint32_t k;
  MurmurHash3_x86_32((void *)(&key), sizeof(gid_t), 14, (void *)&k);
  return(k % table_length);
}

template <typename gid_t,typename lid_t,typename user_t>
unsigned int
Zoltan2_Directory<gid_t,lid_t,user_t>::recommended_hash_size(unsigned int n)
{
    /* Prime number approximately in the middle of the range [2^x..2^(x+1)]
     * is in primes[x-1]. Every prime number stored is approximately two times
     * the previous one, so hash table size doubles every time.
     */
    unsigned int primes[] = {
    3, 7, 13, 23, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593,
    49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469,
    12582917, 25165842, 50331653, 100663319, 201326611, 402653189,
    805306457, 1610612741 } ;

    /* SRSR : err on the side of performance and choose the next largest prime
     * number. One can also choose primes[i-1] below to cut the memory by half.
     */
    int hsize = primes[29] ;
    for (int i = 0 ; i < 30 ; i++)
    {
        if (n <= primes[i])
        {
            /*hsize = (i == 0 ? n : primes[i-1]) ;*/
            hsize = primes[i] ;
            break ;
        }
    }

    return hsize ;
}
#endif

/* TODO: Currently disabled - I need to review the benefits of this and see if
   we need this in the new version. Also how do we best handle this for variable
   sized data. */
/*
template <typename gid_t,typename lid_t,typename user_t>
size_t Zoltan2_Directory<gid_t,lid_t,user_t>::align_size_t(size_t a) const
{
  #define ZOLTAN2_ALIGN_VAL 7U
  return((ZOLTAN2_ALIGN_VAL + a) & ~ZOLTAN2_ALIGN_VAL);
}
*/

} // end namespace Zoltan2

#endif // CONVERT_DIRECTORY_ORIGINAL
