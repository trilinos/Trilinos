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

#ifndef ZOLTAN2_DIRECTORY_IMPL_H_
#define ZOLTAN2_DIRECTORY_IMPL_H_

#include "Zoltan2_Directory.hpp"
#include "Zoltan2_Directory_Comm.hpp"

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

  // get the base size of the user data component
  // for simple mode, this is going to just be the sizeof(user_t)
  // for vector mode, this is sizeof(size_t) for the std::vector length
  // the actual data will be variable
  size_t user_base_size = is_Zoltan2_Directory_Vector() ? sizeof(size_t) :
    size_of_value_type();

  // set up the base size to include the gid and the lid (if used) + user size
  size_t size = sizeof(gid_t) + (use_lid?sizeof(lid_t):0) + user_base_size;

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
  node_map = src.node_map;
  return 0;
}

template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::update(
  size_t count, const gid_t * gid, const lid_t * lid,
  const user_t * user, const int * partition,
  Update_Mode update_mode)
{
  // for conveniece store but maybe this should be a construct property
  // should a directory allow mixed modes (Replace, then Aggregate... for ex)
  this->mode = update_mode;

  const char * yo = "Zoltan2_Directory::update";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN(comm->getRank(), yo, NULL);
  }

  // part of initializing the error checking process
  // for each linked list head, walk its list resetting errcheck
  if(debug_level) {
    for(size_t n = 0; n < node_map.size(); ++n) {
      node_map.value_at(n).errcheck = -1; // not possible processor
    }
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After reset errcheck");
  }

  int err = 0;

  // allocate memory for list of processors to contact
  Teuchos::ArrayRCP<int> procs;
  if(count > 0) {
    procs = Teuchos::arcp(new int[count], 0, count, true);
  }

  // set up the msg sizes for vector mode only
  Teuchos::ArrayRCP<int> msg_sizes;
  int sum_msg_size = 0;
  if(is_Zoltan2_Directory_Vector() && count > 0) {
    msg_sizes = Teuchos::arcp(new int[count], 0, count, true);
    for (size_t i = 0; i < count; i++) {
      size_t msg_size = get_update_msg_size(user[i]);
      sum_msg_size += msg_size;
      msg_sizes[i] = msg_size;
    }
  }
  else {
    sum_msg_size = update_msg_size * count; // simple case
  }

  Teuchos::ArrayRCP<char> sbuff;
  if(sum_msg_size) {
    sbuff = Teuchos::arcp(new char[sum_msg_size], 0, sum_msg_size, true);
  }

  typedef Zoltan2_DD_Update_Msg<gid_t,lid_t> msg_t;

  // build update messages
  int track_offset = 0;
  char * trackptr = sbuff.getRawPtr();
  for (size_t i = 0; i < count; i++) {
    // hash the gid to determine which proc will own it
    procs[i] = hash_proc(gid[i]);

    // this ptr is the beginning of the message which can be different
    // lengths for different gids if using Zoltan2_Directory_Vector
    msg_t *ptr = reinterpret_cast<msg_t*>(trackptr);

    // write in all my info
    ptr->lid_flag       = lid ? 1 : 0;
    ptr->user_flag      = user ? 1 : 0;
    ptr->partition_flag = partition ? 1 : 0;
    ptr->partition      = partition ? partition[i] :  -1;
    ptr->owner          = comm->getRank();

    // now deal with the gid
    gid_t * pgid = ptr->adjData;
    *pgid = gid[i];

    // optionally write the lid if we are using that
    if(use_lid) {
      if(!lid) {
        throw std::logic_error(
          "Did not pass lid values but directory was created to use them!");
      }
      lid_t * plid = reinterpret_cast<lid_t*>(ptr->adjData + 1);
      *plid = lid[i];
    }
    else {
      if(lid) {
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
    if (user) {
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

  // now create efficient communication plan

  // Kokkos mode uses the new refactored C++ communicator class
  Zoltan2_Directory_Comm directoryComm(count, procs, comm,
    ZOLTAN2_DD_UPDATE_MSG_TAG);
  int nrec = directoryComm.getNRec();

  if (err) {
    throw std::logic_error("Zoltan2_Directory::update() Comm_Create error");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Create");
  }

  int sum_recv_sizes = 0;
  if(is_Zoltan2_Directory_Vector()) {
    // Only vector mode has to use the resizing options for getting receive info
    err = directoryComm.resize(msg_sizes,
      ZOLTAN2_DD_RESIZE_MSG_TAG, &sum_recv_sizes);
  }
  else {
    sum_recv_sizes = update_msg_size * nrec;
  }

  if (err) {
    throw std::logic_error("directoryComm.execute_resize error");
  }

  // If dd has no nodes allocated (e.g., first call to DD_Update;
  // create the nodelist and freelist

  // TODO upgrade nrec as size_t and all corresponding changes...
  if(nrec && static_cast<int>(node_map.size()) < nrec) {
    // some things to consider here if we will have multiple update calls in
    // series... how will subsequent calls optimally manage this list since the
    // new gid set may be partially overlapped with the original. Currently the
    // update_local has a mechanism to rehash and increase this size if we run
    // out so skipping this call would be logically ok (but not best performance)
    rehash_node_map(nrec);
  }

  // allocate receive buffer for nrec DD_Update_Msg structures
  Teuchos::ArrayRCP<char> rbuff;
  if(sum_recv_sizes > 0) {
    rbuff = Teuchos::arcp(
      new char[sum_recv_sizes], 0, sum_recv_sizes, true);   // receive buffer
  }

  // send my update messages & receive updates directed to me
  //if resizing we send size 1 because the sizes will be built individually
  const int nbytes = is_Zoltan2_Directory_Vector() ? 1 : update_msg_size;

  err = directoryComm.do_forward(ZOLTAN2_DD_UPDATE_MSG_TAG+1,
    sbuff, nbytes, rbuff);

  if (err) {
    throw std::logic_error("Zoltan2_Directory::update() Comm_Do error");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Do");
  }

  int errcount = 0;

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

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Local update");
  }

  err = 0;

  if (debug_level) {  // overwrite error return if extra checking is on
    err = (errcount) ? 1 : 0;
  }

  if (debug_level)  {
    char str[100];      // used to build message string
    sprintf (str, "Processed %lu GIDs (%d local), %d GID errors", count,
      directoryComm.getNRec(),
      errcount);
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);
  }

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT(comm->getRank(), yo, NULL);
  }

  return err;
}

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
  size_t node_index = node_map.find(*gid);
  if(node_map.valid_at(node_index)) {
    Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node =
      node_map.value_at(node_index);

    // found match, update directory information
    if (lid) {
      node.lid = *lid;
    }
    if (user) {
      update_local_user(user, node.userData);
    }

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
  }

  // gid not found.
  // Create new Zoltan2_Directory_Node<gid_t,lid_t> and fill it in
  Zoltan2_Directory_Node<gid_t,lid_t,user_t> node; // will add to hash at end
  node.free = 0; // TODO - is this necessary - see notes in rehash_node_map
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

  node.partition = partition;
  node.owner = owner;
  node.errcheck = owner;

  if(node_map.insert(*gid, node).failed()) {
    // Need more nodes (that's the assumption here if we have failure)
    // A new update has added more to our local list and we need more capacity.
    // TODO: Decide most efficient scheme. Here we bump to at least 10 or if
    // we're already at the level, increase by 10% increments. I think this will
    // be less efficient for small scale problems, when we probably care less,
    // but more efficient as we scale up.
    size_t new_guess_size = (node_map.size() < 10) ? 10 :
      ( node_map.size() + node_map.size()/10); // adds a minimum of 1
    rehash_node_map(new_guess_size);
    if(node_map.insert(*gid, node).failed()) {
      throw std::logic_error("Hash insert failed. Mem sufficient?....");
    }
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, "Created new directory item");
  }

  if (debug_level > 5) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }

  return 0;
}

template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::find(
 size_t count,
 const gid_t * gid, /* Incoming list of GIDs to get owners proc   */
 lid_t * lid,       /* Outgoing corresponding list of LIDs        */
 user_t * user,     /* Outgoing optional corresponding user data  */
 int * partition,   /* Outgoing optional partition information    */
 int * owner,       /* Outgoing optional list of data owners      */
 bool throw_if_missing)         /* default true - throw if gid is not found.  */
{
  const char * yo = "Zoltan2_Directory::find";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN(comm->getRank(), yo, NULL);
  }

  int err = 0;

  /* allocate memory for processors to contact for directory info */
  Teuchos::ArrayRCP<int> procs;
  if(count > 0) {
    procs = Teuchos::arcp(
      new int[count], 0, count, true); // processors to contact
  }

  /* Setup procs list */
  for (size_t i = 0; i < count; i++) {
    procs[i] = hash_proc(gid[i]); // determines the owner
  }

  // create efficient communication plan
  Zoltan2_Directory_Comm directoryComm(count, procs, comm,
    ZOLTAN2_DD_FIND_MSG_TAG);
  int nrec = directoryComm.getNRec();

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Create");
  }

  if (err) {
    throw std::logic_error("Zoltan2_Directory::find() error");
  }

  Teuchos::ArrayRCP<char> sbuff;
  if(count > 0) {
    sbuff = Teuchos::arcp(new char[find_msg_size*count],
      0, find_msg_size*count, true); // send buffer
  }

  typedef Zoltan2_DD_Find_Msg<gid_t,lid_t> msg_t;

  /* for each GID, fill DD_Find_Msg buffer and contact list */
  char *trackptr = sbuff.getRawPtr();
  for (size_t i = 0; i < count; i++)  {
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

  // allocate receive buffer
  Teuchos::ArrayRCP<char> rbuff;
  if(nrec > 0) {
    rbuff = Teuchos::arcp(new char[nrec * find_msg_size],
      0, nrec * find_msg_size, true);
  }

  const int nbytes = find_msg_size; // just getting length not full vector data

  err = directoryComm.do_forward(ZOLTAN2_DD_FIND_MSG_TAG+1, sbuff,
    nbytes, rbuff);

  if (err) {
    throw std::logic_error("Zoltan2_Directory::find() error");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Do");
  }

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

  Teuchos::ArrayRCP<char>::size_type track_offset_resized = 0;

  Teuchos::ArrayRCP<char> rbuff_resized_build;

  if(is_Zoltan2_Directory_Vector()) {

    if(sum_rmsg_sizes_resized > 0) {
      rbuff_resized_build = Teuchos::arcp(new char[sum_rmsg_sizes_resized],
        0, sum_rmsg_sizes_resized, true);
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

  // This section is handled differently if it's variable sized array
  size_t size_scale = is_Zoltan2_Directory_Vector() ? 1 : find_msg_size;

  Teuchos::ArrayRCP<char> sbuff_resized;
  if(!is_Zoltan2_Directory_Vector()) {
    sbuff_resized = sbuff; // vector mode will set this and fill it
  }

  if(is_Zoltan2_Directory_Vector()) {
    err = directoryComm.do_reverse(ZOLTAN2_DD_FIND_MSG_TAG+2,
      rbuff_resized, size_scale, rmsg_sizes_resized, sbuff_resized);
  }
  else {
    err = directoryComm.do_reverse(ZOLTAN2_DD_FIND_MSG_TAG+2,
      rbuff_resized, size_scale, Teuchos::null, sbuff_resized);
  }

  if (err) {
    throw std::logic_error("Zoltan2_Directory::find() do reverse failed");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Reverse");
  }

  // fill in user supplied lists with returned information
  track_offset_resized = 0;

  // it's not going to be in order... so don't use gid[i]
  char * trackptr_resized = sbuff_resized.getRawPtr();
  for (size_t i = 0; i < count; i++) {

    if(track_offset_resized >= sbuff_resized.size()) {
      printf(
        "%d has gid.size() %d track_offset_resized: %d sbuff_resized: %d\n",
        comm->getRank(),
        (int) count, (int) track_offset_resized, (int) sbuff_resized.size());
      throw std::logic_error("Bad buffer overflow! Internal error.");
    }

    msg_t *ptr = reinterpret_cast<msg_t*>(trackptr_resized);

    if (owner)
      owner[ptr->index] = ptr->proc;
    if (partition)
      partition[ptr->index] = ptr->partition ;
    if (lid) {
      memcpy (&lid[ptr->index], ptr->adjData, sizeof(lid_t));
    }

    user_t * pRead = reinterpret_cast<user_t*>(
      reinterpret_cast<char*>(ptr->adjData) + max_id_size);

    // if find_local failed proc is set to -1. Then we can leave the data
    // untouched - the default behavior is to throw but the unit tests are
    // set up to track and verify each remove id was properly taken out. To do
    // this the test overrides with an optional flag on find() and says do not
    // throw - and expects the data to remain untouched - then validates the
    // data is not changed.
    if(ptr->proc != -1  && user) {
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

  int errcount = 0;
  Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_SUM, 1, &errcount, &err);
  err = (err) ? 1 : 0;

  // if at least one GID was not found, potentially notify caller of error
  if (debug_level > 0) {
    char str[100];      /* diagnostic message string */
    sprintf(str, "Processed %lu GIDs, GIDs not found: %d", count, errcount);
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);
  }

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT(comm->getRank(), yo, NULL);
  }

  return err;
}

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
  // Note performance is better if we first get index, then check valid_at,
  // and use index to call value_at. Alternative is to call exists(*gid) and
  // then node_map.value_at(node_map.find(*gid))) which is slower.
  // TODO: Can this be optimized further?
  size_t node_index = node_map.find(*gid);
  if(node_map.valid_at(node_index))
  {
    const Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node =
      node_map.value_at(node_index);
    /* matching global ID found! Return gid's information */
    if(lid) {
      *lid = node.lid;
    }

    if (user) {
      user_to_raw(node.userData, user);
    }

    if (owner)     *owner     = node.owner;
    if (partition) *partition = node.partition;

    if (debug_level > 5) {
      ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
    }
    return 0; // success point
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

// Print block
template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::print() const
{
  const char * yo = "Zoltan2_Directory::print";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

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
    printf(") ");
    if (use_lid) {
      printf("\tLID (");
      // see above note on the gid and printing
      printf( "TODO: Decide how to print lid of arbitrary type.");
      printf(") ");
    }
    printf ("\tPart %d\n", node.partition);
    printf ("\tOwner %d\n", node.owner);
  }

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

  char str[100];      // used to build message string
  // not much to do here for equivalent to stats
  // TODO: Consider removing stats()
  sprintf(str, "Kokkos unordered map %d nodes.", node_map.size());
  ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);
  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }
}

template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::remove(
  size_t count,
  const gid_t * gid)   /* Incoming list of GIDs to remove */
{
  const char * yo = "Zoltan2_Directory::remove";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

  int err = 0;

  // allocate memory for processor contact list
  Teuchos::ArrayRCP<int> procs;
  Teuchos::ArrayRCP<char> sbuff;
  if(count > 0) {
    procs = Teuchos::arcp( // list of processors to contact
      new int[count], 0, count, true);
    sbuff = Teuchos::arcp( // send buffer
      new char[count*remove_msg_size], 0, count*remove_msg_size, true);
  }

  typedef Zoltan2_DD_Remove_Msg<gid_t,lid_t> msg_t;

  // for each GID, fill in contact list and then message structure
  char * trackptr = sbuff.getRawPtr();
  for (size_t i = 0; i < count; i++)  {
    procs[i] = hash_proc(gid[i]);
    msg_t *ptr = reinterpret_cast<msg_t*>(trackptr);
    ptr->owner = comm->getRank();
    *(ptr->adjData) = gid[i];
    trackptr += remove_msg_size;
  }

  // now create efficient communication plan
  Zoltan2_Directory_Comm directoryComm(count, procs, comm,
    ZOLTAN2_DD_UPDATE_MSG_TAG);
  int nrec = directoryComm.getNRec();

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
  err = directoryComm.do_forward(ZOLTAN2_DD_UPDATE_MSG_TAG+1,
    sbuff, remove_msg_size, rbuff);

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
    node_map.begin_erase();
    for (int i = 0; i < nrec; i++)  {
      msg_t *ptr = reinterpret_cast<msg_t*>(&(rbuff[i*remove_msg_size]));
      err = remove_local(ptr->adjData);
      if (err == 1) { // TODO eliminate warns (1) make all errors
        ++errcount;
      }
    }
    node_map.end_erase();
  } // if nrec > 0

  err = 0;

  if (debug_level) {
    char str[100];      // used to build message string
    sprintf (str, "Processed %zu GIDs (%d local), %d GIDs not found",
      count, nrec, errcount);
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);
    err = (errcount) ? 1 : 0;
  }

  return err;
}

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

  if(node_map.exists(*gid)) {
    node_map.erase(*gid);
    return 0;
  }

  /* We get here only if the global ID has not been found */
  if (debug_level > 5) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }

  return 1;
}

template <typename gid_t,typename lid_t,typename user_t>
unsigned int Zoltan2_Directory<gid_t,lid_t,user_t>::hash_proc(
  const gid_t & gid) const
{
  uint32_t k;

  // Copied from murmurc.3
  // TODO: Will be removed as Kokkos form develops
  #define ZOLTAN2_ROTL32(x,r) (uint32_t) \
    (((uint32_t)(x) << (int8_t)(r)) | ((uint32_t)(x) >> (32 - (int8_t)(r))))

  const void * key = (void *)(&gid);
  int len = sizeof(gid_t);
  uint32_t seed = 14;
  void * out = (void *)&k;

  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 4;
  int i;
  uint32_t h1 = seed;
  uint32_t c1 = 0xcc9e2d51;
  uint32_t c2 = 0x1b873593;
  const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);
  for(i = -nblocks; i; i++)
  {
    uint32_t k1 = blocks[i];
    k1 *= c1;
    k1 = ZOLTAN2_ROTL32(k1,15);
    k1 *= c2;
    h1 ^= k1;
    h1 = ZOLTAN2_ROTL32(h1,13);
    h1 = h1*5+0xe6546b64;
  }
  const uint8_t * tail = (const uint8_t*)(data + nblocks*4);
  uint32_t k1 = 0;
  switch(len & 3)
  {
    case 3: k1 ^= tail[2] << 16;
    case 2: k1 ^= tail[1] << 8;
    case 1: k1 ^= tail[0];
            k1 *= c1; k1 = ZOLTAN2_ROTL32(k1,15); k1 *= c2; h1 ^= k1;
  };
  h1 ^= len;
  h1 ^= h1 >> 16;
  h1 *= 0x85ebca6b;
  h1 ^= h1 >> 13;
  h1 *= 0xc2b2ae35;
  h1 ^= h1 >> 16;
  *(uint32_t*)out = h1;

  return(k % comm->getSize());
}


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

#endif // ZOLTAN2_DIRECTORY_IMPL_H_
