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

#include "Zoltan2_Directory.hpp"

#ifndef CONVERT_DIRECTORY_ORIGINAL

#include <stdlib.h>
#include <stdexcept>

#ifdef CONVERT_DIRECTORY_RELIC
  #include <comm.h>
#endif

#ifdef CONVERT_DIRECTORY_KOKKOS
  #include "Zoltan2_Directory_Comm.hpp"
#endif

// Supplies the current hash - rolled over from zoltan
#include "murmur3.c"

#ifdef CONVERT_DIRECTORY_TPETRA
  #include <unordered_map>
#endif

namespace Zoltan2 {

// To refactor
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

/* Tags for MPI communications.  These need unique values. Arbitrary */
#define ZOLTAN2_DD_FIND_MSG_TAG     29137  /* needs 3 consecutive values */
#define ZOLTAN2_DD_UPDATE_MSG_TAG   29140  /* needs 2 consecutive values */
#define ZOLTAN2_DD_REMOVE_MSG_TAG   29142  /* needs 2 consecutive values */
#define ZOLTAN2_DD_RESIZE_MSG_TAG   29150  /*  */




template <typename gid_t,typename lid_t,typename user_t>
Zoltan2_Directory<gid_t,lid_t,user_t>::Zoltan2_Directory(
  const Zoltan2_Directory & src) : comm(src.comm), use_lid(src.use_lid),
#ifdef CONVERT_DIRECTORY_RELIC
    table_length(src.table_length),
#endif
    debug_level(src.debug_level),contiguous(src.contiguous)
{
  allocate();
  copy(src);
}

template <typename gid_t,typename lid_t,typename user_t>
Zoltan2_Directory<gid_t,lid_t,user_t>::~Zoltan2_Directory()
{
  release();
}

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
  Teuchos::reduceAll<int,size_t>(*comm, Teuchos::REDUCE_MAX, 3, array, max_array);
  Teuchos::reduceAll<int,size_t>(*comm, Teuchos::REDUCE_MIN, 3, array, min_array);
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

  // Note that in the original Zoltan implementaion memory was allocated for the
  // main struct with additional memory for the table - so this step mimics
  // that part of the code where we allocate memory for the table
  table.resize(table_length);
  for(size_t i = 0; i < table.size(); ++i) {
    table[i] = -1;
  }

  nodelistlen = 0;
  nextfreenode = -1;

#endif

  /* frequently used dynamic allocation computed sizes */
  // for user_t int for Zoltan2_Directory_Single, size_of_value_type() is just int
  // for user_t std::vector<int> for Zoltan2_Directory_Vector, size_of_value_type() is int
  size_t size = sizeof(gid_t) + (use_lid?sizeof(lid_t):0) + size_of_value_type();

#ifdef CONVERT_DIRECTORY_RELIC
  nodedata_size = size;
#endif

  update_msg_size = size + sizeof(Zoltan2_DD_Update_Msg<gid_t,lid_t>);

  size = sizeof(gid_t);
  remove_msg_size = size + sizeof(Zoltan2_DD_Remove_Msg<gid_t,lid_t>);

  size = sizeof(lid_t) + size_of_value_type();
  find_msg_size = size + sizeof(Zoltan2_DD_Find_Msg<gid_t,lid_t>);

  /* force alignment */
#ifdef CONVERT_DIRECTORY_RELIC
//  nodedata_size   = align_size_t(nodedata_size);
#endif

//  update_msg_size = align_size_t(update_msg_size);
//  remove_msg_size = align_size_t(remove_msg_size);
//  find_msg_size   = align_size_t(find_msg_size);

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }
}

// Copy Block
template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::copy(
  const Zoltan2_Directory<gid_t,lid_t,user_t> & src) {
#ifdef CONVERT_DIRECTORY_RELIC
  throw std::logic_error("copy not implemented.");
  /*
  if (nodelistlen) {
    nodelist = (Zoltan2_Directory_Node<gid_t,lid_t> *)
      ZOLTAN2_MALLOC(nodelistlen * sizeof(Zoltan2_Directory_Node<gid_t,lid_t>));
    memcpy(nodelist, src.nodelist,
      nodelistlen * sizeof(Zoltan2_Directory_Node<gid_t,lid_t>));

    nodedata = (char *) ZOLTAN2_MALLOC(nodelistlen * nodedata_size);
    memcpy(nodedata, src.nodedata, nodelistlen * nodedata_size);

    for (relice_idx_t i = 0; i < nodelistlen; i++) {
      nodelist[i].gid = static_cast<gid_t*>(nodedata + i*nodedata_size);
    }
  }
  */
#endif

#ifdef CONVERT_DIRECTORY_KOKKOS
  throw std::logic_error("copy not implemented.");
#endif

#ifdef CONVERT_DIRECTORY_TPETRA
  throw std::logic_error("copy not implemented.");
#endif

  return 0;
}

template <typename gid_t,typename lid_t,typename user_t>
void Zoltan2_Directory<gid_t,lid_t,user_t>::release() {
  const char * yo = "Zoltan2_Directory::release";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }
}

template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::update(
  const std::vector<gid_t>& gid, const std::vector<lid_t>& lid,
  const std::vector<user_t>& user, const std::vector<int>& partition,
  Update_Mode mode)
{
  Zoltan2_Directory_Clock clock("update");

  const char * yo = "Zoltan2_Directory::update";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN(comm->getRank(), yo, NULL);
  }

  Zoltan2_Directory_Clock update_setup("update_setup", 1);

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
  update_setup.complete();

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After reset errcheck");
  }

  int err = 0;

  Zoltan2_Directory_Clock guess_hash_size("update_hash_guess", 1);

#ifdef CONVERT_DIRECTORY_KOKKOS
  // TODO decide how to best allocate initial size
  // Ideally something just a little bit bigger than we actually need.
  // This gets more complicated for things like aggregate mode where
  // various procs may all be sending the same id and total count won't
  // be as useful. When we update_local, the node map will check for a failed
  // insert and increase this size, so the decision here is only an issue
  // of performance.
  size_t globalCount = 0;
  Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, gid.size(),
    Teuchos::outArg(globalCount));
  if(node_map.capacity() < globalCount) {
    size_t estimate_local_count = globalCount / comm->getSize();
    estimate_local_count += estimate_local_count / 10; // add 10% buffer
    if(node_map.capacity() < estimate_local_count) {
      // skipping this line would be ok... it would just rehash when it gets
      // the update_local events though efficiency would be less. In a current
      // test of 10 million total nodes, I saw about 2x increase in total time
      // if this line was removed. Also this 'guess' is having a lot of impact
      // on the overall performance, more than I originally realized.
      node_map.rehash(estimate_local_count);
    }
  }
#endif

  guess_hash_size.complete();

#ifdef CONVERT_DIRECTORY_TPETRA
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

  if (contiguous) {
    // For contigous ids, can use Tpetra default Map
    gid_t gnUnique = idMap->getMaxAllGlobalIndex() - indexBase + 1;
    oto_idMap = Teuchos::rcp(new map_t(gnUnique, indexBase, comm));
  }
  else {
    // Since ids may not be contiguous, cannot use default Tpetra Map
    oto_idMap = Tpetra::createOneToOne(idMap);
  }

  oto_idMap = Tpetra::createOneToOne(idMap);

  oto_idVec = Teuchos::rcp(new vector_t(oto_idMap));

  // Create an exporter between the two maps
  rcp_export_t idExporter = Teuchos::rcp(new export_t(idMap, oto_idMap));

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

  // allocate memory for list of processors to contact
  std::vector<int> procs(gid.size());

  std::vector<int> msg_sizes(gid.size());
  int sum_msg_size = 0;
  for (size_t i = 0; i < gid.size(); i++) {
    size_t msg_size = get_update_msg_size(user[i]);
    sum_msg_size += msg_size;
    msg_sizes[i] = msg_size;
  }
  std::vector<char> sbuff(sum_msg_size);

  typedef Zoltan2_DD_Update_Msg<gid_t,lid_t> msg_t;

  Zoltan2_Directory_Clock update_build_raw("update_build_raw", 1);

  int track_offset = 0;
  for (size_t i = 0; i < gid.size(); i++) {
    procs[i] = hash_proc(gid[i]);

    msg_t *ptr = (msg_t*)(&(sbuff[track_offset]));

    ptr->lid_flag       = (lid.size())  ? 1 : 0;
    ptr->user_flag      = (user.size()) ? 1 : 0;
    ptr->partition_flag = (partition.size()) ? 1 : 0;
    ptr->partition      = (partition.size()) ? partition[i] :  -1;
    ptr->owner          = comm->getRank();

    gid_t * pgid = ptr->adjData;

    *pgid = gid[i];

    // TODO: Fix cast
    lid_t * plid = (lid_t*)(ptr->adjData + 1);
    if (lid.size()) {
      if(!use_lid) {
        throw std::logic_error(
          "Passed lid values but directory was created not to use them!");
      }
      *plid = lid[i];
    }
    else {
      *plid = 0;
    }

    // TODO: Fix cast
    user_t * puser = (user_t*)((char*)ptr->adjData + sizeof(gid_t) +
      (use_lid?sizeof(lid_t):0));

    if (user.size()) {
      user_to_raw(user[i], puser);
    }
    else {
      *puser = user_t(); // create a null version... how to handle
    }

    track_offset += get_update_msg_size(user[i]);
  }

  if(track_offset != sum_msg_size) {
    throw std::logic_error("Bad summing!");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After fill contact list");
  }

  update_build_raw.complete();

  // now create efficient communication plan

  Zoltan2_Directory_Clock update_build_plan("update_build_plan", 1);

#ifdef CONVERT_DIRECTORY_KOKKOS
  Zoltan2_Directory_Comm directoryComm(gid.size(), procs, comm,
    ZOLTAN2_DD_UPDATE_MSG_TAG);
  int nrec = directoryComm.getNRec();
#else
  ZOLTAN_COMM_OBJ *plan  = NULL;   // for efficient MPI communication
  int nrec = 0;       // number of receives to expect
  err = Zoltan_Comm_Create (&plan, gid.size(), &(procs[0]),
    Teuchos::getRawMpiComm(*comm), ZOLTAN2_DD_UPDATE_MSG_TAG, &nrec);
#endif

  if (err) {
    throw std::logic_error("Zoltan2_Directory::update() Comm_Create error");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Create");
  }

  update_build_plan.complete();

  Zoltan2_Directory_Clock update_resize("update_resize", 1);

  int sum_recv_sizes;
#ifdef CONVERT_DIRECTORY_KOKKOS
  err = directoryComm.resize(msg_sizes,
   ZOLTAN2_DD_RESIZE_MSG_TAG, &sum_recv_sizes);
#else
  err = Zoltan_Comm_Resize(plan, &(msg_sizes[0]),
    ZOLTAN2_DD_RESIZE_MSG_TAG, &sum_recv_sizes);
#endif
  if (err) {
    throw std::logic_error("directoryComm.execute_resize error");
  }

  update_resize.complete();

  // If dd has no nodes allocated (e.g., first call to DD_Update;
  // create the nodelist and freelist
#ifdef CONVERT_DIRECTORY_RELIC
  if (nrec && nodelistlen == 0) {
    allocate_node_list((relice_idx_t) nrec, 0.);
  }
#endif

  // allocate receive buffer for nrec DD_Update_Msg structures
  std::vector<char> rbuff(sum_recv_sizes);   // receive buffer

  // send my update messages & receive updates directed to me
  const int nbytes = 1;

  Zoltan2_Directory_Clock update_forward("update_forward", 1);

#ifdef CONVERT_DIRECTORY_KOKKOS
  err = directoryComm.do_forward(ZOLTAN2_DD_UPDATE_MSG_TAG+1,
    sbuff, nbytes, rbuff);
#else
  err = Zoltan_Comm_Do(plan, ZOLTAN2_DD_UPDATE_MSG_TAG+1, &(sbuff[0]),
    nbytes, &(rbuff[0]));
#endif

  update_forward.complete();

  if (err) {
    throw std::logic_error("Zoltan2_Directory::update() Comm_Do error");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Do");
  }

  Zoltan2_Directory_Clock update_locals("update_locals", 1);

  // for each message rec'd, update local directory information
  int errcount = 0;
  track_offset = 0;
  for (int i = 0; i < nrec; i++) {
    msg_t *ptr = (msg_t*)(&(rbuff[track_offset]));

    // TODO: Fix casts
    user_t * puser = (ptr->user_flag) ?
      (user_t*)((char*)ptr->adjData + sizeof(gid_t) + (use_lid?sizeof(lid_t):0))
      : NULL;

    err = update_local(ptr->adjData,
      (ptr->lid_flag) ? (lid_t*)(ptr->adjData + 1) : NULL,
      puser,
      (ptr->partition_flag) ? (ptr->partition) : -1,  // illegal partition
      ptr->owner, mode);

    if (err)
      ++errcount;

    track_offset += get_update_msg_size(puser);
  }

  update_locals.complete();

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
 int owner,                  /* gid's current owner (proc number) (in)    */
 Update_Mode mode)
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
  if(node_map.exists(*gid)) {
    Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node =
      node_map.value_at(node_map.find(*gid));
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
        update_local_user(user, node.userData, mode);
      }
#endif

      node.owner = owner;
      if (partition != -1) {
        node.partition = partition;
      }

      // Response to multiple updates to a gid in 1 update cycle
      if (debug_level > 0 && node.errcheck != owner) {
        throw std::logic_error(
          "Zoltan2_Directory::update_local() Multiply defined GID");
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


  // gid not found. Create new Zoltan2_Directory_Node<gid_t,lid_t> and fill it in
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

  user_t * puser = (user_t*)((char*)node.gid + sizeof(gid_t)
    + (use_lid?sizeof(lid_t):0));
  if (user) {
    *puser = *user;
  }
  else {
    *puser = user_t();
  }
#else
  node.lid = lid ? (*lid) : 0;

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

  auto result = node_map.insert(*gid, node); // add to the map
  if(result.failed()) {
    // need more nodes... our initial guess didn't cut it
    // TODO: Decide most efficient scheme. Here we bump to at least 10 or if
    // we're already at the level, increase by 10% increments. I think this will
    // less efficient for small scale problems, when we probably care less,
    // but more efficient as we scale up.
    size_t new_guess_size = (node_map.size() < 10) ? 10 :
      ( node_map.size() + node_map.size()/10); // adds a minimum of 1
    node_map.rehash(new_guess_size);
    result = node_map.insert(*gid, node); // add to the map again
    if(result.failed()) {
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
  std::vector<int> procs(gid.size());  // list of processors to contact

  /* Setup procs list */
  for (size_t i = 0; i < gid.size(); i++) {
    procs[i] = hash_proc(gid[i]);
  }

  // create efficient communication plan
#ifdef CONVERT_DIRECTORY_RELIC
  ZOLTAN_COMM_OBJ *plan  = NULL;     // efficient MPI communication
  int nrec;                          // number of messages to receive
  err = Zoltan_Comm_Create (&plan, gid.size(), &(procs[0]),
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

  std::vector<char> sbuff(find_msg_size*gid.size()); // send buffer

  typedef Zoltan2_DD_Find_Msg<gid_t,lid_t> msg_t;

  /* for each GID, fill DD_Find_Msg buffer and contact list */
  int track_offset = 0;

  for (size_t i = 0; i < gid.size(); i++)  {
    msg_t *ptr = (msg_t*)(&(sbuff[track_offset]));
    ptr->index = i;
    ptr->proc  = procs[i];
    *(ptr->adjData) = gid[i];
    track_offset += find_msg_size;
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After fill");
  }

  if (err) {
    throw std::logic_error("directoryComm.execute_resize error");
  }

  // allocate receive buffer
  std::vector<char> rbuff(nrec * find_msg_size);     // receive buffer

  const int nbytes = find_msg_size; // just getting length not full vector data

#ifdef CONVERT_DIRECTORY_RELIC
  err = Zoltan_Comm_Do (plan, ZOLTAN2_DD_FIND_MSG_TAG+1, &(sbuff[0]),
    nbytes, &(rbuff[0]));
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

  // get find messages directed to me and determine total recv size
  // and a list of individual sizes (rmsg_sizes_resized)
  track_offset = 0;
  std::vector<int> rmsg_sizes_resized(nrec);
  size_t sum_rmsg_sizes_resized = 0;

  for (int i = 0; i < nrec; i++) {
    typedef Zoltan2_DD_Find_Msg<gid_t,lid_t> find_msg_t;
    find_msg_t *msg = (find_msg_t*)(&(rbuff[track_offset]));
    track_offset += find_msg_size;
    size_t find_rmsg_size_resized = get_local_find_msg_size(msg->adjData,
      throw_if_missing);
    rmsg_sizes_resized[i] = find_rmsg_size_resized;
    sum_rmsg_sizes_resized += find_rmsg_size_resized;
  }

  size_t track_offset_resized = 0;

  std::vector<char> rbuff_resized(sum_rmsg_sizes_resized);

  track_offset = 0;
  for (int i = 0; i < nrec; i++) {
    // TODO: Fix cast
    char *ptr = (char*)(&(rbuff[track_offset]));
    char *ptr_resized = (char*)(&(rbuff_resized[track_offset_resized]));
    memcpy(ptr_resized, ptr, find_msg_size);
    track_offset += find_msg_size;
    track_offset_resized += rmsg_sizes_resized[i];
  }

  // Fill it with the true array data
  track_offset_resized = 0;
  for (int i = 0; i < nrec; i++) {
    // TODO: Fix cast
    typedef Zoltan2_DD_Find_Msg<gid_t,lid_t> find_msg_t;
    find_msg_t *ptr = (find_msg_t*)(&(rbuff_resized[track_offset_resized]));
    // TODO: Fix casts
    user_t * puser = (user_t*)(ptr->adjData + 1);

    // In original DD_Find_Local the first two values (gid and lid) are
    // passed as the same, but not sure that is can mean - currently I've
    // got the lid just set NULL and lid is not working
    err = find_local(ptr->adjData, NULL,
      puser, &ptr->partition, &ptr->proc, throw_if_missing);
    track_offset_resized += rmsg_sizes_resized[i];
  }

  if(track_offset_resized != sum_rmsg_sizes_resized) {
    throw std::logic_error("Bad sum!");
  }

  // This section is handled differently if it's variable sized array
  size_t size_scale = is_Zoltan2_Directory_Vector() ? 1 : find_msg_size;

  std::vector<char> sbuff_resized;
  if(!is_Zoltan2_Directory_Vector()) {
    sbuff_resized = sbuff; // vector mode will set this and fill it
  }

  if(is_Zoltan2_Directory_Vector()) {
#ifdef CONVERT_DIRECTORY_RELIC
    err = Zoltan_Comm_Do_Reverse(plan, ZOLTAN2_DD_FIND_MSG_TAG+2,
      &(rbuff_resized[0]), size_scale, &(rmsg_sizes_resized[0]), &(sbuff_resized[0]));
#else
    err = directoryComm.do_reverse(ZOLTAN2_DD_FIND_MSG_TAG+2,
      rbuff_resized, size_scale, rmsg_sizes_resized, sbuff_resized);
#endif
  }
  else {
#ifdef CONVERT_DIRECTORY_RELIC
    err = Zoltan_Comm_Do_Reverse(plan, ZOLTAN2_DD_FIND_MSG_TAG+2,
      &(rbuff_resized[0]), size_scale, NULL, &(sbuff_resized[0]));
#else
    err = directoryComm.do_reverse(ZOLTAN2_DD_FIND_MSG_TAG+2,
      rbuff_resized, size_scale, std::vector<int>(), sbuff_resized);
#endif
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
  for (size_t i = 0; i < gid.size(); i++) {

    if(track_offset_resized >= sbuff_resized.size()) {
    std::cout << "sbuff_resized.size() is now: " << sbuff_resized.size() << std::endl;
    std::cout << "track_offset_resized is now: " << track_offset_resized << std::endl;
      throw std::logic_error("Bad buffer overflow! Internal error.");
    }

    // TODO: Fix cast
    msg_t *ptr = (msg_t*)(&(sbuff_resized[track_offset_resized]));

    if (owner.size())
      owner[ptr->index] = ptr->proc;
    if (partition.size())
      partition[ptr->index] = ptr->partition ;
    if (lid.size()) {
      // TODO - need to redo/fix the lid handling and make it correct
      lid[ptr->index] = 0; // *((lid_t*)ptr->adjData);
    }

    user_t * pRead = (user_t*)(ptr->adjData+1);

    // if find_local failed proc is set to -1. Then we can leave the data
    // untouched - the default behavior is to throw but the unit tests are
    // set up to track and verify each remove id was properly taken out. To do
    // this the test overrides with an optional flag on find() and says do not
    // throw - and expects the data to remain untouched - then validates the
    // data is not changed.
    if(ptr->proc != -1) {
      raw_to_user(pRead, user[ptr->index]);
    }

    // don't use smsg_sizes_resized here because we don't get them back
    // in the same order
    track_offset_resized += get_incoming_find_msg_size(ptr);
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
    sprintf(str, "Processed %lu GIDs, GIDs not found: %d", gid.size(), errcount);
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);
  }

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
  if(node_map.exists(*gid)) {
    Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node =
      node_map.value_at(node_map.find(*gid));
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
          (char*)node.gid + sizeof(gid_t) + (use_lid?sizeof(lid_t):0));
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

      return 0;
#ifdef CONVERT_DIRECTORY_RELIC
    }
#endif
  }

  if (owner != NULL)
    *owner = -1;    /* JDT Added -1 owner not found */

  if (debug_level > 5) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }

  if (debug_level > 0)  {
    throw std::logic_error("GID not found");
  }

  if(throw_if_missing) {
    throw std::logic_error("find_local did not succeed");
  }

  return 0;
}
#endif

#ifdef CONVERT_DIRECTORY_RELIC
template <typename gid_t,typename lid_t,typename user_t>
int Zoltan2_Directory<gid_t,lid_t,user_t>::allocate_node_list(
  relice_idx_t count,         /* Number of GIDs in update list  */
  float overalloc           /* Percentage to extra nodes to
                               allocate (for future dynamic
                               additions).                    */
)
{
  // Allocate node memory and initialize it to be all free.
  // Return error code if memory alloc fails.
  relice_idx_t len = count * (1. + overalloc);
  nodelistlen = len;
  if (len > 0) {
    nodelist = std::vector<Zoltan2_Directory_Node<gid_t,lid_t,user_t>>(len);
    nodedata = std::vector<char>(nodedata_size * len);
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
  throw std::logic_error("UNTESTED CHECKPOINT"); // needs unit testing

  const char * yo = "Zoltan2_Directory::print";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

#ifdef CONVERT_DIRECTORY_TPETRA
  throw std::logic_error("Tpetra not iplemented");
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
  for (int i = 0; i < node_map.size(); i++) {
      Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node =
        node_map.value_at(node_map.find(i));
      printf ("ZOLTAN DD Print(%d): \tList, \tGID ", comm->getRank());
#endif

      printf(") ");
      if (use_lid) {
        printf("\tLID (");
#ifdef CONVERT_DIRECTORY_RELIC
        lid_t * ptr = reinterpret_cast<lid_t*>(node.gid + 1);
        printf( "%zu ", static_cast<size_t>(*ptr));
#endif

#ifdef CONVERT_DIRECTORY_KOKKOS
         printf( "%zu ", static_cast<size_t>(node.lid));
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

// operator=
template <typename gid_t,typename lid_t,typename user_t>
Zoltan2_Directory<gid_t,lid_t,user_t> & Zoltan2_Directory<gid_t,lid_t,user_t>::
  operator=(const Zoltan2_Directory<gid_t,lid_t,user_t> &src)
{
  throw std::logic_error("UNTESTED CHECKPOINT"); // needs unit testing

  release();

  use_lid = src.use_lid;
  debug_level = src.debug_level;

#ifdef CONVERT_DIRECTORY_RELIC
  table_length = src.table_length;
#endif

  allocate();
  copy(src);
  return *this;
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
  std::vector<int> procs(gid.size());   // list of processors to contact

  std::vector<char> sbuff(gid.size()*remove_msg_size);   // send buffer

  typedef Zoltan2_DD_Remove_Msg<gid_t,lid_t> msg_t;

  // for each GID, fill in contact list and then message structure
  for (size_t i = 0; i < gid.size(); i++)  {
    procs[i] = hash_proc(gid[i]);

    // TODO: Fix cast
    msg_t *ptr = (msg_t*)(&(sbuff[i*remove_msg_size]));
    ptr->owner = comm->getRank();
    *(ptr->adjData) = gid[i];
  }

  // now create efficient communication plan
#ifdef CONVERT_DIRECTORY_KOKKOS
  Zoltan2_Directory_Comm directoryComm(gid.size(), procs, comm,
    ZOLTAN2_DD_UPDATE_MSG_TAG);
  int nrec = directoryComm.getNRec();
#else
  ZOLTAN_COMM_OBJ *plan  = NULL;   // for efficient MPI communication
  int nrec;       // number of receives to expect
  err = Zoltan_Comm_Create (&plan, gid.size(), &(procs[0]),
    Teuchos::getRawMpiComm(*comm), ZOLTAN2_DD_REMOVE_MSG_TAG, &nrec);
#endif

  if (err) {
    throw std::logic_error("Zoltan_Comm_Create failed");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO (comm->getRank(), yo, "After Zoltan_Comm_Create");
  }

  // allocate receive buffer for nrec DD_Remove_Msg structures
  std::vector<char> rbuff(nrec*remove_msg_size);   // receive buffer

  // send my update messages & receive updates directed to me
#ifdef CONVERT_DIRECTORY_KOKKOS
  err = directoryComm.do_forward(ZOLTAN2_DD_UPDATE_MSG_TAG+1,
    sbuff, remove_msg_size, rbuff);
#else
  err = Zoltan_Comm_Do (plan, ZOLTAN2_DD_UPDATE_MSG_TAG+1, &(sbuff[0]),
    remove_msg_size, &(rbuff[0]));
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
    // TODO: Fix cast
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
/*
template <typename gid_t,lid_t,typename user_t>
double *Zoltan2_Directory<gid_t,lid_t,user_t>::array_alloc(
  char *file, int lineno, int numdim, ...) {
  const char *yo = "Zoltan2_Array_Alloc";
  int i, j;
  struct dimension {
    long index;  // Number of elements in the dimension
    long total;  // Total number of elements
    long size;   // Size of a single element in bytes
    long off;    // offset from beginning of array
  } dim[4];      // Info about each dimension

  va_list va;           // Current pointer in the argument list
  va_start(va, numdim);

  if (numdim <= 0) {
    fprintf(stderr, "%s (%s: %d) ERROR: number of dimensions, %d, is <=0\n",
            yo, file, lineno, numdim);
    va_end(va);
    return((double *) NULL);
  }
  else if (numdim > 4) {
    fprintf(stderr, "%s (%s: %d) ERROR: number of dimensions, %d, is > 4\n",
            yo, file, lineno, numdim);
    va_end(va);
    return((double *) NULL);
  }

  dim[0].index = va_arg(va, int);

  if (dim[0].index <= 0) {
    va_end(va);
    return((double *) NULL);
  }

  dim[0].total = dim[0].index;
  dim[0].size  = sizeof(void *);
  dim[0].off   = 0;
  for (i = 1; i < numdim; i++) {
    dim[i].index = va_arg(va, int);
    if (dim[i].index <= 0) {
      fprintf(stderr, "WARNING: %s (%s: %d) called with dimension %d <= 0, "
              "%ld; will return NULL\n",
              yo, file, lineno, i+1, dim[i].index);
      va_end(va);
      return((double *) NULL);
    }
    dim[i].total = dim[i-1].total * dim[i].index;
    dim[i].size  = sizeof(void *);
    dim[i].off   = dim[i-1].off + dim[i-1].total * dim[i-1].size;
  }

  dim[numdim-1].size = va_arg(va, int);
  va_end(va);

  // Round up the last offset value so data is properly aligned.

  dim[numdim-1].off = dim[numdim-1].size *
    ((dim[numdim-1].off+dim[numdim-1].size-1)/dim[numdim-1].size);

  long total = dim[numdim-1].off + dim[numdim-1].total * dim[numdim-1].size;

  // TODO check me out - refactors from ptr malloc
  std::vector<double> dfield(total);

  // TODO make more C++ like...
  char *field  = (char *) &(dfield[0]); // The multi-dimensional array
  for (i = 0; i < numdim - 1; i++) {
    char **ptr  = (char **) (field + dim[i].off); // Pointer offset
    char *data = (char *) (field + dim[i+1].off); // Data offset
    for (j = 0; j < dim[i].total; j++) {
      ptr[j] = data + j * dim[i+1].size * dim[i+1].index;
    }
  }
  // TODO - didn't finish refactoring this - make use std::vectors
  return dfield;
}
*/
#endif

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
unsigned int Zoltan2_Directory<gid_t,lid_t,user_t>::hash_proc(const gid_t & key) const
{
  uint32_t k;
  MurmurHash3_x86_32((void *)(&key), sizeof(gid_t), 14, (void *)&k);
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

template <typename gid_t,typename lid_t,typename user_t>
size_t Zoltan2_Directory<gid_t,lid_t,user_t>::align_size_t(size_t a) const
{
  #define ZOLTAN2_ALIGN_VAL 7U
  return((ZOLTAN2_ALIGN_VAL + a) & ~ZOLTAN2_ALIGN_VAL);
}

// Stats block
template <typename gid_t,typename lid_t,typename user_t>
void Zoltan2_Directory<gid_t,lid_t,user_t>::stats() const
{
  throw std::logic_error("UNTESTED CHECKPOINT"); // needs unit testing

  const char    *yo = "Zoltan2_Directory::stats";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

#ifndef CONVERT_DIRECTORY_RELIC
  throw std::logic_error( "stats not implemented yet..." );
#else
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
        char str[100];      // used to build message string
   //     sprintf(str, "GID %zu, Owner %d, Table Index %d.",
   //       (size_t) *nodelist[nodeidx].gid, nodelist[nodeidx].owner, i);
        ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);
      }
      length++;                  /* linked list length */
      node_count++;              /* count of Nodes */
    }
    if (length > maxlength)
      maxlength = length;        /* save length of longest linked list */
  }

  char str[100];      // used to build message string
  sprintf(str, "Hash table size %d, %d nodes on %d lists, max list length %d.",
          table_length, node_count, list_count, maxlength);
  ZOLTAN2_PRINT_INFO (comm->getRank(), yo, str);

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }
#endif
}

} // end namespace Zoltan2

#endif // CONVERT_DIRECTORY_ORIGINAL
