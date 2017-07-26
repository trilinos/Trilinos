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

// This is temporary for prototyping a std::vector
// For now only Kokkos will set this since the others don't do that
// #define TEMP_TRIAL_USER_ARRAY_TYPE

// This is temporary for turning off resizing and using all arrays same length
// #define TEMP_USED_FIXED_SIZE

// This is temporary for making arrays same length - but can still use resize
// #define TEMP_CONSTANT_ARRAY_SIZE

// This is temporary - If arrays are constant length - use this val
#define CONSTANT_ARRAY_SIZE 5 // temp testing

#include "Zoltan2_Directory.hpp"

// This is temporary for timing and debugging - to be deleted
#include "Zoltan2_Directory_Clock.hpp"

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

// TODO: These message structures should become MPI datatypes(KDD)
template <typename gid_t, typename lid_t>
class Zoltan2_DD_Update_Msg {  /* Only used by Zoltan_DD_Update()             */
  public:
   char lid_flag;              /* indicates if LID data are present           */
   char user_flag;             /* indicates if USER data are present          */
   char partition_flag;        /* indicates if optional partition data        */
   int owner;                  /* range [0, nproc-1] or -1                    */
   int partition;
   gid_t adjData[1];           /* TODO: refactor - includes gid & lid & user  */
};

template <typename gid_t, typename lid_t>
class Zoltan2_DD_Find_Msg {  /* Only used by Zoltan_DD_Find()          */
  public:
   int proc;                 /* destination or location                */
   int partition;
   int index;                /* to put things back in order afterward  */
   gid_t adjData[1];         /* TODO: refactor - includes gid and user */
};

template <typename gid_t, typename lid_t>
class Zoltan2_DD_Remove_Msg {  /* Only used by Zoltan_DD_Remove()       */
  public:
   int owner;                  /* range [0, nproc-1] or -1              */
   gid_t adjData[1];           /* TODO: refactor - includes gid         */
};

/* Tags for MPI communications.  These need unique values. Arbitrary */
#define ZOLTAN2_DD_FIND_MSG_TAG     29137  /* needs 3 consecutive values */
#define ZOLTAN2_DD_UPDATE_MSG_TAG   29140  /* needs 2 consecutive values */
#define ZOLTAN2_DD_REMOVE_MSG_TAG   29142  /* needs 2 consecutive values */
#define ZOLTAN2_DD_RESIZE_MSG_TAG   29150  /*  */


template <typename gid_t,typename lid_t,typename user_t>
Zoltan2_Directory<gid_t,lid_t,user_t>::Zoltan2_Directory(
  Teuchos::RCP<const Teuchos::Comm<int> > comm_, bool use_lid_,
#ifdef CONVERT_DIRECTORY_RELIC
    int table_length_,
#endif
    int debug_level_, bool contiguous_) :
    comm(comm_), contiguous(contiguous_), use_lid(use_lid_),
    debug_level(debug_level_)
#ifdef CONVERT_DIRECTORY_RELIC
    , table_length(table_length_)
#endif
{
  Zoltan2_Directory_Clock clock("construct", comm);
  allocate();
}

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
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
  // +1 for vector size - then variable length will include the elements
  size_t user_size = sizeof(user_val_t); // TODO make it write a size_t for array length
#ifdef TEMP_USED_FIXED_SIZE
  user_size += CONSTANT_ARRAY_SIZE * sizeof(user_val_t); // fixed array size
#endif
#else
  size_t user_size = sizeof(user_t); // old way
#endif

  size_t size = sizeof(gid_t) + (use_lid?sizeof(lid_t):0) + user_size;

#ifdef CONVERT_DIRECTORY_RELIC
  nodedata_size = size;
#endif

  update_msg_size = size + sizeof(Zoltan2_DD_Update_Msg<gid_t,lid_t>);

  size = sizeof(gid_t);
  remove_msg_size = size + sizeof(Zoltan2_DD_Remove_Msg<gid_t,lid_t>);

  size = sizeof(lid_t) + user_size;
  find_msg_size = size + sizeof(Zoltan2_DD_Find_Msg<gid_t,lid_t>);

  /* force alignment */
#ifdef CONVERT_DIRECTORY_RELIC
  nodedata_size   = align_size_t(nodedata_size);
#endif

//  update_msg_size = align_size_t(update_msg_size);
//  remove_msg_size = align_size_t(remove_msg_size);
//  find_msg_size   = align_size_t(find_msg_size);

  if (debug_level > 4) {
    ZOLTAN2_TRACE_OUT (comm->getRank(), yo, NULL);
  }
};

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
  Zoltan2_Directory_Clock clock("update", comm);

  const char * yo = "Zoltan2_Directory::update";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN(comm->getRank(), yo, NULL);
  }

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

  int err = 0;

#ifdef CONVERT_DIRECTORY_KOKKOS
  // TODO decide how to best allocate initial size
  size_t globalCount = 0;
  Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, gid.size(),
    Teuchos::outArg(globalCount));
  if(node_map.capacity() < globalCount) {
    node_map.rehash(globalCount);
  }
#endif

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
  // This vector will store number of occurrences of each id across procs
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

  Zoltan2_Directory_Clock buildUpdateMSGClock("build update", comm);

  // allocate memory for list of processors to contact
  std::vector<int> procs(gid.size());

  std::vector<int> msg_sizes(gid.size());
  int sum_msg_size = 0;
  for (size_t i = 0; i < gid.size(); i++) {
    size_t msg_size = update_msg_size;
    if(user.size()) {
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
#ifndef TEMP_USED_FIXED_SIZE
      msg_size += user[i].size() * sizeof(user_val_t);
#endif
#endif
    }
    sum_msg_size += msg_size;
    msg_sizes[i] = msg_size;
  }
  std::vector<char> sbuff(sum_msg_size);

  typedef Zoltan2_DD_Update_Msg<gid_t,lid_t> msg_t;

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
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
      // prototyping vector usage
      // right now occupy the fixed memory space
      user_val_t *pWrite = (user_val_t*)(puser);
      *pWrite = user[i].size();
      //printf( "Write update %lu: ", user[i].size());
      for(size_t n = 0; n < user[i].size(); ++n) {
        ++pWrite;
        *pWrite = user[i][n];
        //printf( "%d ", user[i][n] );
      }
      //printf("\n");
#else
      *puser = user[i];
#endif
    }
    else {
      *puser = user_t(); // create a null version... how to handle
    }

    track_offset += update_msg_size;

#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
#ifndef TEMP_USED_FIXED_SIZE
    track_offset += (user.size() ? (user[i].size() * sizeof(user_val_t)) : 0);
#endif
#endif

  }

  if(track_offset != sum_msg_size) {
    throw std::logic_error("Bad summing!");
  }

  buildUpdateMSGClock.complete();

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After fill contact list");
  }

  // now create efficient communication plan

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

  Zoltan2_Directory_Clock updateSendClock("update send", comm);

#ifdef TEMP_USED_FIXED_SIZE // Not array
  int sum_recv_sizes = nrec * update_msg_size;
#else
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
#endif

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
#ifdef TEMP_USED_FIXED_SIZE // Not array
  const int nbytes = update_msg_size;
#else
  const int nbytes = 1;
#endif

#ifdef CONVERT_DIRECTORY_KOKKOS
  err = directoryComm.do_forward(ZOLTAN2_DD_UPDATE_MSG_TAG+1,
    sbuff, nbytes, rbuff);
#else
  err = Zoltan_Comm_Do(plan, ZOLTAN2_DD_UPDATE_MSG_TAG+1, &(sbuff[0]),
    nbytes, &(rbuff[0]));
#endif

  if (err) {
    throw std::logic_error("Zoltan2_Directory::update() Comm_Do error");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Do");
  }

  updateSendClock.complete();

  Zoltan2_Directory_Clock updateLocalClock("update_local", comm);

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

    track_offset += update_msg_size;

#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
#ifndef TEMP_USED_FIXED_SIZE
    user_val_t * pVal = (user_val_t*)(puser);
    track_offset += (puser ? ((*pVal) * sizeof(user_val_t)) : 0);
#endif
#endif
  }

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

  updateLocalClock.complete();

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
};

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
        #ifndef TEMP_TRIAL_USER_ARRAY_TYPE
        user_t * puser = reinterpret_cast<user_t*>(
          reinterpret_cast<char*>(node.gid) +
          sizeof(gid_t) + (use_lid?sizeof(lid_t):0));
        #endif
        switch(mode) {
          case Update_Mode::Replace:
            #ifdef TEMP_TRIAL_USER_ARRAY_TYPE
              throw std::logic_error("Did not refactor relic mode for Array Replace.");
            #else
              *puser = *user;
            #endif
            break;
          case Update_Mode::Add:
            #ifdef TEMP_TRIAL_USER_ARRAY_TYPE
              throw std::logic_error("Did not refactor relic mode for Array Add.");
            #else
              *puser += *user;
            #endif
            break;
          case Update_Mode::Aggregate:
            #ifdef TEMP_TRIAL_USER_ARRAY_TYPE
              throw std::logic_error("TODO - implement aggregate for Relic.");
            #else
              throw std::logic_error("Aggregate doesn't have meaning for single type.");
            #endif
            break;
        }
      }
#else
      // found match, update directory information
      if (lid) {
        node.lid = *lid;
      }
      if (user) {
        #ifdef TEMP_TRIAL_USER_ARRAY_TYPE
          user_val_t * pRead = (user_val_t*)(user);
          size_t read_array_length = static_cast<size_t>(*pRead);
        #endif
        switch(mode) {
          case Update_Mode::Replace:
            #ifdef TEMP_TRIAL_USER_ARRAY_TYPE
              node.userData.resize(read_array_length); // change to new
              for(size_t i = 0; i < node.userData.size(); ++i) {
                ++pRead;
                node.userData[i] = *pRead;
              }
            #else
              node.userData = *user;
            #endif
            break;
          case Update_Mode::Add:
            #ifdef TEMP_TRIAL_USER_ARRAY_TYPE
              if(node.userData.size() != static_cast<size_t>(read_array_length)) {
                throw std::logic_error("The data lengths are not the same size");
              }
              for(size_t i = 0; i < node.userData.size(); ++i) {
                ++pRead;
                node.userData[i] += *pRead;
              }
            #else
              node.userData += *user;
            #endif
            break;
          case Update_Mode::Aggregate:
            #ifdef TEMP_TRIAL_USER_ARRAY_TYPE
              // Add only unique elements
              // Preserve ordering
              // TODO: Make it faster... optimize
              // First scan the new incoming data
              for(size_t i = 0; i < read_array_length; ++i) {
                ++pRead; // get the next incoming array element (*pRead)
                if(node.userData.size() == 0 ||
                  (*pRead) > node.userData[node.userData.size()-1]) {
                  node.userData.push_back(*pRead); // add first element or at end
                }
                else {
                  for(auto itr = node.userData.begin(); // scan final aggregated
                    itr != node.userData.end(); ++itr) {
                    if((*itr) == (*pRead)) { // do they match
                      break; // break because it's already in there
                    }
                    else if((*itr) > (*pRead)) { // is scanned element larger?
                      node.userData.insert(itr, (*pRead)); // preserve ordering
                      break; // break because once we add it we are done
                    }
                  }
                }
              }
            #else
              throw std::logic_error("Aggregate doesn't have meaning for single type.");
            #endif
            break;
        }
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
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
    user_val_t* pRead = (user_val_t*) user;
    node.userData.resize(static_cast<size_t>(*pRead));
    // printf( "Reading update_local %lu ", node.userData.size());
    for(size_t n = 0; n < node.userData.size(); ++n) {
      ++pRead;
      node.userData[n] = *pRead;
      // printf( "%d ", *pRead);
    }
    // printf("\n");
#else
    node.userData = *user;
#endif
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
     throw std::logic_error(
       "Hash insert failed. Maybe size not setup right....");
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
 std::vector<int>& owner)       /* Outgoing optional list of data owners      */
{
  Zoltan2_Directory_Clock clock("find", comm);

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
  std::vector<int> procs(gid.size());     // list of processors to contact

  /* Setup procs list */
  for (size_t i = 0; i < gid.size(); i++)  {
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

  std::vector<char> sbuff(find_msg_size*gid.size());     // send buffer

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

  // send out find messages across entire system
  Zoltan2_Directory_Clock sendClock("send", comm);

  const int nbytes = find_msg_size; // testing without variable size

#ifdef CONVERT_DIRECTORY_RELIC
  err = Zoltan_Comm_Do (plan, ZOLTAN2_DD_FIND_MSG_TAG+1, &(sbuff[0]),
    nbytes, &(rbuff[0]));
#else
  err = directoryComm.do_forward(ZOLTAN2_DD_FIND_MSG_TAG+1, sbuff,
    nbytes, rbuff);
#endif
  sendClock.complete();

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
  int sum_rmsg_sizes_resized = 0;
  for (int i = 0; i < nrec; i++) {
    size_t find_rmsg_size_resized = find_msg_size;
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
#ifndef TEMP_USED_FIXED_SIZE
    // First read the incoming rbuff
    // TODO: Fix cast
    msg_t *ptr = (msg_t*)(&(rbuff[track_offset]));

    // extract the variable length
    // TODO: Fix casts
    user_t * puser = (user_t*)(ptr->adjData + 1);

    user_val_t * pVal = (user_val_t*)(puser);

    // access our local database and determine the length of the array
    // fill the message and track the new resized values we will use
    // This doesn't collect the array values, just the length at this point.
    err = find_local(ptr->adjData, NULL, // (lid_t*)ptr->adjData,
      puser, &ptr->partition, &ptr->proc, false);
    find_rmsg_size_resized += (*pVal) * sizeof(user_val_t);
#endif
#endif
    track_offset += find_msg_size;
    rmsg_sizes_resized[i] = find_rmsg_size_resized;
    sum_rmsg_sizes_resized += find_rmsg_size_resized;
  }

#ifdef TEMP_USED_FIXED_SIZE // Not array
  std::vector<char> rbuff_resized = rbuff;
#else
  std::vector<char> rbuff_resized(sum_rmsg_sizes_resized);

  int track_offset_resized = 0;
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
    msg_t *ptr = (msg_t*)(&(rbuff_resized[track_offset_resized]));
    // TODO: Fix casts
    user_t * puser = (user_t*)(ptr->adjData + 1);
    err = find_local(ptr->adjData, NULL, // (lid_t*)ptr->adjData,
      puser, &ptr->partition, &ptr->proc, true);
    track_offset_resized += rmsg_sizes_resized[i];
  }
 #endif

  if(track_offset_resized != sum_rmsg_sizes_resized) {
    throw std::logic_error("Bad sum!");
  }

  Zoltan2_Directory_Clock reverseClock("reverse", comm);
#if defined(TEMP_USED_FIXED_SIZE) || not defined(TEMP_TRIAL_USER_ARRAY_TYPE) // Not array
  std::vector<char> sbuff_resized = sbuff;
#ifdef CONVERT_DIRECTORY_RELIC
  err = Zoltan_Comm_Do_Reverse(plan, ZOLTAN2_DD_FIND_MSG_TAG+2, &(rbuff_resized[0]),
    find_msg_size, NULL, &(sbuff_resized[0]));
#else
  err = directoryComm.do_reverse(ZOLTAN2_DD_FIND_MSG_TAG+2, rbuff_resized,
    find_msg_size, std::vector<int>(), sbuff_resized);
#endif

#else
  // send return information back to requester
  // resize the directory to handle the array contents
  std::vector<char> sbuff_resized;  // will be filled in reverse call
#ifdef CONVERT_DIRECTORY_RELIC
  err = Zoltan_Comm_Do_Reverse(plan, ZOLTAN2_DD_FIND_MSG_TAG+2, &(rbuff_resized[0]),
    1, &(rmsg_sizes_resized[0]), &(sbuff_resized[0]));
#else
  err = directoryComm.do_reverse(ZOLTAN2_DD_FIND_MSG_TAG+2, rbuff_resized,
    1, rmsg_sizes_resized, sbuff_resized);
#endif

#endif

  reverseClock.complete();

  if (err) {
    throw std::logic_error("Zoltan2_Directory::find() do reverse failed");
  }

  if (debug_level > 6) {
    ZOLTAN2_PRINT_INFO(comm->getRank(), yo, "After Comm_Reverse");
  }

  // fill in user supplied lists with returned information
  Zoltan2_Directory_Clock fillClock("fill", comm);
  track_offset_resized = 0;

  // it's not going to be in order... so don't use gid[i]
  for (size_t i = 0; i < gid.size(); i++) {
    // TODO: Fix cast
    msg_t *ptr = (msg_t*)(&(sbuff_resized[track_offset_resized]));

    if (owner.size())
      owner[ptr->index] = ptr->proc;
    if (partition.size())
      partition[ptr->index] = ptr->partition ;
    if (lid.size()) {
      lid[ptr->index] = *(ptr->adjData);
    }

    size_t find_msg_size_resized = find_msg_size;

#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
    user_val_t * pRead = ((user_val_t*)(ptr->adjData+1));
    find_msg_size_resized += (*pRead) * sizeof(user_val_t);
    size_t array_length = *pRead;
    user[ptr->index].resize(array_length);
    for(size_t n = 0; n < array_length; ++n) {
      ++pRead;
      user[ptr->index][n] = *pRead;
    }
#else
    if (user.size()) {
      user[ptr->index] = *((user_t*)(ptr->adjData+1));
    }
#endif

    // don't use smsg_sizes_resized here because we don't get them back
    // in the same order
    track_offset_resized += find_msg_size_resized;
  }

  fillClock.complete();

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
 bool bVariableData)         /* data has adjustable length             */
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
    Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node = nodelist[nodeidx];
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
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
        user_val_t * pWrite = (user_val_t*) user;
        *pWrite = static_cast<user_val_t>(node.userData.size());

        if(bVariableData) {
          for(size_t n = 0; n < node.userData.size(); ++n) {
            ++pWrite;
            *pWrite = node.userData[n];
          }
        }
#else
        *user = node.userData;
#endif
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

//  throw std::logic_error("find_local did not succeed");
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
  const char * yo = "Zoltan2_Directory::remove";

  if (debug_level > 4) {
    ZOLTAN2_TRACE_IN (comm->getRank(), yo, NULL);
  }

  int err = 0;

#ifdef CONVERT_DIRECTORY_TPETRA
  if(gid.size() == 0) {
    return 0; // TEMPORARY - just allows this test to pass for now
  }
  throw std::logic_error("Tpetra does not support remove yet.");
#else

  // allocate memory for processor contact list
  std::vector<int> procs(gid.size());   // list of processors to contact

#ifdef CONVERT_DIRECTORY_KOKKOS
  std::vector<size_t> msg_offset(gid.size()+1); // +1 so last can get size
  size_t sum_msg_size = 0;
  for (size_t i = 0; i < msg_offset.size(); i++) {
    size_t this_msg_size = remove_msg_size; // will be varying soon
    msg_offset[i] = sum_msg_size;
    sum_msg_size += this_msg_size;
  }
  // allocate memory for DD_Remove_Msg send buffer
  std::vector<char> sbuff(sum_msg_size);   // send buffer
#else
  std::vector<char> sbuff(gid.size()*remove_msg_size);   // send buffer
#endif

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
  for (int i = 0; i < nrec; i++)  {
    // TODO: Fix cast
    msg_t *ptr = reinterpret_cast<msg_t*>(&(rbuff[i*remove_msg_size]));
    err = remove_local(ptr->adjData);
    if (err == 1) { // TODO eliminate warns (1) make all errors
      ++errcount;
    }
  }

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
    node_map.begin_erase();
    node_map.erase(*gid);
    node_map.end_erase();
    if(node_map.exists(*gid)) {
      throw std::logic_error( "Failed to erase node!" );
    }
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
#endif

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
int Zoltan2_Directory<gid_t,lid_t,user_t>::equal_id(int n, gid_t* a, gid_t* b)
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
unsigned int Zoltan2_Directory<gid_t,lid_t,user_t>::hash_proc(const gid_t & key)
{
  uint32_t k;
  MurmurHash3_x86_32((void *)(&key), sizeof(gid_t), 14, (void *)&k);
  return(k % comm->getSize());
}

#ifdef CONVERT_DIRECTORY_RELIC
template <typename gid_t,typename lid_t,typename user_t>
unsigned int Zoltan2_Directory<gid_t,lid_t,user_t>::hash_table(const gid_t& key)
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
size_t Zoltan2_Directory<gid_t,lid_t,user_t>::align_size_t(size_t a)
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
