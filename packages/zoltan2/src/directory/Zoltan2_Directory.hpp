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
 * 3. Neither the name of the Corporation nor the names of the
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

#ifndef ZOLTAN2_DIRECTORY_H_
#define ZOLTAN2_DIRECTORY_H_

#ifndef CONVERT_DIRECTORY_ORIGINAL

#include <Teuchos_DefaultComm.hpp>
#include <vector>

#ifdef CONVERT_DIRECTORY_TPETRA
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#endif

#ifdef CONVERT_DIRECTORY_KOKKOS
  #include <Kokkos_UnorderedMap.hpp>
  #include <Kokkos_Core.hpp>
#endif

namespace Zoltan2 {

#ifdef CONVERT_DIRECTORY_KOKKOS
template <typename gid_t, typename lid_t, typename user_t>
class Zoltan2_Directory_Node {
  public:
    int owner;       /* processor hosting global ID object    */
    int partition;   /* Optional data                         */
    int errcheck;    /* Error checking(inconsistent updates)  */
    lid_t lid;       /* lid value */
    user_t userData; /* user data */
};
#endif

#ifdef CONVERT_DIRECTORY_RELIC
typedef int relice_idx_t;  /* nodeList index is signed since -1 = NULL */
template <typename gid_t, typename lid_t,typename user_t>
class Zoltan2_Directory_Node {
  public:
    int owner;                 /* processor hosting global ID object    */
    int partition;             /* Optional data                         */
    int errcheck;              /* Error checking(inconsistent updates)  */
    relice_idx_t next;         /* index in nodelist of next DD_Node     */
    gid_t *gid;                /* gid used as key for update & lookup   */
                               /* lid starts at gid + sizeof(gid_t)     */
                               /* user data starts at                   */
                               /*   gid + sizeof(gid_t) + sizeof(lid_t) */
    int free;                  /* flag whether node is free or used     */
};
#endif

template <typename gid_t, typename lid_t, typename user_t>
class Zoltan2_Directory {
#ifdef TEMP_TRIAL_USER_ARRAY_TYPE
  typedef typename user_t::value_type user_val_t;
#endif

#ifdef CONVERT_DIRECTORY_TPETRA
  typedef Tpetra::Map<int, gid_t> map_t;
  typedef Teuchos::RCP<const map_t> rcp_map_t;
  typedef Tpetra::Vector<user_t, int, gid_t> vector_t;
  typedef Teuchos::RCP<vector_t> rcp_vector_t;
  typedef Tpetra::Export<int, gid_t> export_t;
  typedef Teuchos::RCP<const export_t> rcp_export_t;
  typedef Teuchos::ArrayRCP<user_t> vectordata_t;
#endif

public:
  enum Update_Mode {
    Replace = 0,
    Add,
    Aggregate
  };

  Zoltan2_Directory(Teuchos::RCP<const Teuchos::Comm<int> > comm, bool use_lid,
#ifdef CONVERT_DIRECTORY_RELIC
    int table_length,
#endif
    int debug_level, bool contiguous);

  Zoltan2_Directory (const Zoltan2_Directory<gid_t,lid_t,user_t> &dd);

  ~Zoltan2_Directory();

  int update(const std::vector<gid_t>& gid, const std::vector<lid_t>& lid,
    const std::vector<user_t>& user, const std::vector<int>& partition,
    Update_Mode mode);

  int find(const std::vector<gid_t>& gid, std::vector<lid_t>& lid,
    std::vector<user_t>& user, std::vector<int>& partition,
    std::vector<int>& owner);

  int remove(const std::vector<gid_t>& gid);

  int print() const;

  Zoltan2_Directory<gid_t,lid_t,user_t> & operator=
    (const Zoltan2_Directory<gid_t,lid_t,user_t> &src);

  void stats() const;

private:
  void allocate();
  void release(); // release all memory
  int copy(const Zoltan2_Directory<gid_t,lid_t,user_t> &dd);

#ifdef CONVERT_DIRECTORY_RELIC
  int allocate_node_list(relice_idx_t count, float overalloc);
  relice_idx_t allocate_node();
  void free_node(relice_idx_t freenode);
#endif

#ifndef CONVERT_DIRECTORY_TPETRA
  int update_local(gid_t* gid, lid_t* lid, user_t *user,
    int partition, int owner, Update_Mode mode);
  int find_local(gid_t* gid, lid_t* lid, user_t *user,
    int *partition, int *owner, bool bVariableData);
  int remove_local(gid_t* gid);
#endif

  // TODO: finish refactor
  // double* array_alloc(char *file, int lineno, int numdim, ...);

#ifdef CONVERT_DIRECTORY_RELIC
  int equal_id(int n, gid_t* a, gid_t* b);
  unsigned int hash_table(const gid_t& key);
#endif
  size_t align_size_t(size_t a);

  unsigned int hash_proc(const gid_t& key);

  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  bool contiguous;

  bool use_lid;                  /* If false not using lid                 */
  size_t find_msg_size;          /* Total allocation for DD_FIND_MSG       */
  size_t update_msg_size;        /* Total allocation for DD_UPDATE_MSG     */
  size_t remove_msg_size;        /* Total allocation for DD_REMOVE_MSG     */
  int debug_level;               /* Determines actions to multiple updates */

#ifdef CONVERT_DIRECTORY_KOKKOS
  typedef Kokkos::UnorderedMap<gid_t, Zoltan2_Directory_Node<gid_t,lid_t,user_t>>
    node_map_t;
  node_map_t node_map;
#endif

#ifdef CONVERT_DIRECTORY_TPETRA
  rcp_map_t oto_idMap;
  rcp_vector_t oto_idVec;
#endif

#ifdef CONVERT_DIRECTORY_RELIC
  int table_length;                      /* # of heads of linked lists */
  unsigned int recommended_hash_size (unsigned int n);
   /* Memory for storing all nodes in the directory */
  std::vector<Zoltan2_Directory_Node<gid_t,lid_t,user_t>> nodelist;
  relice_idx_t nodelistlen;        /* Length of the nodelist. */
  /* Index of first free node in nodelist; -1 if no nodes are free */
  relice_idx_t nextfreenode;
  std::vector<relice_idx_t> table; /* Hash table heads of link lists */
  std::vector<char> nodedata; /* Memory for storing data in the directory  */
  size_t nodedata_size;          /* Malloc for GID & LID & user storage    */
#endif
};

}; // end namespace Zoltan2

#endif // CONVERT_DIRECTORY_ORIGINAL - should not be using this code at all

#endif
