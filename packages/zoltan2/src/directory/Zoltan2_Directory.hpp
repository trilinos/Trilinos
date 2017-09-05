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

// This is temporary for timing and debugging - to be deleted
#include "Zoltan2_Directory_Clock.hpp"

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

template <typename gid_t, typename lid_t, typename user_t>
class Zoltan2_Directory {
  private:
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

    Zoltan2_Directory(Teuchos::RCP<const Teuchos::Comm<int> > comm_, bool use_lid_,
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
      Zoltan2_Directory_Clock clock("construct");
      allocate();
    }

    Zoltan2_Directory (const Zoltan2_Directory<gid_t,lid_t,user_t> &dd);

    ~Zoltan2_Directory();

    int update(const std::vector<gid_t>& gid, const std::vector<lid_t>& lid,
      const std::vector<user_t>& user, const std::vector<int>& partition,
      Update_Mode mode);

    int find(const std::vector<gid_t>& gid, std::vector<lid_t>& lid,
      std::vector<user_t>& user, std::vector<int>& partition,
      std::vector<int>& owner, bool throw_if_missing = true);

    int remove(const std::vector<gid_t>& gid);

    int print() const;

    Zoltan2_Directory<gid_t,lid_t,user_t> & operator=
      (const Zoltan2_Directory<gid_t,lid_t,user_t> &src);

    void stats() const;

  protected:
    size_t find_msg_size;    /* Total allocation for Zoltan2_DD_FindMsg       */
    size_t update_msg_size;  /* Total allocation for Zoltan2_DD_Update_Msg    */
    size_t remove_msg_size;  /* Total allocation for Zoltan2_DD_Remove_Msg    */

  #ifdef CONVERT_DIRECTORY_KOKKOS
    typedef Kokkos::UnorderedMap<gid_t,
      Zoltan2_Directory_Node<gid_t,lid_t,user_t>> node_map_t;
    node_map_t node_map;
  #endif

  #ifndef CONVERT_DIRECTORY_TPETRA
    int update_local(gid_t* gid, lid_t* lid, user_t *user,
      int partition, int owner, Update_Mode mode);
    int find_local(gid_t* gid, lid_t* lid, user_t *user,
      int *partition, int *owner, bool throw_if_missing = true) const;
    int remove_local(gid_t* gid);
  #endif

    void allocate();
    int copy(const Zoltan2_Directory<gid_t,lid_t,user_t> &dd);

  #ifdef CONVERT_DIRECTORY_RELIC
    int allocate_node_list(relice_idx_t count, float overalloc);
    relice_idx_t allocate_node();
    void free_node(relice_idx_t freenode);
  #endif

    // TODO: finish refactor
    // double* array_alloc(char *file, int lineno, int numdim, ...);

  #ifdef CONVERT_DIRECTORY_RELIC
    int equal_id(int n, gid_t* a, gid_t* b) const;
    unsigned int hash_table(const gid_t& key) const;
  #endif
    size_t align_size_t(size_t a) const;

    unsigned int hash_proc(const gid_t& key) const;

    Teuchos::RCP<const Teuchos::Comm<int> > comm;
    bool contiguous;

    bool use_lid;                  /* If false not using lid                 */
    int debug_level;               /* Determines actions to multiple updates */

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

    virtual void user_to_raw(const user_t & src, user_t * pRaw) const {
      *pRaw = src;
    }

    virtual void raw_to_user(const user_t * pRaw, user_t & dst) const {
      dst = *pRaw;
    }

    virtual size_t size_of_value_type() const { return sizeof(user_t); }

    virtual size_t get_update_msg_size(const user_t & data) const {
      return update_msg_size;
    }

    virtual size_t get_update_msg_size(const user_t * pRaw) const {
      return update_msg_size;
    }

    virtual size_t get_local_find_msg_size(gid_t *gid,
      bool throw_if_missing = true) const {
      return find_msg_size;
    }

    virtual size_t get_incoming_find_msg_size(Zoltan2_DD_Find_Msg<gid_t,lid_t>* msg) const {
      return find_msg_size;
    }

    virtual void update_local_user(const user_t * pRaw, user_t & dst,
      typename Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode mode) = 0;

    virtual bool is_Zoltan2_Directory_Vector() const = 0;
};

template <typename gid_t, typename lid_t, typename user_t>
class Zoltan2_Directory_Simple : public Zoltan2_Directory<gid_t, lid_t, user_t> {
  public:
    typedef user_t user_val_t;

    Zoltan2_Directory_Simple(Teuchos::RCP<const Teuchos::Comm<int> > comm_, bool use_lid_,
    #ifdef CONVERT_DIRECTORY_RELIC
      int table_length_,
    #endif
      int debug_level_, bool contiguous_) :
      Zoltan2_Directory<gid_t, lid_t, user_t>(comm_, use_lid_,
      #ifdef CONVERT_DIRECTORY_RELIC
        table_length_,
      #endif
        debug_level_, contiguous_) {
    }

    virtual void update_local_user(const user_t * pRaw, user_t & dst,
      typename Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode mode) {
      switch(mode) {
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Replace:
          dst = *pRaw;
          break;
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Add:
           dst += *pRaw;
        break;
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Aggregate:
          throw std::logic_error("Aggregate doesn't have meaning for single type.");
          break;
      }
    }

    // awkward to have this at all - so maybe to refactor out with future progress
    virtual bool is_Zoltan2_Directory_Vector() const { return false; }
};

template <typename gid_t, typename lid_t, typename user_t>
class Zoltan2_Directory_Vector : public Zoltan2_Directory<gid_t, lid_t, user_t> {
  public:
    typedef typename user_t::value_type user_val_t;

    Zoltan2_Directory_Vector(Teuchos::RCP<const Teuchos::Comm<int> > comm_, bool use_lid_,
      int debug_level_, bool contiguous_) :
      Zoltan2_Directory<gid_t, lid_t, user_t>(comm_, use_lid_,
        debug_level_, contiguous_) {
    }

  private:
    virtual void user_to_raw(const user_t & src, user_t * pRaw) const {
      user_val_t *pWrite = (user_val_t*)(pRaw);
      *pWrite = src.size();;
      for(size_t n = 0; n < src.size(); ++n) {
        ++pWrite;
        *pWrite = src[n];
      }
    }

    virtual void raw_to_user(const user_t * pRaw, user_t & dst) const {
      user_val_t* pRead = (user_val_t*) pRaw;
      dst.resize(static_cast<size_t>(*pRead));
      for(size_t n = 0; n < dst.size(); ++n) {
        ++pRead;
        dst[n] = *pRead;
      }
    }

    virtual void update_local_user(const user_t * pRaw, user_t & dst,
      typename Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode mode) {
      user_val_t * pRead = (user_val_t*)(pRaw);
      size_t read_array_length = static_cast<size_t>(*pRead);
      switch(mode) {
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Replace: {
          dst.resize(read_array_length); // change to new
          for(size_t i = 0; i < read_array_length; ++i) {
            ++pRead;
            dst[i] = *pRead;
          }
        }
        break;
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Add: {
          if(dst.size() != static_cast<size_t>(read_array_length)) {
            throw std::logic_error("The data lengths are not the same size");
          }
          for(size_t i = 0; i < dst.size(); ++i) {
            ++pRead;
            dst[i] += *pRead;
          }
        }
        break;
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Aggregate: {
          // Add only unique elements
          // Preserve ordering
          // First scan the new incoming data
          for(size_t i = 0; i < read_array_length; ++i) {
            ++pRead; // get the next incoming array element (*pRead)
            if(dst.size() == 0 ||
              (*pRead) > dst[dst.size()-1]) {
              user_val_t theValue = (*pRead);
              dst.push_back(theValue); // add first element or at end
            }
            else {
              for(auto itr = dst.begin(); itr != dst.end(); ++itr) {
                if((*itr) == (*pRead)) { // do they match
                  break; // break because it's already in there
                }
                else if((*itr) > (*pRead)) { // is scanned element larger?
                  dst.insert(itr, (*pRead)); // preserve ordering
                  break; // break because once we add it we are done
                }
              }
            }
          }

        }
        break;
      }
    }

    virtual size_t size_of_value_type() const { return sizeof(typename user_t::value_type); }

    virtual size_t get_update_msg_size(const user_t & data) const {
      return this->update_msg_size + data.size() * size_of_value_type();
    }

    virtual size_t get_update_msg_size(const user_t * pRaw) const {
      user_val_t * pRawValue = (user_val_t*) (pRaw);
      return this->update_msg_size + (pRaw ? ((*pRawValue) * size_of_value_type()) : 0);
    }

    virtual size_t get_local_find_msg_size(gid_t * gid,
      bool throw_if_missing = true) const {
      #ifdef CONVERT_DIRECTORY_KOKKOS
        if(this->node_map.exists(*gid)) {
          Zoltan2_Directory_Node<gid_t,lid_t,user_t> & node =
            this->node_map.value_at(this->node_map.find(*gid));
          return this->find_msg_size + node.userData.size() * sizeof(user_val_t);
        }
        else if(throw_if_missing) {
          throw std::logic_error( "Could not find gid in map." );
        }
        else {
          // not clear yet if we ever want to handle this case or always err out
          // I'm using this right now for the unit testing to validate that the
          // remove command actually works
          return this->find_msg_size; // will not have any data content
        }
      #else
        throw std::logic_error( "Did not implement variable array support for the "
          "non kokkos modes yet." );
      #endif
    }

    virtual size_t get_incoming_find_msg_size(Zoltan2_DD_Find_Msg<gid_t,lid_t>* msg) const {
      user_val_t * pVal = (user_val_t*)(msg->adjData + 1);
      return this->find_msg_size + (*pVal) * sizeof(user_val_t);
    }

    virtual bool is_Zoltan2_Directory_Vector() const { return true; }
};

} // end namespace Zoltan2

#endif // CONVERT_DIRECTORY_ORIGINAL - should not be using this code at all

#endif
