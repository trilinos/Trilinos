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

// Note on file structure
// Currently Zoltan2_Directory_Impl.hpp contains a lot of the method code
// But as this evolved the point of that file is now doubtful - though we may
// like to keep some separation of API and implementation.
//
// The Zoltan2_Directory class is an abstract base class from which we derive:
//    Zoltan2_Directory_Simple - for simple user data of a type (int,long)
//    Zoltan2_Directory_Vector - for std::vector user data of variable length
//
// I have kept all the code specific to the derived classes in this file.


// Block everything out for Original - this is because Original mode points
// back to the original zoltan code and should not be referencing anything here.
#ifndef CONVERT_DIRECTORY_ORIGINAL

#include <Teuchos_DefaultComm.hpp> // currently using Teuchos comm throughout

#ifdef CONVERT_DIRECTORY_TPETRA
  #include <Tpetra_Vector.hpp> // tpetra impleemnts behavior with maps/vectors
#endif

#ifdef CONVERT_DIRECTORY_KOKKOS
  #include <Kokkos_UnorderedMap.hpp> // unordered map stores the local nodes
#endif

// This is temporary for timing and debugging - to be deleted eventually
// The clock times a block of code and also has some debug printing options
// Running test category PERFORMANCE will enable tests which log out time
// data for the modes.
#include "Zoltan2_Directory_Clock.hpp"

namespace Zoltan2 {

#ifdef CONVERT_DIRECTORY_KOKKOS

// The new Kokkos mode maps over the gid using unordered map
// Originally was thinking we don't want ptrs here and want the user data
// to be a natural type so it can cleanly support std::vector. However as this
// evolved I am not sure this is the best decision and can be more costly for
// memory. Now that all the unit testing is in place it would be easier to
// refactor this back to the original packed ptr style if necessary.
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

// relic just implements the node as it was in the orignal zoltan code.
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
// Currently these are implemented as they were in the original zoltan code.
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

/*!  \brief Zoltan2_Directory is an abstract base class.

  The user will implement Zoltan2_Directory_Simple or Zoltan2_Directory_Vector
  in the current form and this class provides all the shared methods.
*/

template <typename gid_t, typename lid_t, typename user_t>
class Zoltan2_Directory {
  public:
    /*! \brief Update_Mode determines how update executes. */
    enum Update_Mode {
      Replace = 0, /*!< \brief The new value replaces the original value. */
      Add, /*!< \brief All values from different procs are summed. */
      Aggregate /*!< \brief For std::vector user data, aggregates all data
                            so for example [1,2,5] and [3,5] becomes [1,2,3,5]*/
    };

    /*! \brief Construct Zoltan2_Directory (abstract class). */
    Zoltan2_Directory(
        Teuchos::RCP<const Teuchos::Comm<int> > comm_,
          /*!< \brief Teuchos comm provided by user. */
        bool use_lid_, /*!< \brief are local IDs being submitted and read. */
    #ifdef CONVERT_DIRECTORY_RELIC
        int table_length_, /*!< \brief table length estimate */
    #endif
        int debug_level_) /*!< \brief debug level controls output */
          : comm(comm_), use_lid(use_lid_),
          debug_level(debug_level_)
    #ifdef CONVERT_DIRECTORY_RELIC
        , table_length(table_length_)
    #endif
    {
      // construct clock tracks total construct time
      Zoltan2_Directory_Clock clock("construct");
    }

    /*! \brief Destructor currently does nothing. */
    virtual ~Zoltan2_Directory() {
    }

    /*! \brief update is called by user to submit new data. */
    int update(
      const std::vector<gid_t>& gid,           /*! \brief gids being updated. */
      const std::vector<lid_t>& lid,                 /*! \brief lids if used. */
      const std::vector<user_t>& user,          /*! \brief user data if used. */
      const std::vector<int>& partition,   /*! \brief partition data if used. */
      Update_Mode update_mode); /*! \brief Can be Replace, Add, or Aggregate. */

    /*! \brief find is called by user to get data back from directory. */
    int find(
      const std::vector<gid_t>& gid,                 /*! \brief gids to find. */
      std::vector<lid_t>& lid,         /*! \brief lids to find if being used. */
      std::vector<user_t>& user,  /*! \brief user data to find if being used. */
      std::vector<int>& partition, /*! \brief partition data to find if used. */
      std::vector<int>& owner,         /*! \brief owner data to find if used. */
      bool throw_if_missing = true); /*! \brief if true will throw if a gid is
        not found. This is used by the unit tests to properly assess if remove
        has worked. */

    /*! \brief remove eliminates these gids from the directory . */
    int remove(
      const std::vector<gid_t>& gid); /*! \brief gids to remove. */

    /*! \brief print output. New Kokkos mode needs further development. */
    int print() const;

    /*! \brief stats. New Kokkos mode needs further development. */
    void stats() const;

    /*! \brief returns true if the directory is handling local ids. */
    bool is_use_lid() const { return use_lid; }

  protected:
  #ifndef CONVERT_DIRECTORY_TPETRA
    // handled updating the local node information when the proc receives
    // a new gid to store or updated data for a preexisting node
    int update_local(gid_t* gid, lid_t* lid, user_t *user,
      int partition, int owner);

    // collect data on the local proc which has been requested by the directory.
    int find_local(gid_t* gid, lid_t* lid, user_t *user,
      int *partition, int *owner, bool throw_if_missing = true) const;

    // remove the locally stored node for this gid
    int remove_local(gid_t* gid);
  #endif

    size_t find_msg_size;    /* Total allocation for Zoltan2_DD_FindMsg       */
    size_t update_msg_size;  /* Total allocation for Zoltan2_DD_Update_Msg    */
    size_t remove_msg_size;  /* Total allocation for Zoltan2_DD_Remove_Msg    */

  #ifdef CONVERT_DIRECTORY_KOKKOS
    // originally the nodes are stored in a hash but in the new Kokkos mode
    // they are stored using Kokkos::UnorderedMap
    typedef Kokkos::UnorderedMap<gid_t,
      Zoltan2_Directory_Node<gid_t,lid_t,user_t>> node_map_t;
    node_map_t node_map;
  #endif

    // this method exists so constructor and copy constructor don't duplicate
    void allocate();

    // this method exists so operator= and copy constructor don't duplicate
    int copy(const Zoltan2_Directory<gid_t,lid_t,user_t> &dd);

  #ifdef CONVERT_DIRECTORY_RELIC
    // original zoltan methods for creating and freeing nodes
    int allocate_node_list(relice_idx_t count, float overalloc);
    relice_idx_t allocate_node();
    void free_node(relice_idx_t freenode);
  #endif

  #ifdef CONVERT_DIRECTORY_RELIC
    int equal_id(int n, gid_t* a, gid_t* b) const;
    unsigned int hash_table(const gid_t& key) const;
  #endif

    // TODO: Decide if this stays and how to incorporate with variable length
    // data. See comments in Zoltan2_Directory_Impl.hpp
    // size_t align_size_t(size_t a) const;

    // take a gid and hash to proc - determines which proc wil own the gid data
    unsigned int hash_proc(const gid_t & gid) const;

    // stores the comm provided by the user
    Teuchos::RCP<const Teuchos::Comm<int> > comm;

    bool use_lid;                  /* If false not using lid                 */
    int debug_level;               /* Determines actions to multiple updates */

    size_t max_id_size;          /* Stores: max(sizeof(gid_t),sizeof(lid_t)) */
    Update_Mode mode; /* Last mode sent using update */

  #ifdef CONVERT_DIRECTORY_RELIC
    int table_length;                      /* # of heads of linked lists */
    unsigned int recommended_hash_size (unsigned int n);
     /* Memory for storing all nodes in the directory */
    Teuchos::ArrayRCP<Zoltan2_Directory_Node<gid_t,lid_t,user_t>> nodelist;
    relice_idx_t nodelistlen;        /* Length of the nodelist. */
    /* Index of first free node in nodelist; -1 if no nodes are free */
    relice_idx_t nextfreenode;
    Teuchos::ArrayRCP<relice_idx_t> table;  /* Hash table heads of link lists */
    Teuchos::ArrayRCP<char> nodedata;     /* Memory for data in the directory */
    size_t nodedata_size;              /* Malloc for GID & LID & user storage */
  #endif

    // abstract methods are implemented below by Zoltan2_Directory_Simple
    // or Zoltan2_Directory_Vector. These methods contain all the places where
    // the code had to be specialized for normal user type (int, long, etc) or
    // std::vector user type of variable length. Possibly we could consider
    // making this all work by templating but exactly how to do that cleanly
    // I am not sure. The class inheritance approach may be easier to understand
    // but it does mean the user has to pick the right class to use.
    virtual bool is_Zoltan2_Directory_Vector() const                      = 0;
    virtual void update_local_user(const user_t * pRaw, user_t & dst)     = 0;
    virtual void user_to_raw(const user_t & src, user_t * pRaw) const     = 0;
    virtual void raw_to_user(const user_t * pRaw, user_t & dst) const     = 0;
    virtual size_t size_of_value_type() const                             = 0;
    virtual size_t get_update_msg_size(const user_t & data) const         = 0;
    virtual size_t get_update_msg_size(const user_t * pRaw) const         = 0;
    virtual size_t get_local_find_msg_size(gid_t *gid,
      bool throw_if_missing = true) const                                 = 0;
    virtual size_t get_incoming_find_msg_size(
      Zoltan2_DD_Find_Msg<gid_t,lid_t>* msg) const                        = 0;

  private:
  #ifdef CONVERT_DIRECTORY_TPETRA
    // Tpetra implemented as an example for comparison - in this case the
    // directory class is more of a wrapper which just implements standard
    // Tpetra calls to export data.
    typedef Tpetra::Map<lid_t, gid_t> map_t;
    typedef Teuchos::RCP<const map_t> rcp_map_t;
    typedef Tpetra::Vector<user_t, lid_t, gid_t> vector_t;
    typedef Teuchos::RCP<vector_t> rcp_vector_t;
    typedef Tpetra::Export<lid_t, gid_t> export_t;
    typedef Teuchos::RCP<const export_t> rcp_export_t;
    typedef Teuchos::ArrayRCP<user_t> vectordata_t;
    rcp_map_t oto_idMap;
    rcp_vector_t oto_idVec;
  #endif
};

template <typename gid_t, typename lid_t, typename user_t>
class Zoltan2_Directory_Simple : public Zoltan2_Directory<gid_t, lid_t, user_t> {
  public:
    typedef user_t user_val_t;

    /*! \brief Constructo directory which handles simple user data types. */
    Zoltan2_Directory_Simple(Teuchos::RCP<const Teuchos::Comm<int> > comm_, bool use_lid_,
    #ifdef CONVERT_DIRECTORY_RELIC
      int table_length_,
    #endif
      int debug_level_) :
      Zoltan2_Directory<gid_t, lid_t, user_t>(comm_, use_lid_,
      #ifdef CONVERT_DIRECTORY_RELIC
        table_length_,
      #endif
        debug_level_) {
      // Note that allocate() must be called in the derived class, not the
      // base class or inheritance of the methods will break
      this->allocate();
    }

    /*! \brief Copy constructor. */
    Zoltan2_Directory_Simple (
      const Zoltan2_Directory_Simple<gid_t,lid_t,user_t> &src) :
        Zoltan2_Directory<gid_t, lid_t, user_t>(src.comm, src.use_lid,
    #ifdef CONVERT_DIRECTORY_RELIC
          src.table_length,
    #endif
          src.debug_level) {
      this->allocate();
      this->copy(src);
    }

    /*! \brief operator= to copy a directory. */
    Zoltan2_Directory_Simple<gid_t,lid_t,user_t> & operator=
      (const Zoltan2_Directory_Simple<gid_t,lid_t,user_t> &src) {
      this->comm = src.comm;
      this->use_lid = src.use_lid;
    #ifdef CONVERT_DIRECTORY_RELIC
      this->table_length = src.table_length;
    #endif
      this->debug_level = src.debug_level;
      this->allocate(); // operator= was setup in derived class so this inherits
      this->copy(src);
      return *this;
    }

  protected:
    // awkward to have this at all - so maybe to refactor out with future progress
    virtual bool is_Zoltan2_Directory_Vector() const { return false; }

    // given raw data from the MPI stream we update user data based on mode
    virtual void update_local_user(const user_t * pRaw, user_t & dst) {
      switch(this->mode) {
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Replace:
          dst = *pRaw;
          break;
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Add:
           dst += *pRaw;
        break;
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Aggregate:
          throw std::logic_error("Aggregate doesn't mean anything for single "
            "types. Must use Zoltan2_Directory_Vector class.");
          break;
      }
    }

    // convert user data to raw data - simple conversion for this class
    virtual void user_to_raw(const user_t & src, user_t * pRaw) const {
      *pRaw = src;
    }

    // convert raw data to user data - simple conversion for this class
    virtual void raw_to_user(const user_t * pRaw, user_t & dst) const {
      dst = *pRaw;
    }

    // get size of the user type which is simply sizeof(user_t) for this class
    virtual size_t size_of_value_type() const { return sizeof(user_t); }

    // for this class, update_msg_size is simple (not variable length)
    virtual size_t get_update_msg_size(const user_t & data) const {
      return this->update_msg_size;
    }

    // for this class, update_msg_size is simple (not variable length)
    virtual size_t get_update_msg_size(const user_t * pRaw) const {
      return this->update_msg_size;
    }

    // for this class, find_msg_size is simple (not variable length)
    virtual size_t get_local_find_msg_size(gid_t *gid,
      bool throw_if_missing = true) const {
      return this->find_msg_size;
    }

    // for this class, find_msg_size is simple (not variable length)
    virtual size_t get_incoming_find_msg_size(
      Zoltan2_DD_Find_Msg<gid_t,lid_t>* msg) const {
      return this->find_msg_size;
    }
};

template <typename gid_t, typename lid_t, typename user_t>
class Zoltan2_Directory_Vector : public Zoltan2_Directory<gid_t, lid_t, user_t> {
  public:
    typedef typename user_t::value_type user_val_t;

    /*! \brief Constructo directory which handles std::vector user data types. */
    Zoltan2_Directory_Vector(Teuchos::RCP<const Teuchos::Comm<int> > comm_, bool use_lid_,
      int debug_level_) :
      Zoltan2_Directory<gid_t, lid_t, user_t>(comm_, use_lid_,
        debug_level_) {
      // Note that allocate() must be called in the derived class, not the
      // base class or inheritance of the methods will break
      this->allocate();
    }

    /*! \brief Copy constructor. */
    Zoltan2_Directory_Vector (
      const Zoltan2_Directory_Vector<gid_t,lid_t,user_t> &src) :
        Zoltan2_Directory<gid_t, lid_t, user_t>(src.comm, src.use_lid,
    #ifdef CONVERT_DIRECTORY_RELIC
          src.table_length,
    #endif
          src.debug_level) {
      this->allocate(); // operator= was setup in derived class so this inherits
      this->copy(src);
    }

    /*! \brief operator= to copy a directory. */
    Zoltan2_Directory_Vector<gid_t,lid_t,user_t> & operator=
      (const Zoltan2_Directory_Vector<gid_t,lid_t,user_t> &src) {
      this->comm = src.comm;
      this->use_lid = src.use_lid;
    #ifdef CONVERT_DIRECTORY_RELIC
      this->table_length = src.table_length;
    #endif
      this->debug_level = src.debug_level;
      this->allocate(); // operator= was setup in derived class so this inherits
      this->copy(src);
      return *this;
    }

  protected:
    // awkward to have this at all - so maybe to refactor out with future progress
    virtual bool is_Zoltan2_Directory_Vector() const { return true; }

    // given raw data from the MPI stream we update user data based on mode
    virtual void update_local_user(const user_t * pRaw, user_t & dst) {
      // we're reading raw data of form: size_t, val, val, val ...
      size_t * pLength = (size_t*)(pRaw);
      size_t read_array_length = *pLength;
      ++pLength; // move up to first element
      user_val_t * pRead = (user_val_t*)(pLength);
      switch(this->mode) {
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Replace: {
          // Note this is the raw_to_user method and we could just call it
          // but Add and Aggregate don't have the equivalent so I've done it
          // this way to keep the pattern.
          dst.resize(read_array_length); // change to new
          for(size_t i = 0; i < read_array_length; ++i) {
            dst[i] = *pRead;
            ++pRead;
          }
        }
        break;
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Add: {
          // ADD currently requires equal length vectors to add each element
          if(dst.size() != static_cast<size_t>(read_array_length)) {
            throw std::logic_error("The data lengths are not the same size");
          }
          // loop through and do the addition
          for(size_t i = 0; i < dst.size(); ++i) {
            dst[i] += *pRead;
            ++pRead;
          }
        }
        break;
        case Zoltan2_Directory<gid_t,lid_t,user_t>::Update_Mode::Aggregate: {
          // Add only unique elements
          // Preserve ordering
          // First scan the new incoming data
          //
          // For example we can have data of:
          //   [1,4,5,7]
          // Then new incoming data is:
          //   [4,5,6,10]
          // Then result would be:
          //   [1,4,5,6,7,10]
          for(size_t i = 0; i < read_array_length; ++i) {
            // handle the cases of dst no size or adding past last element
            if(dst.size() == 0 ||
              (*pRead) > dst[dst.size()-1]) {
              user_val_t theValue = (*pRead);
              dst.push_back(theValue); // add first element or at end
            }
            else {
              // otherwise we are going to insert unless it's not unique
              for(auto itr = dst.begin(); itr != dst.end(); ++itr) {
                if((*itr) == (*pRead)) { // do they match
                  break; // break because it's already in there - do nothing
                }
                else if((*itr) > (*pRead)) { // is scanned element larger?
                  dst.insert(itr, (*pRead)); // preserve ordering
                  break; // break because once we add it we are done
                }
              }
            }
            ++pRead; // get the next incoming array element (*pRead)
          }
        }
        break;
      }
    }

    // write the std::vector as length, x1, x2, x3 ...
    virtual void user_to_raw(const user_t & src, user_t * pRaw) const {
      // we're writing raw data of form: size_t, val, val, val ...
      size_t *pLength = (size_t*)(pRaw);
      *pLength = src.size(); // first write the length
      ++pLength; // move up to first element
      user_val_t *pWrite = (user_val_t*)(pLength);
      for(size_t n = 0; n < src.size(); ++n) {
        *pWrite = src[n]; // now write each element
        ++pWrite;
      }
    }

    // raw comes in as length, x1, x2, x3 ...
    virtual void raw_to_user(const user_t * pRaw, user_t & dst) const {
      // we're reading raw of form: size_t, val, val, val ...
      size_t* pLength = (size_t*) pRaw;
      dst.resize(static_cast<size_t>(*pLength)); // first read the length
      ++pLength; // move up to first element
      user_val_t* pRead = (user_val_t*) pLength;
      for(size_t n = 0; n < dst.size(); ++n) {
        dst[n] = *pRead; // now read each element
        ++pRead;
      }
    }

    // for the std::vector directory, value type is the size of the std::vector
    // template parameter, so for std::vector<int> we want sizeof(int)(
    virtual size_t size_of_value_type() const {
      return sizeof(typename user_t::value_type);
    }

    // the update msg is the base size (includes vector length) plus the size
    // of all the elements.
    virtual size_t get_update_msg_size(const user_t & data) const {
      return this->update_msg_size + data.size() * size_of_value_type();
    }

    // the update msg is the base size (includes vector length) plus the size
    // of all the elements. This is same idea as above method but here we are
    // intepreting raw data (which comes in as length, x1, x2, x3...) so we
    // read the size as the first element. That's all we need to determine the
    // total update_msg_size.
    virtual size_t get_update_msg_size(const user_t * pRaw) const {
      // the first element is size_t (length of the vector)
      size_t * pLength = (size_t*) (pRaw);
      return this->update_msg_size +
        (pRaw ? ((*pLength) * size_of_value_type()) : 0);
    }

    // to get the local find msg size we need to verify the node exists and then
    // read the std::vector size (which is added to base length find_msg_size.
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
          // remove command actually works.
          return this->find_msg_size; // will not have any data content
        }
      #else
        throw std::logic_error( "Did not implement variable array support for "
          "the non kokkos modes yet." );
      #endif
    }

    // here we have a find msg coming in and we need to extract the vector length
    // which is always at the same place in the message.
    virtual size_t get_incoming_find_msg_size(
      Zoltan2_DD_Find_Msg<gid_t,lid_t>* msg) const {
      if(msg->proc == -1) {
        // this happens if we called find for an element which was removed
        // eventually we might just throw on find_local but for the testing,
        // this is allowed, the data is left untouched, and the test validates
        // that the results are as expected based on remove events
        return this->find_msg_size; // no extra data for unfound node
      }
      // the first element of the user data is size_t (length of the vector)
      size_t * pVectorLength =
        (size_t*)(reinterpret_cast<char*>(msg->adjData) + this->max_id_size);
      return this->find_msg_size + (*pVectorLength) * sizeof(user_val_t);
    }
};

} // end namespace Zoltan2

#endif // CONVERT_DIRECTORY_ORIGINAL - should not be using this code at all

#endif
