// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_RCP_NODE_HPP
#define TEUCHOS_RCP_NODE_HPP


/** \file Teuchos_RCPNode.hpp
 *
 * \brief Reference-counted pointer node classes.
 */


#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_any.hpp"
#include "Teuchos_map.hpp"
#include "Teuchos_ENull.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Exceptions.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_toString.hpp"
#include "Teuchos_getBaseObjVoidPtr.hpp"

#if defined(HAVE_TEUCHOSCORE_CXX11) && defined(HAVE_TEUCHOS_THREAD_SAFE) && !defined(DISABLE_ATOMIC_COUNTERS)
#include <atomic>
#define USING_ATOMICS
#endif

namespace Teuchos {

#ifdef USING_ATOMICS
#  define TEUCHOS_RCP_DECL_ATOMIC(VAR, T) std::atomic<T> VAR
#else
#  define TEUCHOS_RCP_DECL_ATOMIC(VAR, T) T VAR
#endif

/** \brief Used to specify a pre or post destruction of extra data
 *
 * \ingroup teuchos_mem_mng_grp
 */
enum EPrePostDestruction { PRE_DESTROY, POST_DESTROY };

/** \brief Used to specify if the pointer is weak or strong.
 *
 * \ingroup teuchos_mem_mng_grp
 */
enum ERCPStrength { RCP_STRONG=0, RCP_WEAK=1 };

/** \brief Used to determine if RCPNode lookup is performed or not.
 *
 * \ingroup teuchos_mem_mng_grp
 */
enum ERCPNodeLookup { RCP_ENABLE_NODE_LOOKUP, RCP_DISABLE_NODE_LOOKUP };

/** \brief . */
inline void debugAssertStrength(ERCPStrength strength)
{
#ifdef TEUCHOS_DEBUG
  switch (strength) {
  case RCP_STRONG:
    // fall through
  case RCP_WEAK:
    return; // Fine
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Teuchos::RCPNode: ERCPStrength enum value "
      << strength << " is invalid (neither RCP_STRONG = " << RCP_STRONG
      << " nor RCP_WEAK = " << RCP_WEAK << ").");
  }
#else
  (void) strength; // Silence "unused variable" compiler warning.
#endif // TEUCHOS_DEBUG
}

/** \brief Traits class specialization for <tt>toString(...)</tt> function for
 * converting from <tt>ERCPStrength</tt> to <tt>std::string</tt>.
 *
 * \ingroup teuchos_mem_mng_grp
 */
template<>
class TEUCHOSCORE_LIB_DLL_EXPORT ToStringTraits<ERCPStrength> {
public:
  static std::string toString( const ERCPStrength &t )
    {
      switch (t) {
        case RCP_STRONG:
          return "RCP_STRONG";
        case RCP_WEAK:
          return "RCP_WEAK";
        default:
          // Should never get here but fall through ...
          break;
      }
      // Should never get here!
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT(true);
#else
      return "";
#endif
    }
};


/** \brief Node class to keep track of address and the reference count for a
 * reference-counted utility class and delete the object.
 *
 * This is not a general user-level class.  This is used in the implementation
 * of all of the reference-counting utility classes.
 *
 * NOTE: The reference counts all start a 0 so the client (i.e. RCPNodeHandle)
 * must increment them from 0 after creation.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class TEUCHOSCORE_LIB_DLL_EXPORT RCPNode {
public:
  /** \brief . */
  RCPNode(bool has_ownership_in)
    : has_ownership_(has_ownership_in), extra_data_map_(NULL)
#ifdef TEUCHOS_DEBUG
    ,insertion_number_(-1)
#endif // TEUCHOS_DEBUG
    {
      count_[RCP_STRONG] = 0;
      count_[RCP_WEAK] = 0;
    }
  /** \brief . */
  virtual ~RCPNode()
    {
      if(extra_data_map_)
        delete extra_data_map_;
    }
  /** \brief attemptIncrementStrongCountFromNonZeroValue() supports weak
   * to strong conversion but this is forward looking code.
   */
  bool attemptIncrementStrongCountFromNonZeroValue()
    {
#ifdef USING_ATOMICS
      // this code follows the boost method
      int strong_count_non_atomic = count_[RCP_STRONG];
      for( ;; ) {
        if (strong_count_non_atomic == 0) {
          return false;
        }
        if (std::atomic_compare_exchange_weak( &count_[RCP_STRONG],
          &strong_count_non_atomic, strong_count_non_atomic + 1)) {
          return true;
        }
      }
#else
      // the non-thread safe version - this fails with threads because
      // strong_count_ can become 0 after the check if it is 0 and we would
      // return true with no valid object
      if (count_[RCP_STRONG] == 0) {
        return false;
      }
      else {
        ++count_[RCP_STRONG];
        return true;
      }
#endif
    }
  /** \brief . */
  int strong_count() const
    {
      return count_[RCP_STRONG];
    }
  /** \brief . */
  int weak_count() const // not atomically safe
    {
      return count_[RCP_WEAK] - (count_[RCP_STRONG] ? 1 : 0 ); // weak is +1 when strong > 0
    }
  /** \brief . */
  void incr_count( const ERCPStrength strength )
    {
      debugAssertStrength(strength);
      if (++count_[strength] == 1) {
        if (strength == RCP_STRONG) {
          ++count_[RCP_WEAK]; // this is the special condition - the first strong creates a weak
        }
      }
    }
  /** \brief . */
  int deincr_count( const ERCPStrength strength )
    {
      debugAssertStrength(strength);
#ifdef BREAK_THREAD_SAFETY_OF_DEINCR_COUNT
      --count_[strength];
      return count_[strength];  // not atomically valid
#else
      return --count_[strength];
#endif
    }
  /** \brief . */
  void has_ownership(bool has_ownership_in)
    {
      has_ownership_ = has_ownership_in;
    }
  /** \brief . */
  bool has_ownership() const
    {
      return has_ownership_;
    }
  /** \brief . */
  void set_extra_data(
    const any &extra_data, const std::string& name,
    EPrePostDestruction destroy_when, bool force_unique );
  /** \brief . */
  any& get_extra_data( const std::string& type_name,
    const std::string& name );
  /** \brief . */
  const any& get_extra_data( const std::string& type_name,
    const std::string& name
    ) const
    {
      return const_cast<RCPNode*>(this)->get_extra_data(type_name, name);
    }
  /** \brief . */
  any* get_optional_extra_data(const std::string& type_name,
    const std::string& name );
  /** \brief . */
  const any* get_optional_extra_data(
    const std::string& type_name, const std::string& name
    ) const
    {
      return const_cast<RCPNode*>(this)->get_optional_extra_data(type_name, name);
    }
  /** \brief . */
  virtual bool is_valid_ptr() const = 0;
  /** \brief . */
  virtual void delete_obj() = 0;
  /** \brief . */
  virtual void throw_invalid_obj_exception(
    const std::string& rcp_type_name,
    const void* rcp_ptr,
    const RCPNode* rcp_node_ptr,
    const void* rcp_obj_ptr
    ) const = 0;
  /** \brief . */
  virtual const std::string get_base_obj_type_name() const = 0;
#ifdef TEUCHOS_DEBUG
  /** \brief . */
  virtual const void* get_base_obj_map_key_void_ptr() const = 0;
#endif
protected:
  /** \brief . */
  void pre_delete_extra_data()
    {
      if(extra_data_map_)
        impl_pre_delete_extra_data();
    }
private:
  struct extra_data_entry_t {
    extra_data_entry_t() : destroy_when(POST_DESTROY) {}
    extra_data_entry_t( const any &_extra_data, EPrePostDestruction _destroy_when )
      : extra_data(_extra_data), destroy_when(_destroy_when)
      {}
    any extra_data;
    EPrePostDestruction destroy_when;
  };
  typedef Teuchos::map<std::string,extra_data_entry_t> extra_data_map_t;

  TEUCHOS_RCP_DECL_ATOMIC(count_[2], int);
  TEUCHOS_RCP_DECL_ATOMIC(has_ownership_, bool);

  extra_data_map_t *extra_data_map_;
  // Above is made a pointer to reduce overhead for the general case when this
  // is not used.  However, this adds just a little bit to the overhead when
  // it is used.
  // Provides the "basic" guarantee!
  void impl_pre_delete_extra_data();
  // Not defined and not to be called
  RCPNode();
  RCPNode(const RCPNode&);
  RCPNode& operator=(const RCPNode&);
#ifdef TEUCHOS_DEBUG
  // removed atomic because mutex handles it - atomic would be redundant
  int insertion_number_;
public:
  void set_insertion_number(int insertion_number_in)
    {
      insertion_number_ = insertion_number_in;
    }
  int insertion_number() const
    {
      return insertion_number_;
    }
#endif // TEUCHOS_DEBUG
};


/** \brief Throw that a pointer passed into an RCP object is null.
 *
 * \relates RCPNode
 */
TEUCHOSCORE_LIB_DLL_EXPORT void throw_null_ptr_error( const std::string &type_name );


#ifdef TEUCHOS_DEBUG
  // to fully implement abort for TEUCHOS_STANDARD_CATCH_STATEMENTS in the cpp
  /** \brief handle std::exception when deleting RCP node. */
  TEUCHOSCORE_LIB_DLL_EXPORT void abort_for_exception_in_destructor(const std::exception &);
  /** \brief handle const int &excpt_code when deleting RCP node. */
  TEUCHOSCORE_LIB_DLL_EXPORT void abort_for_exception_in_destructor(const int &);
  /** \brief handle unknown exception when deleting RCP node. */
  TEUCHOSCORE_LIB_DLL_EXPORT void abort_for_exception_in_destructor();
  // called when RCPNode detects any exception in a destructor
  #define TEUCHOS_CATCH_AND_ABORT                                                   \
  catch(const std::exception &excpt) { abort_for_exception_in_destructor(excpt); }  \
  catch(const int &excpt_code) { abort_for_exception_in_destructor(excpt_code); }   \
  catch(...) { abort_for_exception_in_destructor(); }
#endif

/** \brief Debug-mode RCPNode tracing class.
 *
 * This is a static class that is used to trace all RCP nodes that are created
 * and destroyed and to look-up RCPNodes given an an object's address.  This
 * database is used for several different types of debug-mode runtime checking
 * including a) the detection of cicular references, b) detecting the creation
 * of duplicate owning RCPNode objects for the same reference-counted object,
 * and c) to create weak RCP objects for existing RCPNode objects.
 *
 * This is primarily an internal implementation class but there are a few
 * functions (maked as such below) that can be called by general users to turn
 * on and off node tracing and to print the active RCPNode objects at any
 * time.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class TEUCHOSCORE_LIB_DLL_EXPORT RCPNodeTracer {
public:

  /** \name Public types. */
  //@{

  /** \brief RCP statistics struct. */
  struct RCPNodeStatistics {
    RCPNodeStatistics()
      : maxNumRCPNodes(0), totalNumRCPNodeAllocations(0),
        totalNumRCPNodeDeletions(0)
      {}
    long int maxNumRCPNodes;
    long int totalNumRCPNodeAllocations;
    long int totalNumRCPNodeDeletions;
  };

  //@}

  /** \name General user functions (can be called by any client) */
  //@{

  /** \brief Return if we are tracing active nodes or not.
   *
   * NOTE: This will always return <tt>false</tt> when <tt>TEUCHOS_DEBUG</tt> is
   * not defined.
   */
  static bool isTracingActiveRCPNodes();

#if defined(TEUCHOS_DEBUG) && !defined(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING)
  /** \brief Set if will be tracing active RCP nodes.
   *
   * This will only cause tracing of RCPNode-based objects that are created
   * after this has been called with <tt>true</tt>.  This function can later be
   * called with <tt>false</tt> to turn off tracing RCPNode objects.  This can
   * allow the client to keep track of RCPNode objects that get created in
   * specific blocks of code and can help as a debugging aid.
   *
   * NOTE: This function call will not even compile unless
   * <tt>defined(TEUCHOS_DEBUG) &&
   * !defined(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING)<tt>.  This occurs when
   * Teuchos is configured with <tt>Teuchos_ENABLE_DEBUG=ON</tt> but with
   * <tt>Teuchos_ENABLE_DEBUG_RCP_NODE_TRACING=OFF</tt>.  The motivation for
   * this logic is that if debug support is not enabled, then node tracing can
   * not be enbled.  If node tracing is turned on at configure time then there
   * is no sense in trying to turn it on or off.
   */
  static void setTracingActiveRCPNodes(bool tracingActiveNodes);
#endif

  /** \brief Print the number of active RCPNode objects currently being
   * tracked.
   */
  static int numActiveRCPNodes();

  /** \brief Return the statistics on RCPNode allocations. */
  static RCPNodeStatistics getRCPNodeStatistics() ;

  /** \brief Print the RCPNode allocation statistics. */
  static void printRCPNodeStatistics(
    const RCPNodeStatistics& rcpNodeStatistics, std::ostream &out);

  /** \brief Set if RCPNode usage statistics will be printed when the program
   * ends or not.
   */
  static void setPrintRCPNodeStatisticsOnExit(
    bool printRCPNodeStatisticsOnExit);

  /** \brief Return if RCPNode usage statistics will be printed when the
   * program ends or not.
   */
  static bool getPrintRCPNodeStatisticsOnExit();

  /** \brief Set if printActiveRCPNodes() is called on exit from the
   * program.
   */
  static void setPrintActiveRcpNodesOnExit(bool printActiveRcpNodesOnExit);

  /** \brief Return if printActiveRCPNodes() is called on exit from the
   * program.
   */
  static bool getPrintActiveRcpNodesOnExit();

  /** \brief Print the list of currently active RCP nodes.
   *
   * When the macro <tt>TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE</tt> is
   * defined, this function will print out all of the RCP nodes that are
   * currently active.  This function can be called at any time during a
   * program.
   *
   * When the macro <tt>TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE</tt> is
   * defined this function will get called automatically after the program
   * ends by default and all of the local and global RCP objects have been
   * destroyed.  If any RCP nodes are printed at that time, then this is an
   * indication that there may be some circular references that will caused
   * memory leaks.  You memory checking tool such as valgrind or purify should
   * complain about this!
   */
  static void printActiveRCPNodes(std::ostream &out);

  //@}

  /** \name Internal implementation functions (not to be called by general
   * clients).
   */
  //@{

  /** \brief Add new RCPNode to the global list.
   *
   * Only gets called when RCPNode tracing has been activated.
   */
  static void addNewRCPNode(RCPNode* rcp_node,
    const std::string &info );

  /** \brief Remove an RCPNode from global list.
   *
   * Always gets called in a debug build (<tt>TEUCHOS_DEBUG</tt> defined) when
   * node tracing is enabled.
   */
  static void removeRCPNode( RCPNode* rcp_node );

  /** \brief Get a <tt>const void*</tt> address to be used as the lookup key
   * for an RCPNode given its embedded object's typed pointer.
   *
   * This only returns the base address reliabily for all types if
   * <tt>HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR</tt> is defined.
   *
   * \returns <tt>0</tt> if <tt>p == 0</tt>.
   */
  template<class T>
  static const void* getRCPNodeBaseObjMapKeyVoidPtr(T *p)
    {
#ifdef HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR
      return getBaseObjVoidPtr(p);
#else
      // This will not return the base address for polymorphic types if
      // multiple inheritance and/or virtual bases are used but returning the
      // static_cast should be okay in how it is used.  It is just that the
      // RCPNode tracing support will not always be able to figure out if two
      // pointers of different type are pointing to the same object or not.
      return static_cast<const void*>(p);
#endif
    }

  /** \brief Return a raw pointer to an existing owning RCPNode given its
   * lookup key.
   *
   * \returns <tt>returnVal != 0</tt> if an owning RCPNode exists, <tt>0</tt>
   * otherwsise.
   */
  static RCPNode* getExistingRCPNodeGivenLookupKey(
    const void* lookupKey);

  /** \brief Return a raw pointer to an existing owning RCPNode given the
   * address to the underlying object if it exits.
   *
   * \returns <tt>returnVal != 0</tt> if an owning RCPNode exists, <tt>0</tt>
   * otherwsise.
   */
  template<class T>
  static RCPNode* getExistingRCPNode(T *p)
    {
      return getExistingRCPNodeGivenLookupKey(getRCPNodeBaseObjMapKeyVoidPtr(p));
    }

  /** \brief Header string used in printActiveRCPNodes(). */
  static std::string getActiveRCPNodeHeaderString();

  /** \brief Common error message string on how to debug RCPNode problems. */
  static std::string getCommonDebugNotesString();

  //@}

};


#ifdef TEUCHOS_DEBUG
#  define TEUCHOS_RCP_INSERION_NUMBER_STR() \
      "  insertionNumber:      " << rcp_node_ptr->insertion_number() << "\n"
#else
#  define TEUCHOS_RCP_INSERION_NUMBER_STR()
#endif


/** \brief Templated implementation class of <tt>RCPNode</tt> that has the
 * responsibility for deleting the reference-counted object.
 *
 * \ingroup teuchos_mem_mng_grp
 */
template<class T, class Dealloc_T>
class RCPNodeTmpl : public RCPNode {
public:
  /** \brief For defined types. */
  RCPNodeTmpl(T* p, Dealloc_T dealloc, bool has_ownership_in)
    : RCPNode(has_ownership_in), ptr_(p),
#ifdef TEUCHOS_DEBUG
      base_obj_map_key_void_ptr_(RCPNodeTracer::getRCPNodeBaseObjMapKeyVoidPtr(p)),
      deleted_ptr_(0),
#endif
      dealloc_(dealloc)
    {}
  /** \brief For undefined types . */
  RCPNodeTmpl(T* p, Dealloc_T dealloc, bool has_ownership_in, ENull)
    : RCPNode(has_ownership_in), ptr_(p),
#ifdef TEUCHOS_DEBUG
      base_obj_map_key_void_ptr_(0),
      deleted_ptr_(0),
#endif
      dealloc_(dealloc)
    {}
  /** \brief . */
  Dealloc_T& get_nonconst_dealloc()
    { return dealloc_; }
  /** \brief . */
  const Dealloc_T& get_dealloc() const
    { return dealloc_; }
  /** \brief . */
  ~RCPNodeTmpl()
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_TERMINATION( ptr_!=0,
        "Error, the underlying object must be explicitly deleted before deleting"
        " the node object!" );
#endif
    }
  /** \brief . */
  virtual bool is_valid_ptr() const
    {
      return ptr_ != 0;
    }
  /** \brief Delete the underlying object.
   * Will abort if an exception is detected in the destructor.
   */
  virtual void delete_obj()
    {
      if (ptr_!= 0) {
        this->pre_delete_extra_data(); // Should not throw!
        T* tmp_ptr = ptr_;
#ifdef TEUCHOS_DEBUG
        deleted_ptr_ = tmp_ptr;
#endif
        ptr_ = 0;
        if (has_ownership()) {
#ifdef TEUCHOS_DEBUG
          try {
#endif
            dealloc_.free(tmp_ptr);
#ifdef TEUCHOS_DEBUG
          }
          TEUCHOS_CATCH_AND_ABORT
#endif
        }
      }
    }
  /** \brief . */
  virtual void throw_invalid_obj_exception(
    const std::string& rcp_type_name,
    const void* rcp_ptr,
    const RCPNode* rcp_node_ptr,
    const void* rcp_obj_ptr
    ) const
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG( ptr_!=0, "Internal coding error!" );
      const T* deleted_ptr =
#ifdef TEUCHOS_DEBUG
        deleted_ptr_
#else
        0
#endif
        ;
      TEUCHOS_ASSERT(rcp_node_ptr);
      TEUCHOS_TEST_FOR_EXCEPTION( true, DanglingReferenceError,
        "Error, an attempt has been made to dereference the underlying object\n"
        "from a weak smart pointer object where the underling object has already\n"
        "been deleted since the strong count has already gone to zero.\n"
        "\n"
        "Context information:\n"
        "\n"
        "  RCP type:             " << rcp_type_name << "\n"
        "  RCP address:          " << rcp_ptr << "\n"
        "  RCPNode type:         " << typeName(*this) << "\n"
        "  RCPNode address:      " << rcp_node_ptr << "\n"
        TEUCHOS_RCP_INSERION_NUMBER_STR()
        "  RCP ptr address:      " << rcp_obj_ptr << "\n"
        "  Concrete ptr address: " << deleted_ptr << "\n"
        "\n"
        << RCPNodeTracer::getCommonDebugNotesString()
        );
      // 2008/09/22: rabartl: Above, we do not provide the concreate object
      // type or the concrete object address.  In the case of the concrete
      // object address, in a non-debug build, we don't want to pay a price
      // for extra storage that we strictly don't need.  In the case of the
      // concrete object type name, we don't want to force non-debug built
      // code to have the require that types be fully defined in order to use
      // the memory management software.  This is related to bug 4016.

    }
  /** \brief . */
  const std::string get_base_obj_type_name() const
    {
#ifdef TEUCHOS_DEBUG
      return TypeNameTraits<T>::name();
#else
      return "UnknownType";
#endif
    }
#ifdef TEUCHOS_DEBUG
  /** \brief . */
  const void* get_base_obj_map_key_void_ptr() const
    {
      return base_obj_map_key_void_ptr_;
    }
#endif
private:
  T *ptr_;
#ifdef TEUCHOS_DEBUG
  const void *base_obj_map_key_void_ptr_;
  T *deleted_ptr_;
#endif
  Dealloc_T dealloc_;
  // not defined and not to be called
  RCPNodeTmpl();
  RCPNodeTmpl(const RCPNodeTmpl&);
  RCPNodeTmpl& operator=(const RCPNodeTmpl&);

}; // end class RCPNodeTmpl<T>


/** \brief Sets up node tracing and prints remaining RCPNodes on destruction.
 *
 * This class is used by automataic code that sets up support for RCPNode
 * tracing and for printing of remaining nodes on destruction.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class TEUCHOSCORE_LIB_DLL_EXPORT ActiveRCPNodesSetup {
public:
  /** \brief . */
  ActiveRCPNodesSetup();
  /** \brief . */
  ~ActiveRCPNodesSetup();
  /** \brief . */
  void foo();
private:
  static int count_;
};


} // namespace Teuchos


namespace {
// This static variable is declared before all other static variables that
// depend on RCP or other classes. Therefore, this static variable will be
// deleted *after* all of these other static variables that depend on RCP or
// created classes go away!  This ensures that the node tracing machinery is
// setup and torn down correctly (this is the same trick used by the standard
// stream objects in many compiler implementations).
Teuchos::ActiveRCPNodesSetup local_activeRCPNodesSetup;
} // namespace (anonymous)


namespace Teuchos {

/** \class RCPNodeHandle
 * \brief Handle class that manages the RCPNode's reference counting.
 *
 * \warning This class is <i>not</i> intended for Teuchos users.  It
 *   is an implementation detail of Teuchos' reference-counting
 *   "smart" pointer (RCP) and array (ArrayRCP) classes.
 *
 * NOTE: I (Ross Bartlett) am not generally a big fan of handle classes and
 * greatly prefer smart pointers.  However, this is one case where a handle
 * class makes sense.  First, I want special behavior in some functions when
 * the wrapped RCPNode pointer is null.  Second, I can't use one of the
 * smart-pointer classes because this class is used to implement all of those
 * smart-pointer classes!
 *
 * \ingroup teuchos_mem_mng_grp
 */
class TEUCHOSCORE_LIB_DLL_EXPORT RCPNodeHandle {
public:
  //! Default constructor
  RCPNodeHandle (ENull null_arg = null)
    : node_ (0), strength_ (RCP_STRONG)
  {
    (void) null_arg; // Silence "unused variable" compiler warning.
  }

  //! Constructor that takes a pointer to an RCPNode.
  RCPNodeHandle (RCPNode* node,
                 ERCPStrength strength_in = RCP_STRONG,
                 bool newNode = true)
    : node_ (node), strength_ (strength_in)
  {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT(node);
#endif // TEUCHOS_DEBUG

    bind();

#ifdef TEUCHOS_DEBUG
    // Add the node if this is the first RCPNodeHandle to get it.  We have
    // to add it because unbind() will call the remove_RCPNode(...) function
    // and it needs to match when node tracing is on from the beginning.
    if (RCPNodeTracer::isTracingActiveRCPNodes() && newNode) {
      std::ostringstream os;
      os << "{T=Unknown, ConcreteT=Unknown, p=Unknown,"
         << " has_ownership="<<node_->has_ownership()<<"}";
      RCPNodeTracer::addNewRCPNode(node_, os.str());
    }
#else
    (void) newNode; // Silence "unused variable" compiler warning.
#endif // TEUCHOS_DEBUG
  }

#ifdef TEUCHOS_DEBUG
  /** \brief Only gets called in debug mode. */
  template<typename T>
  RCPNodeHandle (RCPNode* node, T *p, const std::string &T_name,
                 const std::string &ConcreteT_name,
                 const bool has_ownership_in,
                 ERCPStrength strength_in = RCP_STRONG)
    : node_ (node), strength_ (strength_in)
  {
    TEUCHOS_ASSERT(strength_in == RCP_STRONG); // Can't handle weak yet!
    TEUCHOS_ASSERT(node_);
    bind();
    if (RCPNodeTracer::isTracingActiveRCPNodes()) {
      std::ostringstream os;
      os << "{T="<<T_name<<", ConcreteT="<< ConcreteT_name
         <<", p="<<static_cast<const void*>(p)
         <<", has_ownership="<<has_ownership_in<<"}";
      RCPNodeTracer::addNewRCPNode(node_, os.str());
    }
  }
#endif // TEUCHOS_DEBUG

  //! Copy constructor
  RCPNodeHandle (const RCPNodeHandle& node_ref)
    : node_ (node_ref.node_), strength_ (node_ref.strength_)
  {
    bind();
  }

  //! Move constructor
  RCPNodeHandle (RCPNodeHandle&& node_ref)
    : node_ (node_ref.node_), strength_ (node_ref.strength_)
  {
    node_ref.node_ = 0;
    node_ref.strength_ = RCP_STRONG;
  }

  //! Swap the contents of \c node_ref with \c *this.
  void swap (RCPNodeHandle& node_ref) {
    std::swap (node_ref.node_, node_);
    std::swap (node_ref.strength_, strength_);
  }



  /// \brief Null assignment.
  ///
  /// This method satisfies the strong exception guarantee: It either
  /// returns successfully, or throws an exception without modifying
  /// any user-visible state.
  RCPNodeHandle& operator= (ENull) {
    unbind(); // May throw in some cases
    node_ = 0;
    strength_ = RCP_STRONG;
    return *this;
  }

  /// \brief Copy assignment operator.
  ///
  /// This method satisfies the strong exception guarantee: It either
  /// returns successfully, or throws an exception without modifying
  /// any user-visible state.
  RCPNodeHandle& operator= (const RCPNodeHandle& node_ref) {
    // NOTE: Don't need to check assignment to self since user-facing classes
    // do that!
    unbind(); // May throw in some cases
    node_ = node_ref.node_;
    strength_ = node_ref.strength_;
    bind();
    return *this;
  }

  /// \brief Move assignment operator.
  ///
  /// This method satisfies the strong exception guarantee: It either
  /// returns successfully, or throws an exception without modifying
  /// any user-visible state.
  RCPNodeHandle& operator= (RCPNodeHandle&& node_ref) {
    // NOTE: Don't need to check assignment to self since user-facing classes
    // do that!
    unbind(); // May throw in some cases
    node_ = node_ref.node_;
    strength_ = node_ref.strength_;
    node_ref.node_ = 0;
    node_ref.strength_ = RCP_STRONG;
    return *this;
  }

  //! Destructor
  ~RCPNodeHandle() {
    unbind();
  }

  //! Return a strong handle if possible using thread safe atomics
  // otherwise return a null handle
  RCPNodeHandle create_strong_lock() const {
    // make weak handle
    RCPNodeHandle possibleStrongNode(node_, RCP_WEAK, false);
    if (possibleStrongNode.attemptConvertWeakToStrong()) {
      return possibleStrongNode; // success - we have a good strong handle
    }
    return RCPNodeHandle(); // failure - return an empty handle
  }

  //! Return a weak handle.
  RCPNodeHandle create_weak() const {
    if (node_) {
      return RCPNodeHandle(node_, RCP_WEAK, false);
    }
    return RCPNodeHandle();
  }
  //! Return a strong handle.
  RCPNodeHandle create_strong() const {
    if (node_) {
      return RCPNodeHandle(node_, RCP_STRONG, false);
    }
    return RCPNodeHandle();
  }
  //! Return a pointer to the underlying RCPNode.
  RCPNode* node_ptr() const {
    return node_;
  }
  //! Whether the underlying RCPNode is NULL.
  bool is_node_null() const {
    return node_==0;
  }
  /// \brief Whether the underlying pointer is valid.
  ///
  /// \note NULL is a valid pointer; this method returns true in that case.
  bool is_valid_ptr() const {
    if (node_) {
      return node_->is_valid_ptr();
    }
    return true; // Null is a valid ptr!
  }
  /// \brief Whether the RCPNode for which \c node2 is a handle is the
  ///   same RCPNode as this object's RCPNode.
  bool same_node(const RCPNodeHandle &node2) const {
    return node_ == node2.node_;
  }
  //! The strong count for this RCPNode, or 0 if the node is NULL.
  int strong_count() const {
    if (node_) {
      return node_->strong_count();
    }
    return 0;
  }
  //! The weak count for this RCPNode, or 0 if the node is NULL.
  int weak_count() const {
    if (node_) {
      return node_->weak_count(); // Not atomically safe
    }
    return 0;
  }
  //! The sum of the weak and string counts.
  int total_count() const {
    if (node_) {
      return node_->strong_count() + node_->weak_count(); // not atomically safe
    }
    return 0;
  }
  //! The strength of this handle.
  ERCPStrength strength() const {
    return strength_;
  }
  /** \brief . */
  void has_ownership(bool has_ownership_in)
    {
      if (node_)
        node_->has_ownership(has_ownership_in);
    }
  /** \brief . */
  bool has_ownership() const
    {
      if (node_)
        return node_->has_ownership();
      return false;
    }
  /** \brief . */
  void set_extra_data(
    const any &extra_data, const std::string& name,
    EPrePostDestruction destroy_when, bool force_unique
    )
    {
      debug_assert_not_null();
      node_->set_extra_data(extra_data, name, destroy_when, force_unique);
    }
  /** \brief . */
  any& get_extra_data( const std::string& type_name,
    const std::string& name
    )
    {
      debug_assert_not_null();
      return node_->get_extra_data(type_name, name);
    }
  /** \brief . */
  const any& get_extra_data( const std::string& type_name,
    const std::string& name
    ) const
    {
      return const_cast<RCPNodeHandle*>(this)->get_extra_data(type_name, name);
    }
  /** \brief . */
  any* get_optional_extra_data(
    const std::string& type_name, const std::string& name
    )
    {
      debug_assert_not_null();
      return node_->get_optional_extra_data(type_name, name);
    }
  /** \brief . */
  const any* get_optional_extra_data(
    const std::string& type_name, const std::string& name
    ) const
    {
      return const_cast<RCPNodeHandle*>(this)->get_optional_extra_data(type_name, name);
    }
  /** \brief . */
  void debug_assert_not_null() const
    {
#ifdef TEUCHOS_DEBUG
      if (!node_)
        throw_null_ptr_error(typeName(*this));
#endif
    }
  /** \brief . */
  template<class RCPType>
  void assert_valid_ptr(const RCPType& rcp_obj) const
    {
      if (!node_)
        return; // Null is a valid pointer!
      if (!is_valid_ptr()) {
        node_->throw_invalid_obj_exception( typeName(rcp_obj),
          this, node_, rcp_obj.access_private_ptr() );
      }
    }
  /** \brief . */
  template<class RCPType>
  void debug_assert_valid_ptr(const RCPType& rcp_obj) const
    {
#ifdef TEUCHOS_DEBUG
      assert_valid_ptr(rcp_obj);
#endif
    }
#ifdef TEUCHOS_DEBUG
  const void* get_base_obj_map_key_void_ptr() const
    {
      if (node_)
        return node_->get_base_obj_map_key_void_ptr();
      return 0;
    }
#endif
private:
  RCPNode *node_;
  ERCPStrength strength_;
  // atomically safe conversion of a weak handle to a strong handle if
  // possible - if not possible nothing changes
  bool attemptConvertWeakToStrong() {
    if (node_->attemptIncrementStrongCountFromNonZeroValue()) {
      // because we converted strong + 1 we account for this by doing weak - 1
      node_->deincr_count(RCP_WEAK);
      // we have successfully incremented the strong count by one
      strength_ = RCP_STRONG;
      return true;
    }
    return false;
  }
  inline void bind()
    {
      if (node_)
        node_->incr_count(strength_);
    }
  inline void unbind()
    {
      if (node_) {
        if(strength_ == RCP_STRONG) {
          // only strong checks for --strong == 0
          if (node_->deincr_count(RCP_STRONG) == 0) {
            unbindOneStrong();
            // but if strong hits 0 it also decrements weak_count_plus which
            // is weak + (strong != 0)
            if( node_->deincr_count(RCP_WEAK) == 0) {
              unbindOneTotal();
            }
          }
        }
        else if(node_->deincr_count(RCP_WEAK) == 0) {  // weak checks here
          unbindOneTotal();
        }
      }
    }
  void unbindOneStrong();
  void unbindOneTotal();
};


/** \brief Ouput stream operator for RCPNodeHandle.
 *
 * \relates RCPNodeHandle
 */
inline
std::ostream& operator<<(std::ostream& out, const RCPNodeHandle& node)
{
  // mfh 15 Sep 2015: Make sure that NULL pointers print consistently.
  // Clang 3.5 likes to print an empty string in that case, while GCC
  // prints 0.  Thus, we test if the pointer is NULL and print 0 in
  // that case.  This is important for MueLu tests, which compare
  // string print-outs.
  if (node.node_ptr () == NULL) {
    out << "0";
  } else {
    out << node.node_ptr ();
  }
  return out;
}


/** \brief Deletes a (non-owning) RCPNode but not it's underlying object in
 * case of a throw.
 *
 * This class is used in contexts where RCPNodeTracer::addNewRCPNode(...)
 * might thrown an exception for a duplicate node being added.  The assumption
 * is that there must already be an owning (or non-owning) RCP object that
 * will delete the underlying object and therefore this class should *not*
 * call delete_obj()!
 */
class TEUCHOSCORE_LIB_DLL_EXPORT RCPNodeThrowDeleter {
public:
  /** \brief . */
  RCPNodeThrowDeleter(RCPNode *node)
    : node_(node)
    {}
  /** \brief Called with node_!=0 when an exception is thrown.
   *
   * When an exception is not thrown, the client should have called release()
   * before this function is called.
   */
  ~RCPNodeThrowDeleter()
    {
      if (node_) {
        node_->has_ownership(false); // Avoid actually deleting ptr_
        node_->delete_obj(); // Sets the pointer ptr_=0 to allow RCPNode delete
        delete node_;
      }
    }
  /** \brief . */
  RCPNode* get() const
    {
      return node_;
    }
  /** \brief Releaes the RCPNode pointer before the destructor is called. */
  void release()
    {
      node_ = 0;
    }
private:
  RCPNode *node_;
  RCPNodeThrowDeleter(); // Not defined
  RCPNodeThrowDeleter(const RCPNodeThrowDeleter&); // Not defined
  RCPNodeThrowDeleter& operator=(const RCPNodeThrowDeleter&); // Not defined
};


//
// Unit testing support
//


#if defined(TEUCHOS_DEBUG) && !defined(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING)

class SetTracingActiveNodesStack {
public:
  SetTracingActiveNodesStack()
    {RCPNodeTracer::setTracingActiveRCPNodes(true);}
  ~SetTracingActiveNodesStack()
    {RCPNodeTracer::setTracingActiveRCPNodes(false);}
};

#  define SET_RCPNODE_TRACING() Teuchos::SetTracingActiveNodesStack setTracingActiveNodesStack;

#else

#  define SET_RCPNODE_TRACING() (void)0

#endif // defined(TEUCHOS_DEBUG) && !defined(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING)


} // end namespace Teuchos


#endif // TEUCHOS_RCP_NODE_HPP
