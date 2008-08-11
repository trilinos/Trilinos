// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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


namespace Teuchos {


/** \brief Used to specify a pre or post destruction of extra data
 *
 * \ingroup teuchos_mem_mng_grp 
 */
enum EPrePostDestruction { PRE_DESTROY, POST_DESTROY };


/** \brief Node class to keep track of the delete address and the reference
 * count for a reference-counted utility class.
 *
 * This is not a general user-level class.  This is used in the implementation
 * of all of the reference-counting utility classes.
 *
 * \ingroup teuchos_mem_mng_grp 
 */
class RCPNode {
public:
  /** \brief . */
  RCPNode(bool has_ownership_in)
    : count_(1), has_ownership_(has_ownership_in), extra_data_map_(NULL)
    {}
  /** \brief . */
  virtual ~RCPNode()
    {
      if(extra_data_map_)
        delete extra_data_map_;
    }
  /** \brief . */
  int count() const {
    return count_; 
  }
  /** \brief . */
  int incr_count() {
    return ++count_;
  }
  /** \brief . */
  int deincr_count() {
    return --count_;
  }
  /** \brief . */
  void has_ownership(bool has_ownership_in) {
    has_ownership_ = has_ownership_in;
  }
  /** \brief . */
  bool has_ownership() const {
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
    const std::string& name ) const
    {
      return const_cast<RCPNode*>(this)->get_extra_data(type_name,name);
    }
  /** \brief . */
  any* get_optional_extra_data(
    const std::string& type_name, const std::string& name );
  /** \brief . */
  const any* get_optional_extra_data(
    const std::string& type_name, const std::string& name ) const
    {
      return const_cast<RCPNode*>(this)->get_optional_extra_data(type_name,name);
    }
protected:
  /** \brief . */
  void pre_delete_extra_data() {
    if(extra_data_map_)
      impl_pre_delete_extra_data();
  }
private:
  struct extra_data_entry_t {
    extra_data_entry_t() : destroy_when(POST_DESTROY) {}
    extra_data_entry_t( const any &_extra_data, EPrePostDestruction _destroy_when )
      : extra_data(_extra_data), destroy_when(_destroy_when) {}
    any extra_data;
    EPrePostDestruction destroy_when;
  }; 
  typedef Teuchos::map<std::string,extra_data_entry_t> extra_data_map_t;
  int count_;
  bool has_ownership_;
  extra_data_map_t *extra_data_map_;
  // Above is made a pointer to reduce overhead for the general case
  // where this is not used
  void impl_pre_delete_extra_data();
  // Not defined and not to be called
  RCPNode();
  RCPNode(const RCPNode&);
  RCPNode& operator=(const RCPNode&);
}; // end class RCPNode;


/** \brief Implementation class for actually deleting the object.
 *
 * \ingroup teuchos_mem_mng_grp 
 */
template<class T, class Dealloc_T>
class RCPNodeTmpl : public RCPNode {
public:
  /** \brief . */
  RCPNodeTmpl(T* p, Dealloc_T dealloc, bool has_ownership_in)
    : RCPNode(has_ownership_in), ptr_(p), dealloc_(dealloc)
    {}
  /** \brief . */
  Dealloc_T& get_nonconst_dealloc() { return dealloc_; }
  /** \brief . */
  const Dealloc_T& get_dealloc() const { return dealloc_; }
  /** \brief . */
  ~RCPNodeTmpl() {
    this->pre_delete_extra_data();
    if( has_ownership() )
      dealloc_.free(ptr_);
  }
private:
  T *ptr_;
  Dealloc_T dealloc_;
  // not defined and not to be called
  RCPNodeTmpl();
  RCPNodeTmpl(const RCPNodeTmpl&);
  RCPNodeTmpl& operator=(const RCPNodeTmpl&);

}; // end class RCPNodeTmpl<T>


/** \brief Add new RCP to global list.
 *
 * \relates PrintActiveRCPNodes
 */
void add_new_RCPNode( RCPNode* rcp_node, const std::string &info );


/** \brief Remove RCP from global list.
 *
 * \relates PrintActiveRCPNodes
 */
void remove_RCPNode( RCPNode* rcp_node );


/** \brief Print global list on destruction.
 *
 * \ingroup teuchos_mem_mng_grp
 */
class PrintActiveRCPNodes {
public:
  /** \brief . */
  PrintActiveRCPNodes();
  /** \brief . */
  ~PrintActiveRCPNodes();
  /** \brief . */
  void foo();
private:
  static int count_;
};


/** \brief Return if we are tracing active nodes or not.
 *
 * NOTE: This will always return <tt>false</tt> when <tt>TEUCHOS_DEBUG</tt> is
 * not defined.
 *
 * \relates RCPNode
 */
bool isTracingActiveRCPNodes();


#ifdef TEUCHOS_DEBUG

/** \brief Set if we should be tracing active RCP nodes.
 *
 * This will only cause tracing of RCPNode-based objects that are created
 * after this has been called with <tt>true</tt>.  This function can later be
 * called with <tt>false</tt> to turn off tracing RCPNode objects.  This can
 * allow the client to keep track of RCPNode objects that get created in
 * specific blocks of code and can help as a debugging aid.
 *
 * NOTE: This function call will not even compile unless
 * <tt>TEUCHOS_DEBUG</tt> is defined!
 *
 * \relates RCPNode
 */
void setTracingActiveRCPNodes(bool tracingActiveNodes);

#endif // TEUCHOS_DEBUG


/** \brief Print the number of active RCPNode objects being tracked. */
int numActiveRCPNodes();


/** \brief Print the list of currently active RCP nodes.
 *
 * When the macro <tt>TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE</tt> is
 * defined, this function will print out all of the RCP nodes that are
 * currently active.  This function can be called at any time during a
 * program.
 *
 * When the macro <tt>TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE</tt> is
 * defined this function will get called automatically after the program ends
 * and all of the local and global RCP objects have been destroyed.  If any
 * RCP nodes are printed at that time, then this is an indication that there
 * may be some circular references that will caused memory leaks.  You memory
 * checking tool such as valgrind or purify should complain about this!
 *
 * \relates RCPNode
 */
void printActiveRCPNodes(std::ostream &out);


/** \brief Throw that a pointer passed into an RCP object is null. */
void throw_null_ptr_error( const std::string &type_name );


} // end namespace Teuchos


namespace {
// This static variable should be delcared before all other static variables
// that depend on RCP or other classes. Therefore, this static varaible should
// be deleted *after* all of these other static variables that depend on RCP
// or created classes go away!
Teuchos::PrintActiveRCPNodes local_printActiveRCPNodes;
} // namespace


#endif // TEUCHOS_RCP_NODE_HPP
