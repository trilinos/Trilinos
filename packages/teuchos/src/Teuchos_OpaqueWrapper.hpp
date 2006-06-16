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

#ifndef TEUCHOS_OPAQUE_WRAPPER_HPP
#define TEUCHOS_OPAQUE_WRAPPER_HPP

#include "Teuchos_RefCountPtr.hpp"

//#define TEUCHOS_OPAQUE_WRAPPER_ANNOUNCE_FREE

#ifdef TEUCHOS_OPAQUE_WRAPPER_ANNOUNCE_FREE
#  include "Teuchos_VerboseObject.hpp"
#endif // TEUCHOS_OPAQUE_WRAPPER_ANNOUNCE_FREE

namespace Teuchos {

/** \defgroup Teuchos_opaque_wrappers_grp */

/** \brief Base class for wrapped opaque objects.
 *
 * \ingroup Teuchos_opaque_wrappers_grp
 */
template <class Opaque>
class OpaqueWrapper {
public:
  OpaqueWrapper( Opaque opaque )
    : opaque_(opaque)
    {}
  operator Opaque () const
    { return opaque_; }
  Opaque operator()() const
    { return opaque_; }
protected:
  Opaque  opaque_; // Bad in general but ...
private:
  OpaqueWrapper(); // Not defined
  OpaqueWrapper(const OpaqueWrapper&); // Not defined
  OpaqueWrapper& operator=(const OpaqueWrapper&); // Not defined
}; 

/** \brief Subclass for wrapped opaque objects with a free function.
 *
 * \ingroup Teuchos_opaque_wrappers_grp
 */
template <class Opaque, class OpaqueFree>
class OpaqueWrapperWithFree : public OpaqueWrapper<Opaque> {
public:
  OpaqueWrapperWithFree( Opaque opaque, OpaqueFree opaqueFree )
    : OpaqueWrapper<Opaque>(opaque), opaqueFree_(opaqueFree)
    {}
  ~OpaqueWrapperWithFree()
    {
      if(opaqueFree_) {
#ifdef TEUCHOS_OPAQUE_WRAPPER_ANNOUNCE_FREE
        Teuchos::RefCountPtr<Teuchos::FancyOStream>
          out = Teuchos::VerboseObjectBase::getDefaultOStream();
        Teuchos::OSTab tab(out);
        *out << "\nOpaqueWrapperWithFree::~OpaqueWrapperWithFree(): Freeing opaque object"
             << " of type " << typeid(Opaque).name() << "!\n";
#endif // TEUCHOS_OPAQUE_WRAPPER_ANNOUNCE_FREE
        opaqueFree_(&this->opaque_);
      }
    }
private:
  OpaqueFree opaqueFree_;
  OpaqueWrapperWithFree(); // Not defined
  OpaqueWrapperWithFree(const OpaqueWrapperWithFree&); // Not defined
  OpaqueWrapperWithFree& operator=(const OpaqueWrapperWithFree&); // Not defined
}; 

/** \brief Helper function created a new <tt>OpaqueWrapper</tt> object without
 * a free function.
 *
 * \ingroup Teuchos_opaque_wrappers_grp
 */
template <class Opaque>
inline
RefCountPtr<OpaqueWrapper<Opaque> >
opaqueWrapper( Opaque opaque)
{
  return rcp(new OpaqueWrapper<Opaque>(opaque));
}

/** \brief Helper function created a new <tt>OpaqueWrapper</tt> object with a
 * free function.
 *
 * \ingroup Teuchos_opaque_wrappers_grp
 */
template <class Opaque, class OpaqueFree>
inline
RefCountPtr<OpaqueWrapper<Opaque> >
opaqueWrapper( Opaque opaque, OpaqueFree opaqueFree)
{
  return rcp(new OpaqueWrapperWithFree<Opaque,OpaqueFree>(opaque,opaqueFree));
}

} // namespace Teuchos

#endif	// TEUCHOS_OPAQUE_WRAPPER_HPP
