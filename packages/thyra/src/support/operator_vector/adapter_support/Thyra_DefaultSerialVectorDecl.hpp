// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_SERIAL_VECTOR_STD_DECL_HPP
#define THYRA_SERIAL_VECTOR_STD_DECL_HPP

#include "Thyra_SerialVectorBaseDecl.hpp"

namespace Thyra {

/** \brief General concrete <tt>VectorBase</tt> implementation of serial
 * vectors.
 *
 * This class can be used either as a view of a vector data or as a
 * storage for vector data (with any underlying storage type).
 *
 * To create with storage with the dimension of <tt>dim</tt> just call the
 * constructor <tt>%DefaultSerialVector(dim)</tt> or after construction you can
 * call <tt>this->initialize(dim)</tt>.
 *
 * To simply create a view of a vector <tt>v</tt> with stride <tt>vs</tt>,
 * without ownership just call
 * <tt>%DefaultSerialVector(Teuchos::rcp(v,false),vs)</tt> or after construction
 * call <tt>this->initialize(Teuchos::rcp(v,false),vs)</tt>.
 *
 * To use another storage type, such as an <tt>std::vector<Scalar></tt>,
 * construct as:
 *
 \code

 template<class Scalar>
 Teuchos::RefCountPtr<VectorBase<Scalar> > STLVectorSpace<Scalar>::createMember() const
 {
   Teuchos::RefCountPtr<std::vector<Scalar> > stl_v = Teuchos::rcp( new std::vector<Scalar>(dim_) );
   Teuchos::RefCountPtr<Scalar> v = Teuchos::rcp(&(*stl_v)[0],false);
   Teuchos::set_extra_data( stl_v, "stl::vector", &v );
   return Teuchos::rcp( new DefaultSerialVector<Scalar>( v, 1, dim_, Teuchos::rcp(this,false) ) );
 }

 \endcode
 *
 * \ingroup Thyra_Op_Vec_adapters_serial_concrete_std_grp
 */
template<class Scalar>
class DefaultSerialVector : public SerialVectorBase<Scalar> {
public:

  /** \brief . */
  using SerialVectorBase<Scalar>::describe;

  /** @name Constructors/initializers */
  //@{

  /** \brief Calls <tt>this->initialize(vecSpc)</tt>.
   */
  DefaultSerialVector(
    const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vecSpc
    );
  /** \brief Calls <tt>this->initialize(dim)</tt>.
   */
  DefaultSerialVector(
    const Index dim = 0
    );
  /** \brief Calls <tt>this->initialize(v,vs,dim,vecSpc)</tt>.
   */
  DefaultSerialVector(
    const Teuchos::RefCountPtr<Scalar>                          &v
    ,const Index                                                vs
    ,const Index                                                dim
    ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vecSpc = Teuchos::null
    );
  /** \brief Call <tt>this->initialize(v,vs,vecSpc)</tt> with internally dynamically allocated data <tt>v</tt>.
   */
  void initialize(
    const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vecSpc
    );
  /** \brief Call <tt>this->initialize(v,vs,true)</tt> with internally dynamically allocated data <tt>v</tt>.
   */
  void initialize(
    const Index dim
    );
  /** \brief Initialize with storage.
   *
   * @param  v      [in] Smart pointer to array of storage that <tt>*this</tt> will represent.
   * @param  vs     [in] Stride for the storage in <tt>v[]</tt> (see Postconditions).
   * @param  dim    [in] Number of elements in <tt>v[]</tt> this this will represent (see Postconditions).
   * @param  vecSpc
   *                [in] Smart pointer to a <tt>VectorSpaceBase</tt> object that will be used to represent the
   *                vector space for <tt>*this</tt>.  If <tt>vecSpc.get()==NULL</tt> on input, then
   *                a <tt>SerialVectorSpace</tt> object of dimension <tt>dim</tt> is allocated for this
   *                role.  The default is <tt>Teuchos::null</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li> [<tt>vecSpc.get()!=NULL</tt>] <tt>dim == vecSpc->dim()</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>vecSpc.get()!=NULL</tt>] <tt>vecSpc->createMember()</tt> must create vectors that are compatible
   *      with <tt>*this</tt> (i.e. <tt>acquireDetachedView()</tt>, <tt>commitDetachedView()</tt> behave the same as with
   *      this class).
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> [<tt>vecSpc.get()!=NULL</tt>] <tt>vecSpc.get() == this->space().get()</tt>
   * <li> [<tt>vecSpc.get()==NULL</tt>] <tt>dynamic_cast<const DefaultSerialVectorSpace<Scalar>*>(this->space().get()) != NULL</tt>
   * <li> <tt>this->space()->dim() == dim</tt>
   * <li> <tt>this->getRCPtr().get() == v.get()</tt>
   * <li> <tt>this->getPtr() == v.get()</tt>
   * <li> <tt>this->getStride() == vs</tt>
   * </ul>
   *
   * Note that this function is declared virtual so that subclasses
   * can override it to be informed whenever <tt>*this</tt> vector
   * is resized.  An override should call this function as
   * <tt>this->DefaultSerialVector<Scalar>::initialize(...)</tt>.
   */
  virtual void initialize(
    const Teuchos::RefCountPtr<Scalar>                          &v
    ,const Index                                                vs
    ,const Index                                                dim
    ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vecSpc = Teuchos::null
    );

  //@}

  /** @name Accessors (inlined for minimal overhead) */
  //@{

  /** \brief Return non-<tt>const</tt> smart pointer to underlying data */
  Teuchos::RefCountPtr<Scalar> getRCPtr();
  /** \brief Return <tt>const</tt> smart pointer to underlying data */
  Teuchos::RefCountPtr<const Scalar> getRCPtr() const;
  /** \brief Return non-<tt>const</tt> raw pointer to underlying data */
  Scalar* getPtr();
  /** \brief Return <tt>const</tt> raw pointer to underlying data */
  const Scalar* getPtr() const;
  /** \brief Return the stride between entries returned from <tt>getPtr()</tt> */
  Index getStride() const;
  /** \brief Inline call to return dimension */
  Index getDim() const;
  
  //@}

  /** @name Overridden from SerialVectorBase */
  //@{
  /** \brief . */
  void getData( Scalar** values, Index* stride );
  /** \brief . */
  void commitData( Scalar** values );
  /** \brief . */
  void getData( const Scalar** values, Index* stride ) const;
  /** \brief . */
  void freeData( const Scalar** values ) const;
  //@}

  /** @name Overridden from VectorBase */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > space() const;
  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

private:

  // ///////////////////////////////////////
  // Private data members
  
  Teuchos::RefCountPtr<Scalar>                            v_;
  Index                                                   vs_;
  Index                                                   dim_;
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >    space_serial_;

  // ////////////////////////////////
  // Private member functions

  void free_mem();

  // Not defined and not to be called
  DefaultSerialVector(const DefaultSerialVector&);
  DefaultSerialVector& operator=(const DefaultSerialVector&);

}; // end class DefaultSerialVector

// /////////////////////////////////////////////////////
// Inline members

template<class Scalar>
inline
Teuchos::RefCountPtr<Scalar> DefaultSerialVector<Scalar>::getRCPtr()
{
  return v_;
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const Scalar> DefaultSerialVector<Scalar>::getRCPtr() const
{
  return v_;
}

template<class Scalar>
inline
Scalar* DefaultSerialVector<Scalar>::getPtr()
{
  return v_.get();
}

template<class Scalar>
inline
const Scalar* DefaultSerialVector<Scalar>::getPtr() const
{
  return v_.get();
}

template<class Scalar>
inline
Index DefaultSerialVector<Scalar>::getStride() const
{
  return vs_;
}	

template<class Scalar>
inline
Index DefaultSerialVector<Scalar>::getDim() const
{
  return dim_;
}	

} // end namespace Thyra

#endif // THYRA_SERIAL_VECTOR_STD_DECL_HPP
