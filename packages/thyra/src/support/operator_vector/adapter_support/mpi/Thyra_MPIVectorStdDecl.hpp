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

#ifndef THYRA_MPI_VECTOR_STD_DECL_HPP
#define THYRA_MPI_VECTOR_STD_DECL_HPP

#include "Thyra_MPIVectorBaseDecl.hpp"

namespace Thyra {

/** \brief Efficient concrete implementation subclass for SPMD-MPI-based vectors.
 *
 * This subclass provides a very efficient and very general concrete
 * implementation of a <tt>Thyra::VectorBase</tt> object.
 *
 * Objects of this type generally should not be constructed directly
 * by a client but instead by using the concrete vector space subclass
 * <tt>Thyra::MPIVectorSpaceStd</tt> and using the function
 * <tt>Thyra::MPIVectorSpaceStd::createMember()</tt>.
 *
 * The storage type can be anything since a
 * <tt>Teuchos::RefCountPtr</tt> is used to pass in the local values
 * pointer into the constructor and <tt>initialize()</tt>.
 *
 * \ingroup Thyra_Op_Vec_adapters_MPI_concrete_std_grp
 */
template<class Scalar>
class MPIVectorStd : virtual public MPIVectorBase<Scalar> {
public:

  /** @name Constructors/initializers */
  //@{


  /// Construct to uninitialized
  MPIVectorStd();

  /// Calls <tt>initialize()</tt>
  MPIVectorStd(
    const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    &mpiSpace
    ,const Teuchos::RefCountPtr<Scalar>                              &localValues
    ,const Index                                                     stride
    );

  /** \brief Initialize.
   *
   * @param  mpiSpace
   *                   [in] Smart pointer to <tt>MPIVectorSpaceBase</tt> object
   *                   that defines the data distribution for <tt>mpiSpace()</tt> and <tt>space()</tt>.
   * @param  localValues
   *                   [in] Smart pointer to beginning of local strided vector data.
   *                   This array must be at least of dimension <tt>mpiRangeSpace->localDim()*stride</tt>
   *                   and <tt>(&*localValues)[ i*stride ]</tt> gives the local value
   *                   of the zero-based entry <tt>(i)</tt> where <tt>i=0...mpiSpace()->localSubDim()-1</tt>.
   * @param  stride   [in] Stride between local vector elements.
   *
   * Preconditions:<ul>
   * <li><tt>mpiSpace.get()!=NULL</tt>
   * <li><tt>localValues.get()!=NULL</tt>
   * <li><tt>stride != 0</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->getRCptr().get() == localValues.get()</tt>
   * <li> <tt>this->getPtr() == &*localValues</tt>
   * <li> <tt>this->getStride() == stride</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    &mpiSpace
    ,const Teuchos::RefCountPtr<Scalar>                              &localValues
    ,const Index                                                     stride
    );

  /** \brief Set to an uninitialized state.
   *
   * Postconditions:<ul>
   * <li><tt>this->mpiSpace().get() == NULL</tt>.
   */
  void uninitialize(
    Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    *mpiSpace      = NULL
    ,Teuchos::RefCountPtr<Scalar>                              *localValues   = NULL
    ,Index                                                     *stride        = NULL
    );

  //@}

  /** @name Accessors (inlined for minimal overhead) */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<Scalar> getRCPtr();
  /** \brief . */
  Teuchos::RefCountPtr<const Scalar> getRCPtr() const;
  /** \brief . */
  Scalar* getPtr();
  /** \brief . */
  const Scalar* getPtr() const;
  /** \brief . */
  Index getStride() const;
  
  //@}

  /** @name Overridden form Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

  /** @name Overridden from MPIVectorBase */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const;
  /** \brief . */
  void getLocalData( Scalar** localValues, Index* stride );
  /** \brief . */
  void commitLocalData( Scalar* localValues );
  /** \brief . */
  void getLocalData( const Scalar** localValues, Index* stride ) const;
  /** \brief . */
  void freeLocalData( const Scalar* localValues ) const;

  //@}

private:

  // ///////////////////////////////////////
  // Private data members
  
  Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >   mpiSpace_;
  Teuchos::RefCountPtr<Scalar>                              localValues_;
  Index                                                     stride_;

}; // end class MPIVectorStd

// /////////////////////////////////////////////////////
// Inline members

template<class Scalar>
inline
Teuchos::RefCountPtr<Scalar> MPIVectorStd<Scalar>::getRCPtr()
{
  return localValues_;
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const Scalar> MPIVectorStd<Scalar>::getRCPtr() const
{
  return localValues_;
}

template<class Scalar>
inline
Scalar* MPIVectorStd<Scalar>::getPtr()
{
  return localValues_.get();
}

template<class Scalar>
inline
const Scalar* MPIVectorStd<Scalar>::getPtr() const
{
  return localValues_.get();
}

template<class Scalar>
inline
Index MPIVectorStd<Scalar>::getStride() const
{
  return stride_;
}	

} // end namespace Thyra

#endif // THYRA_MPI_VECTOR_STD_DECL_HPP
