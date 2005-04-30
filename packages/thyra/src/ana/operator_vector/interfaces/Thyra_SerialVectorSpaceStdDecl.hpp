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

#ifndef THYRA_SERIAL_VECTOR_SPACE_STD_DECL_HPP
#define THYRA_SERIAL_VECTOR_SPACE_STD_DECL_HPP

#include "Thyra_SerialVectorSpaceBaseDecl.hpp"

namespace Thyra {

/** \brief General concrete <tt>%VectorSpaceBase</tt> subclass for serial vectors and multi-vectors.
 *
 * The default copy constructor and assignment operators are allowed
 * since they have the correct semantics.
 *
 * \ingroup Thyra_Op_Vec_adapters_serial_concrete_std_grp
 */
template<class Scalar>
class SerialVectorSpaceStd : public SerialVectorSpaceBase<Scalar> {
public:

  /** @name Constructors / initializers */
  //@{

  /** \brief Calls <tt>this->initialize()</tt>.
   */
  SerialVectorSpaceStd( int dim = 0 );

  /** \brief Initialize given the dimension of the vector space.
   *
   * @param  dim   [in] The dimension of the vector space.
   */
  void initialize( int dim );

  //@}

  /** @name Overriddend form Teuchos::Describable */
  //@{
  /** \brief . */
  std::string describe() const;
  //@}

  /** @name Overridden from VectorSpece */
  //@{

  /// Returns 0 if uninitialized
  Index dim() const;
  /// Returns a <tt>SerialVector</tt> object.
  Teuchos::RefCountPtr<VectorBase<Scalar> > createMember() const;
  /// Returns a <tt>SerialMultiVector</tt> object.
  Teuchos::RefCountPtr< MultiVectorBase<Scalar> > createMembers(int numMembers) const;
  /// Clones the object as promised
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > clone() const;

  //@}

private:

  int   dim_;

}; // end class SerialVectorSpaceStd

} // end namespace Thyra

#endif // THYRA_SERIAL_VECTOR_SPACE_STD_DECL_HPP
