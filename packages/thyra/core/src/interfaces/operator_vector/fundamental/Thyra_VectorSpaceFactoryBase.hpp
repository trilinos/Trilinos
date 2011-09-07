// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_VECTOR_SPACE_FACTORY_DECL_HPP
#define THYRA_VECTOR_SPACE_FACTORY_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_Describable.hpp"

namespace Thyra {

/** \brief Abstract interface for objects that can create vector spaces of a
 * specified dimension.
 *
 * The primary role that a <tt>%VectorSpaceFactoryBase</tt> object takes is
 * defined in the documentation for the class <tt>VectorSpaceBase</tt> and is
 * related to the domain space of <tt>MultiVectorBase</tt> objects.  However,
 * this is a general factory interface class that can be used to create almost
 * any <tt>VectorSpaceBase</tt> object just given a dimension.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class Scalar>
class VectorSpaceFactoryBase : virtual public Teuchos::Describable {
public:

  /** \brief . */
  virtual ~VectorSpaceFactoryBase() {}

  /** @name Public pure virtual functions that must be overridden */
  //@{

  /** \brief Create a vector space of the given dimension.
   *
   * @param  dim  [in] The dimension of the vector space to create.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>dim > 0</tt> (throw <tt>std::invalid_argument</tt>).
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>return->dim() == dim</tt>
   * </ul>
   *
   * @return Returns a smart reference-counted pointer to a dynamically
   * allocated vector space object that can be used to create vectors and
   * multi-vectors.
   */
  virtual RCP< const VectorSpaceBase<Scalar> > createVecSpc(int dim) const = 0;

  //@}

private:
  
  // Not defined and not to be called
  VectorSpaceFactoryBase<Scalar>&
  operator=(const VectorSpaceFactoryBase<Scalar>&);

};


} // end namespace Thyra


#endif  // THYRA_VECTOR_SPACE_FACTORY_DECL_HPP
