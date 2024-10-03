// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
