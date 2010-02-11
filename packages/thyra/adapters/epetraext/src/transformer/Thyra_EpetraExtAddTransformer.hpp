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

#ifndef THYRA_EPETRAEXT_ADD_TRANSFORMER_HPP
#define THYRA_EPETRAEXT_ADD_TRANSFORMER_HPP

#include "Thyra_LinearOpTransformerBase.hpp"


namespace Thyra {


/** \brief Transformer subclass for adding Epetra/Thyra operators using
 * EpetraExt::MatrixMatrix.
 *
 * \ingroup EpetraExt_Thyra_Op_Vec_adapters_grp
 */
class EpetraExtAddTransformer : public LinearOpTransformerBase<double>
{
public:

  /** \name Overridden from LinearOpTransformerBase. */
  //@{

  /** \brief . */
  virtual bool isCompatible(const LinearOpBase<double> &op_in) const;

  /** \brief . */
  virtual RCP<LinearOpBase<double> > createOutputOp() const;

  /** \brief . */
  virtual void transform(
    const LinearOpBase<double> &op_in,
    const Ptr<LinearOpBase<double> > &op_inout
    ) const;

  //@}

private:
  
};


/** \brief Nonmember constructor.
 *
 * \relates EpetraExtAddTransformer
 */
inline
RCP<EpetraExtAddTransformer>
epetraExtAddTransformer()
{
  return Teuchos::rcp(new EpetraExtAddTransformer());
}


} // namespace Thyra


#endif	// THYRA_EPETRAEXT_ADD_TRANSFORMER_HPP
