//
// @HEADER
// ***********************************************************************
// 
//                           Rythmos Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef Rythmos_INOUT_ARGS_H
#define Rythmos_INOUT_ARGS_H

#include "Thyra_VectorBase.hpp"
//#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"

namespace Rythmos {

//-----------------------------------------------------------------------------
// Class         : InArgs
// Purpose       : Input arguments class for NonlinearModel evaluateModel function
// Special Notes :
// Creator       : Todd Coffey, SNL
// Creation Date : 05/20/05
//-----------------------------------------------------------------------------
template<class Scalar> 
class InArgs
{
  public:
    
    // Destructor
    ~InArgs();

    // Constructor
    InArgs(const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &x, Scalar t)
      { x_ = x; t_ = t; };
    InArgs();

    void set_x(const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &x)
      { x_ = x; };
    const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_x()
      { return(x_); };

    void set_t(Scalar t)
      { t_ = t; };
    Scalar get_t()
      { return(t_); };

  protected:
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_;
    Scalar t_;

};
//
//-----------------------------------------------------------------------------
// Class         : OutArgs
// Purpose       : Output arguments class for NonlinearModel evaluateModel function
// Special Notes :
// Creator       : Todd Coffey, SNL
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
template<class Scalar> 
class OutArgs
{
  public:
    
    // Destructor
    ~OutArgs();

    // Constructor
    OutArgs();

    // Request residual:
    void request_F(const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &F) { F_ = F; };
    const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > get_F()
      { return(F_); };

  protected:
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > F_;
    
};

}// namespace Rythmos 

#endif // Rythmos_INOUT_ARGS_H
