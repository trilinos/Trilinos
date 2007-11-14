// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_PCE_WORKSPACE_HPP
#define SACADO_PCE_WORKSPACE_HPP

#include <vector>
#include "Sacado_PCE_TripleProduct.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Sacado {

  namespace PCE {

    //! Workspace class to store linear system for nonlinear PCE operations
    template <typename BasisT> 
    class Workspace {
    public:

      //! Typename of values
      typedef typename BasisT::value_type value_type;

      //! Ordinal type
      typedef int ordinal_type;

      //! Typename of matrix
      typedef Teuchos::SerialDenseMatrix<ordinal_type,value_type> matrix_type;

      //! Typename of TripleProduct tensor
      typedef TripleProduct<BasisT> tp_type;

      //! Constructor
      Workspace(unsigned int sz);

      //! Destructor
      ~Workspace() {}

      //! Resize workspace
      void resize(unsigned int sz);

      //! Get workspace size
      unsigned int size() const { return sz; }

      //! Get matrix
      matrix_type& getMatrix() { return A; }

      //! Get RHS
      matrix_type& getRHS() { return b; }

      //! Get TripleProduct tensor
      const tp_type& getTripleProduct() const { return Cijk; }

      //! Solve linear system
      ordinal_type solve(ordinal_type s, ordinal_type nrhs);

    protected:

      //! Workspace size
      unsigned int sz;

      //! Matrix
      matrix_type A;

      //! RHS
      matrix_type b;

      //! Pivot array
      std::vector<ordinal_type> piv;

      //! Triple-product tensor
      tp_type Cijk;

      //! LAPACK wrappers
      Teuchos::LAPACK<ordinal_type,value_type> lapack;

    }; // class Workspace

  } // namesspace PCE

} // namespace Sacado

#include "Sacado_PCE_WorkspaceImp.hpp"
      
#endif // SACADO_PCE_WORKSPACE_HPP
