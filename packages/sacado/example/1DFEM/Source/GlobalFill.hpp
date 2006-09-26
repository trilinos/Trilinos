// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef GLOBALFILL_HPP
#define GLOBALFILL_HPP

#include <vector>

#include "Teuchos_RefCountPtr.hpp"

#include "AbstractPDE.hpp"
#include "AbstractQuadrature.hpp"
#include "Mesh.hpp"
#include "AbstractInitPostOp.hpp"

template <typename ScalarT>
class GlobalFill {
public:

  //! Constructor
  GlobalFill(
	  const Teuchos::RefCountPtr<const Mesh>& elementMesh,
	  const Teuchos::RefCountPtr<const AbstractQuadrature>& quadRule,
	  const Teuchos::RefCountPtr< AbstractPDE<ScalarT> >& pdeEquations);
  
  //! Destructor
  ~GlobalFill();

  //! Compute global fill
  void computeGlobalFill(AbstractInitPostOp<ScalarT>& initPostOp);

private:

  //! Private to prohibit copying
  GlobalFill(const GlobalFill&);

  //! Private to prohibit copying
  GlobalFill& operator=(const GlobalFill&);

protected:

  //! Element mesh
  Teuchos::RefCountPtr<const Mesh> mesh;

  //! Quadrature rule
  Teuchos::RefCountPtr<const AbstractQuadrature> quad;

  //! PDE Equations
  Teuchos::RefCountPtr< AbstractPDE<ScalarT> > pde;

  //! Number of nodes per element
  unsigned int nnode;

  //! Number of PDE equations
  unsigned int neqn;

  //! Number of element-level DOF
  unsigned int ndof;

  //! Element solution variables
  std::vector<ScalarT> elem_x;

  //! Element residual variables
  std::vector<ScalarT> elem_f;

};

#include "GlobalFillImpl.hpp"

#endif // GLOBALFILL_HPP
