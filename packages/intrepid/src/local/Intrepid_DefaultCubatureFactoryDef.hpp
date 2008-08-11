// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_DefaultCubatureFactoryDef.hpp
\brief  Definition file for the class Intrepid::DefaultCubatureFactory.
\author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

// create method
template<class Scalar>
Teuchos::RCP<Cubature<Scalar> > DefaultCubatureFactory<Scalar>::create(ECell cell,
                                                                       int degree) {
  // Create generic cubature.
  Teuchos::RCP<Cubature<Scalar> > pickCubature;

  switch (cell) {

    case CELL_EDGE:
    case CELL_TRI:
    case CELL_TET:
    case CELL_PYRAMID:
      pickCubature = Teuchos::rcp( new CubatureDirect<Scalar>(cell, degree) );
    break;

    case CELL_QUAD:
    case CELL_HEX:
    case CELL_TRIPRISM:
      pickCubature = Teuchos::rcp( new CubatureTensor<Scalar>(cell, degree) );
    break;

    default:
      TEST_FOR_EXCEPTION( ( (cell != CELL_EDGE)      && (cell != CELL_TRI)      &&
                            (cell != CELL_TET)       && (cell != CELL_PYRAMID)  &&
                            (cell != CELL_QUAD)      && (cell != CELL_HEX)      &&
                            (cell != CELL_TRIPRISM) ),
                          std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Invalid cell type prevents cubature creation.");
  }

  return pickCubature;
  
}

} // namespace Intrepid
