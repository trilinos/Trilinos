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

/** \file   Intrepid_DefaultBasisFactoryDef.hpp
\brief  Definition file for the class Intrepid::DefaultBasisFactory.
\author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

// create method
template<class Scalar>
Teuchos::RCP<Basis<Scalar> > DefaultBasisFactory<Scalar>::create(EField field,
                                                                 ECell cell,
                                                                 EReconstructionSpace recSpace,
                                                                 int degree,
                                                                 EBasis basisType,
                                                                 ECoordinates coordSys) {
  unsigned long basisCode = field *     100000000 +
                            cell *      1000000   +
                            recSpace *  100000    +
                            degree *    1000      +
                            basisType * 10        +
                            coordSys;

  TEST_FOR_EXCEPTION( (is_null(basisMap_[basisCode])), 
                      std::invalid_argument,
                      ">>> ERROR (DefaultBasisFactory): Invalid set of parameters prevents basis creation.");

  return basisMap_[basisCode];
}

} // namespace Intrepid
