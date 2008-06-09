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

/** \file   Intrepid_DefaultFieldFactoryDef.hpp
\brief  Definition file for the class Intrepid::DefaultFieldFactory.
\author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

// create method
template<class Scalar>
Teuchos::RCP<LocalField<Scalar> > DefaultFieldFactory<Scalar>::create(EField                field,
                                                                      ECell                 cell,
                                                                      EReconstructionSpace  recSpace,
                                                                      int                   polyDegree,
                                                                      EBasis                basisType,
                                                                      ECoordinates          coordSys,
                                                                      int                   cubDegree) {
  // Create basis.
  DefaultBasisFactory<Scalar> BFactory;
  Teuchos::RCP<Basis<Scalar> > basis =  BFactory.create(
    field, cell, recSpace, polyDegree, basisType, coordSys);


  // Create cell and subcell cubatures.

  // cubDegree should be sufficient for all differential operators,
  // including value-value bilinear forms (mass matrices)
  if (cubDegree < 0) {
    cubDegree = polyDegree*2;
  }
  Teuchos::Array<Teuchos::Array< Teuchos::RCP<Cubature<Scalar> > > > allCubs;
  DefaultCubatureFactory<Scalar> CFactory;

  int cDim = MultiCell<Scalar>::getCellDim(cell);
  switch (cDim) {
    case 3: {
      int numFaces = MultiCell<Scalar>::getCellNumSubcells(cell, 2);
      int numEdges = MultiCell<Scalar>::getCellNumSubcells(cell, 1);
      allCubs.resize(3);
      // set cell cubature
      Teuchos::RCP<Cubature<Scalar> > cellCub =  CFactory.create(cell, cubDegree);
      allCubs[0].resize(1);
      allCubs[0][0] = cellCub;
      // set face cubatures
      allCubs[1].resize(numFaces);
      for (int face=0; face<numFaces; face++) {
        Teuchos::RCP<Cubature<Scalar> > faceCub =  CFactory.create(MultiCell<Scalar>::getSubcellType(cell,2,face), cubDegree);
        allCubs[1][face] = faceCub;
      }
      // set edge cubatures
      Teuchos::RCP<Cubature<Scalar> > edgeCub =  CFactory.create(CELL_EDGE, cubDegree);
      allCubs[2].resize(numEdges);
      for (int edge=0; edge<numEdges; edge++) {
        allCubs[2][edge] = edgeCub;
      }
    }
    break;

    case 2: {
      int numEdges = MultiCell<Scalar>::getCellNumSubcells(cell,1);
      allCubs.resize(2);
      // set cell cubature
      Teuchos::RCP<Cubature<Scalar> > cellCub =  CFactory.create(cell, cubDegree);
      allCubs[0].resize(1);
      allCubs[0][0] = cellCub;
      // set edge cubatures
      Teuchos::RCP<Cubature<Scalar> > edgeCub =  CFactory.create(CELL_EDGE, cubDegree);
      allCubs[1].resize(numEdges);
      for (int edge=0; edge<numEdges; edge++) {
        allCubs[1][edge] = edgeCub;
      }
    }
    break;

    case 1: {
      allCubs.resize(1);
      // set cell cubature
      Teuchos::RCP<Cubature<Scalar> > cellCub =  CFactory.create(cell, cubDegree);
      allCubs[0].resize(1);
      allCubs[0][0] = cellCub;
    }
    break;

    default:
      TEST_FOR_EXCEPTION(( (cDim != 3) && (cDim != 2) && (cDim != 1) ),
                         std::invalid_argument,
                         ">>> ERROR (DefaultFieldFactory): Cell dimension not supported!");
  }


  // Create generic field.
  Teuchos::RCP<LocalField<Scalar> > pickField;
  // fill field
  switch (field) {

    case FIELD_FORM_0: {
      pickField = Teuchos::rcp( new LocalForm0<Scalar>(basis, allCubs) );
    }
    break;

    default:
    TEST_FOR_EXCEPTION( (true),
                        std::invalid_argument,
                        ">>> ERROR (DefaultFieldFactory): Invalid field parameter prevents field creation!");
  }

  return pickField;
  
}

} // namespace Intrepid
