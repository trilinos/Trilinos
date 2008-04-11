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

  TEST_FOR_EXCEPTION( (is_null(BMap_[basisCode])),
                      std::invalid_argument,
                      ">>> ERROR (DefaultBasisFactory): Invalid set of parameters prevents basis creation.");

  return BMap_[basisCode];

/*  // Create generic basis.
  Teuchos::RCP<Basis<Scalar> > pickBasis; */

/*  switch(field) {
  
    case FIELD_FORM_0: {

      switch(cell) {
     
        case CELL_TRI: {
        
          switch(recSpace) {

            case RECONSTRUCTION_SPACE_COMPLETE: {

              switch(degree) {

                case 1: {

                  switch(basisType) {

                    case BASIS_FEM_DEFAULT: {

                      switch(coordSys) {

                        case COORDINATES_CARTESIAN: {
                          pickBasis = Teuchos::rcp( new Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>() );
                        }
                        break;

                        default:
                          TEST_FOR_EXCEPTION( (coordSys != COORDINATES_CARTESIAN),
                                              std::invalid_argument,
                                              ">>> ERROR (DefaultBasisFactory): Invalid coordinate system.");
                      }
                    }
                    break; 

                    default:
                      TEST_FOR_EXCEPTION( (basisType != BASIS_FEM_DEFAULT),
                                          std::invalid_argument,
                                          ">>> ERROR (DefaultBasisFactory): Invalid basis type.");
                  }
                }
                break;

                default:
                  TEST_FOR_EXCEPTION( (degree != 1),
                                      std::invalid_argument,
                                      ">>> ERROR (DefaultBasisFactory): Invalid polynomial degree.");
              }
            }
            break;

            default:
              TEST_FOR_EXCEPTION( (recSpace != RECONSTRUCTION_SPACE_COMPLETE),
                                  std::invalid_argument,
                                  ">>> ERROR (DefaultBasisFactory): Invalid reconstruction space.");
          }
        }
        break;

        default:
          TEST_FOR_EXCEPTION( (cell != CELL_TRI),
                              std::invalid_argument,
                              ">>> ERROR (DefaultBasisFactory): Invalid cell type.");
      }
    }
    break;

    default:
      TEST_FOR_EXCEPTION( (field != FIELD_FORM_0),
                          std::invalid_argument,
                          ">>> ERROR (DefaultBasisFactory): Invalid field type.");
  }*/


/*  if (
       (field==FIELD_FORM_0) &&
       (cell==CELL_TRI) &&
       (recSpace==RECONSTRUCTION_SPACE_COMPLETE) &&
       (degree==1) &&
       (basisType==BASIS_FEM_DEFAULT) &&
       (coordSys==COORDINATES_CARTESIAN)
     )
     {
    pickBasis = Teuchos::rcp( new Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>() );
  }
  else {
    TEST_FOR_EXCEPTION( (true),
                        std::invalid_argument,
                        ">>> ERROR (DefaultBasisFactory): Invalid set of parameters prevents basis creation.");
  }*/

  //return pickBasis;
  
}

} // namespace Intrepid
