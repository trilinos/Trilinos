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

/** \file   Intrepid_DefaultBasisFactory.hpp
\brief  Header file for the abstract base class Intrepid::DefaultBasisFactory.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_DEFAULT_BASIS_FACTORY_HPP
#define INTREPID_DEFAULT_BASIS_FACTORY_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Basis.hpp"
#include "Teuchos_RCP.hpp"

/////   list of default basis includes   /////

#include "Intrepid_F0_QUAD_I1_FEM_DEFAULT.hpp"
#include "Intrepid_F0_TRI_C1_FEM_DEFAULT.hpp"
#include "Intrepid_F0_TRI_C2_FEM_DEFAULT.hpp"

///// end of list of default basis includes /////

#ifndef DOXYGEN_SHOULD_SKIP_THIS
///// FIAT-generated element includes here

#include "Intrepid_F0_TRI_C1_FEM_FIAT.hpp"
#include "Intrepid_F0_TRI_C2_FEM_FIAT.hpp"
#include "Intrepid_F0_TRI_C3_FEM_FIAT.hpp"

#include "Intrepid_F0_TET_C1_FEM_FIAT.hpp"
#include "Intrepid_F0_TET_C2_FEM_FIAT.hpp"
#include "Intrepid_F0_TET_C3_FEM_FIAT.hpp"

#include "Intrepid_F2_TRI_I1_FEM_FIAT.hpp"
#include "Intrepid_F2_TRI_I2_FEM_FIAT.hpp"
#include "Intrepid_F2_TRI_I3_FEM_FIAT.hpp"

#include "Intrepid_F2_TET_I1_FEM_FIAT.hpp"
#include "Intrepid_F2_TET_I2_FEM_FIAT.hpp"
#include "Intrepid_F2_TET_I3_FEM_FIAT.hpp"

#include "Intrepid_F1_TRI_I1_FEM_FIAT.hpp"
#include "Intrepid_F1_TRI_I2_FEM_FIAT.hpp"
#include "Intrepid_F1_TRI_I3_FEM_FIAT.hpp"

#include "Intrepid_F1_TET_I1_FEM_FIAT.hpp"
#include "Intrepid_F1_TET_I2_FEM_FIAT.hpp"
#include "Intrepid_F1_TET_I3_FEM_FIAT.hpp"

///// end FIAT-generated element includes
#endif

namespace Intrepid {
  
/** \class Intrepid::DefaultBasisFactory
    \brief A factory class that generates specific instances of bases.
*/
template<class Scalar>
class DefaultBasisFactory {
  private:
    std::map<unsigned long, Teuchos::RCP<Basis<Scalar> > > basisMap_;

  public:
    
  /** \brief Default constructor.
  */
  DefaultBasisFactory() {
    /* List all basis keys with corresponding basis classes.
       Legend:
       F_ - two-digit field code (see EField) - UPPER LIMIT is 41 due to unsigned long integer standard!
       C_ - two-digit cell code (see ECell)
       R  - one-digit reconstruction space code (see EReconstructionSpace)
       D_ - two-digit degree
       B  - two-digit basis code (see EBasis)
       S  - one-digit coordinate system code (see ECoordinates)
       If any of the leading codes are zero (00 or 0), must leave blank!
    */

    /******** F_C_RD_B_S ************************************************************/
    basisMap_[   2001000] = Teuchos::rcp( new Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>() );
    basisMap_[   2002000] = Teuchos::rcp( new Basis_F0_TRI_C2_FEM_DEFAULT<Scalar>() );
    basisMap_[   3101000] = Teuchos::rcp( new Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>() );

    // begin FIAT-generated bases
    basisMap_[   2001020] = Teuchos::rcp( new Basis_F0_TRI_C1_FEM_FIAT<Scalar>() );
    basisMap_[   2002020] = Teuchos::rcp( new Basis_F0_TRI_C2_FEM_FIAT<Scalar>() );
    basisMap_[   2003020] = Teuchos::rcp( new Basis_F0_TRI_C3_FEM_FIAT<Scalar>() );

    basisMap_[   4001020] = Teuchos::rcp( new Basis_F0_TET_C1_FEM_FIAT<Scalar>() );
    basisMap_[   4002020] = Teuchos::rcp( new Basis_F0_TET_C2_FEM_FIAT<Scalar>() );
    basisMap_[   4003020] = Teuchos::rcp( new Basis_F0_TET_C3_FEM_FIAT<Scalar>() );

    basisMap_[ 202101020] = Teuchos::rcp( new Basis_F2_TRI_I1_FEM_FIAT<Scalar>() );
    basisMap_[ 202102020] = Teuchos::rcp( new Basis_F2_TRI_I2_FEM_FIAT<Scalar>() );
    basisMap_[ 202103020] = Teuchos::rcp( new Basis_F2_TRI_I3_FEM_FIAT<Scalar>() );

    basisMap_[ 204101020] = Teuchos::rcp( new Basis_F2_TET_I1_FEM_FIAT<Scalar>() );
    basisMap_[ 204102020] = Teuchos::rcp( new Basis_F2_TET_I2_FEM_FIAT<Scalar>() );
    basisMap_[ 204103020] = Teuchos::rcp( new Basis_F2_TET_I3_FEM_FIAT<Scalar>() );

    basisMap_[ 102101020] = Teuchos::rcp( new Basis_F1_TRI_I1_FEM_FIAT<Scalar>() );
    basisMap_[ 102102020] = Teuchos::rcp( new Basis_F1_TRI_I2_FEM_FIAT<Scalar>() );
    basisMap_[ 102103020] = Teuchos::rcp( new Basis_F1_TRI_I3_FEM_FIAT<Scalar>() );

    basisMap_[ 104101020] = Teuchos::rcp( new Basis_F1_TET_I1_FEM_FIAT<Scalar>() );
    basisMap_[ 104102020] = Teuchos::rcp( new Basis_F1_TET_I2_FEM_FIAT<Scalar>() );
    basisMap_[ 104103020] = Teuchos::rcp( new Basis_F1_TET_I3_FEM_FIAT<Scalar>() );
    // end FIAT-generated bases

  };

  /** \brief Destructor.
  */
  virtual ~DefaultBasisFactory() {};

  /** \brief Factory method.

      \param field       [in]    - Field type (FIELD_FORM_0, etc.).
      \param cell        [in]    - Cell type (CELL_TRI, CELL_QUAD, etc.).
      \param recSpace    [in]    - Reconstruction space type (RECONSTRUCTION_SPACE_COMPLETE, etc.).
      \param degree      [in]    - Polynomial degree.
      \param basisType   [in]    - Basis type (BASIS_FEM_DEFAULT, etc.).
      \param coordSys    [in]    - Coordinate system (COORDINATES_CARTESIAN, etc.).

      \return
              - RCP to basis with given specifications.
  */
  Teuchos::RCP<Basis<Scalar> > create(EField field,
                                      ECell cell,
                                      EReconstructionSpace recSpace,
                                      int degree,
                                      EBasis basisType, 
                                      ECoordinates coordSys);
    
};
  
}// namespace Intrepid

#include "Intrepid_DefaultBasisFactoryDef.hpp"

#endif






















































































































