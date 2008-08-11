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

/** \file   Intrepid_DefaultFieldFactory.hpp
\brief  Header file for the abstract base class Intrepid::DefaultFieldFactory.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_DEFAULT_FIELD_FACTORY_HPP
#define INTREPID_DEFAULT_FIELD_FACTORY_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_LocalField.hpp"
#include "Intrepid_DefaultBasisFactory.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Teuchos_RCP.hpp"

/////   list of default basis includes   /////

#include "Intrepid_LocalForm0.hpp"

///// end of list of default basis includes /////


namespace Intrepid {
  
/** \class Intrepid::DefaultFieldFactory
    \brief A factory class that generates specific instances of forms.
*/
template<class Scalar>
class DefaultFieldFactory {
  private:

  public:
    
  /** \brief Default constructor.
  */
  DefaultFieldFactory() {};

  /** \brief Destructor.
  */
  virtual ~DefaultFieldFactory() {};

  /** \brief Factory method.

      \param field       [in]    - Field type (FIELD_FORM_0, etc.).
      \param cell        [in]    - Cell type (CELL_TRI, CELL_QUAD, etc.).
      \param recSpace    [in]    - Reconstruction space type (RECONSTRUCTION_SPACE_COMPLETE, etc.).
      \param degree      [in]    - Polynomial degree.
      \param basisType   [in]    - Basis type (BASIS_FEM_DEFAULT, etc.).
      \param coordSys    [in]    - Coordinate system (COORDINATES_CARTESIAN, etc.).
      \param cubDegree   [in]    - Cubature accuracy; setting <var>cubDegree</var> to a nonnegative
                                   value overrides Intrepid's default selection of cubature accuracy.
                                   <b>Use with caution!</b>

      \return
              - RCP to field with given specifications.
  */
  Teuchos::RCP<LocalField<Scalar> > create(EField               field,
                                           ECell                cell,
                                           EReconstructionSpace recSpace,
                                           int                  degree,
                                           EBasis               basisType, 
                                           ECoordinates         coordSys,
                                           int                  cubDegree = -1);
    
};
  
}// namespace Intrepid

#include "Intrepid_DefaultFieldFactoryDef.hpp"

#endif
