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

/** \file   Intrepid_DefaultCubatureFactory.hpp
\brief  Header file for the abstract base class Intrepid::DefaultCubatureFactory.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_DEFAULT_CUBATURE_FACTORY_HPP
#define INTREPID_DEFAULT_CUBATURE_FACTORY_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Shards_CellTopology.hpp"
#include "Teuchos_RCP.hpp"

/////   list of default cubature includes   /////

#include "Intrepid_CubatureDirectLineGauss.hpp"
#include "Intrepid_CubatureDirectTriDefault.hpp"
#include "Intrepid_CubatureDirectTetDefault.hpp"
#include "Intrepid_CubatureCompositeTet.hpp"
#include "Intrepid_CubatureTensor.hpp"

///// end of list of default cubature includes /////


namespace Intrepid {
  
/** \class Intrepid::DefaultCubatureFactory
    \brief A factory class that generates specific instances of cubatures.
*/
template<class Scalar, class ArrayPoint=FieldContainer<Scalar>, class ArrayWeight=ArrayPoint >
class DefaultCubatureFactory {
  private:

  public:
    
  /** \brief Default constructor.
  */
  DefaultCubatureFactory() {};

  /** \brief Destructor.
  */
  virtual ~DefaultCubatureFactory() {};

  /** \brief Factory method.

      \param cell        [in]    - Cell topology.
      \param degree      [in]    - Array of polynomial degrees, one for each component cubature.

      \return
              - RCP to cubature with given specifications.
  */
  Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > create(const shards::CellTopology & cellTopology,
                                                                const std::vector<int> & degree);

  /** \brief Factory method.

      \param cell        [in]    - Cell topology.
      \param degree      [in]    - A single polynomial degree, used for all component cubatures.

      \return
              - RCP to cubature with given specifications.
  */
  Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > create(const shards::CellTopology & cellTopology,
                                                                int   degree);
    
};
  
}// namespace Intrepid

#include "Intrepid_DefaultCubatureFactoryDef.hpp"

#endif
