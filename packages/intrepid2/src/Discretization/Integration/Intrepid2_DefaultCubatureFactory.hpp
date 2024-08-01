// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_DefaultCubatureFactory.hpp
    \brief  Header file for the abstract base class Intrepid2::DefaultCubatureFactory.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_DEFAULT_CUBATURE_FACTORY_HPP__
#define __INTREPID2_DEFAULT_CUBATURE_FACTORY_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Utils.hpp"

#include "Shards_CellTopology.hpp"
#include "Teuchos_RCP.hpp"

#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_CubatureDirectLineGauss.hpp"
#include "Intrepid2_CubatureDirectLineGaussJacobi20.hpp"
#include "Intrepid2_CubatureDirectTriDefault.hpp"
#include "Intrepid2_CubatureDirectTriSymmetric.hpp"
#include "Intrepid2_CubatureDirectTetDefault.hpp"
#include "Intrepid2_CubatureDirectTetSymmetric.hpp"
//#include "Intrepid2_CubatureCompositeTet.hpp"
#include "Intrepid2_CubatureTensor.hpp"
#include "Intrepid2_CubatureTensorPyr.hpp"
//#include "Intrepid2_CubaturePolygon.hpp"

#include "Intrepid2_CubaturePolylib.hpp"

namespace Intrepid2 {
  
  /** \class Intrepid2::DefaultCubatureFactory
      \brief A factory class that generates specific instances of cubatures.
  */

  class DefaultCubatureFactory {
  public:

    /** \brief Factory method.

        \param topologyKey [in]    - Key of the cell topology.
        \param degree      [in]    - Array of polynomial degrees, one for each component cubature.

        \return
        - RCP to cubature with given specifications.
    */
    template<typename DeviceType,
             typename pointValueType = double,
             typename weightValueType = double>
    static Teuchos::RCP<Cubature<DeviceType,pointValueType,weightValueType> > 
    create( unsigned                         topologyKey,
            const std::vector<ordinal_type> &degree,
            const EPolyType                  polytype = POLYTYPE_MAX,
            const bool                       symmetric = false );

    /** \brief Factory method.

        \param cell        [in]    - Cell topology.
        \param degree      [in]    - Array of polynomial degrees, one for each component cubature.

        \return
        - RCP to cubature with given specifications.
    */
    template<typename DeviceType,
             typename pointValueType = double,
             typename weightValueType = double>
    static Teuchos::RCP<Cubature<DeviceType,pointValueType,weightValueType> >
    create( const shards::CellTopology       cellTopology,
            const std::vector<ordinal_type> &degree,
            const EPolyType                  polytype = POLYTYPE_MAX,
            const bool                       symmetric = false  );


    /** \brief Factory method.

        \param topologyKey [in]    - Key of the cell topology.
        \param degree      [in]    - A single polynomial degree, used for all component cubatures.

        \return
        - RCP to cubature with given specifications.
    */
    template<typename DeviceType,
             typename pointValueType = double,
             typename weightValueType = double>
    static Teuchos::RCP<Cubature<DeviceType,pointValueType,weightValueType> >
    create( unsigned                    topologyKey,
            const ordinal_type          degree,
            const EPolyType             polytype = POLYTYPE_MAX, 
            const bool                  symmetric = false  );
    
    /** \brief Factory method.

        \param cell        [in]    - Cell topology.
        \param degree      [in]    - A single polynomial degree, used for all component cubatures.

        \return
        - RCP to cubature with given specifications.
    */
    template<typename DeviceType,
             typename pointValueType = double,
             typename weightValueType = double>
    static Teuchos::RCP<Cubature<DeviceType,pointValueType,weightValueType> > 
    create( const shards::CellTopology  cellTopology,
            const ordinal_type          degree,
            const EPolyType             polytype = POLYTYPE_MAX, 
            const bool                  symmetric = false   );


    /** \brief Factory method for polygon cubature.
      
        \param cellTopology   [in]   - Cell topology
        \param cellVertices   [in]   - Vertices of physical cell
        \param degree         [in]   - A single polynomial degree used for all triangles in tessalation
      
        \return 
        - RCP to cubature with given specifications.
    */
    // template<typename DeviceType>
    // template<typename cellVertexValueType, class ...cellVertexProperties>
    // static Teuchos::RCP<Cubature<ExecSpace> > 
    // create( const shards::CellTopology &cellTopology,
    //         const Kokkos::DynRankView<cellVertexValueType,cellVertexProperties> cellVertices,
    //         const ordinal_type degree );
  };
  
}// namespace Intrepid2

#include "Intrepid2_DefaultCubatureFactoryDef.hpp"

#endif
