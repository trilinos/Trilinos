// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
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
#include "Intrepid2_CubatureDirectTetDefault.hpp"
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
    template<typename ExecSpaceType,
             typename pointValueType = double,
             typename weightValueType = double>
    static Teuchos::RCP<Cubature<ExecSpaceType,pointValueType,weightValueType> > 
    create( unsigned                         topologyKey,
            const std::vector<ordinal_type> &degree,
            const EPolyType                  polytype = POLYTYPE_MAX );

    /** \brief Factory method.

        \param cell        [in]    - Cell topology.
        \param degree      [in]    - Array of polynomial degrees, one for each component cubature.

        \return
        - RCP to cubature with given specifications.
    */
    template<typename ExecSpaceType,
             typename pointValueType = double,
             typename weightValueType = double>
    static Teuchos::RCP<Cubature<ExecSpaceType,pointValueType,weightValueType> >
    create( const shards::CellTopology       cellTopology,
            const std::vector<ordinal_type> &degree,
            const EPolyType                  polytype = POLYTYPE_MAX );


    /** \brief Factory method.

        \param topologyKey [in]    - Key of the cell topology.
        \param degree      [in]    - A single polynomial degree, used for all component cubatures.

        \return
        - RCP to cubature with given specifications.
    */
    template<typename ExecSpaceType,
             typename pointValueType = double,
             typename weightValueType = double>
    static Teuchos::RCP<Cubature<ExecSpaceType,pointValueType,weightValueType> >
    create( unsigned                    topologyKey,
            const ordinal_type          degree,
            const EPolyType             polytype = POLYTYPE_MAX );
    
    /** \brief Factory method.

        \param cell        [in]    - Cell topology.
        \param degree      [in]    - A single polynomial degree, used for all component cubatures.

        \return
        - RCP to cubature with given specifications.
    */
    template<typename ExecSpaceType,
             typename pointValueType = double,
             typename weightValueType = double>
    static Teuchos::RCP<Cubature<ExecSpaceType,pointValueType,weightValueType> > 
    create( const shards::CellTopology  cellTopology,
            const ordinal_type          degree,
            const EPolyType             polytype = POLYTYPE_MAX );


    /** \brief Factory method for polygon cubature.
      
        \param cellTopology   [in]   - Cell topology
        \param cellVertices   [in]   - Vertices of physical cell
        \param degree         [in]   - A single polynomial degree used for all triangles in tessalation
      
        \return 
        - RCP to cubature with given specifications.
    */
    // template<typename ExecSpaceType>
    // template<typename cellVertexValueType, class ...cellVertexProperties>
    // static Teuchos::RCP<Cubature<ExecSpace> > 
    // create( const shards::CellTopology &cellTopology,
    //         const Kokkos::DynRankView<cellVertexValueType,cellVertexProperties> cellVertices,
    //         const ordinal_type degree );
  };
  
}// namespace Intrepid2

#include "Intrepid2_DefaultCubatureFactoryDef.hpp"

#endif
