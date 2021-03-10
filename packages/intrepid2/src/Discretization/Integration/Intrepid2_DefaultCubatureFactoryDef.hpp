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

/** \file   Intrepid2_DefaultCubatureFactoryDef.hpp
    \brief  Definition file for the class Intrepid2::DefaultCubatureFactory.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_DEFAULT_CUBATURE_FACTORY_DEF_HPP__
#define __INTREPID2_DEFAULT_CUBATURE_FACTORY_DEF_HPP__

namespace Intrepid2 {

  // first create method
  template<typename DT, typename PT, typename WT>
  Teuchos::RCP<Cubature<DT,PT,WT> > 
  DefaultCubatureFactory::
  create( unsigned topologyKey,
          const std::vector<ordinal_type> &degree,
          const EPolyType                  polytype ) {

    // Create generic cubature.
    Teuchos::RCP<Cubature<DT,PT,WT> > r_val;

    switch (topologyKey) {
    case shards::Line<>::key: {
      INTREPID2_TEST_FOR_EXCEPTION( degree.size() < 1, std::invalid_argument,
                                    ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      if (isValidPolyType(polytype))
        r_val = Teuchos::rcp(new CubaturePolylib<DT,PT,WT>(degree[0], polytype));
      else
        r_val = Teuchos::rcp(new CubatureDirectLineGauss<DT,PT,WT>(degree[0]));
      break;
    }
    case shards::Triangle<>::key: {
      INTREPID2_TEST_FOR_EXCEPTION( degree.size() < 1, std::invalid_argument,
                                    ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      r_val = Teuchos::rcp(new CubatureDirectTriDefault<DT,PT,WT>(degree[0]));
      break;
    }
    case shards::Quadrilateral<>::key: 
    case shards::ShellQuadrilateral<>::key: {
      INTREPID2_TEST_FOR_EXCEPTION( degree.size() < 2, std::invalid_argument,
                                    ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      if (isValidPolyType(polytype)) {
        const auto x_line = CubaturePolylib<DT,PT,WT>(degree[0], polytype);
        const auto y_line = ( degree[1] == degree[0] ? x_line : CubaturePolylib<DT,PT,WT>(degree[1], polytype) );
        r_val = Teuchos::rcp(new CubatureTensor<DT,PT,WT>( x_line, y_line ));
      } else {
        const auto x_line = CubatureDirectLineGauss<DT,PT,WT>(degree[0]);
        const auto y_line = ( degree[1] == degree[0] ? x_line : CubatureDirectLineGauss<DT,PT,WT>(degree[1]) );
        r_val = Teuchos::rcp(new CubatureTensor<DT,PT,WT>( x_line, y_line ));
      }
      break;
    }
    case shards::Tetrahedron<>::key: {
      INTREPID2_TEST_FOR_EXCEPTION( degree.size() < 1, std::invalid_argument,
                                    ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      r_val = Teuchos::rcp(new CubatureDirectTetDefault<DT,PT,WT>(degree[0]));
      break;
    }
    case shards::Hexahedron<>::key: {
      INTREPID2_TEST_FOR_EXCEPTION( degree.size() < 3, std::invalid_argument,
                                    ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      if (isValidPolyType(polytype)) {
        const auto x_line = CubaturePolylib<DT,PT,WT>(degree[0], polytype);
        const auto y_line = ( degree[1] == degree[0] ? x_line : CubaturePolylib<DT,PT,WT>(degree[1], polytype) );
        const auto z_line = ( degree[2] == degree[0] ? x_line : 
                              degree[2] == degree[1] ? y_line : CubaturePolylib<DT,PT,WT>(degree[2], polytype) );

        r_val = Teuchos::rcp(new CubatureTensor<DT,PT,WT>( x_line, y_line, z_line ));
      } else {
        const auto x_line = CubatureDirectLineGauss<DT,PT,WT>(degree[0]);
        const auto y_line = ( degree[1] == degree[0] ? x_line : CubatureDirectLineGauss<DT,PT,WT>(degree[1]) );
        const auto z_line = ( degree[2] == degree[0] ? x_line : 
                              degree[2] == degree[1] ? y_line : CubatureDirectLineGauss<DT,PT,WT>(degree[2]) );

        r_val = Teuchos::rcp(new CubatureTensor<DT,PT,WT>( x_line, y_line, z_line ));
      }
      break;
    }
    case shards::Wedge<>::key: {
      INTREPID2_TEST_FOR_EXCEPTION( degree.size() < 2, std::invalid_argument,
                                    ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.")
        {
          const auto xy_tri = CubatureDirectTriDefault<DT,PT,WT>(degree[0]);
          const auto z_line = CubatureDirectLineGauss<DT,PT,WT>(degree[1]);
          r_val = Teuchos::rcp(new CubatureTensor<DT,PT,WT>( xy_tri, z_line ));
        }
      break;
    }
    case shards::Pyramid<>::key: {
      INTREPID2_TEST_FOR_EXCEPTION( degree.size() < 1, std::invalid_argument,
                                    ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      {
        // if direct gauss is used for pyramid 
        // need over-integration to account for the additional transformation
        const auto xy_line = CubatureDirectLineGauss<DT,PT,WT>(degree[0]);
        const auto z_line  = CubatureDirectLineGaussJacobi20<DT,PT,WT>(degree[0]);
        r_val = Teuchos::rcp(new CubatureTensorPyr<DT,PT,WT>( xy_line, xy_line, z_line ));
      }
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( ( (topologyKey != shards::Line<>::key)               &&
                                      (topologyKey != shards::Triangle<>::key)           &&
                                      (topologyKey != shards::Quadrilateral<>::key)      &&
                                      (topologyKey != shards::ShellQuadrilateral<>::key) &&
                                      (topologyKey != shards::Tetrahedron<>::key)        &&
                                      (topologyKey != shards::Hexahedron<>::key)         &&
                                      (topologyKey != shards::Pyramid<>::key)            &&
                                      (topologyKey != shards::Wedge<>::key) ),
                                    std::invalid_argument,
                                    ">>> ERROR (DefaultCubatureFactory): Invalid cell topology prevents cubature creation.");
    }
    }
    return r_val;
  }

  template<typename DT, typename PT, typename WT>
  Teuchos::RCP<Cubature<DT,PT,WT> >
  DefaultCubatureFactory::
  create( const shards::CellTopology       cellTopology,
          const std::vector<ordinal_type> &degree,
          const EPolyType                  polytype) {
       return create<DT,PT,WT>(cellTopology.getBaseKey(), degree, polytype);
  }


  template<typename DT, typename PT, typename WT>
  Teuchos::RCP<Cubature<DT,PT,WT> > 
  DefaultCubatureFactory::
  create( unsigned topologyKey,
          const ordinal_type         degree,
          const EPolyType            polytype ) {
    // uniform order for 3 axes
    const std::vector<ordinal_type> degreeArray(3, degree);
    return create<DT,PT,WT>(topologyKey, degreeArray, polytype);
  }

  template<typename DT, typename PT, typename WT>
  Teuchos::RCP<Cubature<DT,PT,WT> >
  DefaultCubatureFactory::
  create( const shards::CellTopology cellTopology,
          const ordinal_type         degree,
          const EPolyType            polytype ) {
    // uniform order for 3 axes
    const std::vector<ordinal_type> degreeArray(3, degree);
    return create<DT,PT,WT>(cellTopology.getBaseKey(), degreeArray, polytype);
  }


  // template<typename DT>
  // template<typename cellVertexValueType, class ...cellVertexProperties>
  // Teuchos::RCP<Cubature<DT,PT,WT> > 
  // DefaultCubatureFactory::
  // create<DT,PT,WT>(const shards::CellTopology& cellTopology,
  //                       const Kokkos::DynRankView<cellVertexValueType,cellVertexProperties> cellVertices,
  //                       ordinal_type degree){
  //   return Teuchos::rcp(new CubaturePolygon<DT,PT,WT>(cellTopology,cellVertices, degree));
  // }

} // namespace Intrepid2

#endif
