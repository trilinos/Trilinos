// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_DefaultCubatureFactoryDef.hpp
\brief  Definition file for the class Intrepid::DefaultCubatureFactory.
\author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

// first create method
template<class Scalar, class ArrayPoint, class ArrayWeight>
Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > DefaultCubatureFactory<Scalar,ArrayPoint,ArrayWeight>::create(
  const shards::CellTopology & cellTopology,
  const std::vector<int> & degree) {

  // Create generic cubature.
  Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > pickCubature;

  switch (cellTopology.getBaseCellTopologyData()->key) {

    case shards::Line<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 1), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      pickCubature = Teuchos::rcp(new CubatureDirectLineGauss<Scalar,ArrayPoint,ArrayWeight>(degree[0]));
      break;

    case shards::Triangle<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 1), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      pickCubature = Teuchos::rcp(new CubatureDirectTriDefault<Scalar,ArrayPoint,ArrayWeight>(degree[0]));
      break;

    case shards::Quadrilateral<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 2), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      {
      std::vector< Teuchos::RCP< Cubature<Scalar,ArrayPoint,ArrayWeight> > > lineCubs(2);
      lineCubs[0]  = Teuchos::rcp(new CubatureDirectLineGauss<Scalar,ArrayPoint,ArrayWeight>(degree[0]));
      lineCubs[1]  = Teuchos::rcp(new CubatureDirectLineGauss<Scalar,ArrayPoint,ArrayWeight>(degree[1]));
      pickCubature = Teuchos::rcp(new CubatureTensor<Scalar,ArrayPoint,ArrayWeight>(lineCubs));
      }
      break;

    case shards::Tetrahedron<11>::key:
          TEST_FOR_EXCEPTION( (degree.size() < 1), std::invalid_argument,
                              ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
          pickCubature = Teuchos::rcp(new CubatureCompositeTet<Scalar,ArrayPoint,ArrayWeight>(degree[0]));
          break;

    case shards::Tetrahedron<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 1), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      pickCubature = Teuchos::rcp(new CubatureDirectTetDefault<Scalar,ArrayPoint,ArrayWeight>(degree[0]));
      break;

    case shards::Hexahedron<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 3), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      {
      std::vector< Teuchos::RCP< Cubature<Scalar,ArrayPoint,ArrayWeight> > > lineCubs(3);
      lineCubs[0]  = Teuchos::rcp(new CubatureDirectLineGauss<Scalar,ArrayPoint,ArrayWeight>(degree[0]));
      lineCubs[1]  = Teuchos::rcp(new CubatureDirectLineGauss<Scalar,ArrayPoint,ArrayWeight>(degree[1]));
      lineCubs[2]  = Teuchos::rcp(new CubatureDirectLineGauss<Scalar,ArrayPoint,ArrayWeight>(degree[2]));
      pickCubature = Teuchos::rcp(new CubatureTensor<Scalar,ArrayPoint,ArrayWeight>(lineCubs));
      }
      break;

    case shards::Wedge<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 2), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.")
      {
      std::vector< Teuchos::RCP< Cubature<Scalar,ArrayPoint,ArrayWeight> > > miscCubs(2);
      miscCubs[0]  = Teuchos::rcp(new CubatureDirectTriDefault<Scalar,ArrayPoint,ArrayWeight>(degree[0]));
      miscCubs[1]  = Teuchos::rcp(new CubatureDirectLineGauss<Scalar,ArrayPoint,ArrayWeight>(degree[1]));
      pickCubature = Teuchos::rcp(new CubatureTensor<Scalar,ArrayPoint,ArrayWeight>(miscCubs));
      }
      break;

    default:
      TEST_FOR_EXCEPTION( ( (cellTopology.getBaseCellTopologyData()->key != shards::Line<>::key)             &&
                            (cellTopology.getBaseCellTopologyData()->key != shards::Triangle<>::key)         &&
                            (cellTopology.getBaseCellTopologyData()->key != shards::Quadrilateral<>::key)    &&
                            (cellTopology.getBaseCellTopologyData()->key != shards::Tetrahedron<>::key)      &&
                            (cellTopology.getBaseCellTopologyData()->key != shards::Hexahedron<>::key)       &&
                            (cellTopology.getBaseCellTopologyData()->key != shards::Wedge<>::key) ),
                          std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Invalid cell topology prevents cubature creation.");
  }

  return pickCubature;
}



// second create method
template<class Scalar, class ArrayPoint, class ArrayWeight>
Teuchos::RCP<Cubature<Scalar,ArrayPoint,ArrayWeight> > DefaultCubatureFactory<Scalar,ArrayPoint,ArrayWeight>::create(
  const shards::CellTopology & cellTopology, int degree) {
  std::vector<int> d(3);
  d[0] = degree; d[1] = degree; d[2] = degree;
  return create(cellTopology, d);
}


} // namespace Intrepid
