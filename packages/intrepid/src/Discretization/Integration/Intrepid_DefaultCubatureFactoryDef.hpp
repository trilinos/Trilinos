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

/** \file   Intrepid_DefaultCubatureFactoryDef.hpp
\brief  Definition file for the class Intrepid::DefaultCubatureFactory.
\author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

// first create method
template<class Scalar, class ArrayType>
Teuchos::RCP<Cubature<Scalar,ArrayType> > DefaultCubatureFactory<Scalar,ArrayType>::create(
  const shards::CellTopology & cellTopology,
  const std::vector<int> & degree) {

  // Create generic cubature.
  Teuchos::RCP<Cubature<Scalar,ArrayType> > pickCubature;

  switch (cellTopology.getBaseTopology()->key) {

    case shards::Line<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 1), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      pickCubature = Teuchos::rcp(new CubatureDirectLineGauss<Scalar,ArrayType>(degree[0]));
      break;

    case shards::Triangle<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 1), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      pickCubature = Teuchos::rcp(new CubatureDirectTriDefault<Scalar,ArrayType>(degree[0]));
      break;

    case shards::Quadrilateral<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 2), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      {
      std::vector< Teuchos::RCP< Cubature<Scalar,ArrayType> > > lineCubs(2);
      lineCubs[0]  = Teuchos::rcp(new CubatureDirectLineGauss<double>(degree[0]));
      lineCubs[1]  = Teuchos::rcp(new CubatureDirectLineGauss<double>(degree[1]));
      pickCubature = Teuchos::rcp(new CubatureTensor<Scalar,ArrayType>(lineCubs));
      }
      break;

    case shards::Tetrahedron<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 1), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      pickCubature = Teuchos::rcp(new CubatureDirectTetDefault<Scalar,ArrayType>(degree[0]));
      break;

    case shards::Hexahedron<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 3), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.");
      {
      std::vector< Teuchos::RCP< Cubature<Scalar,ArrayType> > > lineCubs(3);
      lineCubs[0]  = Teuchos::rcp(new CubatureDirectLineGauss<double>(degree[0]));
      lineCubs[1]  = Teuchos::rcp(new CubatureDirectLineGauss<double>(degree[1]));
      lineCubs[2]  = Teuchos::rcp(new CubatureDirectLineGauss<double>(degree[2]));
      pickCubature = Teuchos::rcp(new CubatureTensor<Scalar,ArrayType>(lineCubs));
      }
      break;

    case shards::Wedge<>::key:
      TEST_FOR_EXCEPTION( (degree.size() < 2), std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Provided degree array is of insufficient length.")
      {
      std::vector< Teuchos::RCP< Cubature<Scalar,ArrayType> > > miscCubs(2);
      miscCubs[0]  = Teuchos::rcp(new CubatureDirectTriDefault<double>(degree[0]));
      miscCubs[1]  = Teuchos::rcp(new CubatureDirectLineGauss<double>(degree[1]));
      pickCubature = Teuchos::rcp(new CubatureTensor<Scalar,ArrayType>(miscCubs));
      }
      break;

    default:
      TEST_FOR_EXCEPTION( ( (cellTopology.getBaseTopology()->key != shards::Line<>::key)             &&
                            (cellTopology.getBaseTopology()->key != shards::Triangle<>::key)         &&
                            (cellTopology.getBaseTopology()->key != shards::Quadrilateral<>::key)    &&
                            (cellTopology.getBaseTopology()->key != shards::Tetrahedron<>::key)      &&
                            (cellTopology.getBaseTopology()->key != shards::Hexahedron<>::key)       &&
                            (cellTopology.getBaseTopology()->key != shards::Wedge<>::key) ),
                          std::invalid_argument,
                          ">>> ERROR (DefaultCubatureFactory): Invalid cell topology prevents cubature creation.");
  }

  return pickCubature;
}



// second create method
template<class Scalar, class ArrayType>
Teuchos::RCP<Cubature<Scalar,ArrayType> > DefaultCubatureFactory<Scalar,ArrayType>::create(
  const shards::CellTopology & cellTopology, int degree) {
  std::vector<int> d(3);
  d[0] = degree; d[1] = degree; d[2] = degree;
  return create(cellTopology, d);
}


} // namespace Intrepid
