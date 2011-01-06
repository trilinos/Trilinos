
/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <Intrepid_FieldContainer.hpp>
#define REDS1 1


#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"

#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include <Intrepid_Basis.hpp>
#include "Shards_CellTopology.hpp"


#include "Intrepid_CellTools.hpp"

int mymainpgi(int argc, char **argv) { 
#if 0
  if( new Intrepid::Basis_HGRAD_HEX_C1_FEM<double, Intrepid::FieldContainer<double> >())
    {
    }
#endif

  if (new Intrepid::DefaultCubatureFactory<double>)
    {
    }
#if 1
  //Intrepid::Basis<double, MDArray > BasisType;
  {
    using namespace Intrepid;
    using namespace shards;

    CellTopology cell_topo(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    unsigned numCells = 1;
    unsigned numNodes = cell_topo.getNodeCount();
    unsigned spaceDim = cell_topo.getDimension();

    // Rank-3 array with dimensions (C,N,D) for the node coordinates of 3 traingle cells
    FieldContainer<double> cellNodes(numCells, numNodes, spaceDim);

    DefaultCubatureFactory<double> cubFactory;                                              // create cubature factory
    unsigned cubDegree = 2;                                                                      // set cubature degree, e.g. 2
    Teuchos::RCP<Cubature<double> > myCub = cubFactory.create(cell_topo, cubDegree);         // create default cubature

    unsigned numCubPoints = myCub->getNumPoints();                                               // retrieve number of cubature points

    FieldContainer<double> cub_points(numCubPoints, spaceDim);
    FieldContainer<double> cub_weights(numCubPoints);

    // Rank-4 array (C,P,D,D) for the Jacobian and its inverse and Rank-2 array (C,P) for its determinant
    FieldContainer<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> jacobian_det(numCells, numCubPoints);

    myCub->getCubature(cub_points, cub_weights);                                          // retrieve cubature points and weights
    // Methods to compute cell Jacobians, their inverses and their determinants

    CellTools<double>::setJacobian(jacobian, cub_points, cellNodes, cell_topo);           // compute cell Jacobians
#if 0

    CellTools<double>::setJacobianInv(jacobian_inv, jacobian);                            // compute inverses of cell Jacobians
    CellTools<double>::setJacobianDet(jacobian_det, jacobian);                            // compute determinants of cell Jacobians
#endif
  }
#endif
  return 0;
}
