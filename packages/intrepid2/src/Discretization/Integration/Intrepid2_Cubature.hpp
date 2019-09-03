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

/** \file   Intrepid2_Cubature.hpp
    \brief  Header file for the Intrepid2::Cubature class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_HPP__
#define __INTREPID2_CUBATURE_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2 {

  /* \struct Intrepid2::CubatureTemplate
      \brief  Template for the cubature rules used by Intrepid. Cubature template consists of
              cubature points and cubature weights. Intrepid provides a collection of cubature
              templates for most standard cell topologies. The templates are defined in reference
              coordinates using a standard reference cell for each canonical cell type. Cubature
              points are always specified by a triple of (X,Y,Z) coordinates even if the cell
              dimension is less than 3. The unused dimensions should be padded by zeroes.

              For example, a set of Gauss rules on [-1,1] looks as the following array of CubatureTemplate structs:

  \verbatim
  cubature_rule[4] =
  {                                         // Collection of Gauss rules on [-1,1]
    {
      1,                                      ----> number of points in the rule
      {{0.0,0.0,0.0}},                        ----> X,Y,Z coordinates of the cubature points
      {0.5}                                   ----> the cubature weight
    },
    {
      2,
      {{-sqrt(1.0/3.0),0.0,0.0},
       {+sqrt(1.0/3.0),0.0,0.0}},
      {1.0,1.0}
    },
    {
      3,
      {{-sqrt(3.0/5.0),0.0,0.0},
       {0.0,0.0,0.0},
       {+sqrt(3.0/5.0),0.0,0.0}},
      {5.0/9.0, 8.0/9.0,5.0/9.0}
    },
    {
      4,
      {{-sqrt((3.0+4.0*sqrt(0.3))/7.0),0.0,0.0},
       {-sqrt((3.0-4.0*sqrt(0.3))/7.0),0.0,0.0},
       {+sqrt((3.0-4.0*sqrt(0.3))/7.0),0.0,0.0},
       {+sqrt((3.0+4.0*sqrt(0.3))/7.0),0.0,0.0}},
      //
      {0.5-sqrt(10.0/3.0)/12.0,
       0.5+sqrt(10.0/3.0)/12.0,
       0.5+sqrt(10.0/3.0)/12.0,
       0.5-sqrt(10.0/3.0)/12.0}
    }
  };  // end Gauss
  \endverbatim

      Also see data member documentation.
  */


  /** \class Intrepid2::Cubature
      \brief Defines the base class for cubature (integration) rules in Intrepid.

      Cubature template (rule) consists of cubature points and cubature weights.
      Intrepid provides a small collection of frequently used cubature rule templates
      for FEM reconstructions on simplices (edge, tri, tet) and the pyramid cell,
      defined in the derived classes of CubatureDirect.

      For quad, hex, and triprism cells cubature templates are tensor products of CubatureDirect
      templates. The tensor-product cubatures are defined in the derived class CubatureTensor.
  */
  template<typename ExecSpaceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class Cubature {
  public:

    typedef Kokkos::DynRankView<pointValueType,Kokkos::LayoutStride,ExecSpaceType>   PointViewType;
    typedef Kokkos::DynRankView<weightValueType,Kokkos::LayoutStride,ExecSpaceType>  weightViewType;

    /** \brief Returns cubature points and weights
        (return arrays must be pre-sized/pre-allocated).

        \param cubPoints       [out]     - Array containing the cubature points.
        \param cubWeights      [out]     - Array of corresponding cubature weights.
    */
    virtual
    void
    getCubature( PointViewType  /* cubPoints */,
                 weightViewType /* cubWeights */ ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Cubature::getCubature): this method should be over-riden by derived classes.");
    }

    /** \brief Returns cubature points and weights on physical cells
        (return arrays must be pre-sized/pre-allocated).

        \param cubPoints       [out]     - Array containing the cubature points.
        \param cubWeights      [out]     - Array of corresponding cubature weights.
        \param cellVertices    [in]      - Array containing the cell vertices.
    */
    virtual
    void
    getCubature( PointViewType  /* cubPoints */,
                 weightViewType /* cubWeights */,
                 PointViewType  /* cellVertices */) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Cubature::getCubature): this method should be over-riden by derived classes.");
    }

    /** \brief Returns the number of cubature points.
     */
    virtual
    ordinal_type
    getNumPoints() const {
      INTREPID2_TEST_FOR_WARNING( true, 
                                  ">>> ERROR (Cubature::getNumPoints): this method should be over-riden by derived classes.");
      return 0;
    }


    /** \brief Returns dimension of the integration domain.
     */
    virtual
    ordinal_type
    getDimension() const {
      INTREPID2_TEST_FOR_WARNING( true, 
                                  ">>> ERROR (Cubature::getDimension): this method should be over-riden by derived classes.");
      return 0;
    }

    /** \brief Returns dimension of the integration domain.
     */
    virtual
    ordinal_type
    getAccuracy() const {
      INTREPID2_TEST_FOR_WARNING( true, 
                                  ">>> ERROR (Cubature::getDimension): this method should be over-riden by derived classes.");
      return 0;
    }

    /** \brief Returns cubature name.
     */
    virtual
    const char*
    getName() const  {
      return "Cubature";
    }

    Cubature() = default;
    virtual ~Cubature() {}

  };

}// end namespace Intrepid2


/*! \mainpage INTREPID2 Documentation (Development Version)

  \section intro_sec Introduction

  Intrepid2 is an extension of <a href="https://trilinos.org/packages/intrepid/">Intrepid</a>, 
  a library of interoperable tools for compatible discretizations of Partial
  Differential Equations (PDEs). Intrepid2 utilizes <a href="https://github.com/kokkos/">Kokkos</a>
  dynamic rank views as the default multidimensional array type, which enables the
  use of Intrepid2 on heterogeneous architectures. 

  \section overview_sec Overview

  Current release of %Intrepid2 includes the following features:
  \li Default finite element basis functions for <em>H(grad)</em>, <em>H(curl)</em>,
       <em>H(div)</em> and <em>H(vol)</em> spaces of orders up to 2 on standard cell
       topologies in 1D, 2D and 3D</li>
  \li High-order (up to 10) basis functions for <em>H(grad)</em>, <em>H(curl)</em>,
       <em>H(div)</em> and <em>H(vol)</em> spaces on select cell topologies</li>
  \li Pullbacks (transformations) from reference coordinate frame of <em>H(grad)</em>,
       <em>H(curl)</em>, <em>H(div)</em> and <em>H(vol)</em> fields</li>
  \li Pullbacks of gradient, curl and divergence of <em>H(grad)</em>, <em>H(curl)</em>,
       <em>H(div)</em> fields</li>
  \li Cubature rules of orders up to 20 on most standard 1D, 2D and 3D cell topologies</li>

  \section quickstart_sec Quick Start

  Familiarity with with the following concepts, objects, and tools is required:
  \li <a href="https://trilinos.org/packages/shards/">Shards</a> cell topologies,
  \li numerical integration / Intrepid2::Cubature,
  \li discrete (e.g. finite element) bases / Intrepid2::Basis / \ref basis_page,
  \li cell mappings and transformations / Intrepid2::CellTools / \ref cell_tools_page, and
  \li function mappings (pullbacks) / Intrepid2::FunctionSpaceTools / \ref function_space_tools_page.

  The following example demonstrates, in 7 steps, the computation of finite element
  stiffness matrices on a set of tetrahedral cells using a piecewise linear basis
  and an appropriate integration rule.

  \subsection topo_qs_sec Step 1: Select a cell topology

  \code
      shards::CellTopology cellType = shards::getCellTopologyData< shards::Tetrahedron<> >(); // cell type: tetrahedron
      ordinal_type spaceDim = cellType->getDimension();                                       // retrieve spatial dimension
      ordinal_type numNodes = cellType->getNodeCount();                                       // retrieve number of 0-cells (nodes)
  \endcode

  We additionally set the number of computational cells \c numCells.


  \subsection integration_qs_sec Step 2: Select integration (cubature) rule

  \code
      DefaultCubatureFactory<double> cubFactory;                                              // create cubature factory
      ordinal_type cubDegree = 2;                                                             // set cubature degree, e.g. 2
      Teuchos::RCP<Cubature<double> > myCub = cubFactory.create(cellType, cubDegree);         // create default cubature
      ordinal_type numCubPoints = myCub->getNumPoints();                                      // retrieve number of cubature points
  \endcode


  \subsection bases_qs_sec Step 3: Select discrete basis

  \code
      Basis_HGRAD_TET_C1_FEM<double, DynRankView<double> > tetBasis;                          // create tet basis
      ordinal_type numFields = tetBasis.getCardinality();                                     // get basis cardinality
  \endcode


  \subsection mdarray_qs_sec Step 4: Format multi-dimensional arrays

  \code
      DynRankView<double> cub_points(numCubPoints, spaceDim);
      DynRankView<double> cub_weights(numCubPoints);

      DynRankView<double> cell_nodes(numCells, numNodes, spaceDim);

      DynRankView<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
      DynRankView<double> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
      DynRankView<double> jacobian_det(numCells, numCubPoints);
      DynRankView<double> weighted_measure(numCells, numCubPoints);

      DynRankView<double> grad_at_cub_points(numFields, numCubPoints, spaceDim);
      DynRankView<double> transformed_grad_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
      DynRankView<double> weighted_transformed_grad_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
      DynRankView<double> stiffness_matrices(numCells, numFields, numFields);
  \endcode

  We assume that the array \c cell_nodes is filled with nodes defining a set of computational (physical) cells.

  \subsection tabulate_qs_sec Step 5: Evaluate differential operator applied to basis at cubature points

  \code
      myCub->getCubature(cub_points, cub_weights);                                          // retrieve cubature points and weights
      tetBasis.getValues(grad_at_cub_points, cub_points, OPERATOR_GRAD);                    // evaluate grad operator at cubature points
  \endcode


  \subsection ct_qs_sec Step 6: Apply cell tools

  \code
      CellTools<double>::setJacobian(jacobian, cub_points, cell_nodes, cellType);           // compute cell Jacobians
      CellTools<double>::setJacobianInv(jacobian_inv, jacobian);                            // compute inverses of cell Jacobians
      CellTools<double>::setJacobianDet(jacobian_det, jacobian);                            // compute determinants of cell Jacobians
  \endcode


  \subsection fst_qs_sec Step 7: Apply function space tools

  \code
      FunctionSpaceTools::computeCellMeasure<double>(weighted_measure,                      // compute weighted cell measure
                                                     jacobian_det,
                                                     cub_weights);
      FunctionSpaceTools::HGRADtransformGRAD<double>(transformed_grad_at_cub_points,        // transform reference gradients into physical space
                                                     jacobian_inv,
                                                     grad_at_cub_points);
      FunctionSpaceTools::multiplyMeasure<double>(weighted_transformed_grad_at_cub_points,  // multiply with weighted measure
                                                  weighted_measure,
                                                  transformed_grad_at_cub_points);
      FunctionSpaceTools::integrate<double>(stiffness_matrices,                             // compute stiffness matrices
                                            transformed_grad_at_cub_points,
                                            weighted_transformed_grad_at_cub_points,
                                            COMP_CPP);
  \endcode

  The computed (local) stiffness matrices can now be used in the assembly of a (global)
  discrete differential operator, e.g. a discrete Laplacian.

  \subsection doen_qs_sec Done!


*/

#endif
