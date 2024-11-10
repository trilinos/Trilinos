// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_02.hpp
\brief  Patch test for the Intrepid2::Basis_HGRAD_HEX_C1_FEM class.
\author Created by P. Bochev, D. Ridzal, K. Peterson.
*/

#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

  namespace Test {


    template<typename viewOut, typename viewIn>
    void rhsFunc(viewOut, const viewIn, int, int, int);

    template<typename ValueType, typename DeviceSpaceType>
    void neumann(Kokkos::DynRankView<ValueType,DeviceSpaceType>  result,
                 const Kokkos::DynRankView<ValueType,DeviceSpaceType> points,
                 const Kokkos::DynRankView<ValueType,DeviceSpaceType> jacs,
                 const shards::CellTopology   & ,
                 int, int, int, int);

    template<typename viewOut, typename viewIn>
    void u_exact(viewOut, const viewIn, int, int, int);

    /// right-hand side function
    template<typename viewOut, typename viewIn>
    void rhsFunc(viewOut result,
                 const viewIn points,
                 int xd,
                 int yd,
                 int zd) {

      int x = 0, y = 1, z = 2;

      // second x-derivatives of u
      if (xd > 1) {
        for (size_type cell=0; cell<result.extent(0); cell++) {
          for (size_type pt=0; pt<result.extent(1); pt++) {
            result(cell,pt) = - xd*(xd-1)*std::pow(points(cell,pt,x), xd-2) *
                                std::pow(points(cell,pt,y), yd) * std::pow(points(cell,pt,z), zd);
          }
        }
      }

      // second y-derivatives of u
      if (yd > 1) {
        for (size_type cell=0; cell<result.extent(0); cell++) {
          for (size_type pt=0; pt<result.extent(1); pt++) {
            result(cell,pt) -=  yd*(yd-1)*std::pow(points(cell,pt,y), yd-2) *
                                std::pow(points(cell,pt,x), xd) * std::pow(points(cell,pt,z), zd);
          }
        }
      }

      // second z-derivatives of u
      if (zd > 1) {
        for (size_type cell=0; cell<result.extent(0); cell++) {
          for (size_type pt=0; pt<result.extent(1); pt++) {
            result(cell,pt) -=  zd*(zd-1)*std::pow(points(cell,pt,z), zd-2) *
                                std::pow(points(cell,pt,x), xd) * std::pow(points(cell,pt,y), yd);
          }
        }
      }

      // add u
      for (size_type cell=0; cell<result.extent(0); cell++) {
        for (size_type pt=0; pt<result.extent(1); pt++) {
          result(cell,pt) +=  std::pow(points(cell,pt,x), xd) * std::pow(points(cell,pt,y), yd) * std::pow(points(cell,pt,z), zd);
        }
      }

    }


    /// neumann boundary conditions
    template<typename ValueType, typename DeviceSpaceType>
    void neumann(Kokkos::DynRankView<ValueType,DeviceSpaceType>  result,
                 const Kokkos::DynRankView<ValueType,DeviceSpaceType> points,
                 const Kokkos::DynRankView<ValueType,DeviceSpaceType> jacs,
                 const shards::CellTopology   & parentCell,
                 int sideOrdinal, int xd, int yd, int zd) {

      int x = 0, y = 1, z = 2;

      int numCells  = result.extent(0);
      int numPoints = result.extent(1);

      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
      DynRankView grad_u("grad_u", numCells, numPoints, 3);
      DynRankView side_normals("side_normals", numCells, numPoints, 3);
      DynRankView normal_lengths("normal_lengths", numCells, numPoints);
      DynRankView side_normals_out("side_normals_out", numCells, numPoints, 3);  

      // first x-derivatives of u
      if (xd > 0) {
        for (int cell=0; cell<numCells; cell++) {
          for (int pt=0; pt<numPoints; pt++) {
            grad_u(cell,pt,x) = xd*std::pow(points(cell,pt,x), xd-1) *
                                std::pow(points(cell,pt,y), yd) * std::pow(points(cell,pt,z), zd);
          }
        }
      }

      // first y-derivatives of u
      if (yd > 0) {
        for (int cell=0; cell<numCells; cell++) {
          for (int pt=0; pt<numPoints; pt++) {
            grad_u(cell,pt,y) = yd*std::pow(points(cell,pt,y), yd-1) *
                                std::pow(points(cell,pt,x), xd) * std::pow(points(cell,pt,z), zd);
          }
        }
      }

      // first z-derivatives of u
      if (zd > 0) {
        for (int cell=0; cell<numCells; cell++) {
          for (int pt=0; pt<numPoints; pt++) {
            grad_u(cell,pt,z) = zd*std::pow(points(cell,pt,z), zd-1) *
                                std::pow(points(cell,pt,x), xd) * std::pow(points(cell,pt,y), yd);
          }
        }
      }
      
      CellTools<DeviceSpaceType>::getPhysicalSideNormals(side_normals, jacs, sideOrdinal, parentCell);

      // scale normals
      RealSpaceTools<DeviceSpaceType>::vectorNorm(normal_lengths, side_normals, NORM_TWO);
      FunctionSpaceTools<DeviceSpaceType>::scalarMultiplyDataData(side_normals_out, normal_lengths, side_normals, true); 

      FunctionSpaceTools<DeviceSpaceType>::dotMultiplyDataData(result, grad_u, side_normals_out);

    }

    /// exact solution
    template<typename viewOut, typename viewIn>
    void u_exact(viewOut result, const viewIn points, int xd, int yd, int zd) {
      int x = 0, y = 1, z = 2;
      for (size_type cell=0; cell<result.extent(0); cell++) {
        for (size_type pt=0; pt<result.extent(1); pt++) {
          result(cell,pt) = std::pow(points(pt,x), xd)*std::pow(points(pt,y), yd)*std::pow(points(pt,z), zd);
        }
      }
    }



    template<typename ValueType, typename DeviceSpaceType>
    int HGRAD_HEX_C1_FEM_Test02(const bool verbose) {
  
      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      *outStream 
        << "===============================================================================\n" \
        << "|                                                                             |\n" \
        << "|                    Unit Test (Basis_HGRAD_HEX_C1_FEM)                       |\n" \
        << "|                                                                             |\n" \
        << "|     1) Patch test involving mass and stiffness matrices,                    |\n" \
        << "|        for the Neumann problem on a physical parallelepiped                 |\n" \
        << "|        AND a reference hex Omega with boundary Gamma.                       |\n" \
        << "|                                                                             |\n" \
        << "|        - div (grad u) + u = f  in Omega,  (grad u) . n = g  on Gamma        |\n" \
        << "|                                                                             |\n" \
        << "|        For a generic parallelepiped, the basis recovers a complete          |\n" \
        << "|        polynomial space of order 1. On a (scaled and/or translated)         |\n" \
        << "|        reference hex, the basis recovers a complete tensor product          |\n" \
        << "|        space of order 1 (i.e. incl. xy, xz, yz, xyz term).                  |\n" \
        << "|                                                                             |\n" \
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
        << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
        << "|                                                                             |\n" \
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
        << "|                                                                             |\n" \
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;

      const ValueType tol = tolerence();
      int errorFlag = 0;

      typedef ValueType pointValueType;
      typedef ValueType weightValueType;

      typedef CubatureDirectLineGauss<DeviceSpaceType,pointValueType,weightValueType> CubatureLineType;
      typedef CubatureTensor<DeviceSpaceType,pointValueType,weightValueType> CubatureTensorType;
      
      typedef RealSpaceTools<DeviceSpaceType> rst;
      typedef FunctionSpaceTools<DeviceSpaceType> fst;
      typedef CellTools<DeviceSpaceType> cet;
  

      outStream -> precision(16);

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Patch test                                                          |\n"
        << "===============================================================================\n";



      try {

        int max_order = 1;                                                                    // max total order of polynomial solution
        //DefaultCubatureFactory<double>  cubFactory;                                           // create factory
        shards::CellTopology cell(shards::getCellTopologyData< shards::Hexahedron<> >());     // create parent cell topology
        shards::CellTopology side(shards::getCellTopologyData< shards::Quadrilateral<> >());  // create relevant subcell (side) topology
        int cellDim = cell.getDimension();
        int sideDim = side.getDimension();
        unsigned numSides = 6;

        // Define array containing points at which the solution is evaluated, on the reference tet.
        int numIntervals = 10;
        int numInterpPoints = (numIntervals + 1)*(numIntervals + 1)*(numIntervals + 1);
        DynRankView ConstructWithLabel(interp_points_ref, numInterpPoints, 3);
        int counter = 0;
        for (int k=0; k<=numIntervals; k++) {
          for (int j=0; j<=numIntervals; j++) {
            for (int i=0; i<=numIntervals; i++) {
              interp_points_ref(counter,0) = i*(1.0/numIntervals)-1.0;
              interp_points_ref(counter,1) = j*(1.0/numIntervals)-1.0;
              interp_points_ref(counter,2) = k*(1.0/numIntervals)-1.0;
              counter++;
            }
          }
        }

        /* Parent cell definition. */
        DynRankView ConstructWithLabel(cell_nodes, 2, 1, 8, cellDim);

        // Generic parallelepiped.
        cell_nodes(0, 0, 0, 0) = -5.0;
        cell_nodes(0, 0, 0, 1) = -1.0;
        cell_nodes(0, 0, 0, 2) = 0.0;
        cell_nodes(0, 0, 1, 0) = 4.0;
        cell_nodes(0, 0, 1, 1) = 1.0;
        cell_nodes(0, 0, 1, 2) = 1.0;
        cell_nodes(0, 0, 2, 0) = 8.0;
        cell_nodes(0, 0, 2, 1) = 3.0;
        cell_nodes(0, 0, 2, 2) = 1.0;
        cell_nodes(0, 0, 3, 0) = -1.0;
        cell_nodes(0, 0, 3, 1) = 1.0;
        cell_nodes(0, 0, 3, 2) = 0.0;
        cell_nodes(0, 0, 4, 0) = 5.0;
        cell_nodes(0, 0, 4, 1) = 9.0;
        cell_nodes(0, 0, 4, 2) = 1.0;
        cell_nodes(0, 0, 5, 0) = 14.0;
        cell_nodes(0, 0, 5, 1) = 11.0;
        cell_nodes(0, 0, 5, 2) = 2.0;
        cell_nodes(0, 0, 6, 0) = 18.0;
        cell_nodes(0, 0, 6, 1) = 13.0;
        cell_nodes(0, 0, 6, 2) = 2.0;
        cell_nodes(0, 0, 7, 0) = 9.0;
        cell_nodes(0, 0, 7, 1) = 11.0;
        cell_nodes(0, 0, 7, 2) = 1.0;
        // Reference hex.
        cell_nodes(1, 0, 0, 0) = -1.0;
        cell_nodes(1, 0, 0, 1) = -1.0;
        cell_nodes(1, 0, 0, 2) = -1.0;
        cell_nodes(1, 0, 1, 0) = 1.0;
        cell_nodes(1, 0, 1, 1) = -1.0;
        cell_nodes(1, 0, 1, 2) = -1.0;
        cell_nodes(1, 0, 2, 0) = 1.0;
        cell_nodes(1, 0, 2, 1) = 1.0;
        cell_nodes(1, 0, 2, 2) = -1.0;
        cell_nodes(1, 0, 3, 0) = -1.0;
        cell_nodes(1, 0, 3, 1) = 1.0;
        cell_nodes(1, 0, 3, 2) = -1.0;
        cell_nodes(1, 0, 4, 0) = -1.0;
        cell_nodes(1, 0, 4, 1) = -1.0;
        cell_nodes(1, 0, 4, 2) = 1.0;
        cell_nodes(1, 0, 5, 0) = 1.0;
        cell_nodes(1, 0, 5, 1) = -1.0;
        cell_nodes(1, 0, 5, 2) = 1.0;
        cell_nodes(1, 0, 6, 0) = 1.0;
        cell_nodes(1, 0, 6, 1) = 1.0;
        cell_nodes(1, 0, 6, 2) = 1.0;
        cell_nodes(1, 0, 7, 0) = -1.0;
        cell_nodes(1, 0, 7, 1) = 1.0;
        cell_nodes(1, 0, 7, 2) = 1.0;

        std::stringstream mystream[2];
        mystream[0].str("\n>> Now testing basis on a generic parallelepiped ...\n");
        mystream[1].str("\n>> Now testing basis on the reference hex ...\n");


        for (int pcell = 0; pcell < 2; pcell++) {
          *outStream << mystream[pcell].str();
          
          auto acell_nodes = Kokkos::subdynrankview(cell_nodes, pcell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
          DynRankView ConstructWithLabel(interp_points, 1, numInterpPoints, cellDim);
          cet::mapToPhysicalFrame(interp_points, interp_points_ref, acell_nodes, cell); //!!!!!!! cell_nodes[pcells]
          auto interp_points_sub = Kokkos::subdynrankview(interp_points,0, Kokkos::ALL(), Kokkos::ALL());

          for (int x_order=0; x_order <= max_order; x_order++) {
            int max_y_order = max_order;
            if (pcell == 0) {
              max_y_order -= x_order;
            }
            for (int y_order=0; y_order <= max_y_order; y_order++) {
              int max_z_order = max_order;
              if (pcell == 0) {
                max_z_order -= x_order;
                max_z_order -= y_order;
              }
              for (int z_order=0; z_order <= max_z_order; z_order++) {

                // evaluate exact solution
                DynRankView ConstructWithLabel(exact_solution, 1, numInterpPoints);
                u_exact(exact_solution, interp_points_sub, x_order, y_order, z_order);

                int basis_order = 1;

                // set test tolerance;
                double zero = basis_order*basis_order*basis_order*100*tol;

                //create basis
                // for virtual function, value and point types are declared in the class
                typedef ValueType outputValueType;

                Basis_HGRAD_HEX_C1_FEM<DeviceSpaceType,outputValueType,pointValueType> basis;
                //Teuchos::RCP<Basis<double,FieldContainer<double> > > basis =
                //  Teuchos::rcp(new Basis_HGRAD_HEX_C1_FEM<double,FieldContainer<double> >() );
                int numFields = basis.getCardinality();

                // create cubatures
                const auto x_line = CubatureLineType(2*basis_order);
                const auto y_line = CubatureLineType(2*basis_order);
                const auto z_line = CubatureLineType(2*basis_order);
                
                CubatureTensorType cellCub( x_line, y_line, z_line );
                CubatureTensorType sideCub( x_line, y_line );
                
                //Teuchos::RCP<Cubature<double> > cellCub = cubFactory.create(cell, 2*basis_order);
                //Teuchos::RCP<Cubature<double> > sideCub = cubFactory.create(side, 2*basis_order);
                int numCubPointsCell = cellCub.getNumPoints();
                int numCubPointsSide = sideCub.getNumPoints();

                /* Computational arrays. */
                /* Section 1: Related to parent cell integration. */
                DynRankView ConstructWithLabel(cub_points_cell, numCubPointsCell, cellDim);
                DynRankView ConstructWithLabel(cub_points_cell_physical, 1, numCubPointsCell, cellDim);
                DynRankView ConstructWithLabel(cub_weights_cell, numCubPointsCell);
                DynRankView ConstructWithLabel(jacobian_cell, 1, numCubPointsCell, cellDim, cellDim);
                DynRankView ConstructWithLabel(jacobian_inv_cell, 1, numCubPointsCell, cellDim, cellDim);
                DynRankView ConstructWithLabel(jacobian_det_cell, 1, numCubPointsCell);
                DynRankView ConstructWithLabel(weighted_measure_cell, 1, numCubPointsCell);

                DynRankView ConstructWithLabel(value_of_basis_at_cub_points_cell, numFields, numCubPointsCell);
                DynRankView ConstructWithLabel(transformed_value_of_basis_at_cub_points_cell, 1, numFields, numCubPointsCell);
                DynRankView ConstructWithLabel(weighted_transformed_value_of_basis_at_cub_points_cell, 1, numFields, numCubPointsCell);
                DynRankView ConstructWithLabel(grad_of_basis_at_cub_points_cell, numFields, numCubPointsCell, cellDim);
                DynRankView ConstructWithLabel(transformed_grad_of_basis_at_cub_points_cell, 1, numFields, numCubPointsCell, cellDim);
                DynRankView ConstructWithLabel(weighted_transformed_grad_of_basis_at_cub_points_cell, 1, numFields, numCubPointsCell, cellDim);
                DynRankView ConstructWithLabel(fe_matrix, 1, numFields, numFields);

                DynRankView ConstructWithLabel(rhs_at_cub_points_cell_physical, 1, numCubPointsCell);
                DynRankView ConstructWithLabel(rhs_and_soln_vector, 1, numFields);

                /* Section 2: Related to subcell (side) integration. */
                DynRankView ConstructWithLabel(cub_points_side, numCubPointsSide, sideDim);
                DynRankView ConstructWithLabel(cub_weights_side, numCubPointsSide);
                DynRankView ConstructWithLabel(cub_points_side_refcell, numCubPointsSide, cellDim);
                DynRankView ConstructWithLabel(cub_points_side_physical, 1, numCubPointsSide, cellDim);
                DynRankView ConstructWithLabel(jacobian_side_refcell, 1, numCubPointsSide, cellDim, cellDim);
                DynRankView ConstructWithLabel(scratch_side, numCubPointsSide*cellDim*cellDim);
                DynRankView ConstructWithLabel(jacobian_det_side_refcell, 1, numCubPointsSide);
                DynRankView ConstructWithLabel(weighted_measure_side_refcell, 1, numCubPointsSide);

                DynRankView ConstructWithLabel(value_of_basis_at_cub_points_side_refcell, numFields, numCubPointsSide);
                DynRankView ConstructWithLabel(transformed_value_of_basis_at_cub_points_side_refcell, 1, numFields, numCubPointsSide);
                DynRankView ConstructWithLabel(weighted_transformed_value_of_basis_at_cub_points_side_refcell, 1, numFields, numCubPointsSide);
                DynRankView ConstructWithLabel(neumann_data_at_cub_points_side_physical, 1, numCubPointsSide);
                DynRankView ConstructWithLabel(neumann_fields_per_side, 1, numFields);

                /* Section 3: Related to global interpolant. */
                DynRankView ConstructWithLabel(value_of_basis_at_interp_points_ref, numFields, numInterpPoints);
                DynRankView ConstructWithLabel(transformed_value_of_basis_at_interp_points_ref, 1, numFields, numInterpPoints);
                DynRankView ConstructWithLabel(interpolant, 1, numInterpPoints);

                



                /******************* START COMPUTATION ***********************/

                // get cubature points and weights
                cellCub.getCubature(cub_points_cell, cub_weights_cell);

                // compute geometric cell information
                cet::setJacobian(jacobian_cell, cub_points_cell, acell_nodes, cell);
                cet::setJacobianInv(jacobian_inv_cell, jacobian_cell);
                cet::setJacobianDet(jacobian_det_cell, jacobian_cell);

                // compute weighted measure
                fst::computeCellMeasure(weighted_measure_cell, jacobian_det_cell, cub_weights_cell);
     
                ///////////////////////////
                // Computing mass matrices:
                // tabulate values of basis functions at (reference) cubature points
                basis.getValues(value_of_basis_at_cub_points_cell, cub_points_cell, OPERATOR_VALUE);

                // transform values of basis functions 
                fst::HGRADtransformVALUE(transformed_value_of_basis_at_cub_points_cell,
                                                                value_of_basis_at_cub_points_cell);

                // multiply with weighted measure
                fst::multiplyMeasure(weighted_transformed_value_of_basis_at_cub_points_cell,
                                                            weighted_measure_cell,
                                                            transformed_value_of_basis_at_cub_points_cell);

                // compute mass matrices
                fst::integrate(fe_matrix,
                                                      transformed_value_of_basis_at_cub_points_cell,
                                                      weighted_transformed_value_of_basis_at_cub_points_cell);
                                                      
                ///////////////////////////

                ////////////////////////////////
                // Computing stiffness matrices:
                // tabulate gradients of basis functions at (reference) cubature points
                basis.getValues(grad_of_basis_at_cub_points_cell, cub_points_cell, OPERATOR_GRAD);

                // transform gradients of basis functions 
                fst::HGRADtransformGRAD(transformed_grad_of_basis_at_cub_points_cell,
                                                               jacobian_inv_cell,
                                                               grad_of_basis_at_cub_points_cell);

                // multiply with weighted measure
                fst::multiplyMeasure(weighted_transformed_grad_of_basis_at_cub_points_cell,
                                                            weighted_measure_cell,
                                                            transformed_grad_of_basis_at_cub_points_cell);

                // compute stiffness matrices and sum into fe_matrix
                fst::integrate(fe_matrix,
                                                      transformed_grad_of_basis_at_cub_points_cell,
                                                      weighted_transformed_grad_of_basis_at_cub_points_cell,
                                                      true);
                ////////////////////////////////

                ///////////////////////////////
                // Computing RHS contributions:
                // map cell (reference) cubature points to physical space
                cet::mapToPhysicalFrame(cub_points_cell_physical, cub_points_cell, acell_nodes, cell);

                // evaluate rhs function
                rhsFunc(rhs_at_cub_points_cell_physical, cub_points_cell_physical, x_order, y_order, z_order);

                // compute rhs
                fst::integrate(rhs_and_soln_vector,
                                                      rhs_at_cub_points_cell_physical,
                                                      weighted_transformed_value_of_basis_at_cub_points_cell);

                // compute neumann b.c. contributions and adjust rhs
                sideCub.getCubature(cub_points_side, cub_weights_side);
                for (unsigned i=0; i<numSides; i++) {
                  // compute geometric cell information
                  cet::mapToReferenceSubcell(cub_points_side_refcell, cub_points_side, sideDim, (int)i, cell);
                  cet::setJacobian(jacobian_side_refcell, cub_points_side_refcell, acell_nodes, cell);
                  cet::setJacobianDet(jacobian_det_side_refcell, jacobian_side_refcell);

                  // compute weighted face measure
                  fst::computeFaceMeasure(weighted_measure_side_refcell,
                                                                 jacobian_side_refcell,
                                                                 cub_weights_side,
                                                                 i, cell, scratch_side);

                  // tabulate values of basis functions at side cubature points, in the reference parent cell domain
                  basis.getValues(value_of_basis_at_cub_points_side_refcell, cub_points_side_refcell, OPERATOR_VALUE);
                  // transform 
                  fst::HGRADtransformVALUE(transformed_value_of_basis_at_cub_points_side_refcell,
                                                                  value_of_basis_at_cub_points_side_refcell);

                  // multiply with weighted measure
                  fst::multiplyMeasure(weighted_transformed_value_of_basis_at_cub_points_side_refcell,
                                                              weighted_measure_side_refcell,
                                                              transformed_value_of_basis_at_cub_points_side_refcell);

                  // compute Neumann data
                  // map side cubature points in reference parent cell domain to physical space
                  cet::mapToPhysicalFrame(cub_points_side_physical, cub_points_side_refcell, acell_nodes, cell);
                  // now compute data
                  
                  neumann(neumann_data_at_cub_points_side_physical, cub_points_side_physical, jacobian_side_refcell,
                          cell, (int)i, x_order, y_order, z_order);

                  fst::integrate(neumann_fields_per_side,
                                                        neumann_data_at_cub_points_side_physical,
                                                        weighted_transformed_value_of_basis_at_cub_points_side_refcell);

                  // adjust RHS
                  rst::add(rhs_and_soln_vector, neumann_fields_per_side);;
             
                }
                ///////////////////////////////
                
                /////////////////////////////
                // Solution of linear system:
                
                //std::vector<std::vector<std::vector<ValueType> > > _fe_matrix(1, std::vector<std::vector<ValueType> >(numFields, std::vector<ValueType>(numFields, 0)));
                std::vector<ValueType> _fe_matrix(numFields*numFields); 
                
                
                std::vector<int> ipiv(numFields);
                //std::vector<std::vector<ValueType> > _rhs_and_soln_vector(1, std::vector<ValueType>(numFields, 0));
                std::vector<ValueType> _rhs_and_soln_vector(numFields);
                
                
                for (auto i=0; i< numFields; i++) {
                  _rhs_and_soln_vector[i] = rhs_and_soln_vector(0,i);
                  for (auto j=0; j< numFields; j++)
                    _fe_matrix[i*numFields+j] = fe_matrix(0,i,j);
                }
                
                
                int info = 0;
                Teuchos::LAPACK<int, double> solver;
                solver.GESV(numFields, 1, &_fe_matrix[0], numFields, &ipiv[0], &_rhs_and_soln_vector[0], numFields, &info);
                
                for (auto i=0; i< numFields; i++)
                  rhs_and_soln_vector(0,i) = _rhs_and_soln_vector[i];
                
                /////////////////////////////

                ////////////////////////
                // Building interpolant:
                // evaluate basis at interpolation points
                basis.getValues(value_of_basis_at_interp_points_ref, interp_points_ref, OPERATOR_VALUE);
                // transform values of basis functions 
                fst::HGRADtransformVALUE(transformed_value_of_basis_at_interp_points_ref,
                                                                value_of_basis_at_interp_points_ref);
                fst::evaluate(interpolant, rhs_and_soln_vector, transformed_value_of_basis_at_interp_points_ref);
                ////////////////////////

                /******************* END COMPUTATION ***********************/
            
                rst::subtract(interpolant, exact_solution);

                auto interpolant_h = Kokkos::create_mirror_view(interpolant);
                Kokkos::deep_copy(interpolant_h, interpolant);
                ValueType relNorm = rst::Serial::vectorNorm(interpolant_h, NORM_TWO);
                              
                *outStream << "\nRelative norm-2 error between exact solution polynomial of order ("
                           << x_order << ", " << y_order << ", " << z_order
                           << ") and finite element interpolant of order " << basis_order << ": "
                           << relNorm << "\n";

                if (std::isnan(relNorm) || (relNorm > zero)) {
                  *outStream << "\n\nPatch test failed for solution polynomial order ("
                             << x_order << ", " << y_order << ", " << z_order << ") and basis order " << basis_order << "\n\n";
                  errorFlag++;
                }
              } // end for z_order
            } // end for y_order
          } // end for x_order
        } // end for pcell

      }
      // Catch unexpected errors
      catch (std::logic_error &err) {
            *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << "\n\n";
            *outStream << "-------------------------------------------------------------------------------" << "\n";
        errorFlag = -1000;
      };

      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";

      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      return errorFlag;
    }
  }
}
