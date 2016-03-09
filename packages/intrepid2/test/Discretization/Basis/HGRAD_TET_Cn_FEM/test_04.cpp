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

/**
   \file test_04.cpp
   \brief  Patch test for the Intrepid2::Basis_HGRAD_TRI_Cn_FEM class assembly.
   \author Kyungjoo Kim
*/

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Intrepid2_FieldContainer.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )
#include "Intrepid2_BasisSet.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_ToyMesh.hpp"
#endif

using namespace Intrepid2;

typedef double value_type;
typedef double magnitude_type;

// ----------------------------------------------------------------------
// Command line input variable
//
//
static int  maxp = 5; // INTREPID2_MAX_ORDER;
static int  nelement = 2;
static bool verbose = false;
static bool apply_orientation = true;

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )

// ----------------------------------------------------------------------
// Function declaration
//
//
#undef MatVal
#define MatVal(Mat, row, col) Mat.getData().get()[(row) + (col)*(Mat.dimension(1))]

void eval_exact(FieldContainer<value_type> &       val,
                const FieldContainer<value_type> & points,
                const int                          nx,
                const int                          ny,
                const int                          nz);

void eval_rhs(FieldContainer<value_type> &       f,
              const FieldContainer<value_type> & points,
              const int                          nx,
              const int                          ny,
              const int                          nz);

void eval_neumann(FieldContainer<value_type>       & neumann,
                  const FieldContainer<value_type> & points,
                  const FieldContainer<value_type> & jacs,
                  const shards::CellTopology   &     parent_cell,
                  const int                          side_ordinal,
                  const int                          nx,
                  const int                          ny,
                  const int                          nz);

void fill_cell_nodes(FieldContainer<value_type> & cell_nodes,
                     const FieldContainer<value_type> & nodes,
                     const int *element,
                     const int nvert,
                     const int ndim);

void build_element_matrix_and_rhs(FieldContainer<value_type> & A,
                                  FieldContainer<value_type> & b,
                                  DefaultCubatureFactory<value_type> & cubature_factory,
                                  const BasisSet_HGRAD_TRI_Cn_FEM<value_type,FieldContainer<value_type> > &basis_set,
                                  const int *element,
                                  const int *boundary,
                                  const FieldContainer<value_type> & cell_nodes,
                                  const Orientation ort,
                                  const int nx,
                                  const int ny,
                                  const int nz);

void assemble_element_matrix_and_rhs(FieldContainer<value_type> &A_asm,
                                     FieldContainer<value_type> &b_asm,
                                     const FieldContainer<value_type> &A,
                                     const FieldContainer<value_type> &b,
                                     const int (*l2g)[2],
                                     const int nnodes);

void retrieve_element_solution(FieldContainer<value_type> &sol,
                               const FieldContainer<value_type> &b_asm,
                               const int (*l2g)[2],
                               const int nnodes);

void compute_element_error(magnitude_type & interpolation_error,
                           magnitude_type & solution_norm,
                           const int *element,
                           const FieldContainer<value_type> & cell_nodes,
                           const BasisSet_HGRAD_TRI_Cn_FEM<value_type,FieldContainer<value_type> > &basis_set,
                           const FieldContainer<value_type> &sol,
                           const Orientation ort,
                           const int nx,
                           const int ny,
                           const int nz);

// ----------------------------------------------------------------------
// Function implementation
//
//
void fill_cell_nodes(FieldContainer<value_type> & cell_nodes,
                     const FieldContainer<value_type> & nodes,
                     const int *element,
                     const int nvert,
                     const int ndim) {
  for (int i=0;i<nvert;++i)
    for (int j=0;j<ndim;++j)
      cell_nodes(0, i, j) = nodes(0, element[i], j);
}

void eval_exact(FieldContainer<value_type> &       val,
                const FieldContainer<value_type> & points,
                const int                          nx,     // order x
                const int                          ny,     // order y
                const int                          nz) {   // order z
  // sol = x^nx + y^ny + z^nz
  //
  // result is a multi-vector (nsol x npts)
  // points are stored (npts x dim)
  //
  const int npts = val.dimension(1);
  for (int i=0;i<npts;++i) {
    const value_type x = points(i, 0), y = points(i, 1), z = points(i, 2);
    val(0, i) = std::pow(x, nx)*std::pow(y, ny)*std::pow(z, nz);
  }
}

void eval_rhs(FieldContainer<value_type> &       f,
              const FieldContainer<value_type> & points,
              const int                          nx,
              const int                          ny,
              const int                          nz) {
  // f = - \laplace u + u, where u is from exact solution
  //
  // f is a vector (1 x npts)
  // points is a multi-vector (1 x npts x dim)

  const int npts = f.dimension(1);

  // laplacian
  {
    for (int j=0;j<npts;++j)
      f(0, j) = 0.0;

    // second x-derivatives of u
    if (nx > 1) {
      for (int j=0;j<npts;++j) {
        const value_type x = points(0, j, 0), y = points(0, j, 1), z = points(0, j, 2);
        f(0, j) -= nx*(nx-1)*std::pow(x, nx-2) * std::pow(y, ny) * std::pow(z, nz);
      }
    }

    // second y-derivatives of u
    if (ny > 1) {
      for (int j=0;j<npts;++j) {
        const value_type x = points(0, j, 0), y = points(0, j, 1), z = points(0, j, 2);
        f(0, j) -=  ny*(ny-1)*std::pow(y, ny-2) * std::pow(z, nz) * std::pow(x, nx);
      }
    }

    // second z-derivatives of u
    if (nz > 1) {
      for (int j=0;j<npts;++j) {
        const value_type x = points(0, j, 0), y = points(0, j, 1), z = points(0, j, 2);
        f(0, j) -=  nz*(nz-1)*std::pow(z, nz-2) * std::pow(x, nx) * std::pow(y, ny);
      }
    }
  }

  // u
  for (int j=0;j<npts;++j) {
    const value_type x = points(0, j, 0), y = points(0, j, 1), z = points(0, j, 2);
    f(0, j) +=  std::pow(x, nx) * std::pow(y, ny) * std::pow(z, nz);
  }
}

void eval_neumann(FieldContainer<value_type>       & neumann,
                  const FieldContainer<value_type> & points,
                  const FieldContainer<value_type> & jacs,
                  const shards::CellTopology   &     cell,
                  const int                          side_ordinal,
                  const int                          nx,
                  const int                          ny,
                  const int                          nz) {
  // neumann = \grad u \dot normal
  const int npts = neumann.dimension(1);
  const int dim = 3;

  FieldContainer<value_type> grad_u(1, npts, dim);
  {
    // first x-derivatives of u
    if (nx > 0) {
      for (int j=0;j<npts;++j) {
        const value_type x = points(0, j, 0), y = points(0, j, 1), z = points(0, j, 2);
        grad_u(0, j, 0) = nx*std::pow(x, nx-1) * std::pow(y, ny) * std::pow(z, nz);
      }
    }

    // first y-derivatives of u
    if (ny > 0) {
      for (int j=0;j<npts;++j) {
        const value_type x = points(0, j, 0), y = points(0, j, 1), z = points(0, j, 2);
        grad_u(0, j, 1) = ny*std::pow(y, ny-1) * std::pow(z, nz) * std::pow(x, nx);
      }
    }

    // first z-derivatives of u
    if (nz > 0) {
      for (int j=0;j<npts;++j) {
        const value_type x = points(0, j, 0), y = points(0, j, 1), z = points(0, j, 2);
        grad_u(0, j, 2) = nz*std::pow(z, nz-1) * std::pow(x, nx) * std::pow(y, ny);
      }
    }
  }

  FieldContainer<value_type> side_normals(1, npts, dim);
  CellTools<value_type>::getPhysicalSideNormals(side_normals, jacs, side_ordinal, cell);

  // scale normals
  FieldContainer<value_type> normal_lengths(1, npts);
  
  // for minus jacobian set the direction opposite
  RealSpaceTools<value_type>::vectorNorm(normal_lengths, side_normals, NORM_TWO);
  FunctionSpaceTools::scalarMultiplyDataData<value_type>(side_normals, normal_lengths, side_normals, true);
  FunctionSpaceTools::dotMultiplyDataData<value_type>(neumann, grad_u, side_normals);
}

void build_element_matrix_and_rhs(FieldContainer<value_type> & A,
                                  FieldContainer<value_type> & b,
                                  DefaultCubatureFactory<value_type> & cubature_factory,
                                  const BasisSet_HGRAD_TET_Cn_FEM<value_type,FieldContainer<value_type> > &basis_set,
                                  const int *element,
                                  const int *boundary,
                                  const FieldContainer<value_type> & cell_nodes,
                                  const Orientation ort,
                                  const int nx,
                                  const int ny,
                                  const int nz) {

  // Step 0: initilization
  const auto &cell_basis = basis_set.getCellBasis();
  const auto &side_basis = basis_set.getTriangleBasis();

  const shards::CellTopology cell_topo = cell_basis.getBaseCellTopology();
  const shards::CellTopology side_topo = side_basis.getBaseCellTopology();

  const int nbf_cell = cell_basis.getCardinality();
  //const int nbf_side = side_basis.getCardinality();

  const int ndim_cell = cell_topo.getDimension();
  const int ndim_side = side_topo.getDimension();

  //const int nside = cell_topo.getEdgeCount();

  const int p = cell_basis.getDegree();

  // Step 1: create cubature data for integration
  Teuchos::RCP<Cubature<value_type> > cell_cub = cubature_factory.create(cell_topo, 2*p);
  Teuchos::RCP<Cubature<value_type> > side_cub = cubature_factory.create(side_topo, 2*p);

  const int npts_cell_cub = cell_cub->getNumPoints();
  const int npts_side_cub = side_cub->getNumPoints();

  // - cell related containers
  FieldContainer<value_type> cub_points_cell(npts_cell_cub, ndim_cell);
  FieldContainer<value_type> cub_points_cell_physical(1, npts_cell_cub, ndim_cell);
  FieldContainer<value_type> cub_weights_cell(npts_cell_cub);

  FieldContainer<value_type> jacobian_cell(1, npts_cell_cub, ndim_cell, ndim_cell);
  FieldContainer<value_type> jacobian_inv_cell(1, npts_cell_cub, ndim_cell, ndim_cell);
  FieldContainer<value_type> jacobian_det_cell(1, npts_cell_cub);

  FieldContainer<value_type> weighted_measure_cell(1, npts_cell_cub);

  FieldContainer<value_type> value_of_basis_at_cub_points_cell(nbf_cell, npts_cell_cub);
  FieldContainer<value_type> value_of_reordered_basis_at_cub_points_cell(nbf_cell, npts_cell_cub);

  FieldContainer<value_type> transformed_value_of_basis_at_cub_points_cell(1, nbf_cell, npts_cell_cub);
  FieldContainer<value_type> weighted_transformed_value_of_basis_at_cub_points_cell(1, nbf_cell, npts_cell_cub);

  FieldContainer<value_type> grad_of_basis_at_cub_points_cell(nbf_cell, npts_cell_cub, ndim_cell);
  FieldContainer<value_type> grad_of_reordered_basis_at_cub_points_cell(nbf_cell, npts_cell_cub, ndim_cell);

  FieldContainer<value_type> transformed_grad_of_basis_at_cub_points_cell(1, nbf_cell, npts_cell_cub, ndim_cell);
  FieldContainer<value_type> weighted_transformed_grad_of_basis_at_cub_points_cell(1, nbf_cell, npts_cell_cub, ndim_cell);

  FieldContainer<value_type> rhs_at_cub_points_cell_physical(1, npts_cell_cub);
  FieldContainer<value_type> rhs_and_soln_vector(1, nbf_cell);

  // - subcell related containders
  FieldContainer<value_type> cub_points_side(npts_side_cub, ndim_side);
  FieldContainer<value_type> cub_weights_side(npts_side_cub);
  FieldContainer<value_type> cub_points_side_refcell(npts_side_cub, ndim_cell);
  FieldContainer<value_type> cub_points_side_physical(1, npts_side_cub, ndim_cell);

  FieldContainer<value_type> jacobian_side_refcell(1, npts_side_cub, ndim_cell, ndim_cell);
  FieldContainer<value_type> jacobian_det_side_refcell(1, npts_side_cub);

  FieldContainer<value_type> weighted_measure_side_refcell(1, npts_side_cub);

  FieldContainer<value_type> value_of_basis_at_cub_points_side_refcell(nbf_cell, npts_side_cub);
  FieldContainer<value_type> value_of_reordered_basis_at_cub_points_side_refcell(nbf_cell, npts_side_cub);

  FieldContainer<value_type> transformed_value_of_basis_at_cub_points_side_refcell(1, nbf_cell, npts_side_cub);
  FieldContainer<value_type> weighted_transformed_value_of_basis_at_cub_points_side_refcell(1, nbf_cell, npts_side_cub);

  FieldContainer<value_type> neumann_data_at_cub_points_side_physical(1, npts_side_cub);
  FieldContainer<value_type> neumann_fields_per_side(1, nbf_cell);

  // get cubature points and weights
  cell_cub->getCubature(cub_points_cell, cub_weights_cell);

  CellTools<value_type>::setJacobian   (jacobian_cell, cub_points_cell, cell_nodes, cell_topo);
  CellTools<value_type>::setJacobianInv(jacobian_inv_cell, jacobian_cell);
  CellTools<value_type>::setJacobianDet(jacobian_det_cell, jacobian_cell);

  // compute weighted measure
  FunctionSpaceTools::computeCellMeasure<value_type>(weighted_measure_cell,
                                                     jacobian_det_cell, cub_weights_cell);

  // Step 1: mass matrix: tabulate values of basis functions at cubature points
  cell_basis.getValues(value_of_basis_at_cub_points_cell, cub_points_cell, OPERATOR_VALUE);

  if (apply_orientation) {
    OrientationTools<value_type>::verbose = false;
    OrientationTools<value_type>::getBasisFunctionsByTopology(value_of_reordered_basis_at_cub_points_cell,
                                                              value_of_basis_at_cub_points_cell,
                                                              cell_basis);
    OrientationTools<value_type>::getModifiedBasisFunctions(value_of_basis_at_cub_points_cell,
                                                            value_of_reordered_basis_at_cub_points_cell,
                                                            basis_set,
                                                            ort);
    OrientationTools<value_type>::verbose = false;
  }

  // transform values of basis functions
  FunctionSpaceTools::HGRADtransformVALUE<value_type>(transformed_value_of_basis_at_cub_points_cell,
                                                      value_of_basis_at_cub_points_cell);


  // multiply with weighted measure
  FunctionSpaceTools::multiplyMeasure<value_type>(weighted_transformed_value_of_basis_at_cub_points_cell,
                                                  weighted_measure_cell,
                                                  transformed_value_of_basis_at_cub_points_cell);

  // integrate
  FunctionSpaceTools::integrate<value_type>(A,
                                            transformed_value_of_basis_at_cub_points_cell,
                                            weighted_transformed_value_of_basis_at_cub_points_cell,
                                            COMP_CPP);

  // Step 2: stiffness matrix: tabulate grad values of basis functions at cubature points
  cell_basis.getValues(grad_of_basis_at_cub_points_cell, cub_points_cell, OPERATOR_GRAD);
  if (apply_orientation) {
    OrientationTools<value_type>::getBasisFunctionsByTopology(grad_of_reordered_basis_at_cub_points_cell,
                                                              grad_of_basis_at_cub_points_cell,
                                                              cell_basis);
    OrientationTools<value_type>::getModifiedBasisFunctions(grad_of_basis_at_cub_points_cell,
                                                            grad_of_reordered_basis_at_cub_points_cell,
                                                            basis_set,
                                                            ort);
  }

  // transform gradients of basis functions
  FunctionSpaceTools::HGRADtransformGRAD<value_type>(transformed_grad_of_basis_at_cub_points_cell,
                                                     jacobian_inv_cell,
                                                     grad_of_basis_at_cub_points_cell);

  // multiply with weighted measure
  FunctionSpaceTools::multiplyMeasure<value_type>(weighted_transformed_grad_of_basis_at_cub_points_cell,
                                                  weighted_measure_cell,
                                                  transformed_grad_of_basis_at_cub_points_cell);

  // compute stiffness matrices and sum into fe_matrix
  FunctionSpaceTools::integrate<value_type>(A,
                                            transformed_grad_of_basis_at_cub_points_cell,
                                            weighted_transformed_grad_of_basis_at_cub_points_cell,
                                            COMP_CPP,
                                            true);

  // Step 3: compute rhs function
  CellTools<value_type>::mapToPhysicalFrame(cub_points_cell_physical, cub_points_cell, cell_nodes, cell_topo);

  // evaluate rhs function
  eval_rhs(rhs_at_cub_points_cell_physical,
           cub_points_cell_physical,
           nx, ny, nz);

  // compute rhs
  FunctionSpaceTools::integrate<value_type>(b,
                                            rhs_at_cub_points_cell_physical,
                                            weighted_transformed_value_of_basis_at_cub_points_cell,
                                            COMP_CPP);

  // Step 4: compute boundary condition
  side_cub->getCubature(cub_points_side, cub_weights_side);
  const int nside = cell_topo.getSideCount();
  for (int i=0;i<nside;++i) {
    if (boundary[i]) {
      // compute geometric cell information
      CellTools<value_type>::mapToReferenceSubcell(cub_points_side_refcell, cub_points_side, ndim_side, i, cell_topo);
      CellTools<value_type>::setJacobian   (jacobian_side_refcell, cub_points_side_refcell, cell_nodes, cell_topo);
      CellTools<value_type>::setJacobianDet(jacobian_det_side_refcell, jacobian_side_refcell);

      // compute weighted edge measure
      FunctionSpaceTools::computeFaceMeasure<value_type>(weighted_measure_side_refcell,
                                                         jacobian_side_refcell,
                                                         cub_weights_side,
                                                         i,
                                                         cell_topo);

      // tabulate values of basis functions at side cubature points, in the reference parent cell domain
      cell_basis.getValues(value_of_basis_at_cub_points_side_refcell, cub_points_side_refcell, OPERATOR_VALUE);
      if (apply_orientation) {
        OrientationTools<value_type>::getBasisFunctionsByTopology(value_of_reordered_basis_at_cub_points_side_refcell,
                                                                  value_of_basis_at_cub_points_side_refcell,
                                                                  cell_basis);
        OrientationTools<value_type>::getModifiedBasisFunctions(value_of_basis_at_cub_points_side_refcell,
                                                                value_of_reordered_basis_at_cub_points_side_refcell,
                                                                basis_set,
                                                                ort);
      }

      // transform
      FunctionSpaceTools::HGRADtransformVALUE<value_type>(transformed_value_of_basis_at_cub_points_side_refcell,
                                                          value_of_basis_at_cub_points_side_refcell);

      // multiply with weighted measure
      FunctionSpaceTools::multiplyMeasure<value_type>(weighted_transformed_value_of_basis_at_cub_points_side_refcell,
                                                      weighted_measure_side_refcell,
                                                      transformed_value_of_basis_at_cub_points_side_refcell);

      // compute neumann boundary

      // map side cubature points in reference parent cell domain to physical space
      CellTools<value_type>::mapToPhysicalFrame(cub_points_side_physical, cub_points_side_refcell, cell_nodes, cell_topo);

      // now compute data
      eval_neumann(neumann_data_at_cub_points_side_physical,
                   cub_points_side_physical,
                   jacobian_side_refcell,
                   cell_topo,
                   i,
                   nx, ny, nz);

      FunctionSpaceTools::integrate<value_type>(neumann_fields_per_side,
                                                neumann_data_at_cub_points_side_physical,
                                                weighted_transformed_value_of_basis_at_cub_points_side_refcell,
                                                COMP_CPP);

      // adjust rhs
      RealSpaceTools<value_type>::add(b, neumann_fields_per_side);;
    }
  }
}

void assemble_element_matrix_and_rhs(FieldContainer<value_type> &A_asm,
                                     FieldContainer<value_type> &b_asm,
                                     const FieldContainer<value_type> &A,
                                     const FieldContainer<value_type> &b,
                                     const int (*l2g)[2],
                                     const int nnodes) {
  // nnodes represent the number of topological nodes
  const int local = 0, global = 1;
  for (int i=0;i<nnodes;++i) {
    const int
      l_offm = l2g[i][local],
      g_offm = l2g[i][global],
      ndofm  = l2g[i+1][local] - l_offm;

    for (int km=0;km<ndofm;++km)
      MatVal(b_asm, g_offm+km, 0) += MatVal(b, l_offm+km, 0);

    for (int j=0;j<nnodes;++j) {
      const int
        l_offn = l2g[j][local],
        g_offn = l2g[j][global],
        ndofn  = l2g[j+1][local] - l_offn;

      for (int km=0;km<ndofm;++km)
        for (int kn=0;kn<ndofn;++kn)
          MatVal(A_asm, g_offm+km, g_offn+kn) += MatVal(A, l_offm+km, l_offn+kn);
    }
  }
}

void retrieve_element_solution(FieldContainer<value_type> &b,
                               const FieldContainer<value_type> &b_asm,
                               const int (*l2g)[2],
                               const int nnodes) {
  // nnodes represent the number of topological nodes
  const int local = 0, global = 1;
  for (int i=0;i<nnodes;++i) {
    const int
      l_offm = l2g[i][local],
      g_offm = l2g[i][global],
      ndofm  = l2g[i+1][local] - l_offm;

    for (int km=0;km<ndofm;++km)
      MatVal(b, l_offm+km, 0) = MatVal(b_asm, g_offm+km, 0);
  }
}


void compute_element_error(magnitude_type & interpolation_error,
                           magnitude_type & solution_norm,
                           const int *element,
                           const FieldContainer<value_type> & cell_nodes,
                           const BasisSet_HGRAD_TET_Cn_FEM<value_type,FieldContainer<value_type> > &basis_set,
                           const FieldContainer<value_type> &sol,
                           const Orientation ort,
                           const int nx,
                           const int ny,
                           const int nz) {
  // initialize return values
  interpolation_error = 0.0;
  solution_norm = 0.0;

  // general environment
  const auto &cell_basis = basis_set.getCellBasis();
  const shards::CellTopology cell_topo = cell_basis.getBaseCellTopology();

  const int nbf_cell = cell_basis.getCardinality();
  const int ndim_cell = cell_topo.getDimension();

  // create points to evaluate in the reference cell
  const int order = 10;
  const int npts = PointTools::getLatticeSize(cell_topo, order, 1);

  FieldContainer<value_type> ref_cell_pts(npts, ndim_cell);
  PointTools::getLattice<value_type>(ref_cell_pts,
                                     cell_topo,
                                     order, 1);

  // map the points to physical frame
  FieldContainer<value_type> phy_cell_pts(1, npts, ndim_cell);
  CellTools<value_type>::mapToPhysicalFrame(phy_cell_pts, ref_cell_pts, cell_nodes, cell_topo);
  phy_cell_pts.resize(npts, ndim_cell);

  // Step 1: compute L2 error

  // evaluate exact solution
  FieldContainer<double> exact_solution_val(1, npts);
  eval_exact(exact_solution_val,
             phy_cell_pts,
             nx, ny, nz);

  // evaluate basis at interpolation points
  FieldContainer<value_type> value_of_basis_at_ref_cell_pts(nbf_cell, npts);
  FieldContainer<value_type> value_of_reordered_basis_at_ref_cell_pts(nbf_cell, npts);

  cell_basis.getValues(value_of_basis_at_ref_cell_pts, ref_cell_pts, OPERATOR_VALUE);
  if (apply_orientation) {
    OrientationTools<value_type>::getBasisFunctionsByTopology(value_of_reordered_basis_at_ref_cell_pts,
                                                              value_of_basis_at_ref_cell_pts,
                                                              cell_basis);
    OrientationTools<value_type>::getModifiedBasisFunctions(value_of_basis_at_ref_cell_pts,
                                                            value_of_reordered_basis_at_ref_cell_pts,
                                                            basis_set,
                                                            ort);
  }

  // transform values of basis functions
  FieldContainer<double> transformed_value_of_basis_at_ref_cell_pts(1, nbf_cell, npts);
  FunctionSpaceTools::HGRADtransformVALUE<value_type>(transformed_value_of_basis_at_ref_cell_pts,
                                                      value_of_basis_at_ref_cell_pts);

  FieldContainer<double> interpolant(1, npts);
  FunctionSpaceTools::evaluate<value_type>(interpolant,
                                           sol,
                                           transformed_value_of_basis_at_ref_cell_pts);

  // compute error and magnitude of solution
  RealSpaceTools<value_type>::subtract(interpolant, exact_solution_val);

  interpolation_error += RealSpaceTools<value_type>::vectorNorm(interpolant.getData().get(),
                                                                interpolant.dimension(1),
                                                                NORM_TWO);

  solution_norm += RealSpaceTools<value_type>::vectorNorm(exact_solution_val.getData().get(),
                                                          exact_solution_val.dimension(1),
                                                          NORM_TWO);

  // Step 2: compute H1 error
  // skip for now, not meaningful for this unit test
}

#endif

// ----------------------------------------------------------------------
// Main
//
//
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Kokkos::initialize();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;

  for (int i=0;i<argc;++i) {
    if ((strcmp(argv[i],"--nelement")          == 0)) { nelement = atoi(argv[++i]); continue;}
    if ((strcmp(argv[i],"--apply-orientation") == 0)) { apply_orientation = atoi(argv[++i]); continue;}
    if ((strcmp(argv[i],"--verbose")           == 0)) { verbose  = atoi(argv[++i]); continue;}
    if ((strcmp(argv[i],"--maxp")              == 0)) { maxp     = atoi(argv[++i]); continue;}
  }

  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  *outStream << std::scientific;
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                    Unit Test (Basis_HGRAD_TET_Cn_FEM)                       |\n" \
    << "|                                                                             |\n" \
    << "|     1) Patch test involving mass and stiffness matrices,                    |\n" \
    << "|        for the Neumann problem on a triangular patch                        |\n" \
    << "|        Omega with boundary Gamma.                                           |\n" \
    << "|                                                                             |\n" \
    << "|        - div (grad u) + u = f  in Omega,  (grad u) . n = g  on Gamma        |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                      Kyungjoo Kim  (kyukim@sandia.gov).                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n" \
    << "| TEST 4: Patch test for high order assembly                                  |\n" \
    << "===============================================================================\n";

  int r_val = 0;

  // precision control
  outStream->precision(3);

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )

  try {
    // test setup
    const int ndim = 3;
    FieldContainer<value_type> base_nodes(1, 5, ndim);
    base_nodes(0, 0, 0) = 0.0;
    base_nodes(0, 0, 1) = 0.0;
    base_nodes(0, 0, 2) = 0.0;
    
    base_nodes(0, 1, 0) = 1.0;
    base_nodes(0, 1, 1) = 0.0;
    base_nodes(0, 1, 2) = 0.0;
    
    base_nodes(0, 2, 0) = 0.0;
    base_nodes(0, 2, 1) = 1.0;
    base_nodes(0, 2, 2) = 0.0;
    
    base_nodes(0, 3, 0) = 0.0;
    base_nodes(0, 3, 1) = 0.0;
    base_nodes(0, 3, 2) = 1.0;

    base_nodes(0, 4, 0) = 1.0;
    base_nodes(0, 4, 1) = 1.0;
    base_nodes(0, 4, 2) = 1.0;

    // element 0 has globally permuted face node
    const int elt_0[6][4] = { { 0, 1, 2, 3 },
                              { 0, 3, 1, 2 }, 
                              { 0, 2, 3, 1 }, 
                              { 0, 1, 3, 2 }, 
                              { 0, 3, 2, 1 },
                              { 0, 2, 1, 3 } };
    
    // element 1 is locally permuted 
    int elt_1[4] = { 1, 2, 3 ,4 };

    DefaultCubatureFactory<value_type> cubature_factory;

    // for all test orders
    for (int nx=0;nx<=maxp;++nx) {
      for (int ny=0;ny<=maxp-nx;++ny) {
        for (int nz=0;nz<=maxp-nx-ny;++nz) {
          // polynomial order of approximation
          const int minp = std::max(nx+ny+nz, 1);
          
          // test for all basis above order p
          const EPointType pointtype[] = { POINTTYPE_EQUISPACED, POINTTYPE_WARPBLEND };
          for (int ptype=0;ptype<1;++ptype) {
            for (int p=minp;p<=maxp;++p) {
              *outStream << "\n"                                        \
                         << "===============================================================================\n" \
                         << "  Order (nx,ny,nz,p) = " << nx << ", " << ny << ", " << nz << ", " << p << " , PointType = " << EPointTypeToString(pointtype[ptype]) << "\n" \
                         << "===============================================================================\n";
              
              BasisSet_HGRAD_TET_Cn_FEM<value_type,FieldContainer<value_type> > basis_set(p, pointtype[ptype]);
              const auto& basis = basis_set.getCellBasis();
              const shards::CellTopology cell = basis.getBaseCellTopology();

              const int nbf = basis.getCardinality();

              const int nvert = cell.getVertexCount();
              const int nedge = cell.getEdgeCount();
              const int nface = cell.getFaceCount();
              

              FieldContainer<value_type> nodes(1, 5, ndim); 
              FieldContainer<value_type> cell_nodes(1, nvert, ndim);

              // ignore the subdimension; the matrix is always considered as 1D array
              FieldContainer<value_type> A(1, nbf, nbf), b(1, nbf);
              
              // ***** Test for different orientations *****
              for (int conf0=3;conf0<6;++conf0) { 
                for (int ino=0;ino<4;++ino) {                                                                   
                  nodes(0, elt_0[conf0][ino], 0) = base_nodes(0, ino, 0);                                       
                  nodes(0, elt_0[conf0][ino], 1) = base_nodes(0, ino, 1);                                       
                  nodes(0, elt_0[conf0][ino], 2) = base_nodes(0, ino, 2);                                       
                }                                                                                               
                nodes(0, 4, 0) = base_nodes(0, 4, 0);                                                           
                nodes(0, 4, 1) = base_nodes(0, 4, 1);   
                nodes(0, 4, 2) = base_nodes(0, 4, 2);   

                elt_1[0] = 1;
                elt_1[1] = 2;
                elt_1[2] = 3;
                elt_1[3] = 4;

                // for all permutations of element 1
                for (int conf1=0;conf1<24;++conf1) {
                  // filter out left handed element                                                             
                  fill_cell_nodes(cell_nodes,                                                                   
                                  nodes,                                                                        
                                  elt_1,                                                                        
                                  nvert, ndim);                                                                 
                  if (OrientationTools<value_type>::isLeftHandedCell(cell_nodes)) {      
                    // skip left handed
                  } else {
                    const int *element[2] = { elt_0[conf0], elt_1 };
                    *outStream << "\n"                      \
                               << "  Element 0 is configured " << conf0 << " "
                               << "(" << element[0][0] 
                               << "," << element[0][1] 
                               << "," << element[0][2] 
                               << "," << element[0][3] 
                               << ")"
                               << "  Element 1 is configured " << conf1 << " "
                               << "(" << element[1][0] 
                               << "," << element[1][1] 
                               << "," << element[1][2] 
                               << "," << element[1][3] 
                               << ")"
                               << "\n";
              
                    if (verbose) {
                      *outStream << " - Element nodal connectivity - \n";
                      for (int iel=0;iel<nelement;++iel)
                        *outStream << " iel = " << std::setw(4) << iel
                                   << ", nodes = "
                                   << std::setw(4) << element[iel][0]
                                   << std::setw(4) << element[iel][1]
                                   << std::setw(4) << element[iel][2]
                                   << std::setw(4) << element[iel][3]
                                   << "\n";
                    }

                    // Step 0: count one-to-one mapping between high order nodes and dofs
                    Example::ToyMesh mesh;
                    int local2global[2][16][2], boundary[2][4], off_global = 0;

                    const int nnodes_per_element
                      = cell.getVertexCount()
                      + cell.getEdgeCount()
                      + cell.getFaceCount()
                      + 1;

                    for (int iel=0;iel<nelement;++iel) 
                      mesh.getLocalToGlobalMap(local2global[iel], off_global, basis, element[iel]);

                    for (int iel=0;iel<nelement;++iel)
                      mesh.getBoundaryFlags(boundary[iel], cell, element[iel]);

                    if (verbose) {
                      *outStream << " - Element one-to-one local2global map -\n";
                      for (int iel=0;iel<nelement;++iel) {
                        *outStream << " iel = " << std::setw(4) << iel << "\n";
                        for (int i=0;i<(nnodes_per_element+1);++i) {
                          *outStream << "   local = " << std::setw(4) << local2global[iel][i][0]
                                     << "   global = " << std::setw(4) << local2global[iel][i][1]
                                     << "\n";
                        }
                      }
                      *outStream << " - Element boundary flags -\n";
                      const int nside = cell.getSideCount();
                      for (int iel=0;iel<nelement;++iel) {
                        *outStream << " iel = " << std::setw(4) << iel << "\n";
                        for (int i=0;i<nside;++i) {
                          *outStream << "   side = " << std::setw(4) << i
                                     << "   boundary = " << std::setw(4) << boundary[iel][i]
                                     << "\n";
                        }
                      }
                    }

                    // Step 1: assembly
                    const int ndofs = off_global;
                    FieldContainer<value_type> A_asm(1, ndofs, ndofs), b_asm(1, ndofs);
              
                    for (int iel=0;iel<nelement;++iel) {
                      // Step 1.1: create element matrices
                      Orientation ort = Orientation::getOrientation(cell, element[iel]);
                      fill_cell_nodes(cell_nodes,
                                      nodes,
                                      element[iel],
                                      nvert,
                                      ndim);

                      build_element_matrix_and_rhs(A, b,
                                                   cubature_factory,
                                                   basis_set,
                                                   element[iel],
                                                   boundary[iel],
                                                   cell_nodes,
                                                   ort,
                                                   nx, ny, nz);

                      // if p is bigger than 4, not worth to look at the matrix
                      if (verbose && p < 5) {
                        *outStream << " - Element matrix and rhs, iel = " << iel << "\n";
                        *outStream << std::showpos;
                        for (int i=0;i<nbf;++i) {
                          for (int j=0;j<nbf;++j)
                            *outStream << MatVal(A, i, j) << " ";
                          *outStream << ":: " << MatVal(b, i, 0) << "\n";
                        }
                        *outStream << std::noshowpos;
                      }

                      // Step 1.2: assemble high order elements
                      assemble_element_matrix_and_rhs(A_asm, b_asm,
                                                      A, b,
                                                      local2global[iel],
                                                      nnodes_per_element);
                    }

                    if (verbose && p < 5) {
                      *outStream << " - Assembled element matrix and rhs -\n";
                      *outStream << std::showpos;
                      for (int i=0;i<ndofs;++i) {
                        for (int j=0;j<ndofs;++j)
                          *outStream << MatVal(A_asm, i, j) << " ";
                        *outStream << ":: " << MatVal(b_asm, i, 0) << "\n";
                      }
                      *outStream << std::noshowpos;
                    }

                    // Step 2: solve the system of equations
                    int info = 0;
                    Teuchos::LAPACK<int,value_type> lapack;
                    FieldContainer<int> ipiv(ndofs);
                    lapack.GESV(ndofs, 1, &A_asm(0,0,0), ndofs, &ipiv(0,0), &b_asm(0,0), ndofs, &info);
                    TEUCHOS_TEST_FOR_EXCEPTION( info != 0, std::runtime_error,
                                                ">>> ERROR (Intrepid::HGRAD_TRI_Cn::Test 04): " \
                                                "LAPACK solve fails");

                    // Step 3: construct interpolant and check solutions
                    magnitude_type interpolation_error = 0, solution_norm =0;
                    for (int iel=0;iel<nelement;++iel) {
                      retrieve_element_solution(b,
                                                b_asm,
                                                local2global[iel],
                                                nnodes_per_element);

                      if (verbose && p < 5) {
                        *outStream << " - Element solution, iel = " << iel << "\n";
                        *outStream << std::showpos;
                        for (int i=0;i<nbf;++i) {
                          *outStream << MatVal(b, i, 0) << "\n";
                        }
                        *outStream << std::noshowpos;
                      }

                      magnitude_type
                        element_interpolation_error = 0,
                        element_solution_norm = 0;

                      Orientation ort = Orientation::getOrientation(cell, element[iel]);
                      fill_cell_nodes(cell_nodes,
                                      nodes,
                                      element[iel],
                                      nvert,
                                      ndim);

                      compute_element_error(element_interpolation_error,
                                            element_solution_norm,
                                            element[iel],
                                            cell_nodes,
                                            basis_set,
                                            b,
                                            ort,
                                            nx, ny, nz);

                      interpolation_error += element_interpolation_error;
                      solution_norm       += element_solution_norm;

                      {
                        int edge_orts[6], face_orts[4];
                        ort.getEdgeOrientation(edge_orts, nedge);
                        ort.getFaceOrientation(face_orts, nface);
                        *outStream << "   iel = " << std::setw(4) << iel
                                   << ", edge orientation = "
                                   << edge_orts[0]
                                   << edge_orts[1]
                                   << edge_orts[2]
                                   << edge_orts[3]
                                   << edge_orts[4]
                                   << edge_orts[5]
                                   << ", face orientation = "
                                   << face_orts[0]
                                   << face_orts[1]
                                   << face_orts[2]
                                   << face_orts[3]
                                   << " , error = " << element_interpolation_error
                                   << " , solution norm = " << element_solution_norm
                                   << " , relative error = " << (element_interpolation_error/element_solution_norm)
                                   << "\n";
                      }
                      const magnitude_type relative_error = interpolation_error/solution_norm;
                      const magnitude_type tol = p*p*p*200*INTREPID_TOL;
                
                      if (relative_error > tol) {
                        ++r_val;
                        *outStream << "\n\nPatch test failed: \n"
                                   << "    exact polynomial (nx, ny) = " << std::setw(4) << nx << ", " << std::setw(4) << ny << ", " << std::setw(4) << nz << "\n"
                                   << "    basis order               = " << std::setw(4) << p << "\n"
                                   << "    orientation configuration = " << std::setw(4) << conf0 << std::setw(4) << conf1 << "\n"
                                   << "    relative error            = " << std::setw(4) << relative_error << "\n"
                                   << "    tolerance                 = " << std::setw(4) << tol << "\n";
                      }
                    }
                  }
                  // for next iteration                                                                         
                  std::next_permutation(elt_1, elt_1+4); 
                } // end of conf1
              } // end of conf0
            } // end of p
          } // end of point type
        } // end of nz
      } // end of ny
    } // end of nx
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    r_val = -1000;
  };
#else
  *outStream << "\t This test is for high order element assembly. \n"
             << "\t Use -D INTREPID_USING_EXPERIMENTAL_HIGH_ORDER in CMAKE_CXX_FLAGS \n";
#endif
  
  if (r_val != 0)
    std::cout << "End Result: TEST FAILED  :: r_val = " << r_val << "\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  Kokkos::finalize();
  
  return r_val;
}
