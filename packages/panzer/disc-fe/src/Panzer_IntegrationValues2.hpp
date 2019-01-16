// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER


#ifndef PANZER_INTEGRATION_VALUES2_HPP
#define PANZER_INTEGRATION_VALUES2_HPP

#include "Teuchos_RCP.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_ArrayTraits.hpp"
#include "Panzer_Dimension.hpp"
#include "Phalanx_MDField.hpp"
#include "Intrepid2_Cubature.hpp"

namespace panzer {

  template <typename Scalar>
  class IntegrationValues2 {
  public:
    typedef typename ArrayTraits<Scalar,PHX::MDField<Scalar> >::size_type size_type;

    typedef PHX::MDField<Scalar,IP>               Array_IP;
    typedef PHX::MDField<Scalar,IP,Dim>           Array_IPDim;
    typedef PHX::MDField<Scalar,Point>            Array_Point;
    typedef PHX::MDField<Scalar,Cell,IP>          Array_CellIP;
    typedef PHX::MDField<Scalar,Cell,IP,Dim>      Array_CellIPDim;
    typedef PHX::MDField<Scalar,Cell,IP,Dim,Dim>  Array_CellIPDimDim;
    typedef PHX::MDField<Scalar,Cell,BASIS,Dim>   Array_CellBASISDim;

    /**
     * \brief Base constructor
     *
     * \param[in] pre Prefix to apply to all generated arrays
     * \param[in] allocArrays Should the arrays be allocated?
     *
     */
    IntegrationValues2(const std::string & pre="",bool allocArrays=false)
        : alloc_arrays(allocArrays), prefix(pre), ddims_(1,0) {}

    //! Sizes/allocates memory for arrays
    // Note: setupArrays will call setupArraysForNodeRule if rule is 1D side integrator
    void setupArrays(const Teuchos::RCP<const panzer::IntegrationRule>& ir);
    void setupArraysForNodeRule(const Teuchos::RCP<const panzer::IntegrationRule>& ir);

    /** \brief Evaluate basis values.
     *
     * Must be called after setupArrays
     *
     * \param[in] vertex_coordinates Cell vertices (not basis coordinates)
     * \param[in] num_cells Subset of cells to apply to (starting at 0) defaults to vertex_coordinates.extent(0)
     *
     */
    void evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
                        const int num_cells = -1);

    /**
     * \brief Make the ordering of the integration points unique
     *
     * This is required when we need side/surface integration points to line up
     *
     * Must be called after getCubature and evaluateDefaultValues
     *
     * \param[in] num_cells Subset of cells to apply to (starting at 0) defaults to vertex_coordinates.extent(0)
     */
    void makeOrderingUnique(const int num_cells = -1);

    /** \brief Match IP.

       Optionally provide IP coordinates for an element 'other' that
       shares the same side. If provided, a permutation of the
       cubature points is calculated so that the integration values
       are ordered according to the other element's. This permutation
       is then applied so that all fields are ordered accordingly in
       their IP dimension.

        @param num_cells [in] (optional) number of cells in the
        workset. This can be less than the workset size. If set to
        zero, extent(0) of the evaluated array is used which equates
        to the workset size.
    */
//    void evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
//                        const PHX::MDField<Scalar,Cell,IP,Dim> & other_ip_coordinates,
//                        const int num_cells = -1);

    // Universal cubature points/weights arrays - never valid
    // FIXME: These should never be used (use ref_ip_coordinates instead). Some tests use these objects so they remains public.
    Array_IPDim cub_points;
    Array_IP cub_weights;

    // Node coordinates defining mesh
    Array_CellBASISDim node_coordinates;

    // Jacobian values at integration points - valid for VOLUME, SIDE, SURFACE, CV
    Array_CellIPDimDim jac;

    // Inverse Jacobian values at integration points - valid for VOLUME, SIDE, SURFACE, CV
    Array_CellIPDimDim jac_inv;

    // Jacobian determinant values at integration points - valid for VOLUME, SIDE, SURFACE, CV
    Array_CellIP jac_det;

    // Weighted measure (cubature weights * jacobian) - valid for VOLUME, SIDE, SURFACE, CV
    Array_CellIP weighted_measure;

    // Weighted normals (?)
    Array_CellIPDim weighted_normals;    // <Cell,IP,Dim>

    // Surface normals at cubature points - only valid for SURFACE integration
    Array_CellIPDim surface_normals;

    // Surface rotation matrices at cubature points - only valid for SURFACE integration
    // Note: This an orthogonal matrix where the first row is the "normal" direction and the remaining two rows lie in the hyperplane
    // it is used to rotate into and out of the face normal reference frame
    Array_CellIPDimDim surface_rotation_matrices;

    // Integration rule defining values
    Teuchos::RCP<const panzer::IntegrationRule> int_rule;

    // Cubature used to derive values
    Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double>> intrepid_cubature;

    // for Shakib stabilization - only valid for VOLUME, SIDE, SURFACE
    Array_CellIPDimDim covarient;
    Array_CellIPDimDim contravarient;
    Array_CellIP norm_contravarient;

    // Integration points (mesh space) - valid for VOLUME, SIDE, SURFACE, CV
    Array_CellIPDim ip_coordinates;

    // Integration points (reference space) - only valid for SURFACE, CV
    Array_CellIPDim ref_ip_coordinates;

    /**
     * \brief Define a unique coordinate ordering for a given cell in coords
     *
     * \note This call is public only for testing purposes
     *
     * Used for side/surface integration points. Compute a unique ordering in a cell and
     * point offset.
     *
     * \param[in] coords Coordinates array (cell,IP,Dim)
     * \param[in] cell   Cell index
     * \param[in] offset Offset into the points
     * \param[out] order Ordering array on output, correctly sized on input
     *                   (offset + order.size() <= coords.extent(1))
     */
    static void uniqueCoordOrdering(Array_CellIPDim & coords,
                                    int cell,
                                    int offset,
                                    std::vector<int> & order);

    /**
     * \brief Swap the ordering of quadrature points in a specified cell.
     *
     * \note This call is public only for testing purposes - it should only be used internally
     *
     * \param[in] cell   Cell index
     * \param[in] a      Quadrature point a
     * \param[in] b      Quadrature point b
     */
    void swapQuadraturePoints(int cell,int a,int b);


  private:

    typedef PHX::MDField<Scalar> ArrayDynamic;
    typedef PHX::MDField<double> DblArrayDynamic;

    // Scratch space - we use 'dyn_' to state that these are dynamic rank views - required by some Intrepid2 calls
    DblArrayDynamic dyn_cub_points, dyn_side_cub_points, dyn_cub_weights;
    DblArrayDynamic dyn_phys_cub_points, dyn_phys_cub_weights, dyn_phys_cub_norms, dyn_node_coordinates;
    Array_Point scratch_for_compute_side_measure; // <Point> size: span() == jac.span()

    // TODO: Make this a utility function that only exists in source file
    Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double>> getIntrepidCubature(const panzer::IntegrationRule & ir) const;

    bool alloc_arrays;
    std::string prefix;
    std::vector<PHX::index_size_type> ddims_;

    /**
     * \brief Setup integration for VOLUME, SIDE
     *
     * Must be called after setupArrays
     *
     * Fills:
     *   ip_coordinates
     *   ref_ip_coordinates
     *   jac
     *   jac_inv
     *   jac_det
     *
     *
     * This call also fills
     *   dyn_cub_points
     *   dyn_cub_weights
     *   dyn_side_cub_points
     *
     * \note node_coordinates must be set!!
     *
     * \param[in] num_cells Subset of cells to apply to (starting at 0) defaults to vertex_coordinates.extent(0)
     */
    void getCubature(const int in_num_cells);

    /**
     * \brief Setup integration for CV
     *
     * Must be called after setupArrays
     *
     * Fills:
     *   ip_coordinates
     *   ref_ip_coordinates
     *   weighted_measure
     *   weighted_normals
     *   node_coordinates
     *
     * This call also fills
     *   dyn_phys_cub_points
     *   dyn_phys_cub_norms
     *   dyn_phys_cub_weights
     *   dyn_node_coordinates
     *
     * \note node_coordinates must be set!!
     *
     * \param[in] num_cells Subset of cells to apply to (starting at 0) defaults to vertex_coordinates.extent(0)
     */
    void getCubatureCV(const int in_num_cells);

    /**
     * \brief Setup integration for SURFACE
     *
     * Must be called after setupArrays
     *
     * Fills:
     *   ip_coordinates
     *   ref_ip_coordinates
     *   weighted_measure
     *   surface_normals
     *   surface_rotation_matrices
     *   jac
     *   jac_inv
     *   jac_det
     *
     * \param[in] num_cells Subset of cells to apply to (starting at 0) defaults to vertex_coordinates.extent(0)
     */
    void getCubatureSurface(const int in_num_cells);

    /**
     * \brief Evaluate default values for VOLUME and SIDE integration
     *
     * Must be called after getCubature<...>
     *
     * Fills:
     *   covariant
     *   contravarient
     *   norm_contravarient
     *
     * \param[in] num_cells Subset of cells to apply to (starting at 0) defaults to vertex_coordinates.extent(0)
     */
    void evaluateDefaultValues(const int in_num_cells);

  };

} // namespace panzer

#endif
