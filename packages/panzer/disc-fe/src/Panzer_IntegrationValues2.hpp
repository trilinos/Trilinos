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

  class SubcellConnectivity;

  /**
   * \brief Swap the ordering of quadrature points in a specified cell.
   *
   * \param[in] cell   Cell index
   * \param[in] a      Quadrature point a
   * \param[in] b      Quadrature point b
   *
   * NOTE: this should be a member function of IntegrationValues but
   * lambda capture generates a ton of cuda compiler warnings. Making
   * it a stand alone function fixes this.
   */
  template<typename Scalar,
           typename T0,typename T1,typename T2,typename T3,
           typename T4,typename T5,typename T6,typename T7>
  KOKKOS_INLINE_FUNCTION
  void swapQuadraturePoints(int cell,int a,int b,
                            T0& ref_ip_coordinates,
                            T1& ip_coordinates,
                            T2& weighted_measure,
                            T3& jac,
                            T4& jac_det,
                            T5& jac_inv,
                            T6& surface_normals,
                            T7& surface_rotation_matrices);

  template <typename Scalar>
  class IntegrationValues2 {
  public:
    typedef typename ArrayTraits<Scalar,PHX::MDField<Scalar> >::size_type size_type;

    typedef PHX::MDField<Scalar> ArrayDynamic;
    typedef PHX::MDField<double> DblArrayDynamic;

    typedef PHX::MDField<Scalar,IP> Array_IP;
    typedef PHX::MDField<Scalar,IP,Dim> Array_IPDim;

    typedef PHX::MDField<Scalar,Point> Array_Point;
    typedef PHX::MDField<Scalar,Cell,IP> Array_CellIP;
    typedef PHX::MDField<Scalar,Cell,IP,Dim> Array_CellIPDim;
    typedef PHX::MDField<Scalar,Cell,IP,Dim,Dim> Array_CellIPDimDim;

    typedef PHX::MDField<Scalar,Cell,BASIS,Dim> Array_CellBASISDim;

    IntegrationValues2(const std::string & pre="",bool allocArrays=false)
        : alloc_arrays(allocArrays), prefix(pre), ddims_(1,0) {}

    //! Sizes/allocates memory for arrays
    void setupArrays(const Teuchos::RCP<const panzer::IntegrationRule>& ir);

    void setupArraysForNodeRule(const Teuchos::RCP<const panzer::IntegrationRule>& ir);

    /** \brief Evaluate basis values.

        @param vertex_coordinates [in] Cell vertex coordinates, not
        basis coordinates.
        @param num_cells [in] (optional) number of cells in the
        workset. This can be less than the workset size. If set to
        zero, extent(0) of the evaluated array is used which equates
        to the workset size.
        @param face_connectivity [in] (optional) connectivity used to
        enforce quadrature alignment for surface integration.
     */
    void evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
                        const int num_cells = -1,
                        const Teuchos::RCP<const SubcellConnectivity> & face_connectivity = Teuchos::null);

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
    void evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
                        const PHX::MDField<Scalar,Cell,IP,Dim> & other_ip_coordinates,
                        const int num_cells = -1);

    Array_IPDim cub_points;              // <IP,Dim>
    Array_IPDim side_cub_points;         // <IP,Dim> points on face topology (dim-1)
    Array_IP cub_weights;                // <IP>
    Array_CellBASISDim node_coordinates; // <Cell,BASIS,Dim>
    Array_CellIPDimDim jac;              // <Cell,IP,Dim,Dim>
    Array_CellIPDimDim jac_inv;          // <Cell,IP,Dim,Dim>
    Array_CellIP jac_det;                // <Cell,IP>
    Array_CellIP weighted_measure;       // <Cell,IP>
    Array_CellIPDim weighted_normals;    // <Cell,IP,Dim>

    Array_CellIPDim surface_normals;    // <Cell,IP,Dim>
    Array_CellIPDimDim surface_rotation_matrices;    // <Cell,IP,Dim,Dim>
      // this (appears) is a matrix where the first row is the "normal" direction
      // and the remaining two rows lie in the hyperplane

    Teuchos::RCP<const panzer::IntegrationRule> int_rule;

    Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double>> intrepid_cubature;

    // for Shakib stabilization <Cell,IP,Dim,Dim>
    Array_CellIPDimDim covarient;
    Array_CellIPDimDim contravarient;
    Array_CellIP norm_contravarient;

    // integration points
    Array_CellIPDim ip_coordinates;      // <Cell,IP,Dim>
    Array_CellIPDim ref_ip_coordinates;  // <Cell,IP,Dim> for Control Volumes or Surface integrals

    DblArrayDynamic dyn_cub_points, dyn_side_cub_points, dyn_cub_weights;
    DblArrayDynamic dyn_phys_cub_points, dyn_phys_cub_weights, dyn_phys_cub_norms, dyn_node_coordinates;

    Array_Point scratch_for_compute_side_measure; // <Point> size: span() == jac.span()

    /// This should be a private method, but using lambdas on cuda forces this to be public.
    void evaluateRemainingValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & in_node_coordinates, const int in_num_cells);

    /// This should be a private method, but using lambdas on cuda forces this to be public.
    void getCubatureCV(const PHX::MDField<Scalar,Cell,NODE,Dim> & in_node_coordinates, const int in_num_cells);

    /// This should be a private method, but using lambdas on cuda forces this to be public.
    void generateSurfaceCubatureValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & in_node_coordinates, const int in_num_cells,const SubcellConnectivity & face_connectivity);

  protected:

    // TODO: Make this a utility function that only exists in source file
    Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double>> getIntrepidCubature(const panzer::IntegrationRule & ir) const;

  private:
    bool alloc_arrays;
    std::string prefix;
    std::vector<PHX::index_size_type> ddims_;

    void getCubature(const PHX::MDField<Scalar,Cell,NODE,Dim> & in_node_coordinates, const int in_num_cells);
    void evaluateValuesCV(const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,const int in_num_cells);
  };

} // namespace panzer

#endif
