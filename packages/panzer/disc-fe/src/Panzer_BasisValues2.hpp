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

#ifndef PANZER_BASIS_VALUES2_HPP
#define PANZER_BASIS_VALUES2_HPP

#include "Teuchos_RCP.hpp"

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_Orientation.hpp"

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_ArrayTraits.hpp"

namespace panzer {

  /** Data structure that holds all evaluated fields associated
    * with a basis fucntion and integration rule. This class will
    * allocate the memory and evaluate the basis functions. The 
    * orientations must be applied using the  
    * <code>applyOrientations</code> method.
    */
  template <typename Scalar>
  class BasisValues2 { 
  public:
    typedef typename ArrayTraits<Scalar,PHX::MDField<Scalar,void,void,void,void,void,void,void,void> >::size_type size_type;

    typedef PHX::MDField<Scalar> ArrayDynamic;
    typedef PHX::MDField<Scalar,BASIS,IP,void,void,void,void,void,void> Array_BasisIP;
    typedef PHX::MDField<Scalar,Cell,BASIS,IP,void,void,void,void,void> Array_CellBasisIP;
    typedef PHX::MDField<Scalar,BASIS,IP,Dim,void,void,void,void,void> Array_BasisIPDim;
    typedef PHX::MDField<Scalar,Cell,BASIS,IP,Dim,void,void,void,void> Array_CellBasisIPDim;
    typedef PHX::MDField<Scalar,BASIS,Dim,void,void,void,void,void,void> Array_BasisDim;
    typedef PHX::MDField<Scalar,Cell,BASIS,Dim,void,void,void,void,void> Array_CellBasisDim;

    BasisValues2(const std::string & pre="",bool allocArrays=false,bool buildWeighted=false) 
        : build_weighted(buildWeighted), alloc_arrays(allocArrays), prefix(pre)
        , ddims_(1,0), references_evaluated(false) {}
 
    //! Sizes/allocates memory for arrays
    void setupArrays(const Teuchos::RCP<const panzer::BasisIRLayout>& basis,
                     bool computeDerivatives=true);

    void evaluateValues(const PHX::MDField<Scalar,IP,Dim,void,void,void,void,void,void> & cub_points,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
                        const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                        const int in_num_cells = -1);

    void evaluateValues(const PHX::MDField<Scalar,IP,Dim,void,void,void,void,void,void> & cub_points,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
                        const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                        const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                        const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
                        bool use_vertex_coordinates=true,
                        const int in_num_cells = -1);

    void evaluateValuesCV(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cell_cub_points,
                          const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
                          const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
                          const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv);

    void evaluateValues(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
                        const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
                        const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                        const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                        const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
                        bool use_vertex_coordinates=true,
                        const int in_num_cells = -1);


    //! Method to apply orientations to a basis values container.
    // some evaluators use this apply orientation (this will be deprecated)
    void applyOrientations(const PHX::MDField<const Scalar,Cell,BASIS> & orientations);

    // this is used in workset factory
    void applyOrientations(const std::vector<Intrepid2::Orientation> & orientations,
                           const int in_num_cells = -1);

    void setExtendedDimensions(const std::vector<PHX::index_size_type> & ddims)
    { ddims_ = ddims; }

    PureBasis::EElementSpace getElementSpace() const;

    Array_BasisIP     basis_ref_scalar;           // <BASIS,IP> 
    Array_CellBasisIP basis_scalar;               // <Cell,BASIS,IP> 

    Array_BasisIPDim     basis_ref_vector;           // <BASIS,IP,Dim> 
    Array_CellBasisIPDim basis_vector;               // <Cell,BASIS,IP,Dim>

    Array_BasisIPDim     grad_basis_ref;             // <BASIS,IP,Dim>
    Array_CellBasisIPDim grad_basis;                 // <Cell,BASIS,IP,Dim>

    Array_BasisIP     curl_basis_ref_scalar;         // <BASIS,IP> - 2D!
    Array_CellBasisIP curl_basis_scalar;             // <Cell,BASIS,IP> - 2D!

    Array_BasisIPDim     curl_basis_ref_vector;      // <BASIS,IP,Dim> - 3D!
    Array_CellBasisIPDim curl_basis_vector;          // <Cell,BASIS,IP,Dim> - 3D!

    Array_BasisIP     div_basis_ref;           // <BASIS,IP> 
    Array_CellBasisIP div_basis;               // <Cell,BASIS,IP> 

    Array_CellBasisIP weighted_basis_scalar;                  // <Cell,BASIS,IP> 
    Array_CellBasisIPDim weighted_basis_vector;               // <Cell,BASIS,IP,Dim> 
    Array_CellBasisIPDim weighted_grad_basis;                 // <Cell,BASIS,IP,Dim>
    Array_CellBasisIP weighted_curl_basis_scalar;             // <Cell,BASIS,IP>
    Array_CellBasisIPDim weighted_curl_basis_vector;          // <Cell,BASIS,IP,Dim>
    Array_CellBasisIP weighted_div_basis;                     // <Cell,BASIS,IP> 

    /** Carterisan coordinates for basis coefficients

        NOTE: This quantity is not always available.  Certain bases
        may not have a corresponding coordiante value
    */
    Array_BasisDim basis_coordinates_ref;     // <BASIS,Dim>

    /** Carterisan coordinates for basis coefficients

        NOTE: This quantity is not always available.  Certain bases
        may not have a corresponding coordiante value
    */
    Array_CellBasisDim basis_coordinates;     // <Cell,BASIS,Dim>

    Teuchos::RCP<const panzer::BasisIRLayout> basis_layout;
    
    Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space,Scalar,Scalar> > intrepid_basis;

    bool compute_derivatives;
    bool build_weighted;
    bool alloc_arrays;
    std::string prefix;
    std::vector<PHX::index_size_type> ddims_;

  protected:


    void evaluateValues_Const(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
                              const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                              const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                              const int in_num_cells);

    void evaluateValues_HVol(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
                             const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
                             const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                             const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                             const int in_num_cells);

    void evaluateValues_HGrad(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
                              const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                              const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                              const int in_num_cells);

    void evaluateValues_HCurl(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
                              const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
                              const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
                              const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                              const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                              const int in_num_cells);

    void evaluateValues_HDiv(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
                             const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
                             const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
                             const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                             const int in_num_cells);
 
  private:

    void evaluateBasisCoordinates(const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
                                  const int in_num_cells = -1);

    /** Evaluate the reference values for the basis functions needed
      *
      * Node after this call references_evaluated=true 
      */
    void evaluateReferenceValues(const PHX::MDField<Scalar,IP,Dim> & cub_points,bool compute_derivatives,bool use_vertex_coordinates);

    bool references_evaluated;
  };

} // namespace panzer

#endif
