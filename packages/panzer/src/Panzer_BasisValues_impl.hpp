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

#ifndef PANZER_BASIS_VALUES_IMPL_HPP
#define PANZER_BASIS_VALUES_IMPL_HPP

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_Basis.hpp"

#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_Dimension.hpp"

// ***********************************************************
// * Specializations of setupArrays() for different array types
// ***********************************************************

namespace panzer {

  //! Used for default parameters
  template <typename Scalar,typename Array,typename ArrayOrientation>
  const Array BasisValues<Scalar,Array,ArrayOrientation>::dummyArray;    

  // ***********************************************************
  // * Evaluation of values without weithed measure - NOT specialized
  // ***********************************************************

  template<typename Scalar, typename Array,typename ArrayOrientation>
  inline
  void panzer::BasisValues<Scalar,Array,ArrayOrientation>::
  evaluateValues(const Array& cub_points,
                 const Array& jac,
                   const Array& jac_det,
                 const Array& jac_inv)
  { 
     // substitute dummy array for weighted measure
     evaluateValues(cub_points,jac,jac_det,jac_inv,dummyArray,dummyArray);
  }

  template<typename Scalar, typename Array,typename ArrayOrientation>
  inline
  void panzer::BasisValues<Scalar,Array,ArrayOrientation>::
  evaluateValues(const Array& cub_points,
                 const Array& jac,
                   const Array& jac_det,
                 const Array& jac_inv,
                 const Array& node_coordinates)
  { 
     // substitute dummy array for weighted measure
     evaluateValues(cub_points,jac,jac_det,jac_inv,dummyArray,node_coordinates);
  }
  
  // ***********************************************************
  // * Evaluation of values - NOT specialized
  // ***********************************************************
  template<typename Scalar, typename Array,typename ArrayOrientation>
  inline
  void panzer::BasisValues<Scalar,Array,ArrayOrientation>::
  evaluateValues(const Array& cub_points,
                 const Array& jac,
                   const Array& jac_det,
                 const Array& jac_inv,
                 const Array& weighted_measure,
                 const Array& node_coordinates)
  {    
    bool buildWeighted = ((&weighted_measure)!= &dummyArray);
    bool buildBasisPoints = ((&node_coordinates)!= &dummyArray);

    // first grab basis descriptor
    PureBasis::EElementSpace elmtspace = getElementSpace();
    int spaceDim = basis_layout->dimension();

    intrepid_basis->getValues(basis_ref, cub_points, 
                              Intrepid::OPERATOR_VALUE);

    if(elmtspace==PureBasis::CONST) {

       Intrepid::FunctionSpaceTools::
         HGRADtransformVALUE<Scalar>(basis,
                                     basis_ref);

       if(buildWeighted) {
          Intrepid::FunctionSpaceTools::
            multiplyMeasure<Scalar>(weighted_basis, 
                                       weighted_measure, 
                                    basis);
       }
    }
    else if(elmtspace==PureBasis::HGRAD) {

       Intrepid::FunctionSpaceTools::
         HGRADtransformVALUE<Scalar>(basis,
                                     basis_ref);
        
       if(compute_derivatives) {
         intrepid_basis->getValues(grad_basis_ref, cub_points, 
                                        Intrepid::OPERATOR_GRAD);
       
         Intrepid::FunctionSpaceTools::
           HGRADtransformGRAD<Scalar>(grad_basis,
                                      jac_inv,
                                        grad_basis_ref);
       }

       if(buildWeighted) {
          Intrepid::FunctionSpaceTools::
            multiplyMeasure<Scalar>(weighted_basis, 
                                       weighted_measure, 
                                    basis);
       
          if(compute_derivatives) {
             Intrepid::FunctionSpaceTools::
               multiplyMeasure<Scalar>(weighted_grad_basis, 
                                       weighted_measure, 
                                       grad_basis);
          }
       }
    }
    else if(elmtspace==PureBasis::HCURL) {
       Intrepid::FunctionSpaceTools::
         HCURLtransformVALUE<Scalar>(basis,
                                     jac_inv,
                                     basis_ref);

       if(compute_derivatives) {
          intrepid_basis->getValues(curl_basis_ref, cub_points, 
                                       Intrepid::OPERATOR_CURL);
   
          // the CURL transform differs depending on spatial dimension
          // "In 2D the de Rham complex collapses." But I look at it like this:
          // In 2D the curl is simply  $-\partial_x u_1+\partial_y u_0$ which
          // has the same derivative structure as the divergence in 2D 
          // $\partial_x u_0 + \partial_y u_1$ which means the transformation
          // is the same.
          if(spaceDim==2) 
             Intrepid::FunctionSpaceTools::
               HDIVtransformDIV<Scalar>(curl_basis,
                                           jac_det,   // note only volume deformation is needed!
                                                   // this relates directly to this being in
                                                   // the divergence space in 2D!
                                           curl_basis_ref);
          else if(spaceDim==3)
             Intrepid::FunctionSpaceTools::
               HCURLtransformCURL<Scalar>(curl_basis,
                                             jac,
                                             jac_det,
                                             curl_basis_ref);
          else
             TEUCHOS_ASSERT(false); // what you doin?
       }

       if(buildWeighted) {
          Intrepid::FunctionSpaceTools::
            multiplyMeasure<Scalar>(weighted_basis, 
                                    weighted_measure, 
                                    basis);
   
          if(compute_derivatives) {
             Intrepid::FunctionSpaceTools::
               multiplyMeasure<Scalar>(weighted_curl_basis, 
                                       weighted_measure, 
                                       curl_basis);
          }
       }
    }
    else if(elmtspace==PureBasis::HDIV) {
       Intrepid::FunctionSpaceTools::
         HDIVtransformVALUE<Scalar>(basis,
                                    jac,
                                    jac_det,
                                    basis_ref);

       if(compute_derivatives) {
          intrepid_basis->getValues(div_basis_ref, cub_points, 
                                       Intrepid::OPERATOR_DIV);
   
          Intrepid::FunctionSpaceTools::
             HDIVtransformDIV<Scalar>(div_basis,
                                      jac_det,   
                                      div_basis_ref);
       }

       if(buildWeighted) {
          Intrepid::FunctionSpaceTools::
            multiplyMeasure<Scalar>(weighted_basis, 
                                    weighted_measure, 
                                    basis);
   
          if(compute_derivatives) {
             Intrepid::FunctionSpaceTools::
               multiplyMeasure<Scalar>(weighted_div_basis, 
                                       weighted_measure, 
                                       div_basis);
          }
       }
    }
    

    // If basis supports coordinate values at basis points, then
    // compute these values
    if(buildBasisPoints) {
      Teuchos::RCP<Intrepid::DofCoordsInterface<Array> > coords;
      coords = Teuchos::rcp_dynamic_cast<Intrepid::DofCoordsInterface<Array> >(intrepid_basis);
      if (!Teuchos::is_null(coords)) {
        coords->getDofCoords(basis_coordinates_ref);
        Intrepid::CellTools<Scalar> cell_tools;
        cell_tools.mapToPhysicalFrame(basis_coordinates, 
                                      basis_coordinates_ref,
                                      node_coordinates,
                                      intrepid_basis->getBaseCellTopology());
      }
    }

  }

  template<typename Scalar, typename Array,typename ArrayOrientation>
  PureBasis::EElementSpace BasisValues<Scalar,Array,ArrayOrientation>::getElementSpace() const
  { return basis_layout->getBasis()->getElementSpace(); }

  // **************************************************************

  // method for applying orientations
  template<typename Scalar,typename Array,typename ArrayOrientation>
  void BasisValues<Scalar,Array,ArrayOrientation>::
  applyOrientations(const ArrayOrientation& orientations)
  {
    PureBasis::EElementSpace elmtspace = getElementSpace();

    if(elmtspace==PureBasis::HCURL) {
       bool buildWeighted = (weighted_basis.size()!=0);

       // setup the orientations for the trial space
       Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(basis,orientations);

       if(compute_derivatives)
          Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(curl_basis,orientations);

       // setup the orientations for the test space
       if(buildWeighted) {
         Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(weighted_basis,orientations);

         if(compute_derivatives)
            Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(weighted_curl_basis,orientations);
       }
    }
    else if(elmtspace==PureBasis::HDIV) {
       bool buildWeighted = (weighted_basis.size()!=0);

       // setup the orientations for the trial space
       Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(basis,orientations);

       if(compute_derivatives)
          Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(div_basis,orientations);

       // setup the orientations for the test space
       if(buildWeighted) {
         Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(weighted_basis,orientations);

         if(compute_derivatives)
            Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(weighted_div_basis,orientations);
       }
    }
  }

} // end namespace panzer

#endif
