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

#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_BasisValues2.hpp"

#include "Panzer_CommonArrayFactories.hpp"

#include "Intrepid_FunctionSpaceTools.hpp"


namespace panzer {

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,IP,Dim,void,void,void,void,void,void> & cub_points,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
               const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv)
{
  PHX::MDField<Scalar,Cell,IP> weighted_measure;
  PHX::MDField<Scalar,Cell,NODE,Dim> vertex_coordinates;
  build_weighted = false; 
  evaluateValues(cub_points,jac,jac_det,jac_inv,weighted_measure,vertex_coordinates,false);
}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,IP,Dim,void,void,void,void,void,void> & cub_points,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
               const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
               const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
               const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
               bool use_vertex_coordinates)
{
  MDFieldArrayFactory af("",ddims_,true);
 
  int num_dim   = basis_layout->dimension();

  // currently this just copies on the basis objects are converted
  // in intrepid
  evaluateReferenceValues(cub_points,compute_derivatives,use_vertex_coordinates);

  PureBasis::EElementSpace elmtspace = getElementSpace();
  if(elmtspace==PureBasis::CONST ||
     elmtspace==PureBasis::HGRAD) {
    Intrepid::FunctionSpaceTools::
      HGRADtransformVALUE<Scalar>(basis_scalar,
                                  basis_ref_scalar);

    if(build_weighted) {
      Intrepid::FunctionSpaceTools::
        multiplyMeasure<Scalar>(weighted_basis_scalar, 
                                    weighted_measure, 
                                    basis_scalar);
    }
  }
  else if(elmtspace==PureBasis::HCURL) {
    Intrepid::FunctionSpaceTools::
      HCURLtransformVALUE<Scalar>(basis_vector,
                                     jac_inv,
                                     basis_ref_vector);

    if(build_weighted) {
      Intrepid::FunctionSpaceTools::
        multiplyMeasure<Scalar>(weighted_basis_vector, 
                                    weighted_measure, 
                                    basis_vector);
    }
  }
  else if(elmtspace==PureBasis::HDIV)
  {
    Intrepid::FunctionSpaceTools::
      HDIVtransformVALUE<Scalar>(basis_vector,
                                      jac,
                                      jac_det,
                                      basis_ref_vector);

    if(build_weighted) {
      Intrepid::FunctionSpaceTools::
        multiplyMeasure<Scalar>(weighted_basis_vector, 
                                    weighted_measure, 
                                    basis_vector);
    }
  }
  else { TEUCHOS_ASSERT(false); }

  if(elmtspace==PureBasis::HGRAD && compute_derivatives) {
    Intrepid::FunctionSpaceTools::
      HGRADtransformGRAD<Scalar>(grad_basis,
                                     jac_inv,
                                     grad_basis_ref);

    if(build_weighted) {
      Intrepid::FunctionSpaceTools::
                 multiplyMeasure<Scalar>(weighted_grad_basis, 
                                             weighted_measure, 
                                             grad_basis);
    }
  }
  else if(elmtspace==PureBasis::HCURL && num_dim==2 && compute_derivatives) {
    Intrepid::FunctionSpaceTools::
      HDIVtransformDIV<Scalar>(curl_basis_scalar,
                               jac_det,   // note only volume deformation is needed!
                                          // this relates directly to this being in
                                          // the divergence space in 2D!
                               curl_basis_ref_scalar);

    if(build_weighted) {
      Intrepid::FunctionSpaceTools::
                 multiplyMeasure<Scalar>(weighted_curl_basis_scalar, 
                                         weighted_measure, 
                                         curl_basis_scalar);
    }
  }
  else if(elmtspace==PureBasis::HCURL && num_dim==3 && compute_derivatives) {
    Intrepid::FunctionSpaceTools::
      HCURLtransformCURL<Scalar>(curl_basis_vector,
                                     jac,
                                     jac_det,
                                     curl_basis_ref_vector);

    if(build_weighted) {
      Intrepid::FunctionSpaceTools::
                 multiplyMeasure<Scalar>(weighted_curl_basis_vector, 
                                             weighted_measure, 
                                             curl_basis_vector);
    }
  }
  else if(elmtspace==PureBasis::HDIV && compute_derivatives) {
    Intrepid::FunctionSpaceTools::
      HDIVtransformDIV<Scalar>(div_basis,
                                   jac_det,
                                   div_basis_ref);

    if(build_weighted) {
      Intrepid::FunctionSpaceTools::
                 multiplyMeasure<Scalar>(weighted_div_basis, 
                                             weighted_measure, 
                                             div_basis);
    }
  }

  // If basis supports coordinate values at basis points, then
  // compute these values
  if(use_vertex_coordinates) {
    Teuchos::RCP<Intrepid::DofCoordsInterface<ArrayDynamic> > coords
        = Teuchos::rcp_dynamic_cast<Intrepid::DofCoordsInterface<ArrayDynamic> >(intrepid_basis);
    if (!Teuchos::is_null(coords)) {
/*
      ArrayDynamic dyn_basis_coordinates_ref = af.buildArray<Scalar,BASIS,Dim>("basis_coordinates_ref",basis_coordinates_ref.dimension(0),basis_coordinates_ref.dimension(1));
      coords->getDofCoords(dyn_basis_coordinates_ref);

      // fill in basis coordinates
      for (size_type i = 0; i < basis_coordinates_ref.dimension(0); ++i)
        for (size_type j = 0; j < basis_coordinates_ref.dimension(1); ++j)
           basis_coordinates_ref(i,j) = dyn_basis_coordinates_ref(i,j); 
*/

      Intrepid::CellTools<Scalar> cell_tools;
      cell_tools.mapToPhysicalFrame(basis_coordinates, 
                                    basis_coordinates_ref,
                                    vertex_coordinates,
                                    intrepid_basis->getBaseCellTopology());
    }
  }
}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateReferenceValues(const PHX::MDField<Scalar,IP,Dim> & cub_points,bool compute_derivatives,bool use_vertex_coordinates)
{
  MDFieldArrayFactory af("",ddims_,true);

  int num_quad    = basis_layout->numPoints();
  int num_dim   = basis_layout->dimension();
  int num_card  = basis_layout->cardinality();

  ArrayDynamic dyn_cub_points = af.buildArray<Scalar,IP,Dim>("dyn_cub_points",  num_quad,num_dim);

  for (size_type ip = 0; ip < num_quad; ++ip)
    for (size_type d = 0; d < num_dim; ++d)
      dyn_cub_points(ip,d) = cub_points(ip,d);

  PureBasis::EElementSpace elmtspace = getElementSpace();
  if(elmtspace==PureBasis::HGRAD || elmtspace==PureBasis::CONST) {
    ArrayDynamic dyn_basis_ref_scalar = af.buildArray<Scalar,BASIS,IP>("dyn_basis_ref_scalar",num_card,num_quad);

    intrepid_basis->getValues(dyn_basis_ref_scalar, dyn_cub_points, 
                              Intrepid::OPERATOR_VALUE);

    for (size_type b = 0; b < num_card; ++b)
      for (size_type ip = 0; ip < num_quad; ++ip) 
        basis_ref_scalar(b,ip) = dyn_basis_ref_scalar(b,ip);
  }
  else if(elmtspace==PureBasis::HDIV || elmtspace==PureBasis::HCURL) {
    ArrayDynamic dyn_basis_ref_vector = af.buildArray<Scalar,BASIS,IP,Dim>("dyn_basis_ref_vector",num_card,num_quad,num_dim);

    intrepid_basis->getValues(dyn_basis_ref_vector, dyn_cub_points, 
                              Intrepid::OPERATOR_VALUE);

    for (size_type b = 0; b < num_card; ++b)
      for (size_type ip = 0; ip < num_quad; ++ip) 
        for (size_type d = 0; d < num_dim; ++d) 
           basis_ref_vector(b,ip,d) = dyn_basis_ref_vector(b,ip,d);
  }
  else { TEUCHOS_ASSERT(false); }

  if(elmtspace==PureBasis::HGRAD && compute_derivatives) {
    ArrayDynamic dyn_grad_basis_ref = af.buildArray<Scalar,BASIS,IP,Dim>("dyn_basis_ref_vector",num_card,num_quad,num_dim);

    intrepid_basis->getValues(dyn_grad_basis_ref, dyn_cub_points, 
                              Intrepid::OPERATOR_GRAD);

    for (size_type b = 0; b < num_card; ++b)
      for (size_type ip = 0; ip < num_quad; ++ip) 
        for (size_type d = 0; d < num_dim; ++d) 
           grad_basis_ref(b,ip,d) = dyn_grad_basis_ref(b,ip,d);
  }
  else if(elmtspace==PureBasis::HCURL && compute_derivatives && num_dim==2) {
    ArrayDynamic dyn_curl_basis_ref = af.buildArray<Scalar,BASIS,IP>("dyn_curl_basis_ref_scalar",num_card,num_quad);

    intrepid_basis->getValues(dyn_curl_basis_ref, dyn_cub_points, 
                              Intrepid::OPERATOR_CURL);

    for (size_type b = 0; b < num_card; ++b)
      for (size_type ip = 0; ip < num_quad; ++ip) 
        curl_basis_ref_scalar(b,ip) = dyn_curl_basis_ref(b,ip);
  }
  else if(elmtspace==PureBasis::HCURL && compute_derivatives && num_dim==3) {
    ArrayDynamic dyn_curl_basis_ref = af.buildArray<Scalar,BASIS,IP,Dim>("dyn_curl_basis_ref_vector",num_card,num_quad,num_dim);

    intrepid_basis->getValues(dyn_curl_basis_ref, dyn_cub_points, 
                              Intrepid::OPERATOR_CURL);

    for (size_type b = 0; b < num_card; ++b)
      for (size_type ip = 0; ip < num_quad; ++ip) 
        for (size_type d = 0; d < num_dim; ++d) 
           curl_basis_ref_vector(b,ip,d) = dyn_curl_basis_ref(b,ip,d);
  }
  else if(elmtspace==PureBasis::HDIV && compute_derivatives) {
    ArrayDynamic dyn_div_basis_ref = af.buildArray<Scalar,BASIS,IP>("dyn_div_basis_ref_scalar",num_card,num_quad);

    intrepid_basis->getValues(dyn_div_basis_ref, dyn_cub_points, 
                              Intrepid::OPERATOR_DIV);

    for (size_type b = 0; b < num_card; ++b)
      for (size_type ip = 0; ip < num_quad; ++ip) 
        div_basis_ref(b,ip) = dyn_div_basis_ref(b,ip);
  }


  if(use_vertex_coordinates) {
    Teuchos::RCP<Intrepid::DofCoordsInterface<ArrayDynamic> > coords
        = Teuchos::rcp_dynamic_cast<Intrepid::DofCoordsInterface<ArrayDynamic> >(intrepid_basis);
    if (!Teuchos::is_null(coords)) {
      ArrayDynamic dyn_basis_coordinates_ref = af.buildArray<Scalar,BASIS,Dim>("basis_coordinates_ref",basis_coordinates_ref.dimension(0),basis_coordinates_ref.dimension(1));
      coords->getDofCoords(dyn_basis_coordinates_ref);

      // fill in basis coordinates
      for (size_type i = 0; i < basis_coordinates_ref.dimension(0); ++i)
        for (size_type j = 0; j < basis_coordinates_ref.dimension(1); ++j)
          basis_coordinates_ref(i,j) = dyn_basis_coordinates_ref(i,j); 
    }
  }

  references_evaluated = true;
}

// method for applying orientations
template <typename Scalar>
void BasisValues2<Scalar>::
applyOrientations(const PHX::MDField<Scalar,Cell,BASIS> & orientations)
{
  int num_cell  = orientations.dimension(0);
  int num_basis = orientations.dimension(1);
  int num_dim   = basis_layout->dimension();
  int num_ip    = basis_layout->numPoints();
  PureBasis::EElementSpace elmtspace = getElementSpace();

  if(elmtspace==PureBasis::HCURL && num_dim==2) {

    // setup the orientations for the trial space
    // Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(basis_vector,orientations);

    for (int c=0; c<num_cell; c++)
      for (int b=0; b<num_basis; b++)
        for (int p=0; p<num_ip; p++)
          for (int d=0; d<num_dim; d++)
           basis_vector(c, b, p, d) *= orientations(c, b);

    if(compute_derivatives) {
      // Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(curl_basis_scalar,orientations);
      for (int c=0; c<num_cell; c++)
        for (int b=0; b<num_basis; b++)
          for (int p=0; p<num_ip; p++)
            curl_basis_scalar(c, b, p) *= orientations(c, b);
    }

    // setup the orientations for the test space
    if(build_weighted) {
      Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(weighted_basis_vector,orientations);

      if(compute_derivatives)
        Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(weighted_curl_basis_scalar,orientations);
    }
  }
  else if(elmtspace==PureBasis::HCURL && num_dim==3) {

    // setup the orientations for the trial space
    // Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(basis_vector,orientations);

    for (int c=0; c<num_cell; c++)
      for (int b=0; b<num_basis; b++)
        for (int p=0; p<num_ip; p++)
          for (int d=0; d<num_dim; d++)
           basis_vector(c, b, p, d) *= orientations(c, b);

    if(compute_derivatives) {
      // Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(curl_basis_vector,orientations);
      for (int c=0; c<num_cell; c++)
        for (int b=0; b<num_basis; b++)
          for (int p=0; p<num_ip; p++)
            for (int d=0; d<num_dim; d++)
              curl_basis_vector(c, b, p,d) *= orientations(c, b);
    }

    // setup the orientations for the test space
    if(build_weighted) {
      Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(weighted_basis_vector,orientations);

      if(compute_derivatives)
        Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(weighted_curl_basis_vector,orientations);
    }
  }
  else if(elmtspace==PureBasis::HDIV) {
    // setup the orientations for the trial space
    // Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(basis_vector,orientations);

    for (int c=0; c<num_cell; c++)
      for (int b=0; b<num_basis; b++)
        for (int p=0; p<num_ip; p++)
          for (int d=0; d<num_dim; d++)
           basis_vector(c, b, p, d) *= orientations(c, b);

    if(compute_derivatives) {
      // Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(div_basis,orientations);

      for (int c=0; c<num_cell; c++)
        for (int b=0; b<num_basis; b++)
          for (int p=0; p<num_ip; p++)
            div_basis(c, b, p) *= orientations(c, b);
    }

    // setup the orientations for the test space
    if(build_weighted) {
      Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(weighted_basis_vector,orientations);

      if(compute_derivatives)
        Intrepid::FunctionSpaceTools::applyFieldSigns<Scalar>(weighted_div_basis,orientations);
    }
  }
}

template <typename Scalar>
PureBasis::EElementSpace BasisValues2<Scalar>::getElementSpace() const
{ return basis_layout->getBasis()->getElementSpace(); }

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
setupArrays(const Teuchos::RCP<const panzer::BasisIRLayout>& layout,
            bool computeDerivatives)
{
  MDFieldArrayFactory af(prefix,alloc_arrays);

  compute_derivatives = computeDerivatives;
  basis_layout = layout;
  Teuchos::RCP<const panzer::PureBasis> basisDesc = layout->getBasis();

  // for convience pull out basis and quadrature information
  int num_quad = layout->numPoints();
  int dim      = basisDesc->dimension();
  int card     = basisDesc->cardinality();
  int numcells = basisDesc->numCells();
  panzer::PureBasis::EElementSpace elmtspace = basisDesc->getElementSpace();
  Teuchos::RCP<const shards::CellTopology> cellTopo = basisDesc->getCellTopology();
  
  intrepid_basis = basisDesc->getIntrepidBasis<Scalar,ArrayDynamic>();
  
  // allocate field containers
  // field sizes defined by http://trilinos.sandia.gov/packages/docs/dev/packages/intrepid/doc/html/basis_page.html#basis_md_array_sec

  // compute basis fields
  if(elmtspace==panzer::PureBasis::HGRAD) {
     // HGRAD is a nodal field

     // build values
     ///////////////////////////////////////////////
     basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("basis_ref",card,num_quad); // F, P
     basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("basis",numcells,card,num_quad);

     if(build_weighted)
       weighted_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("weighted_basis",numcells,card,num_quad);

     // build gradients
     ///////////////////////////////////////////////

     if(compute_derivatives) {
       grad_basis_ref = af.buildStaticArray<Scalar,BASIS,IP,Dim>("grad_basis_ref",card,num_quad,dim); // F, P, D
       grad_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("grad_basis",numcells,card,num_quad,dim);

       if(build_weighted)
         weighted_grad_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("weighted_grad_basis",numcells,card,num_quad,dim);
     }

     // build curl
     ///////////////////////////////////////////////

     // nothing - HGRAD does not support CURL operation
  }
  else if(elmtspace==panzer::PureBasis::HCURL) {
     // HCURL is a vector field

     // build values
     ///////////////////////////////////////////////

     basis_ref_vector = af.buildStaticArray<Scalar,BASIS,IP,Dim>("basis_ref",card,num_quad,dim); // F, P, D
     basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("basis",numcells,card,num_quad,dim);

     if(build_weighted)
       weighted_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("weighted_basis",numcells,card,num_quad,dim);

     // build gradients
     ///////////////////////////////////////////////

     // nothing - HCURL does not support GRAD operation

     // build curl
     ///////////////////////////////////////////////

     if(compute_derivatives) {
       if(dim==2) {
          // curl of HCURL basis is not dimension dependent
          curl_basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("curl_basis_ref",card,num_quad); // F, P
          curl_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("curl_basis",numcells,card,num_quad);
  
          if(build_weighted)
            weighted_curl_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("weighted_curl_basis",numcells,card,num_quad);
       }
       else if(dim==3){
          curl_basis_ref_vector = af.buildStaticArray<Scalar,BASIS,IP,Dim>("curl_basis_ref",card,num_quad,dim); // F, P, D
          curl_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("curl_basis",numcells,card,num_quad,dim);
  
          if(build_weighted)
            weighted_curl_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("weighted_curl_basis",numcells,card,num_quad,dim);
       }
       else { TEUCHOS_ASSERT(false); } // what do I do with 1D?
     }
  }
  else if(elmtspace==panzer::PureBasis::HDIV) {
     // HDIV is a vector field

     // build values
     ///////////////////////////////////////////////

     basis_ref_vector = af.buildStaticArray<Scalar,BASIS,IP,Dim>("basis_ref",card,num_quad,dim); // F, P, D
     basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("basis",numcells,card,num_quad,dim);

     if(build_weighted)
       weighted_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("weighted_basis",numcells,card,num_quad,dim);

     // build gradients
     ///////////////////////////////////////////////

     // nothing - HCURL does not support GRAD operation

     // build curl
     ///////////////////////////////////////////////

     // nothing - HDIV does not support CURL operation

     // build div
     ///////////////////////////////////////////////

     if(compute_derivatives) {
       div_basis_ref = af.buildStaticArray<Scalar,BASIS,IP>("div_basis_ref",card,num_quad); // F, P
       div_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP>("div_basis",numcells,card,num_quad);
  
       if(build_weighted)
         weighted_div_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP>("weighted_div_basis",numcells,card,num_quad);
     }
  }
  else if(elmtspace==panzer::PureBasis::CONST) {
     // CONST is a nodal field

     // build values
     ///////////////////////////////////////////////
     basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("basis_ref",card,num_quad); // F, P
     basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("basis",numcells,card,num_quad);

     if(build_weighted)
       weighted_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("weighted_basis",numcells,card,num_quad);

     // build gradients
     ///////////////////////////////////////////////

     // nothing - CONST does not support GRAD operation

     // build curl
     ///////////////////////////////////////////////

     // nothing - CONST does not support CURL operation

     // build div 
     ///////////////////////////////////////////////

     // nothing - CONST does not support DIV operation
  }
  else { TEUCHOS_ASSERT(false); }

  basis_coordinates_ref = af.buildStaticArray<Scalar,BASIS,Dim>("basis_coordinates_ref",card,dim);
  basis_coordinates = af.buildStaticArray<Scalar,Cell,BASIS,Dim>("basis_coordinates",numcells,card,dim);
}

// do some explicit instantiation so things compile faster.

#define BASIS_VALUES_INSTANTIATION(SCALAR) \
template class BasisValues2<SCALAR>;

BASIS_VALUES_INSTANTIATION(panzer::Traits::RealType)
BASIS_VALUES_INSTANTIATION(panzer::Traits::FadType)
#ifdef HAVE_STOKHOS
  BASIS_VALUES_INSTANTIATION(panzer::Traits::SGType)
  BASIS_VALUES_INSTANTIATION(panzer::Traits::SGFadType)
#endif

} // namespace panzer
