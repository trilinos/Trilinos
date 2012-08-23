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
  template <typename Scalar,typename Array>
  const Array BasisValues<Scalar,Array>::dummyArray;    

  template<typename Scalar,typename Array>
  template <typename ArrayFactory>
  void panzer::BasisValues<Scalar,Array>::
  setupArrays(const Teuchos::RCP<const panzer::BasisIRLayout>& layout,
              const ArrayFactory & af)
  {
    basis_layout = layout;
    Teuchos::RCP<const panzer::PureBasis> basisDesc = layout->getBasis();

    // for convience pull out basis and quadrature information
    int num_quad = layout->getNumPoints();
    int dim      = basisDesc->getDimension();
    int card     = basisDesc->getCardinality();
    int numcells = basisDesc->getNumCells();
    panzer::PureBasis::EElementSpace elmtspace = basisDesc->getElementSpace();
    Teuchos::RCP<const shards::CellTopology> cellTopo = basisDesc->getCellTopology();
    
    intrepid_basis = basisDesc->getIntrepidBasis<Scalar,Array>();
    
    // allocate field containers
    // field sizes defined by http://trilinos.sandia.gov/packages/docs/dev/packages/intrepid/doc/html/basis_page.html#basis_md_array_sec
  
    // compute basis fields
    if(elmtspace==panzer::PureBasis::HGRAD) {
       // HGRAD is a nodal field

       // build values
       ///////////////////////////////////////////////

       basis_ref = af.template buildArray<Scalar,BASIS,IP>("basis_ref",card,num_quad); // F, P
       basis = af.template buildArray<Scalar,Cell,BASIS,IP>("basis",numcells,card,num_quad);

       weighted_basis = af.template buildArray<Scalar,Cell,BASIS,IP>("weighted_basis",numcells,card,num_quad);

       // build gradients
       ///////////////////////////////////////////////

       grad_basis_ref = af.template buildArray<Scalar,BASIS,IP,Dim>("grad_basis_ref",card,num_quad,dim); // F, P, D
       grad_basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("grad_basis",numcells,card,num_quad,dim);

       weighted_grad_basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("weighted_grad_basis",numcells,card,num_quad,dim);

       // build curl
       ///////////////////////////////////////////////

       // nothing - HGRAD does not support CURL operation
    }
    else if(elmtspace==panzer::PureBasis::HCURL) {
       // HCURL is a vector field

       // build values
       ///////////////////////////////////////////////

       basis_ref = af.template buildArray<Scalar,BASIS,IP,Dim>("basis_ref",card,num_quad,dim); // F, P, D
       basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("basis",numcells,card,num_quad,dim);

       weighted_basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("weighted_basis",numcells,card,num_quad,dim);

       // build gradients
       ///////////////////////////////////////////////

       // nothing - HCURL does not support GRAD operation

       // build curl
       ///////////////////////////////////////////////

       if(dim==2) {
          // curl of HCURL basis is not dimension dependent
          curl_basis_ref = af.template buildArray<Scalar,BASIS,IP>("curl_basis_ref",card,num_quad); // F, P
          curl_basis = af.template buildArray<Scalar,Cell,BASIS,IP>("curl_basis",numcells,card,num_quad);

          weighted_curl_basis = af.template buildArray<Scalar,Cell,BASIS,IP>("weighted_curl_basis",numcells,card,num_quad);
       }
       else if(dim==3){
          curl_basis_ref = af.template buildArray<Scalar,BASIS,IP,Dim>("curl_basis_ref",card,num_quad,dim); // F, P, D
          curl_basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("curl_basis",numcells,card,num_quad,dim);

          weighted_curl_basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("weighted_curl_basis",numcells,card,num_quad,dim);
       }
       else { TEUCHOS_ASSERT(false); } // what do I do with 1D?
    }
    else { TEUCHOS_ASSERT(false); }

    basis_coordinates_ref = af.template buildArray<Scalar,BASIS,Dim>("basis_coordinates_ref",card,dim);
    basis_coordinates = af.template buildArray<Scalar,Cell,BASIS,Dim>("basis_coordinates",numcells,card,dim);
  }
  
  // ***********************************************************
  // * Evaluation of values without weithed measure - NOT specialized
  // ***********************************************************

  template<typename Scalar, typename Array>
  inline
  void panzer::BasisValues<Scalar,Array>::
  evaluateValues(const Array& cub_points,
		 const Array& jac,
 	 	 const Array& jac_det,
		 const Array& jac_inv)
  { 
     // substitute dummy array for weighted measure
     evaluateValues(cub_points,jac,jac_det,jac_inv,dummyArray,dummyArray);
  }

  template<typename Scalar, typename Array>
  inline
  void panzer::BasisValues<Scalar,Array>::
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
  template<typename Scalar, typename Array>
  inline
  void panzer::BasisValues<Scalar,Array>::
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
    int spaceDim = basis_layout->getDimension();

    intrepid_basis->getValues(basis_ref, cub_points, 
			      Intrepid::OPERATOR_VALUE);

    if(elmtspace==PureBasis::HGRAD) {
       intrepid_basis->getValues(grad_basis_ref, cub_points, 
   			         Intrepid::OPERATOR_GRAD);
    
       Intrepid::FunctionSpaceTools::
         HGRADtransformVALUE<Scalar>(basis,
   				  basis_ref);
       
       Intrepid::FunctionSpaceTools::
         HGRADtransformGRAD<Scalar>(grad_basis,
   				 jac_inv,
   				 grad_basis_ref);

       if(buildWeighted) {
          Intrepid::FunctionSpaceTools::
            multiplyMeasure<Scalar>(weighted_basis, 
      			         weighted_measure, 
   			         basis);
       
          Intrepid::FunctionSpaceTools::
            multiplyMeasure<Scalar>(weighted_grad_basis, 
   	   		         weighted_measure, 
   			         grad_basis);
       }
    }
    else if(elmtspace==PureBasis::HCURL) {
       intrepid_basis->getValues(curl_basis_ref, cub_points, 
   			         Intrepid::OPERATOR_CURL);

       Intrepid::FunctionSpaceTools::
         HCURLtransformVALUE<Scalar>(basis,
                                     jac_inv,
   		  		     basis_ref);

       // the CURL transform differs depending on spactial dimension
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

       if(buildWeighted) {
          Intrepid::FunctionSpaceTools::
            multiplyMeasure<Scalar>(weighted_basis, 
   	   		         weighted_measure, 
   			         basis);
   
          Intrepid::FunctionSpaceTools::
            multiplyMeasure<Scalar>(weighted_curl_basis, 
   	   		         weighted_measure, 
   			         curl_basis);
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

  template<typename Scalar, typename Array>
  PureBasis::EElementSpace BasisValues<Scalar,Array>::getElementSpace() const
  { return basis_layout->getBasis()->getElementSpace(); }

}

#endif
