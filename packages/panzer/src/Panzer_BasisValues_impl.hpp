#ifndef PANZER_BASIS_VALUES_IMPL_HPP
#define PANZER_BASIS_VALUES_IMPL_HPP

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_Basis.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"

// ***********************************************************
// * Specializations of setupArrays() for different array types
// ***********************************************************


// * Specialization for Intrepid::FieldContainer<double>

namespace panzer {
  
  template<>
  inline
  void panzer::BasisValues<double,Intrepid::FieldContainer<double> >::
  setupArrays(const Teuchos::RCP<panzer::BasisIRLayout>& b)
  {
    basis_layout = b;
    Teuchos::RCP<const panzer::PureBasis> basisDesc = b->getBasis();

    // for convience pull out basis and quadrature information
    int num_quad = b->getNumIntPoints();
    int dim      = basisDesc->getDimension();
    int card     = basisDesc->getCardinality();
    int numcells = basisDesc->getNumCells();
    panzer::PureBasis::EElementSpace elmtspace = basisDesc->getElementSpace();
    
    intrepid_basis = basisDesc->getIntrepidBasis();
    
    // allocate field containers
    // field sizes defined by http://trilinos.sandia.gov/packages/docs/dev/packages/intrepid/doc/html/basis_page.html#basis_md_array_sec
  
    // compute basis fields
    if(elmtspace==panzer::PureBasis::HGRAD) {
       // HGRAD is a nodal field

       basis_ref = Intrepid::FieldContainer<double>(card,num_quad); // F, P
       basis = Intrepid::FieldContainer<double>(numcells,card,num_quad);

       weighted_basis = Intrepid::FieldContainer<double>(numcells,card,num_quad);
    }
    else if(elmtspace==panzer::PureBasis::HCURL) {
       // HCURL is a vector field

       basis_ref = Intrepid::FieldContainer<double>(card,num_quad,dim); // F, P, D
       basis = Intrepid::FieldContainer<double>(numcells,card,num_quad,dim);

       weighted_basis = Intrepid::FieldContainer<double>(numcells,card,num_quad,dim);
    }
    else { TEUCHOS_ASSERT(false); }

    // Gradient operator
    if(elmtspace==panzer::PureBasis::HGRAD) {
       grad_basis_ref = Intrepid::FieldContainer<double>(card,num_quad,dim); // F, P, D
       grad_basis = Intrepid::FieldContainer<double>(numcells,card,num_quad,dim);

       weighted_grad_basis = Intrepid::FieldContainer<double>(numcells,card,num_quad,dim);
    }
    else if(elmtspace==panzer::PureBasis::HCURL) {
       // nothing
    }
    else { TEUCHOS_ASSERT(false); }

    // Curl operator
    if(elmtspace==panzer::PureBasis::HGRAD) {
       // nothing
    }
    else if(elmtspace==panzer::PureBasis::HCURL) {
       if(dim==2) {
          // curl of HCURL basis is not dimension dependent
          curl_basis_ref = Intrepid::FieldContainer<double>(card,num_quad); // F, P
          curl_basis = Intrepid::FieldContainer<double>(numcells,card,num_quad);

          weighted_curl_basis = Intrepid::FieldContainer<double>(numcells,card,num_quad);
       }
       else if(dim==3){
          curl_basis_ref = Intrepid::FieldContainer<double>(card,num_quad,dim); // F, P, D
          curl_basis = Intrepid::FieldContainer<double>(numcells,card,num_quad,dim);

          weighted_curl_basis = Intrepid::FieldContainer<double>(numcells,card,num_quad,dim);
       }
       else { TEUCHOS_ASSERT(false); } // what do I do with 1D?
    }
    else { TEUCHOS_ASSERT(false); }
    
    basis_coordinates_ref = Intrepid::FieldContainer<double>(card,dim);
    basis_coordinates = Intrepid::FieldContainer<double>(numcells,card,dim);
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
    // first grab basis descriptor
    PureBasis::EElementSpace elmtspace = getElementSpace();

    intrepid_basis->getValues(basis_ref, cub_points, 
			      Intrepid::OPERATOR_VALUE);
    
    if(elmtspace==PureBasis::HGRAD) {
       intrepid_basis->getValues(grad_basis_ref, cub_points, 
   			         Intrepid::OPERATOR_GRAD);
    
       Intrepid::FunctionSpaceTools::
         HGRADtransformVALUE<double>(basis,
   				  basis_ref);
       
       Intrepid::FunctionSpaceTools::
         HGRADtransformGRAD<double>(grad_basis,
   				 jac_inv,
   				 grad_basis_ref);

       Intrepid::FunctionSpaceTools::
         multiplyMeasure<double>(weighted_basis, 
   			         weighted_measure, 
			         basis);
    
       Intrepid::FunctionSpaceTools::
         multiplyMeasure<double>(weighted_grad_basis, 
	   		         weighted_measure, 
			         grad_basis);
    }
    else if(elmtspace==PureBasis::HCURL) {
       intrepid_basis->getValues(curl_basis_ref, cub_points, 
   			         Intrepid::OPERATOR_CURL);

       Intrepid::FunctionSpaceTools::
         HCURLtransformVALUE<double>(basis,
                                     jac_inv,
   		  		     basis_ref);
       
       Intrepid::FunctionSpaceTools::
         HCURLtransformCURL<double>(curl_basis,
   				    jac,
   				    jac_det,
   				    curl_basis_ref);

       Intrepid::FunctionSpaceTools::
         multiplyMeasure<double>(weighted_basis, 
	   		         weighted_measure, 
			         basis);

       Intrepid::FunctionSpaceTools::
         multiplyMeasure<double>(weighted_curl_basis, 
	   		         weighted_measure, 
			         curl_basis);
   }
    

    // If basis supports coordinate values at basis points, then
    // compute these values
    {
      Teuchos::RCP<Intrepid::DofCoordsInterface<Array> > coords;
      coords = Teuchos::rcp_dynamic_cast<Intrepid::DofCoordsInterface<Array> >(intrepid_basis);
      if (!Teuchos::is_null(coords)) {
	coords->getDofCoords(basis_coordinates_ref);
	Intrepid::CellTools<Scalar> cell_tools;
	cell_tools.mapToPhysicalFrame(basis_coordinates, 
				      basis_coordinates_ref,
				      node_coordinates,
				      basis_layout->getIntrepidBasis()->getBaseCellTopology());
      }
    }

  }

  template<typename Scalar, typename Array>
  PureBasis::EElementSpace BasisValues<Scalar,Array>::getElementSpace() const
  { return basis_layout->getBasis()->getElementSpace(); }

}

#endif
