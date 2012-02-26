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

  //! Used for default parameters
  template <typename Scalar,typename Array>
  const Array BasisValues<Scalar,Array>::dummyArray;    
  
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
    Teuchos::RCP<const shards::CellTopology> cellTopo = basisDesc->getCellTopology();
    
    intrepid_basis = basisDesc->getIntrepidBasis();
    
    // allocate field containers
    // field sizes defined by http://trilinos.sandia.gov/packages/docs/dev/packages/intrepid/doc/html/basis_page.html#basis_md_array_sec

    // make sure edge orientations are included
    if(elmtspace==panzer::PureBasis::HCURL) {
       // allocate orientations
       subcell_orientation = Intrepid::FieldContainer<double>(numcells,cellTopo->getEdgeCount());
       basis_orientation = Intrepid::FieldContainer<double>(numcells,intrepid_basis->getCardinality());
    }
  
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
		 const Array& node_coordinates,
                 const Array& orientation)
  {    
    // first grab basis descriptor
    PureBasis::EElementSpace elmtspace = getElementSpace();
    int spaceDim = basis_layout->getDimension();

    intrepid_basis->getValues(basis_ref, cub_points, 
			      Intrepid::OPERATOR_VALUE);


    // make sure edge orientations are included for HCURL bases
    if(elmtspace==panzer::PureBasis::HCURL) {
       TEUCHOS_TEST_FOR_EXCEPTION(&orientation==&dummyArray,std::invalid_argument,
                                  "BasisValues::evaluateValues using HCURL bases requires edge orientations");
       TEUCHOS_TEST_FOR_EXCEPTION(subcell_orientation.size()!=orientation.size(),std::invalid_argument,
                                  "BasisValues::evaluateValues orientation is the wrong size for HCURL bases: " 
                                 << subcell_orientation.size() << " " << orientation.size());

       // simply copy orientations
       for(int i=0;i<subcell_orientation.size();i++)
          subcell_orientation[i] = orientation[i];

       extendOrientationToBasis(PureBasis::HCURL,*intrepid_basis,orientation,basis_orientation);
    }
    
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
         applyFieldSigns<double>(basis,basis_orientation);

       // the CURL transform differs depending on spactial dimension
       // "In 2D the de Rham complex collapses." But I look at it like this:
       // In 2D the curl is simply  $-\partial_x u_1+\partial_y u_0$ which
       // has the same derivative structure as the divergence in 2D 
       // $\partial_x u_0 + \partial_y u_1$ which means the transformation
       // is the same.
       if(spaceDim==2) 
          Intrepid::FunctionSpaceTools::
            HDIVtransformDIV<double>(curl_basis,
   			             jac_det,   // note only volume deformation is needed!
                                                // this relates directly to this being in
                                                // the divergence space in 2D!
   			             curl_basis_ref);
       else if(spaceDim==3)
          Intrepid::FunctionSpaceTools::
            HCURLtransformCURL<double>(curl_basis,
   	                               jac,
   				       jac_det,
   				       curl_basis_ref);
       else
          TEUCHOS_ASSERT(false); // what you doin?

       Intrepid::FunctionSpaceTools::
         applyFieldSigns<double>(curl_basis,basis_orientation);

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

  template<typename Scalar, typename Array>
  void BasisValues<Scalar,Array>::
  extendOrientationToBasis(panzer::PureBasis::EElementSpace space,
                           Intrepid::Basis<double,Array> & intrBasis,
                           const Array & inOrientation,
                           Array & outOrientation) const
  {
     TEUCHOS_ASSERT(space==panzer::PureBasis::HCURL || space==panzer::PureBasis::HDIV);

     shards::CellTopology cellTopo = intrBasis.getBaseCellTopology();

     TEUCHOS_TEST_FOR_EXCEPTION(inOrientation.dimension(0)!=outOrientation.dimension(0),std::invalid_argument,
                                "BasisValues::extendOrientationToBasis workset size for in/out orientations do not match");
     TEUCHOS_TEST_FOR_EXCEPTION(outOrientation.dimension(1)!=intrBasis.getCardinality(),std::invalid_argument,
                                "BasisValues::extendOrientationToBasis dimension(1) of output array is incorrect");
     TEUCHOS_TEST_FOR_EXCEPTION(inOrientation.dimension(1)!=cellTopo.getEdgeCount(),std::invalid_argument,
                                "BasisValues::extendOrientationToBasis dimension(1) of input array is incorrect");

    
     int numcells = inOrientation.dimension(0);
     unsigned dim = cellTopo.getDimension();
     unsigned subcellInd = panzer::PureBasis::HCURL ? 1 : dim-1; // HCURL lives on edges, HDIV lives on dim-1 (edges in 2D, faces in 3D)
 
     const std::vector<std::vector<int> > & ordinalData = intrBasis.getDofOrdinalData()[subcellInd];
     TEUCHOS_TEST_FOR_EXCEPTION(ordinalData.size()!=cellTopo.getSubcellCount(subcellInd),std::logic_error,
                                "BasisValues::extendOrientationToBasis ordinal data size does not match subcell count");

     // initialize with all ones
     for(int i=0;i<outOrientation.size();i++)
        outOrientation[i] = 1.0;

     // loop over subcells of appropriate dimension
     for(std::size_t d=0;d<ordinalData.size();d++) {
        // basis on sub cell index
        for(std::size_t b=0;b<ordinalData[d].size();b++) {
           int fieldIndex = ordinalData[d][b];

           // loop over cells
           for(int c=0;c<numcells;c++) 
              outOrientation(c,fieldIndex) = inOrientation(c,d); 
        } // end for b
     } // end for d
  }

}

#endif
