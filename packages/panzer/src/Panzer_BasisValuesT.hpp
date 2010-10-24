
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
  setupArrays(const Teuchos::RCP<panzer::Basis>& b)
  {
    panzer_basis = b;
    
    intrepid_basis = panzer::createIntrepidBasis<double,Intrepid::FieldContainer<double> >(b->name(), b->getDimension());
    
    basis_ref = Intrepid::FieldContainer<double>(b->basis_ref->dimension(0),
						 b->basis_ref->dimension(1));
    
    basis = Intrepid::FieldContainer<double>(b->basis->dimension(0),
					     b->basis->dimension(1),
					     b->basis->dimension(2));
    
    grad_basis_ref = 
      Intrepid::FieldContainer<double>(b->basis_grad_ref->dimension(0),
				       b->basis_grad_ref->dimension(1),
				       b->basis_grad_ref->dimension(2));
    
    grad_basis = 
      Intrepid::FieldContainer<double>(b->basis_grad->dimension(0),
				       b->basis_grad->dimension(1),
				       b->basis_grad->dimension(2),
				       b->basis_grad->dimension(3));
    
    weighted_basis = 
      Intrepid::FieldContainer<double>(b->basis->dimension(0),
				       b->basis->dimension(1),
				       b->basis->dimension(2));
    
    weighted_grad_basis = 
      Intrepid::FieldContainer<double>(b->basis_grad->dimension(0),
				       b->basis_grad->dimension(1),
				       b->basis_grad->dimension(2),
				       b->basis_grad->dimension(3));
    
  }
  
  // ***********************************************************
  // * Evaluation of values - NOT specialized
  // ***********************************************************
  template<typename Scalar, typename Array>
  inline
  void panzer::BasisValues<Scalar,Array>::
  evaluateValues(const Array& cub_points,
		 const Array& jac_inv,
		 const Array& weighted_measure)
  {    
    intrepid_basis->getValues(basis_ref, cub_points, 
			      Intrepid::OPERATOR_VALUE);
    
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

}
