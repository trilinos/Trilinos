
#ifndef PANZER_BASIS_VALUES_UTILITIES_HPP
#define PANZER_BASIS_VALUES_UTILITIES_HPP

#include "Teuchos_RCP.hpp"

#include "Intrepid_Basis.hpp"
#include "Panzer_Basis.hpp"

namespace panzer {

  template<typename Scalar,typename Array>
  struct BasisValues { 
    
    //! Sizes/allocates memory for arrays
    void setupArrays(const Teuchos::RCP<panzer::Basis>& basis);
   
    void evaluateValues(const Array& cub_points,
			const Array& jac_inv,
			const Array& weighted_measure);

    Array basis_ref;           // <BASIS,IP>
    Array basis;               // <Cell,BASIS,IP>
    Array grad_basis_ref;      // <BASIS,IP,Dim>
    Array grad_basis;          // <Cell,BASIS,IP,Dim>
    Array weighted_basis;      // <Cell,BASIS,IP>
    Array weighted_grad_basis; // <Cell,BASIS,IP,Dim>

    Teuchos::RCP<panzer::Basis> panzer_basis;
    
    Teuchos::RCP<Intrepid::Basis<double,Array> > intrepid_basis;
  };

} // namespace panzer

#include "Panzer_BasisValuesT.hpp"

#endif

