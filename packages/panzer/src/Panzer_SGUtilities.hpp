#ifndef __Panzer_SGUtilities_hpp__

#include "Panzer_config.hpp"

#ifdef HAVE_STOKHOS

#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Panzer_Traits.hpp"

#include <vector>

namespace panzer {
namespace sg_utils {

   //! Converts a vector value into a scalar depending on the scalar type
   template <typename ScalarT>
   inline void vectorToValue(const std::vector<double> & in_vector,
                             const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion,
                             ScalarT & value)
   { value = in_vector[0]; }


   //! Converts a vector value into a scalar depending on the scalar type (specific for SGResidual)
   inline void vectorToValue(const std::vector<double> & in_vector,
                             const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion,
                             panzer::Traits::SGResidual::ScalarT & value)
   { 
      TEUCHOS_ASSERT(!Teuchos::is_null(expansion)); // fail if expansion is null

      value = 0.0;
      value.reset(expansion); // jamb in expansion here
      value.copyForWrite();

      TEUCHOS_ASSERT((int) in_vector.size()<=value.size()); // fail if too many coefficients

      // now blindly use fastAccessCoeff
      for(std::size_t i=0;i<in_vector.size();i++)
         value.fastAccessCoeff(i) = in_vector[i];
      for(int i=in_vector.size();i<value.size();i++) // make sure to zero out unused values
         value.fastAccessCoeff(i) = 0.0;
   }

   //! Converts a vector value into a scalar depending on the scalar type (specific for SGJacobian)
   inline void vectorToValue(const std::vector<double> & in_vector,
                             const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion,
                             panzer::Traits::SGJacobian::ScalarT & value)
   { 
      TEUCHOS_ASSERT(!Teuchos::is_null(expansion)); // fail if expansion is null

      value = 0.0;
      value.val().reset(expansion); // jamb in expansion here
      value.val().copyForWrite();

      TEUCHOS_ASSERT((int) in_vector.size()<=value.val().size()); // fail if too many coefficients

      // now blindly use fastAccessCoeff
      for(std::size_t i=0;i<in_vector.size();i++)
         value.val().fastAccessCoeff(i) = in_vector[i];
      for(int i=in_vector.size();i<value.val().size();i++)
         value.val().fastAccessCoeff(i) = 0.0;
   }
}
}

#endif
#endif 
