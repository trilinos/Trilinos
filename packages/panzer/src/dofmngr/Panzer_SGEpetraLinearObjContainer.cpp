#include "Panzer_config.hpp"

#ifdef HAVE_STOKHOS

#include "Panzer_SGEpetraLinearObjContainer.hpp"

namespace panzer {

SGEpetraLinearObjContainer::SGEpetraLinearObjContainer(const SGEpetraLinearObjContainer::CoeffVector & coeffs,
                                                       const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion)
   : coeffs_(coeffs), expansion_(expansion)
{ }

void SGEpetraLinearObjContainer::initialize()
{
   for(CoeffVector::iterator itr=begin();itr!=end();itr++) 
      (*itr)->initialize();
}

}

#endif
