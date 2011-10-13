#include "Panzer_config.hpp"

#ifdef HAVE_STOKHOS

#include "Panzer_Traits.hpp"
#include "Panzer_SGEpetraLinearObjFactory.hpp"

#ifdef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION

#include "Panzer_SGEpetraLinearObjFactoryT.hpp"

template class panzer::SGEpetraLinearObjFactory<panzer::Traits,int>;
template class panzer::SGEpetraLinearObjFactory<panzer::Traits,short>;

#endif

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
