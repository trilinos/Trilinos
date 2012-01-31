#include "Panzer_config.hpp"

#include "Panzer_DOFManager_decl.hpp"
#include "Panzer_DOFManager_impl.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

template class panzer::DOFManager<int,long int>;
template class panzer::DOFManager<int,int>;
template class panzer::DOFManager<short,int>;
template class panzer::DOFManager<char,long int>;

#endif

// FEI includes
#include "fei_Factory_Trilinos.hpp"


using Teuchos::RCP;

// needed for faster implementation
///////////////////////////////////////////////
namespace panzer {

// Function is "helpers" for DOFManager::getOwnedIndices
///////////////////////////////////////////////////////////////////////////

template < >
void getOwnedIndices_T<int>(const fei::SharedPtr<fei::VectorSpace> & vs,std::vector<int> & indices) 
{
   int numIndices, ni;
   numIndices = vs->getNumIndices_Owned();
   indices.resize(numIndices);

   // directly write to int indices
   vs->getIndices_Owned(numIndices,&indices[0],ni);
}

///////////////////////////////////////////////////////////////////////////

// Function is "helper" for DOFManager::getOwnedAndSharedIndices
///////////////////////////////////////////////////////////////////////////

template < >
void getOwnedAndSharedIndices_T<int>(const fei::SharedPtr<fei::VectorSpace> & vs,std::vector<int> & indices) 
{
   // get the global indices
   vs->getIndices_SharedAndOwned(indices);
}

///////////////////////////////////////////////////////////////////////////

}
