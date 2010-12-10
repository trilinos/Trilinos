#include "Panzer_DOFManager.hpp"

// FEI includes
#include "fei_Factory_Trilinos.hpp"

// build it just for fun
#ifndef NO_EXPLICIT_TEMPLATE_INSTANTIATION

template class panzer::DOFManager<int,long int>;
template class panzer::DOFManager<int,int>;
template class panzer::DOFManager<char,long int>;

#endif

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
