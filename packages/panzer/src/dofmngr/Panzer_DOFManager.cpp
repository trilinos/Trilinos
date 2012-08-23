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
