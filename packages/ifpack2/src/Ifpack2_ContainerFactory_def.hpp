/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_CONTAINERFACTORY_DEF_H
#define IFPACK2_CONTAINERFACTORY_DEF_H

#include "Ifpack2_ContainerFactory_decl.hpp"
#include "Ifpack2_TriDiContainer.hpp"
#include "Ifpack2_DenseContainer.hpp"
#include "Ifpack2_SparseContainer.hpp"
#include "Ifpack2_BandedContainer.hpp"
#include "Ifpack2_BlockTriDiContainer.hpp"
#include "Ifpack2_ILUT.hpp"
#include "Teuchos_ArrayView.hpp"

#include <sstream>

namespace Ifpack2 {

template<typename MatrixType>
void ContainerFactory<MatrixType>::
registerDefaults()
{
  registerContainer<Ifpack2::TriDiContainer<MatrixType, scalar_type>>("TriDi");
  registerContainer<Ifpack2::DenseContainer<MatrixType, scalar_type>>("Dense");
  registerContainer<Ifpack2::BandedContainer<MatrixType, scalar_type>>("Banded");
  registerContainer<SparseContainer<MatrixType, ILUT<MatrixType>>>("SparseILUT");
#ifdef HAVE_IFPACK2_AMESOS2
  registerContainer<SparseContainer<MatrixType, Details::Amesos2Wrapper<MatrixType>>>("SparseAmesos");
  registerContainer<SparseContainer<MatrixType, Details::Amesos2Wrapper<MatrixType>>>("SparseAmesos2");
#endif
#ifdef HAVE_IFPACK2_EXPERIMENTAL_KOKKOSKERNELS_FEATURES
  registerContainer<Ifpack2::BlockTriDiContainer<MatrixType>>("BlockTriDi");
#endif
  registeredDefaults = true;
}

template<typename MatrixType>
template<typename ContainerType>
void ContainerFactory<MatrixType>::
registerContainer(std::string containerType)
{
  //overwrite any existing registration with the same name
  table[containerType] = Teuchos::rcp(new Details::ContainerFactoryEntry<MatrixType, ContainerType>());
}

template<typename MatrixType>
Teuchos::RCP<typename ContainerFactory<MatrixType>::BaseContainer>
ContainerFactory<MatrixType>::
build(std::string containerType, const Teuchos::RCP<const MatrixType>& A,
    const Teuchos::Array<Teuchos::Array<local_ordinal_type>>& localRows, const Teuchos::RCP<const import_type> importer,
    int OverlapLevel, scalar_type DampingFactor)
{
  if(!registeredDefaults)
  {
    registerDefaults();
  }
  //In the case that Amesos2 isn't enabled, provide a better error message than the generic one
  #ifndef HAVE_IFPACK2_AMESOS2
  if(containerType == "SparseAmesos" || containerType == "SparseAmesos2")
  {
    throw std::invalid_argument("Container type SparseAmesos (aka SparseAmesos2) was requested but Amesos2 isn't enabled.\n"
                                "Add the CMake option \"-D Trilinos_ENABLE_Amesos2=ON\" to enable it.");
  }
  #endif
  auto it = table.find(containerType);
  if(it == table.end())
  {
    std::ostringstream oss;
    oss << "Container type \"" << containerType << "\" not registered.\n";
    oss << "Call ContainerFactory<MatrixType>::registerContainer<ContainerType>(containerName) first.\n";
    oss << "Currently registered Container types: ";
    for(auto r : table)
    {
      oss << '\"' << r.first << "\" ";
    }
    //remove the single trailing space from final message
    auto str = oss.str();
    str = str.substr(0, str.length() - 1);
    throw std::invalid_argument(str);
  }
  return it->second->build(A, localRows, importer, OverlapLevel, DampingFactor);
}

template<typename MatrixType>
void ContainerFactory<MatrixType>::
deregisterContainer(std::string containerType)
{
  auto it = table.find(containerType);
  if(it != table.end())
  {
    table.erase(it);
  }
}

// Definitions of static data

template<typename MatrixType>
std::map<std::string, Teuchos::RCP<Details::ContainerFactoryEntryBase<MatrixType>>> ContainerFactory<MatrixType>::table;

template<typename MatrixType>
bool ContainerFactory<MatrixType>::registeredDefaults;     //this will initially be false


} // namespace Ifpack2

#define IFPACK2_CONTAINERFACTORY_INSTANT(S,LO,GO,N) \
template struct Ifpack2::ContainerFactory<Tpetra::RowMatrix<S, LO, GO, N>>;

#endif // IFPACK2_DETAILS_CONTAINERFACTORY_H
