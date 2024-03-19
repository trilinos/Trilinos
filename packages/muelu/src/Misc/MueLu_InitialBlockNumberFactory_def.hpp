// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_INITIALBLOCKNUMBER_FACTORY_DEF_HPP
#define MUELU_INITIALBLOCKNUMBER_FACTORY_DEF_HPP

#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_IO.hpp"

#include "MueLu_InitialBlockNumberFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> InitialBlockNumberFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: block diagonal: interleaved blocksize");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InitialBlockNumberFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InitialBlockNumberFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);
  const ParameterList& pL = GetParameterList();

  RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, "A");
  LO blocksize  = as<LO>(pL.get<int>("aggregation: block diagonal: interleaved blocksize"));

  GetOStream(Statistics1) << "Generating interleaved blocking with " << blocksize << " equations" << std::endl;
  RCP<LocalOrdinalVector> BlockNumber = LocalOrdinalVectorFactory::Build(A->getRowMap(), false);
  Teuchos::ArrayRCP<LO> bn_data       = BlockNumber->getDataNonConst(0);
  for (LO i = 0; i < (LO)A->getRowMap()->getLocalNumElements(); i++)
    bn_data[i] = i % blocksize;

  Set(currentLevel, "BlockNumber", BlockNumber);
}

}  // namespace MueLu

#endif  // MUELU_INITIALBLOCKNUMBER_FACTORY_DEF_HPP
