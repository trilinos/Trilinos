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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_PermutationFactory_def.hpp
 *
 *  Created on: Nov 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_PERMUTATIONFACTORY_DEF_HPP_
#define MUELU_PERMUTATIONFACTORY_DEF_HPP_

#include <vector>
#include <queue>

#include "MueLu_PermutationFactory_decl.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_StridedMap.hpp>    // for nDofsPerNode...
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_AlgebraicPermutationStrategy.hpp"

#undef DEBUG_OUTPUT

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PermutationFactory()
  { }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~PermutationFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<const ParameterList> PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set< RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A to be permuted.");

  validParamList->set< std::string >           ("PermutationRowMapName", "", "Name of input row map for which rows the permutation shall be done. (default='')");
  validParamList->set< RCP<const FactoryBase> >("PermutationRowMapFactory", Teuchos::null, "Generating factory of the input row map for the permutation.");

  validParamList->set< std::string >           ("PermutationStrategy", "Algebraic", "Permutation strategy (default = 'Algebraic', 'Local'");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "A");

  const ParameterList & pL = GetParameterList();
  std::string mapName                        = pL.get<std::string> ("PermutationRowMapName");
  Teuchos::RCP<const FactoryBase> mapFactory = GetFactory          ("PermutationRowMapFactory");

  if(mapName.length() > 0 ) {
    currentLevel.DeclareInput(mapName,mapFactory.get(),this);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
  FactoryMonitor m(*this, "Permutation Factory ", currentLevel);

  Teuchos::RCP<Matrix> A = Get< Teuchos::RCP<Matrix> > (currentLevel, "A");

  const ParameterList & pL = GetParameterList();
  std::string mapName                        = pL.get<std::string> ("PermutationRowMapName");
  Teuchos::RCP<const FactoryBase> mapFactory = GetFactory          ("PermutationRowMapFactory");

  Teuchos::RCP<const Map> permRowMap = Teuchos::null;
  if(mapName.length() > 0 ) {
    permRowMap = currentLevel.Get<RCP<const Map> >(mapName,mapFactory.get());
  } else {
    permRowMap = A->getRowMap(); // use full row map of A
  }

  /*TEUCHOS_TEST_FOR_EXCEPTION(true,
                              std::logic_error,
                              "`MatrixType' has incorrect value (" << MatrixType << ") in input to function CreateCrsMatrix()."
                              << "Check the documentation for a list of valid choices");
*/

  std::string strStrategy = pL.get<std::string> ("PermutationStrategy");
  if( strStrategy == "Algebraic" ) {
  // TODO switch between different permutation strategies
  Teuchos::RCP<AlgebraicPermutationStrategy> permStrat = Teuchos::rcp(new AlgebraicPermutationStrategy());
  permStrat->BuildPermutation(A,permRowMap,currentLevel,this);
  }


}

} // namespace MueLu


#endif /* MUELU_PERMUTATIONFACTORY_DEF_HPP_ */
