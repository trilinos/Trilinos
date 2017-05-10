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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_MUELU_SRC_GRAPH_MUELU_UNSMOOSHFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_GRAPH_MUELU_UNSMOOSHFACTORY_DEF_HPP_


#include "MueLu_Monitor.hpp"

#include "MueLu_UnsmooshFactory_decl.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::UnsmooshFactory() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("P",                  Teuchos::null, "Generating factory of the prolongator P");
    validParamList->set< RCP<const FactoryBase> >("DofStatus",          Teuchos::null, "Generating factory for dofStatus array (usually the VariableDofLaplacdianFactory)");

    validParamList->set< int  >                  ("maxDofPerNode", 1,     "Maximum number of DOFs per node");
    validParamList->set< bool >                  ("fineIsPadded" , false, "true if finest level input matrix is padded");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    //const ParameterList& pL = GetParameterList();
    Input(currentLevel, "P");
    Input(currentLevel, "DofStatus");
    //Input(currentLevel, "Coordinates");

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);
    typedef Teuchos::ScalarTraits<SC> STS;

    const ParameterList  & pL = GetParameterList();

    RCP<Matrix> amalgP = Get< RCP<Matrix> >(currentLevel, "P");

    //Teuchos::RCP< const Teuchos::Comm< int > > comm = amalgP->getRowMap()->getComm();
    //Xpetra::UnderlyingLib lib = amalgP->getRowMap()->lib();

    Teuchos::Array<char> dofStatus = Get<Teuchos::Array<char> >(currentLevel, "DofStatus");

    /*for (size_t i = 0; i < dofStatus.size(); i++) {
      std::cout << i << " " << dofStatus[i] << std::endl;
    }*/

    int maxDofPerNode = pL.get<int> ("maxDofPerNode");
    bool fineIsPadded = pL.get<bool>("fineIsPadded");

    // extract CRS information from amalgamated prolongation operator
    Teuchos::ArrayRCP<const size_t> amalgRowPtr(amalgP->getNodeNumRows());
    Teuchos::ArrayRCP<const LocalOrdinal> amalgCols(amalgP->getNodeNumEntries());
    Teuchos::ArrayRCP<const Scalar> amalgVals(amalgP->getNodeNumEntries());
    Teuchos::RCP<CrsMatrixWrap> amalgPwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(amalgP);
    Teuchos::RCP<CrsMatrix> amalgPcrs = amalgPwrap->getCrsMatrix();
    amalgPcrs->getAllValues(amalgRowPtr, amalgCols, amalgVals);

    // calculate number of dof rows for new prolongator
    size_t paddedNrows = amalgP->getRowMap()->getNodeNumElements() * Teuchos::as<size_t>(maxDofPerNode);

    // reserve CSR arrays for new prolongation operator
    Teuchos::ArrayRCP<size_t> newPRowPtr(paddedNrows+1);
    Teuchos::ArrayRCP<LocalOrdinal> newPCols(amalgP->getNodeNumEntries() * maxDofPerNode);
    Teuchos::ArrayRCP<Scalar> newPVals(amalgP->getNodeNumEntries() * maxDofPerNode);

    size_t rowCount = 0; // actual number of (local) in unamalgamated prolongator
    if(fineIsPadded == true) {
      // build prolongation operator for padded fine level matrices.
      // Note: padded fine level dofs are transfered by injection.
      // That is, these interpolation stencils do not take averages of
      // coarse level variables. Further, fine level Dirichlet points
      // also use injection.

      size_t cnt = 0; // local id counter
      for (size_t i = 0; i < amalgRowPtr.size() - 1; i++) {
        // determine number of entries in amalgamated dof row i
        size_t rowLength = amalgRowPtr[i+1] - amalgRowPtr[i];

        // loop over dofs per node (unamalgamation)
        for(int j = 0; j < maxDofPerNode; j++) {
          newPRowPtr[i*maxDofPerNode+j] = cnt;
          if (dofStatus[i*maxDofPerNode+j] == 's') { // add only "standard" dofs to unamalgamated prolongator
            // loop over column entries in amalgamated P
            for (size_t k = 0; k < rowLength; k++) {
              newPCols[cnt  ] = amalgCols[k+amalgRowPtr[i]] * maxDofPerNode + j;
              newPVals[cnt++] = amalgVals[k+amalgRowPtr[i]];
            }

          }
        }
      }

      newPRowPtr[paddedNrows] = cnt; // close row CSR array
      rowCount = paddedNrows;
    } else {
      // Build prolongation operator for non-padded fine level matrices.
      // Need to map from non-padded dofs to padded dofs. For this, look
      // at the status array and skip padded dofs.

      size_t cnt = 0; // local id counter
      for (size_t i = 0; i < amalgRowPtr.size() - 1; i++) {
        // determine number of entries in amalgamated dof row i
        size_t rowLength = amalgRowPtr[i+1] - amalgRowPtr[i];

        // loop over dofs per node (unamalgamation)
        for(int j = 0; j < maxDofPerNode; j++) {
          // no interpolation for padded fine dofs as they do not exist

          if (dofStatus[i*maxDofPerNode+j] == 's') { // add only "standard" dofs to unamalgamated prolongator
            newPRowPtr[rowCount++] = cnt;
            // loop over column entries in amalgamated P
            for (size_t k = 0; k < rowLength; k++) {
              newPCols[cnt  ] = amalgCols[k+amalgRowPtr[i]] * maxDofPerNode + j;
              newPVals[cnt++] = amalgVals[k+amalgRowPtr[i]];
            }

          }
          if (dofStatus[i*maxDofPerNode+j] == 'd') { // Dirichlet handling
            newPRowPtr[rowCount++] = cnt;
          }
        }
      }

      newPRowPtr[rowCount] = cnt; // close row CSR array
    } // fineIsPadded == false

    // TODO assemble new P



    Set(currentLevel,"P",amalgP);
  }


} /* MueLu */


#endif /* PACKAGES_MUELU_SRC_GRAPH_MUELU_UNSMOOSHFACTORY_DEF_HPP_ */
