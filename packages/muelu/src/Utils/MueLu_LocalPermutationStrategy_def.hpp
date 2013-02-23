/*
 * MueLu_LocalPermutationStrategy_def.hpp
 *
 *  Created on: Feb 19, 2013
 *      Author: tobias
 */

#ifndef MUELU_LOCALPERMUTATIONSTRATEGY_DEF_HPP_
#define MUELU_LOCALPERMUTATIONSTRATEGY_DEF_HPP_

#include <algorithm>

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_LocalPermutationStrategy_decl.hpp"

namespace MueLu {

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LocalPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildPermutation(const Teuchos::RCP<Matrix> & A, const Teuchos::RCP<const Map> permRowMap, Level & currentLevel, const FactoryBase* genFactory) const {

    size_t nDofsPerNode = 1;
    if (A->IsView("stridedMaps")) {
      Teuchos::RCP<const Map> permRowMapStrided = A->getRowMap("stridedMaps");
      nDofsPerNode = Teuchos::rcp_dynamic_cast<const StridedMap>(permRowMapStrided)->getFixedBlockSize();
    }

    //////////////////
    std::vector<std::pair<GlobalOrdinal,GlobalOrdinal> > RowColPairs;

    // loop over local nodes
    // TODO what about nOffset?
    for ( LocalOrdinal node = 0; node < A->getRowMap()->getNodeNumElements()/nDofsPerNode; node++) {

      Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar> subBlockMatrix(nDofsPerNode, nDofsPerNode, true);

      std::vector<GlobalOrdinal> growIds(nDofsPerNode);

      for ( LocalOrdinal lrdof = 0; lrdof < nDofsPerNode; lrdof++) { // TODO more complicated for variable dofs per node
        GlobalOrdinal grow = getGlobalDofId(A, node, lrdof);
        growIds[lrdof] = grow;

        //if(permRowMap->isNodeGlobalElement(grow) == true) continue;

        // extract local row information from matrix
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        A->getLocalRowView(A->getRowMap()->getLocalElement(grow), indices, vals);

        // find column entry with max absolute value
        //GlobalOrdinal gMaxValIdx = 0;
        //Scalar norm1 = 0.0;
        Scalar maxVal = 0.0;
        for (size_t j = 0; j < Teuchos::as<size_t>(indices.size()); j++) {
          //norm1 += std::abs(vals[j]);
          if(std::abs(vals[j]) > maxVal) {
            maxVal = std::abs(vals[j]);
            //gMaxValIdx = A->getColMap()->getGlobalElement(indices[j]);
          }
        }

        GlobalOrdinal grnodeid = globalDofId2globalNodeId(A,grow);

        for (size_t j = 0; j < Teuchos::as<size_t>(indices.size()); j++) {
          GlobalOrdinal gcol = A->getColMap()->getGlobalElement(indices[j]);
          GlobalOrdinal gcnodeid = globalDofId2globalNodeId(A,gcol); // -> global node id
          if (grnodeid == gcnodeid) {
            if(maxVal != 0.0) {
              subBlockMatrix(lrdof, gcol % nDofsPerNode) = vals[j]/maxVal;
            } else
            {
              subBlockMatrix(lrdof, gcol % nDofsPerNode) = vals[j]; // there is a problem
              std::cout << "maxVal never should be zero!!!!" << std::endl;
            }
          }
        }
      }

      // now we have the sub block matrix

      // build permutation string
      std::stringstream ss;
      for(size_t t = 0; t<nDofsPerNode; t++)
        ss << t;
      std::string cs = ss.str();
      std::vector<std::string> result_perms;
      do {
        result_perms.push_back(cs);
        //std::cout << result_perms.back() << std::endl;
      } while (std::next_permutation(cs.begin(),cs.end()));

      std::vector<Scalar> performance_vector = std::vector<Scalar>(result_perms.size());
      for(size_t t = 0; t < result_perms.size(); t++) {
        std::string s = result_perms[t];
        Scalar value = 1.0;
        for(size_t len=0; len<s.length(); len++) {
          int col = static_cast<int>(s[len]-'0');
          value = value * subBlockMatrix(len,col);
        }
        performance_vector[t] = value;
      }

      // find permutation with maximum performance value
      Scalar maxVal = 0.0;
      size_t maxPerformancePermutationIdx = 0;
      for (size_t j = 0; j < Teuchos::as<size_t>(performance_vector.size()); j++) {
        if(std::abs(performance_vector[j]) > maxVal) {
          maxVal = std::abs(performance_vector[j]);
          maxPerformancePermutationIdx = j;
        }
      }

      // the best permutation is
      //std::cout << "The best permutation is " << result_perms[maxPerformancePermutationIdx] << std::endl;

      std::string bestPerformancePermutation = result_perms[maxPerformancePermutationIdx] ;
      for(size_t t = 0; t<nDofsPerNode; t++) {
        int col = static_cast<int>(bestPerformancePermutation[t]-'0');
        RowColPairs.push_back(std::make_pair(growIds[t],growIds[col]));
      }

    } // end loop over local nodes

    //
    /*typename std::vector<std::pair<GlobalOrdinal,GlobalOrdinal> >::iterator pl = RowColPairs.begin();
    while(pl != RowColPairs.end()) {
      std::cout << (*pl).first << " , " << (*pl).second << std::endl;
      pl++;
    }*/

    // build Pperm and Qperm vectors
    Teuchos::RCP<Vector> Pperm = VectorFactory::Build(A->getRowMap());
    Teuchos::RCP<Vector> Qperm = VectorFactory::Build(A->getDomainMap());

    Pperm->putScalar(0.0);
    Qperm->putScalar(0.0);

    Teuchos::ArrayRCP<Scalar> PpermData = Pperm->getDataNonConst(0);
    Teuchos::ArrayRCP<Scalar> QpermData = Qperm->getDataNonConst(0);

    typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::iterator p = RowColPairs.begin();
    while(p != RowColPairs.end() ) {
      GlobalOrdinal ik = (*p).first;
      GlobalOrdinal jk = (*p).second;

      LocalOrdinal lik = A->getRowMap()->getLocalElement(ik);
      LocalOrdinal ljk = A->getDomainMap()->getLocalElement(jk);

      Pperm->replaceLocalValue(lik,ik);
      Qperm->replaceLocalValue(ljk,ik);

      p = RowColPairs.erase(p);
    }

    if(RowColPairs.size()>0) GetOStream(Warnings0,0) << "MueLu::LocalPermutationStrategy: There are Row/col pairs left!" << std::endl;

    // Qperm should be fine
    // build matrices

    // create new empty Matrix
    Teuchos::RCP<CrsMatrixWrap> permPTmatrix = Teuchos::rcp(new CrsMatrixWrap(A->getRowMap(),1,Xpetra::StaticProfile));
    Teuchos::RCP<CrsMatrixWrap> permQTmatrix = Teuchos::rcp(new CrsMatrixWrap(A->getRowMap(),1,Xpetra::StaticProfile));

    for(size_t row=0; row<A->getNodeNumRows(); row++) {
      Teuchos::ArrayRCP<GlobalOrdinal> indoutP(1,Teuchos::as<GO>(PpermData[row])); // column idx for Perm^T
      Teuchos::ArrayRCP<GlobalOrdinal> indoutQ(1,Teuchos::as<GO>(QpermData[row])); // column idx for Qperm
      Teuchos::ArrayRCP<Scalar> valout(1,1.0);
      permPTmatrix->insertGlobalValues(A->getRowMap()->getGlobalElement(row), indoutP.view(0,indoutP.size()), valout.view(0,valout.size()));
      permQTmatrix->insertGlobalValues (A->getRowMap()->getGlobalElement(row), indoutQ.view(0,indoutQ.size()), valout.view(0,valout.size()));
    }

    permPTmatrix->fillComplete();
    permQTmatrix->fillComplete();

    Teuchos::RCP<Matrix> permPmatrix = Utils2::Transpose(permPTmatrix,true);

    for(size_t row=0; row<permPTmatrix->getNodeNumRows(); row++) {
      if(permPTmatrix->getNumEntriesInLocalRow(row) != 1)
        GetOStream(Warnings0,0) <<"#entries in row " << row << " of permPTmatrix is " << permPTmatrix->getNumEntriesInLocalRow(row) << std::endl;
      if(permPmatrix->getNumEntriesInLocalRow(row) != 1)
        GetOStream(Warnings0,0) <<"#entries in row " << row << " of permPmatrix is " << permPmatrix->getNumEntriesInLocalRow(row) << std::endl;
      if(permQTmatrix->getNumEntriesInLocalRow(row) != 1)
        GetOStream(Warnings0,0) <<"#entries in row " << row << " of permQmatrix is " << permQTmatrix->getNumEntriesInLocalRow(row) << std::endl;
    }

    // build permP * A * permQT
    Teuchos::RCP<Matrix> ApermQt = Utils::Multiply(*A, false, *permQTmatrix, false);
    Teuchos::RCP<Matrix> permPApermQt = Utils::Multiply(*permPmatrix, false, *ApermQt, false);

    /*
    MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("A.mat", *A);
    MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("permP.mat", *permPmatrix);
    MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("permQt.mat", *permQTmatrix);
    MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("permPApermQt.mat", *permPApermQt);
    */
    // build scaling matrix
    Teuchos::RCP<Vector> diagVec = VectorFactory::Build(permPApermQt->getRowMap(),true);
    Teuchos::RCP<Vector> invDiagVec = VectorFactory::Build(permPApermQt->getRowMap(),true);
    Teuchos::ArrayRCP< const Scalar > diagVecData = diagVec->getData(0);
    Teuchos::ArrayRCP< Scalar > invDiagVecData = invDiagVec->getDataNonConst(0);

    permPApermQt->getLocalDiagCopy(*diagVec);
    for(size_t i = 0; i<diagVec->getMap()->getNodeNumElements(); ++i) {
      if(diagVecData[i] != 0.0)
        invDiagVecData[i] = 1/diagVecData[i];
      else {
        invDiagVecData[i] = 1.0;
        GetOStream(Statistics0,0) << "MueLu::LocalPermutationStrategy: found zero on diagonal in row " << i << std::endl;
      }
    }

    Teuchos::RCP<CrsMatrixWrap> diagScalingOp = Teuchos::rcp(new CrsMatrixWrap(permPApermQt->getRowMap(),1,Xpetra::StaticProfile));

    for(size_t row=0; row<A->getNodeNumRows(); row++) {
      Teuchos::ArrayRCP<GlobalOrdinal> indout(1,permPApermQt->getRowMap()->getGlobalElement(row)); // column idx for Perm^T
      Teuchos::ArrayRCP<Scalar> valout(1,invDiagVecData[row]);
      diagScalingOp->insertGlobalValues(A->getRowMap()->getGlobalElement(row), indout.view(0,indout.size()), valout.view(0,valout.size()));
    }
    diagScalingOp->fillComplete();

    Teuchos::RCP<Matrix> scaledA = Utils::Multiply(*diagScalingOp, false, *permPApermQt, false);
    currentLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(scaledA), genFactory/*this*/);

    currentLevel.Set("permA", Teuchos::rcp_dynamic_cast<Matrix>(permPApermQt), genFactory/*this*/);  // TODO careful with this!!!
    currentLevel.Set("permP", Teuchos::rcp_dynamic_cast<Matrix>(permPmatrix), genFactory/*this*/);
    currentLevel.Set("permQT", Teuchos::rcp_dynamic_cast<Matrix>(permQTmatrix), genFactory/*this*/);
    currentLevel.Set("permScaling", Teuchos::rcp_dynamic_cast<Matrix>(diagScalingOp), genFactory/*this*/);

    //// count row permutations
    // count zeros on diagonal in P -> number of row permutations
    Teuchos::RCP<Vector> diagPVec = VectorFactory::Build(permPmatrix->getRowMap(),true);
    permPmatrix->getLocalDiagCopy(*diagPVec);
    Teuchos::ArrayRCP< const Scalar > diagPVecData = diagPVec->getData(0);
    LocalOrdinal lNumRowPermutations = 0;
    GlobalOrdinal gNumRowPermutations = 0;
    for(size_t i = 0; i<diagPVec->getMap()->getNodeNumElements(); ++i) {
      if(diagPVecData[i] == 0.0) {
        lNumRowPermutations++;
      }
    }

    // sum up all entries in multipleColRequests over all processors
    sumAll(diagPVec->getMap()->getComm(), (LocalOrdinal)lNumRowPermutations, gNumRowPermutations);

    //// count column permutations
    // count zeros on diagonal in Q^T -> number of column permutations
    Teuchos::RCP<Vector> diagQTVec = VectorFactory::Build(permQTmatrix->getRowMap(),true);
    permQTmatrix->getLocalDiagCopy(*diagQTVec);
    Teuchos::ArrayRCP< const Scalar > diagQTVecData = diagQTVec->getData(0);
    LocalOrdinal lNumColPermutations = 0;
    GlobalOrdinal gNumColPermutations = 0;
    for(size_t i = 0; i<diagQTVec->getMap()->getNodeNumElements(); ++i) {
      if(diagQTVecData[i] == 0.0) {
        lNumColPermutations++;
      }
    }

    // sum up all entries in multipleColRequests over all processors
    sumAll(diagQTVec->getMap()->getComm(), (LocalOrdinal)lNumColPermutations, gNumColPermutations);

    currentLevel.Set("#RowPermutations", gNumRowPermutations, genFactory/*this*/);
    currentLevel.Set("#ColPermutations", gNumColPermutations, genFactory/*this*/);
    currentLevel.Set("#WideRangeRowPermutations", 0, genFactory/*this*/);
    currentLevel.Set("#WideRangeColPermutations", 0, genFactory/*this*/);

    GetOStream(Statistics0, 0) << "#Row    permutations/max possible permutations: " << gNumRowPermutations << "/" << diagPVec->getMap()->getGlobalNumElements() << std::endl;
    GetOStream(Statistics0, 0) << "#Column permutations/max possible permutations: " << gNumColPermutations << "/" << diagQTVec->getMap()->getGlobalNumElements() << std::endl;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GlobalOrdinal LocalPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::getGlobalDofId(const Teuchos::RCP<Matrix> & A, LocalOrdinal localNodeId, LocalOrdinal localDof) const {
    size_t nDofsPerNode = 1;
    if (A->IsView("stridedMaps")) {
      Teuchos::RCP<const Map> permRowMapStrided = A->getRowMap("stridedMaps");
      nDofsPerNode = Teuchos::rcp_dynamic_cast<const StridedMap>(permRowMapStrided)->getFixedBlockSize();
    }

    LocalOrdinal localDofId = localNodeId * nDofsPerNode + localDof;

    return A->getRowMap()->getGlobalElement(localDofId);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GlobalOrdinal LocalPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::globalDofId2globalNodeId( const Teuchos::RCP<Matrix> & A, GlobalOrdinal grid ) const {
    size_t nDofsPerNode = 1;
    if (A->IsView("stridedMaps")) {
      Teuchos::RCP<const Map> permRowMapStrided = A->getRowMap("stridedMaps");
      nDofsPerNode = Teuchos::rcp_dynamic_cast<const StridedMap>(permRowMapStrided)->getFixedBlockSize();
    }

    return (GlobalOrdinal) grid / (GlobalOrdinal)nDofsPerNode; // TODO what about nOffset???
  }

} // namespace MueLu


#endif /* MUELU_LOCALPERMUTATIONSTRATEGY_DEF_HPP_ */
