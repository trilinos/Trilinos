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
  void
  LocalPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  BuildPermutation (const Teuchos::RCP<Matrix> & A,
                    const Teuchos::RCP<const Map> permRowMap,
                    Level & currentLevel,
                    const FactoryBase* genFactory) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;
    typedef Teuchos::ScalarTraits<MT> STM;

#ifndef HAVE_MUELU_INST_COMPLEX_INT_INT // TODO remove this -> check scalar = std::complex
    size_t nDofsPerNode = 1;
    if (A->IsView("stridedMaps")) {
      RCP<const Map> permRowMapStrided = A->getRowMap("stridedMaps");
      nDofsPerNode = rcp_dynamic_cast<const StridedMap>(permRowMapStrided)->getFixedBlockSize();
    }

    //////////////////
    std::vector<std::pair<GO, GO> > RowColPairs;

    // mfh 14 Mar 2013: Amortize the overhead of accessing the row and
    // column Maps by hoisting access to them outside the loops.
    // There are a lot of places in this method where you could do
    // that for loop bounds; see the loop below, for example.
    const Map& rowMap = * (A->getRowMap ());
    const Map& colMap = * (A->getColMap ());
    const Map& domMap = * (A->getDomainMap ());

    // loop over local nodes
    // TODO what about nOffset?
    for (LO node = 0; node < as<LO> (A->getRowMap()->getNodeNumElements()/nDofsPerNode); ++node) {
      // FIXME (mfh 14 Mar 2013) SerialDenseMatrix only works (where
      // "works" means "with the BLAS and LAPACK") for OrdinalType=int
      // (its first template parameter).
      Teuchos::SerialDenseMatrix<LO, Scalar> subBlockMatrix (nDofsPerNode, nDofsPerNode, true);

      std::vector<GO> growIds (nDofsPerNode);

      for (LO lrdof = 0; lrdof < as<GO> (nDofsPerNode); ++lrdof) { // TODO more complicated for variable dofs per node
        const GO grow = getGlobalDofId (A, node, lrdof);
        growIds[lrdof] = grow;

        //if(permRowMap->isNodeGlobalElement(grow) == true) continue;

        // extract local row information from matrix
        ArrayView<const LO> indices;
        ArrayView<const Scalar> vals;
        A->getLocalRowView (rowMap.getLocalElement (grow), indices, vals);

        // Find column entry with max absolute value
        MT maxVal = STM::zero ();
        for (size_t j = 0; j < as<size_t> (indices.size ()); ++j) {
          const MT curMag = STS::magnitude (vals[j]);
          if (curMag > maxVal) {
            maxVal = curMag;
          }
        }

        const GO grnodeid = globalDofId2globalNodeId (A, grow);

        for (size_t j = 0; j < as<size_t> (indices.size ()); ++j) {
          const GO gcol = colMap.getGlobalElement (indices[j]);
          const GO gcnodeid = globalDofId2globalNodeId (A, gcol); // -> global node id
          if (grnodeid == gcnodeid) {
            if (maxVal != STM::zero ()) {
              subBlockMatrix (lrdof, gcol % nDofsPerNode) = vals[j] / maxVal;
            } else {
              subBlockMatrix (lrdof, gcol % nDofsPerNode) = vals[j]; // there is a problem
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
        Scalar value = STS::one ();
        for (size_t len = 0; len < s.length (); ++len) {
          const int col = static_cast<int> (s[len] - '0');
          value = value * subBlockMatrix (len, col);
        }
        performance_vector[t] = value;
      }

      // find permutation with maximum performance value
      MT maxVal = STM::one ();
      size_t maxPerformancePermutationIdx = 0;
      for (size_t j = 0; j < Teuchos::as<size_t>(performance_vector.size()); j++) {
        const MT curMag = STS::magnitude (performance_vector[j]);
        if (curMag > maxVal) {
          maxVal = curMag;
          maxPerformancePermutationIdx = j;
        }
      }

      // the best permutation is
      //std::cout << "The best permutation is " << result_perms[maxPerformancePermutationIdx] << std::endl;

      std::string bestPerformancePermutation = result_perms[maxPerformancePermutationIdx] ;
      for(size_t t = 0; t<nDofsPerNode; t++) {
        const int col = static_cast<int> (bestPerformancePermutation[t] - '0');
        RowColPairs.push_back (std::make_pair (growIds[t], growIds[col]));
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

    Pperm->putScalar (STS::zero ());
    Qperm->putScalar (STS::zero ());

    Teuchos::ArrayRCP<Scalar> PpermData = Pperm->getDataNonConst(0);
    Teuchos::ArrayRCP<Scalar> QpermData = Qperm->getDataNonConst(0);

    typename std::vector<std::pair<GO, GO> >::iterator p = RowColPairs.begin();
    while (p != RowColPairs.end()) {
      const GO ik = p->first;
      const GO jk = p->second;

      const LO lik = rowMap.getLocalElement (ik);
      const LO ljk = domMap.getLocalElement (jk);

      Pperm->replaceLocalValue (lik, ik);
      Qperm->replaceLocalValue (ljk, ik);

      p = RowColPairs.erase (p);
    }

    if (RowColPairs.size() > 0) {
      GetOStream(Warnings0, 0) << "MueLu::LocalPermutationStrategy: "
        "There are Row/col pairs left!" << std::endl;
    }

    // Qperm should be fine
    // build matrices

    // create new empty Matrix
    RCP<CrsMatrixWrap> permPTmatrix =
      rcp (new CrsMatrixWrap (A->getRowMap (), 1, Xpetra::StaticProfile));
    RCP<CrsMatrixWrap> permQTmatrix =
      rcp (new CrsMatrixWrap (A->getRowMap (), 1, Xpetra::StaticProfile));

    for (size_t row = 0; row < A->getNodeNumRows (); ++row) {
      ArrayRCP<GO> indoutP (1, Teuchos::as<GO> (PpermData[row])); // column idx for Perm^T
      ArrayRCP<GO> indoutQ (1, Teuchos::as<GO> (QpermData[row])); // column idx for Qperm
      ArrayRCP<Scalar> valout (1, STS::one ());
      permPTmatrix->insertGlobalValues (rowMap.getGlobalElement (row),
                                        indoutP.view (0, indoutP.size ()),
                                        valout.view (0, valout.size ()));
      permQTmatrix->insertGlobalValues (rowMap.getGlobalElement (row),
                                        indoutQ.view (0,indoutQ.size ()),
                                        valout.view (0,valout.size ()));
    }

    // FIXME (mfh 14 Mar 2013) Do we need to be passing in the domain
    // and range Maps here?
    permPTmatrix->fillComplete();
    permQTmatrix->fillComplete();

    RCP<Matrix> permPmatrix = Utils2::Transpose (permPTmatrix, true);

    for (size_t row=0; row < permPTmatrix->getNodeNumRows (); ++row) {
      if (permPTmatrix->getNumEntriesInLocalRow(row) != 1)
        GetOStream(Warnings0,0) <<"#entries in row " << row << " of permPTmatrix is " << permPTmatrix->getNumEntriesInLocalRow(row) << std::endl;
      if (permPmatrix->getNumEntriesInLocalRow(row) != 1)
        GetOStream(Warnings0,0) <<"#entries in row " << row << " of permPmatrix is " << permPmatrix->getNumEntriesInLocalRow(row) << std::endl;
      if (permQTmatrix->getNumEntriesInLocalRow(row) != 1)
        GetOStream(Warnings0,0) <<"#entries in row " << row << " of permQmatrix is " << permQTmatrix->getNumEntriesInLocalRow(row) << std::endl;
    }

    // build permP * A * permQT
    RCP<Matrix> ApermQt = Utils::Multiply(*A, false, *permQTmatrix, false);
    RCP<Matrix> permPApermQt = Utils::Multiply(*permPmatrix, false, *ApermQt, false);

    /*
    MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("A.mat", *A);
    MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("permP.mat", *permPmatrix);
    MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("permQt.mat", *permQTmatrix);
    MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("permPApermQt.mat", *permPApermQt);
    */
    // build scaling matrix
    RCP<Vector> diagVec = VectorFactory::Build(permPApermQt->getRowMap(),true);
    RCP<Vector> invDiagVec = VectorFactory::Build(permPApermQt->getRowMap(),true);
    ArrayRCP<const Scalar> diagVecData = diagVec->getData(0);
    ArrayRCP<Scalar> invDiagVecData = invDiagVec->getDataNonConst(0);

    permPApermQt->getLocalDiagCopy(*diagVec);
    for (size_t i = 0; i<diagVec->getMap()->getNodeNumElements(); ++i) {
      if(diagVecData[i] != STS::zero ())
        invDiagVecData[i] = STS::one () / diagVecData[i];
      else {
        invDiagVecData[i] = STS::one ();
        GetOStream(Statistics0,0) << "MueLu::LocalPermutationStrategy: found zero on diagonal in row " << i << std::endl;
      }
    }

    RCP<CrsMatrixWrap> diagScalingOp = rcp (new CrsMatrixWrap (permPApermQt->getRowMap(),1,Xpetra::StaticProfile));

    for (size_t row=0; row<A->getNodeNumRows(); ++row) {
      ArrayRCP<GO> indout(1,permPApermQt->getRowMap()->getGlobalElement(row)); // column idx for Perm^T
      ArrayRCP<Scalar> valout(1,invDiagVecData[row]);
      diagScalingOp->insertGlobalValues (rowMap.getGlobalElement (row),
                                         indout.view(0,indout.size()),
                                         valout.view(0,valout.size()));
    }

    // FIXME (mfh 14 Mar 2013) Do we need to be passing in the domain
    // and range Maps here?
    diagScalingOp->fillComplete();

    RCP<Matrix> scaledA = Utils::Multiply(*diagScalingOp, false, *permPApermQt, false);
    currentLevel.Set("A", rcp_dynamic_cast<Matrix>(scaledA), genFactory/*this*/);

    currentLevel.Set("permA", rcp_dynamic_cast<Matrix>(permPApermQt), genFactory/*this*/);  // TODO careful with this!!!
    currentLevel.Set("permP", rcp_dynamic_cast<Matrix>(permPmatrix), genFactory/*this*/);
    currentLevel.Set("permQT", rcp_dynamic_cast<Matrix>(permQTmatrix), genFactory/*this*/);
    currentLevel.Set("permScaling", rcp_dynamic_cast<Matrix>(diagScalingOp), genFactory/*this*/);

    //// count row permutations
    // count zeros on diagonal in P -> number of row permutations
    RCP<Vector> diagPVec = VectorFactory::Build(permPmatrix->getRowMap(),true);
    permPmatrix->getLocalDiagCopy(*diagPVec);
    ArrayRCP<const Scalar> diagPVecData = diagPVec->getData(0);
    LO lNumRowPermutations = 0;
    GO gNumRowPermutations = 0;
    for (size_t i = 0; i < diagPVec->getMap()->getNodeNumElements(); ++i) {
      if (diagPVecData[i] == STS::zero ()) {
        lNumRowPermutations++;
      }
    }

    // sum up all entries in multipleColRequests over all processors
    sumAll (diagPVec->getMap()->getComm(), as<LO> (lNumRowPermutations),
            gNumRowPermutations);

    //// count column permutations
    // count zeros on diagonal in Q^T -> number of column permutations
    RCP<Vector> diagQTVec = VectorFactory::Build(permQTmatrix->getRowMap(),true);
    permQTmatrix->getLocalDiagCopy(*diagQTVec);
    ArrayRCP<const Scalar> diagQTVecData = diagQTVec->getData(0);
    LO lNumColPermutations = 0;
    GO gNumColPermutations = 0;
    for(size_t i = 0; i<diagQTVec->getMap()->getNodeNumElements(); ++i) {
      if(diagQTVecData[i] == STS::zero ()) {
        lNumColPermutations++;
      }
    }

    // sum up all entries in multipleColRequests over all processors
    sumAll (diagQTVec->getMap()->getComm(), as<LO> (lNumColPermutations),
            gNumColPermutations);

    currentLevel.Set("#RowPermutations", gNumRowPermutations, genFactory/*this*/);
    currentLevel.Set("#ColPermutations", gNumColPermutations, genFactory/*this*/);
    currentLevel.Set("#WideRangeRowPermutations", 0, genFactory/*this*/);
    currentLevel.Set("#WideRangeColPermutations", 0, genFactory/*this*/);

    GetOStream(Statistics0, 0) << "#Row    permutations/max possible permutations: " << gNumRowPermutations << "/" << diagPVec->getMap()->getGlobalNumElements() << std::endl;
    GetOStream(Statistics0, 0) << "#Column permutations/max possible permutations: " << gNumColPermutations << "/" << diagQTVec->getMap()->getGlobalNumElements() << std::endl;

#else
#warning PermutationFactory not compiling/working for Scalar==complex.
#endif // #ifndef HAVE_MUELU_INST_COMPLEX_INT_INT
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
