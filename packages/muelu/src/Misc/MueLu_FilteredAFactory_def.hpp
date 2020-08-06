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
#ifndef MUELU_FILTEREDAFACTORY_DEF_HPP
#define MUELU_FILTEREDAFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_FilteredAFactory_decl.hpp"

#include "MueLu_FactoryManager.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class T>
  void sort_and_unique(T & array) {
    std::sort(array.begin(),array.end());
    std::unique(array.begin(),array.end());
  }
      


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("filtered matrix: use lumping");
    SET_VALID_ENTRY("filtered matrix: reuse graph");
    SET_VALID_ENTRY("filtered matrix: reuse eigenvalue");
    SET_VALID_ENTRY("filtered matrix: use root stencil");
    SET_VALID_ENTRY("filtered matrix: Dirichlet threshold");
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used for filtering");
    validParamList->set< RCP<const FactoryBase> >("Graph",          Teuchos::null, "Generating factory for coalesced filtered graph");
    validParamList->set< RCP<const FactoryBase> >("Filtering",      Teuchos::null, "Generating factory for filtering boolean");

    
    // Only need these for the "use root stencil" option
    validParamList->set< RCP<const FactoryBase> >("Aggregates",         Teuchos::null, "Generating factory of the aggregates");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "Filtering");
    Input(currentLevel, "Graph");
    const ParameterList& pL = GetParameterList();
    if(pL.isParameter("filtered matrix: use root stencil") && pL.get<bool>("filtered matrix: use root stencil") == true){
      Input(currentLevel, "Aggregates");
      Input(currentLevel, "UnAmalgamationInfo");
    }    
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Matrix filtering", currentLevel);

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
    if (Get<bool>(currentLevel, "Filtering") == false) {
      GetOStream(Runtime0) << "Filtered matrix is not being constructed as no filtering is being done" << std::endl;
      Set(currentLevel, "A", A);
      return;
    }

    const ParameterList& pL = GetParameterList();
    bool lumping = pL.get<bool>("filtered matrix: use lumping");
    if (lumping)
      GetOStream(Runtime0) << "Lumping dropped entries" << std::endl;
    bool use_root_stencil = lumping && pL.get<bool>("filtered matrix: use root stencil");
    if (use_root_stencil)
      GetOStream(Runtime0) << "Using root stencil for dropping" << std::endl;
    double dirichlet_threshold = pL.get<double>("filtered matrix: Dirichlet threshold");    
    if(dirichlet_threshold >= 0.0)
      GetOStream(Runtime0) << "Filtering Dirichlet threshold of "<<dirichlet_threshold << std::endl;

    if(pL.get<bool>("filtered matrix: reuse graph"))
      GetOStream(Runtime0) << "Reusing graph"<<std::endl;
    else
      GetOStream(Runtime0) << "Generating new graph"<<std::endl;
      

    RCP<GraphBase> G = Get< RCP<GraphBase> >(currentLevel, "Graph");

    RCP<ParameterList> fillCompleteParams(new ParameterList);
    fillCompleteParams->set("No Nonlocal Changes", true);

    RCP<Matrix> filteredA;
    if (pL.get<bool>("filtered matrix: reuse graph")) {
      filteredA = MatrixFactory::Build(A->getCrsGraph());
      filteredA->resumeFill();

      BuildReuse(*A, *G, lumping, dirichlet_threshold,*filteredA);

      filteredA->fillComplete(fillCompleteParams);

    } else {
      filteredA = MatrixFactory::Build(A->getRowMap(), A->getColMap(), A->getNodeMaxNumRowEntries());

      if(use_root_stencil)
	BuildNewUsingRootStencil(*A, *G, dirichlet_threshold, currentLevel,*filteredA);
      else
	BuildNew(*A, *G, lumping, dirichlet_threshold,*filteredA);

      filteredA->fillComplete(A->getDomainMap(), A->getRangeMap(), fillCompleteParams);
    }

    filteredA->SetFixedBlockSize(A->GetFixedBlockSize());

    if (pL.get<bool>("filtered matrix: reuse eigenvalue")) {
      // Reuse max eigenvalue from A
      // It is unclear what eigenvalue is the best for the smoothing, but we already may have
      // the D^{-1}A estimate in A, may as well use it.
      // NOTE: ML does that too
      filteredA->SetMaxEigenvalueEstimate(A->GetMaxEigenvalueEstimate());
    }

    Set(currentLevel, "A", filteredA);
  }

// Epetra's API allows direct access to row array.
// Tpetra's API does not, providing only ArrayView<const .>
// But in most situations we are currently interested in, it is safe to assume
// that the view is to the actual data. So this macro directs the code to do
// const_cast, and modify the entries directly. This allows us to avoid
// replaceLocalValues() call which is quite expensive due to all the searches.
#define ASSUME_DIRECT_ACCESS_TO_ROW

  // Both Epetra and Tpetra matrix-matrix multiply use the following trick:
  // if an entry of the left matrix is zero, it does not compute or store the
  // zero value.
  //
  // This trick allows us to bypass constructing a new matrix. Instead, we
  // make a deep copy of the original one, and fill it in with zeros, which
  // are ignored during the prolongator smoothing.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildReuse(const Matrix& A, const GraphBase& G, const bool lumping, double dirichletThresh, Matrix& filteredA) const {
    using TST = typename Teuchos::ScalarTraits<SC>;
    SC zero = TST::zero();


    size_t blkSize = A.GetFixedBlockSize();

    ArrayView<const LO> inds;
    ArrayView<const SC> valsA;
#ifdef ASSUME_DIRECT_ACCESS_TO_ROW
    ArrayView<SC>       vals;
#else
    Array<SC>           vals;
#endif

    Array<char> filter( A.getColMap()->getNodeNumElements(), 0);

    size_t numGRows = G.GetNodeNumVertices();
    for (size_t i = 0; i < numGRows; i++) {
      // Set up filtering array
      ArrayView<const LO> indsG = G.getNeighborVertices(i);
      for (size_t j = 0; j < as<size_t>(indsG.size()); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter[indsG[j]*blkSize+k] = 1;

      for (size_t k = 0; k < blkSize; k++) {
        LO row = i*blkSize + k;

        A.getLocalRowView(row, inds, valsA);

        size_t nnz = inds.size();
        if (nnz == 0)
          continue;

#ifdef ASSUME_DIRECT_ACCESS_TO_ROW
        // Transform ArrayView<const SC> into ArrayView<SC>
        ArrayView<const SC> vals1;
        filteredA.getLocalRowView(row, inds, vals1);
        vals = ArrayView<SC>(const_cast<SC*>(vals1.getRawPtr()), nnz);

        memcpy(vals.getRawPtr(), valsA.getRawPtr(), nnz*sizeof(SC));
#else
        vals = Array<SC>(valsA);
#endif

	SC ZERO = Teuchos::ScalarTraits<SC>::zero();
	SC ONE = Teuchos::ScalarTraits<SC>::one();
	SC A_rowsum = ZERO, F_rowsum = ZERO;
	for(LO l = 0; l < (LO)inds.size(); l++) 
	  A_rowsum += valsA[l];

        if (lumping == false) {
          for (size_t j = 0; j < nnz; j++)
            if (!filter[inds[j]])
              vals[j] = zero;

        } else {
          LO diagIndex = -1;
          SC diagExtra = zero;

          for (size_t j = 0; j < nnz; j++) {
            if (filter[inds[j]]) {
              if (inds[j] == row) {
                // Remember diagonal position
                diagIndex = j;
              }
              continue;
            }

            diagExtra += vals[j];

            vals[j] = zero;
          }

          // Lump dropped entries
          // NOTE
          //  * Does it make sense to lump for elasticity?
          //  * Is it different for diffusion and elasticity?
	  SC diagA = ZERO;
	  if (diagIndex != -1) {
	    diagA = vals[diagIndex];	
	    vals[diagIndex] += diagExtra;
	    if(dirichletThresh >= 0.0 && TST::real(vals[diagIndex]) <= dirichletThresh) {
	      printf("WARNING: row %d diag(Afiltered) = %8.2e diag(A)=%8.2e\n",row,vals[diagIndex],diagA);
	      for(LO l = 0; l < (LO)nnz; l++) 
		F_rowsum += vals[l];
	      printf("       : A rowsum = %8.2e F rowsum = %8.2e\n",A_rowsum,F_rowsum);	    	    
	      vals[diagIndex] = TST::one();
	    }
	  }
	  	  
        }

#ifndef ASSUME_DIRECT_ACCESS_TO_ROW
        // Because we used a column map in the construction of the matrix
        // we can just use insertLocalValues here instead of insertGlobalValues
        filteredA.replaceLocalValues(row, inds, vals);
#endif
      }

      // Reset filtering array
      for (size_t j = 0; j < as<size_t> (indsG.size()); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter[indsG[j]*blkSize+k] = 0;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildNew(const Matrix& A, const GraphBase& G, const bool lumping, double dirichletThresh, Matrix& filteredA) const {
    using TST = typename Teuchos::ScalarTraits<SC>;
    SC zero = Teuchos::ScalarTraits<SC>::zero();

    size_t blkSize = A.GetFixedBlockSize();

    ArrayView<const LO> indsA;
    ArrayView<const SC> valsA;
    Array<LO>           inds;
    Array<SC>           vals;

    Array<char> filter(blkSize * G.GetImportMap()->getNodeNumElements(), 0);

    size_t numGRows = G.GetNodeNumVertices();
    for (size_t i = 0; i < numGRows; i++) {
      // Set up filtering array
      ArrayView<const LO> indsG = G.getNeighborVertices(i);
      for (size_t j = 0; j < as<size_t>(indsG.size()); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter[indsG[j]*blkSize+k] = 1;

      for (size_t k = 0; k < blkSize; k++) {
        LO row = i*blkSize + k;

        A.getLocalRowView(row, indsA, valsA);

        size_t nnz = indsA.size();
        if (nnz == 0)
          continue;

        inds.resize(indsA.size());
        vals.resize(valsA.size());

        size_t numInds = 0;
        if (lumping == false) {
          for (size_t j = 0; j < nnz; j++)
            if (filter[indsA[j]]) {
              inds[numInds] = indsA[j];
              vals[numInds] = valsA[j];
              numInds++;
            }

        } else {
          LO diagIndex = -1;
          SC diagExtra = zero;

          for (size_t j = 0; j < nnz; j++) {
            if (filter[indsA[j]]) {
              inds[numInds] = indsA[j];
              vals[numInds] = valsA[j];

              // Remember diagonal position
              if (inds[numInds] == row)
                diagIndex = numInds;

              numInds++;

            } else {
              diagExtra += valsA[j];
            }
          }

          // Lump dropped entries
          // NOTE
          //  * Does it make sense to lump for elasticity?
          //  * Is it different for diffusion and elasticity?
          if (diagIndex != -1) {
            vals[diagIndex] += diagExtra;
	    if(dirichletThresh >= 0.0 && TST::real(vals[diagIndex]) <= dirichletThresh) {
	      //	      SC A_rowsum = ZERO, F_rowsum = ZERO;
	      //	      printf("WARNING: row %d diag(Afiltered) = %8.2e diag(A)=%8.2e\n",row,vals[diagIndex],diagA);
	      //	      for(LO l = 0; l < (LO)nnz; l++) 
	      //		F_rowsum += vals[l];
	      //	      printf("       : A rowsum = %8.2e F rowsum = %8.2e\n",A_rowsum,F_rowsum);	    	    
	      vals[diagIndex] = TST::one();
	    }
	  }
	      	
        }
        inds.resize(numInds);
        vals.resize(numInds);



        // Because we used a column map in the construction of the matrix
        // we can just use insertLocalValues here instead of insertGlobalValues
        filteredA.insertLocalValues(row, inds, vals);
      }

      // Reset filtering array
      for (size_t j = 0; j < as<size_t> (indsG.size()); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter[indsG[j]*blkSize+k] = 0;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildNewUsingRootStencil(const Matrix& A, const GraphBase& G, double dirichletThresh, Level& currentLevel,  Matrix& filteredA) const {
    using TST = typename Teuchos::ScalarTraits<SC>;
    SC ZERO = Teuchos::ScalarTraits<SC>::zero();
    SC ONE = Teuchos::ScalarTraits<SC>::one();
    LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();

    size_t numNodes = G.GetNodeNumVertices();
    size_t blkSize  = A.GetFixedBlockSize();
    size_t numRows  = A.getMap()->getNodeNumElements();
    ArrayView<const LO> indsA;
    ArrayView<const SC> valsA;
    Array<LO>           inds(A.getNodeMaxNumRowEntries());
    Array<SC>           vals(A.getNodeMaxNumRowEntries());

    RCP<Aggregates>            aggregates     = Get< RCP<Aggregates> >           (currentLevel, "Aggregates");
    RCP<AmalgamationInfo>      amalgInfo      = Get< RCP<AmalgamationInfo> >     (currentLevel, "UnAmalgamationInfo");
    LO                          numAggs       = aggregates->GetNumAggregates();
    Teuchos::ArrayRCP<const LO> vertex2AggId  = aggregates->GetVertex2AggId()->getData(0);


    
    // HAX
    RCP<const Map> rowMap = A.getRowMap();
    RCP<const Map> colMap = A.getColMap();
    bool goodMap = MueLu::Utilities<SC,LO,GO,NO>::MapsAreNested(*rowMap, *colMap);
    if(goodMap) printf("CMS: Maps are good\n");
    else printf("CMS: Maps are NOT good\n");


    // HAX
    //    BuildNew(A,G,true,filteredA);
    //    filteredA = A;
#ifdef JUST_USE_A
    std::cout<<"CMS: Hacking the diagonal"<<std::endl;
    for(LO row=0; row<(LO) filteredA.getRowMap()->getNodeNumElements(); row++) {
      A.getLocalRowView(row, indsA, valsA);
      vals.resize(valsA.size());
      for(LO j=0; j<(LO) indsA.size(); j++) {
	vals[j] = ((row==indsA[j]) && (valsA[j] < 1e-6)) ? 1e-6 : valsA[j];
      }
      filteredA.insertLocalValues(row, indsA, vals);
    }
    return;
#endif
    // END HAX


    // Lists of nodes in each aggregate
    struct {
      Array<LO> ptr,nodes;
    } nodesInAgg;
    aggregates->ComputeNodesInAggregate(nodesInAgg.ptr, nodesInAgg.nodes);

    // Unamalgamate

    //    bool goodMap = MueLu::Utilities<SC,LO,GO,NO>::MapsAreNested(*rowMap, *colMap);
    //    TEUCHOS_TEST_FOR_EXCEPTION(goodMap,
    //			       Exceptions::RuntimeError,"MueLu::FilteredAFactory::BuildNewUsingRootStencil: Requires nested row/col maps");


    // Find the biggest aggregate size in *nodes*
    LO maxAggSize=0;
    for(LO i=0; i<numAggs; i++) 
      maxAggSize = std::max(maxAggSize,nodesInAgg.ptr[i+1] - nodesInAgg.ptr[i]);


    // Loop over all the aggregates
    std::vector<LO> goodAggNeighbors(G.getNodeMaxNumRowEntries());
    std::vector<LO> badAggNeighbors(std::min(G.getNodeMaxNumRowEntries()*maxAggSize,numNodes));

    size_t numNewDrops=0;
    size_t numOldDrops=0;
    size_t numFixedDiags=0;

    for(LO i=0; i<numAggs; i++) {
      LO numNodesInAggregate = nodesInAgg.ptr[i+1] - nodesInAgg.ptr[i];
      if(numNodesInAggregate == 0) continue;

      // Find the root *node*
      LO root_node = INVALID;
      for (LO k=nodesInAgg.ptr[i]; k < nodesInAgg.ptr[i+1]; k++) {
	if(aggregates->IsRoot(nodesInAgg.nodes[k])) {
	  root_node = k; break;
	}
      }
      
      TEUCHOS_TEST_FOR_EXCEPTION(root_node == INVALID,
				 Exceptions::RuntimeError,"MueLu::FilteredAFactory::BuildNewUsingRootStencil: Cannot find root node");

      // Find the list of "good" node neighbors (aka nodes which border the root node in the Graph G)      
      ArrayView<const LO> goodNodeNeighbors  = G.getNeighborVertices(root_node);

      // Now find the list of "good" aggregate neighbors (aka the aggregates neighbor the root node in the Graph G)
      goodAggNeighbors.resize(0);
      for(LO k=0; k<(LO) goodNodeNeighbors.size(); k++) {
	goodAggNeighbors.push_back(vertex2AggId[goodNodeNeighbors[k]]);
      }
      sort_and_unique(goodAggNeighbors);
      
      // Now we get the list of "bad" aggregate neighbors (aka aggregates which border the 
      // root node in the original matrix A, which are not goodNodeNeighbors).  Since we
      // don't have an amalgamated version of the original matrix, we use the matrix directly
      badAggNeighbors.resize(0);
      for(LO j = 0; j < (LO)blkSize; j++) {
	LO row = amalgInfo->ComputeLocalDOF(root_node,j);
	A.getLocalRowView(row, indsA, valsA);
	for(LO k=0; k<(LO)indsA.size(); k++) {
	  if(indsA[k] < (LO)numRows) {
	    int node = amalgInfo->ComputeLocalNode(indsA[k]);
	    int agg = vertex2AggId[node];
	    if(!std::binary_search(goodAggNeighbors.begin(),goodAggNeighbors.end(),agg))
	      badAggNeighbors.push_back(agg);
	  }
	}	
      }
      sort_and_unique(badAggNeighbors);
  
      
      // Now lets fill the rows in this aggregate and figure out the diagonal lumping
      // We loop over each node in the aggregate and then over the neighbors of that node

      for(LO k=nodesInAgg.ptr[i]; k<nodesInAgg.ptr[i+1]; k++) {
	LO node = nodesInAgg.nodes[k];
	//	ArrayView<const LO> indsG = G.getNeighborVertices(nodesInAgg.nodes[k]);
	//	for (LO j = 0; j < (LO)indsG.size(); j++) {
	for (LO m = 0; m < (LO)blkSize; m++) {
	  LO row = amalgInfo->ComputeLocalDOF(node,m);
	  //	  if(row ==7651)
	    //	    printf(" agg = %d blkSize = %d l = %d m = %d\n",
	  A.getLocalRowView(row, indsA, valsA);
	  
	  LO diagIndex = INVALID;
	  SC diagExtra = ZERO;
	  LO numInds = 0;
	  inds.resize(indsA.size());
	  vals.resize(valsA.size());
	  
	  for(LO l = 0; l < (LO)indsA.size(); l++) {
	    bool is_good = true; // Assume off-processor entries are good. FIXME?
	    if (indsA[l] == row) {
	      diagIndex = l;
	      inds[numInds] = indsA[l];
	      vals[numInds] = valsA[l];
	      numInds++;
	      continue;
	    }
	    
	    if(indsA[l] < (LO)numRows)  {
	      int Anode = amalgInfo->ComputeLocalNode(indsA[l]);
	      int agg = vertex2AggId[Anode];
	      if(std::binary_search(badAggNeighbors.begin(),badAggNeighbors.end(),agg))
		is_good = false;
	    }
	    
	    if(is_good){
	      inds[numInds] = indsA[l];
	      vals[numInds] = valsA[l];
	      numInds++;
	    }
	    else {
	      if(valsA[l] == ZERO) numOldDrops++;
	      else numNewDrops++;
	      diagExtra += valsA[l];
	      inds[numInds] = indsA[l];
	      vals[numInds]=0;
	      numInds++;
	    }
	  } //end for l "indsA.size()" loop
	    
	  if (diagIndex != INVALID) {
	    SC diagA = valsA[diagIndex];
	    vals[diagIndex] += diagExtra;
	    SC A_rowsum=ZERO, A_absrowsum = ZERO, F_rowsum = ZERO;
	    
	    if(dirichletThresh >= 0.0 && TST::real(vals[diagIndex]) <= dirichletThresh) {
	      
	      printf("WARNING: row %d (diagIndex=%d) diag(Afiltered) = %8.2e diag(A)=%8.2e\n",row,diagIndex,vals[diagIndex],diagA);
	      for(LO l = 0; l < (LO)indsA.size(); l++) {		  
		A_rowsum += valsA[l];
		A_absrowsum+=std::abs(valsA[l]);
	      }
	      for(LO l = 0; l < (LO)numInds; l++) 
		F_rowsum += vals[l];
	      printf("       : A rowsum = %8.2e |A| rowsum = %8.2e rowsum = %8.2e\n",A_rowsum,A_absrowsum,F_rowsum);
	      if(row==7444) {
		printf("       vals =");
		for(LO l = 0; l < (LO)indsA.size(); l++)
		  printf("%d(%8.2e) ",indsA[l],valsA[l]);
		printf("\n");
	      }
	      
	      
	      vals[diagIndex] = TST::one();
	      numFixedDiags++;
	    }
	  }
	  else {
	    GetOStream(Runtime0)<<"WARNING: Row "<<row<<" has no diagonal "<<std::endl;
	  }
	  
	  inds.resize(numInds);
	  vals.resize(numInds);
	  filteredA.insertLocalValues(row, inds, vals);
	  
	  if(numInds == 0) 
	    GetOStream(Runtime0)<<"WARNING: Row "<<row<<" has no entries in Afiltered. A has "<<indsA.size()<<" entries"<<std::endl;
	}//end m "blkSize" loop
      }// end k loop over number of nodes in this agg
    }//end i loop over numAggs

    size_t g_newDrops = 0, g_oldDrops = 0, g_fixedDiags=0;
    
    MueLu_sumAll(A.getRowMap()->getComm(), numNewDrops, g_newDrops);
    MueLu_sumAll(A.getRowMap()->getComm(), numOldDrops, g_oldDrops);
    MueLu_sumAll(A.getRowMap()->getComm(), numFixedDiags, g_fixedDiags);
    GetOStream(Runtime0)<< "Filtering out "<<g_newDrops<<" edges, in addition to the "<<g_oldDrops<<" edges dropped earlier"<<std::endl;
    GetOStream(Runtime0)<< "Fixing "<< g_fixedDiags<<" zero diagonal values" <<std::endl;
  }


} //namespace MueLu

#endif // MUELU_FILTEREDAFACTORY_DEF_HPP
