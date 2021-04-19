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
#ifndef MUELU_CLASSICALPFACTORY_DEF_HPP
#define MUELU_CLASSICALPFACTORY_DEF_HPP

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_IO.hpp>

#include <Teuchos_OrdinalTraits.hpp>

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_ClassicalPFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_AmalgamationInfo.hpp"

#include "KokkosGraph_Distance1ColorHandle.hpp"
#include "KokkosGraph_Distance1Color.hpp"


//#define CMS_DEBUG

namespace { 

template<class SC>
int Sign(SC val) {
  using STS = typename Teuchos::ScalarTraits<SC>;
  typename STS::magnitudeType SC_ZERO = STS::zero();
  if(STS::real(val) > SC_ZERO) return 1;
  else if(STS::real(val) < SC_ZERO) return -1;
  else return 0;
}

}// anonymous namepsace

namespace MueLu {




  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: deterministic");
    SET_VALID_ENTRY("aggregation: coloring algorithm");
    SET_VALID_ENTRY("aggregation: classical scheme");
    {
      typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
      validParamList->getEntry("aggregation: classical scheme").setValidator(rcp(new validatorType(Teuchos::tuple<std::string>("direct","ext+i"), "aggregation: classical scheme")));
                                                                        
    }

#undef SET_VALID_ENTRY
    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
    validParamList->set< RCP<const FactoryBase> >("Graph",       null, "Generating factory of the graph");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");

    //    validParamList->set< RCP<const FactoryBase> >("Nullspace",      Teuchos::null, "Generating factory of the nullspace");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {
    Input(fineLevel, "A");
    Input(fineLevel, "Graph");
    Input(fineLevel, "DofsPerNode");    
    Input(fineLevel, "UnAmalgamationInfo");

    //    Input(fineLevel, "Nullspace");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);
    using STS = Teuchos::ScalarTraits<SC>;

    // We start by assuming that someone did a reasonable strength of connection
    // algorithm before we start to get our Graph, DofsPerNode and UnAmalgamationInfo

    // We begin by getting a MIS (from a graph coloring) and then at that point we need
    // to start generating entries for the prolongator.   
    RCP<Matrix>      A              = Get< RCP<Matrix> >      (fineLevel, "A");
    RCP<const GraphBase> graph      = Get< RCP<GraphBase> >(fineLevel, "Graph");
    LO nDofsPerNode                 = Get<LO>(fineLevel, "DofsPerNode");
    RCP<AmalgamationInfo> amalgInfo = Get< RCP<AmalgamationInfo> >     (fineLevel, "UnAmalgamationInfo");
    //    RCP<MultiVector> fineNullspace = Get< RCP<MultiVector> > (fineLevel, "Nullspace");
    RCP<Matrix> P;
    SC SC_ZERO = STS::zero();
    LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
    const ParameterList& pL = GetParameterList();

    // FIXME: This guy doesn't work right now for NumPDEs != 1
    TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() != 1, Exceptions::RuntimeError,"ClassicalPFactory: Multiple PDEs per node not supported yet");

    // FIXME: This does not work in parallel yet
    TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->getComm()->getSize() !=  1,Exceptions::RuntimeError,"ClassicalPFactory: MPI Ranks > 1 not supported yet");


    /* ============================================================= */
    /* Phase 1 : Compute an initial MIS                              */
    /* ============================================================= */
    colors_view_type myColors_d;
    LO numColors=0;
    // FIXME:  This is not going to respect coloring at or near processor
    // boundaries, so that'll need to get cleaned up lated
    {
      SubFactoryMonitor sfm(*this,"GraphColoring",coarseLevel);
      DoGraphColoring(*graph,myColors_d,numColors);
    }
    auto myColors_h = Kokkos::create_mirror_view(myColors_d);
    Kokkos::deep_copy(myColors_h ,myColors_d);

    // FIXME: This coloring will either need to be done MPI parallel, or
    // there needs to be a cleanup phase to fix mistakes

    /* ============================================================= */
    /* Phase 2 : Mark the C-Points                                   */
    /* ============================================================= */
    auto boundaryNodes = graph->GetBoundaryNodeMap();
    Teuchos::Array<point_type> myPointType(myColors_h.size());
    LO num_c_points = 0, num_d_points=0;
    for(LO i=0; i<(LO)myColors_h.size(); i++) {
      if(boundaryNodes[i]) {
        myPointType[i] = DIRICHLET_PT;
        num_d_points++;
      }
      else if ((LO)myColors_h[i] == 1) {
        myPointType[i] = C_PT;
        num_c_points++;
      }
      else
        myPointType[i] = F_PT;
    }
    LO num_f_points = (LO)myColors_h.size() - num_d_points - num_c_points;
    // FIXME: This array will need to be ghosted so we can get the point_types
    // of the neighbors

    // FIXME:  These stats will need to be reduced
    GetOStream(Statistics1) << "ClassicalPFactory: C/F/D = "<<num_c_points<<"/"<<num_f_points<<"/"<<num_d_points<<std::endl;

    /* Generate the Coarse map */
    RCP<const Map> coarseMap;
    Array<LO> cpoint2pcol, pcol2cpoint;
    {
      SubFactoryMonitor sfm(*this,"Coarse Map",coarseLevel);
      GenerateCoarseMap(*A->getRowMap(),myPointType,coarseMap,cpoint2pcol,pcol2cpoint);
      Set(fineLevel, "CoarseMap", coarseMap);
    }
    // FIXME: cpoint2ccol needs to get ghosted

    // Generate edge strength flags (this will make everything easier later)
    Teuchos::Array<size_t> eis_rowptr;
    Teuchos::Array<bool> edgeIsStrong;
    {
      SubFactoryMonitor sfm(*this,"Strength Flags",coarseLevel);
      GenerateStrengthFlags(*A,*graph,eis_rowptr,edgeIsStrong);
    }

    // Phase 3: Generate the P matrix
    // FIXME: coarseColMap needs to ghosted
    RCP<const Map> coarseColMap = coarseMap;
    RCP<const Map> coarseDomainMap = coarseMap;

    std::string scheme = pL.get<std::string>("aggregation: classical scheme");
    if(scheme == "ext+i" || scheme == "classical")
      Coarsen_Ext_Plus_I(*A,*graph,coarseColMap,coarseDomainMap,num_c_points,num_f_points,myPointType,cpoint2pcol,pcol2cpoint,eis_rowptr,edgeIsStrong,P);
    else if(scheme == "direct")
      Coarsen_Direct(*A,*graph,coarseColMap,coarseDomainMap,num_c_points,num_f_points,myPointType,cpoint2pcol,pcol2cpoint,eis_rowptr,edgeIsStrong,P);

    //    Xpetra::IO<SC,LO,GO,NO>::Write("classical_p.mat", *P);
    Set(coarseLevel, "P", P);

    //RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(coarseMap, fineNullspace->getNumVectors());
    //    P->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS, Teuchos::ScalarTraits<SC>::one(), Teuchos::ScalarTraits<SC>::zero());
    //    Set(coarseLevel, "Nullspace", coarseNullspace);

    if (IsPrint(Statistics1)) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*P, "P", params);
    }
  }

/* ************************************************************************* */
template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Coarsen_Direct(const Matrix & A,const GraphBase & graph,  RCP<const Map> & coarseColMap, RCP<const Map> & coarseDomainMap, LO num_c_points, LO num_f_points, Teuchos::Array<point_type> & myPointType, const Teuchos::Array<LO> & cpoint2pcol, const Teuchos::Array<LO> & pcol2cpoint, Teuchos::Array<size_t> & eis_rowptr, Teuchos::Array<bool> & edgeIsStrong, RCP<Matrix> & P) const {

    /* ============================================================= */
    /* Phase 3 : Direct Interpolation                                */
    /* We do not use De Sterck, Falgout, Nolting and Yang (2008)     */
    /* here.  Instead we follow:                                     */
    /* Trottenberg, Oosterlee and Schueller, Multigrid, 2001.        */
    /* ============================================================= */    
    /* Definitions:                                                        */
    /* F = F-points                                                        */
    /* C = C-points                                                        */
    /* N_i = non-zero neighbors of node i                                  */
    /* S_i = {j\in N_i | j strongly influences i } [strong neighbors of i] */
    /* F_i^s = F \cap S_i [strong F-neighbors of i]                        */
    /* C_i^s = C \cap S_i [strong C-neighbors of i]                        */
    /* P_i = Set of interpolatory variables for row i [here = C_i^s]       */

    /* (A.2.17) from p. 426                                                */ 
    /* a_ij^- = {  a_ij,  if a_ij < 0                                      */ 
    /*          {     0,  otherwise                                        */ 
    /* a_ij^+ = {  a_ij,  if a_ij > 0                                      */ 
    /*          {     0,  otherwise                                        */ 
    /* P_i^- =  P_i \cap {k | a_ij^- != 0 and a_ij^- = a_ij}               */
    /*          [strong C-neighbors with negative edges]                   */
    /* P_i^+ =  P_i \cap {k | a_ij^+ != 0 and a_ij^+ = a_ij}               */
    /*          [strong C-neighbors with positive edges]                   */


    /* de Sterck et al., gives us this:                                                      */
    /* Rewritten Equation (6) on p. 119                                                      */
    /* w_ij = - a_ji / a_ii \frac{\sum_{k\in N_i} a_ik} {\sum k\inC_i^s} a_ik},   j\in C_i^s */

    /* Trottenberg et al. (A.7.6) and (A.7.7) on p. 479 gives this:                          */
    /* alpha_i = \frac{ \sum_{j\in N_i} a_ij^- }{ \sum_{k\in P_i} a_ik^- }                   */
    /* beta_i  = \frac{ \sum_{j\in N_i} a_ij^+ }{ \sum_{k\in P_i} a_ik^+ }                   */
    /* w_ik    = { - alpha_i (a_ik / a_ii),   if k\in P_i^-                                  */
    /*           { -  beta_i (a_ik / a_ii),   if k\in P_i^+                                  */  
    /* NOTE: The text says to modify, if  P_i^+ is zero but it isn't entirely clear how that */
    /* works.  We'll follow the PyAMG implementation in a few important ways.                */
    
   

    // Initial (estimated) allocation
    // NOTE: If we only used Tpetra, then we could use these guys as is, but because Epetra, we can't, so there
    // needs to be a copy below.
    using STS = typename Teuchos::ScalarTraits<SC>;
    size_t Nrows = A.getNodeNumRows();
    double c_point_density = (double)num_c_points / (num_c_points+num_f_points);
    double mean_strong_neighbors_per_row = (double) graph.GetNodeNumEdges() / graph.GetNodeNumVertices();
    //    double mean_neighbors_per_row = (double)A.getNodeNumEntries() / Nrows;
    double nnz_per_row_est = c_point_density*mean_strong_neighbors_per_row;

    size_t nnz_est = std::max(Nrows,std::min((size_t)A.getNodeNumEntries(),(size_t)(nnz_per_row_est*Nrows)));
    SC SC_ZERO = Teuchos::ScalarTraits<SC>::zero();
    Array<size_t> tmp_rowptr(Nrows+1);
    Array<LO> tmp_colind(nnz_est);

    // Algorithm (count+realloc)
    // For each row, i, 
    // - Count the number of elements in \hat{C}_j, aka [C-neighbors and C-neighbors of strong F-neighbors of i]   
    size_t ct=0;
    for(LO row=0; row < (LO) Nrows; row++) {
      size_t row_start = eis_rowptr[row];
      ArrayView<const LO> indices;
      ArrayView<const SC> vals;
      std::set<LO> C_hat;

      if(myPointType[row] == C_PT) {
        // C-Points get a single 1 in their row
        C_hat.insert(cpoint2pcol[row]);
      }
      else {
        // C-neighbors of row 
        A.getLocalRowView(row, indices, vals);
        for(LO j=0; j<indices.size(); j++)
          if(myPointType[indices[j]] == C_PT && edgeIsStrong[row_start + j])
            C_hat.insert(cpoint2pcol[indices[j]]);
      }// end else 
      
      // Realloc if needed
      if(ct + (size_t)C_hat.size() > (size_t)tmp_colind.size()) {
        tmp_colind.resize(std::max(ct+(size_t)C_hat.size(),(size_t)2*tmp_colind.size()));
      }
      
      // Copy
      std::copy(C_hat.begin(), C_hat.end(),tmp_colind.begin()+ct);
      ct+=C_hat.size();
      tmp_rowptr[row+1] = tmp_rowptr[row] + C_hat.size();
    }
    // Resize down
    tmp_colind.resize(tmp_rowptr[Nrows]);  

    // Allocate memory & copy
    P = rcp(new CrsMatrixWrap(A.getRowMap(), coarseColMap, 0));
    RCP<CrsMatrix> PCrs   = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();
    ArrayRCP<size_t>  P_rowptr;
    ArrayRCP<LO>      P_colind;
    ArrayRCP<SC>      P_values;
    PCrs->allocateAllValues(tmp_rowptr[Nrows], P_rowptr, P_colind, P_values);
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_rowptr.size() !=P_rowptr.size(), Exceptions::RuntimeError,"ClassicalPFactory: Allocation size error (rowptr)");
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_colind.size() !=P_colind.size(), Exceptions::RuntimeError,"ClassicalPFactory: Allocation size error (colind)");

    // FIXME:  This can be short-circuited for Tpetra, if we decide we care
    std::copy(tmp_rowptr.begin(),tmp_rowptr.end(), P_rowptr.begin());
    std::copy(tmp_colind.begin(),tmp_colind.end(), P_colind.begin());     


    // Algorithm (numeric)
    for(LO i=0; i < (LO)Nrows; i++) {
      if (myPointType[i] == C_PT) {
        // C Points get a single 1 in their row
        P_values[P_rowptr[i]] = Teuchos::ScalarTraits<SC>::one();  
#ifdef CMS_DEBUG        
        // DEBUG
        printf("** A(%d,:) is a C-Point.\n",i);
#endif
      }
      else {
        /* Trottenberg et al. (A.7.6) and (A.7.7) on p. 479 gives this:                          */
        /* alpha_i = \frac{ \sum_{j\in N_i} a_ij^- }{ \sum_{k\in P_i} a_ik^- }                   */
        /* beta_i  = \frac{ \sum_{j\in N_i} a_ij^+ }{ \sum_{k\in P_i} a_ik^+ }                   */
        /* w_ik    = { - alpha_i (a_ik / a_ii),   if k\in P_i^-                                  */
        /*           { -  beta_i (a_ik / a_ii),   if k\in P_i^+                                  */  
        ArrayView<const LO> A_indices_i, A_incides_k;
        ArrayView<const SC> A_vals_i, A_indices_k;
        A.getLocalRowView(i, A_indices_i, A_vals_i);
        size_t row_start = eis_rowptr[i];
        
        ArrayView<LO> P_indices_i  = P_colind.view(P_rowptr[i],P_rowptr[i+1] - P_rowptr[i]);
        ArrayView<SC> P_vals_i     = P_values.view(P_rowptr[i],P_rowptr[i+1] - P_rowptr[i]);
        
#ifdef CMS_DEBUG          
        // DEBUG
        {
          char mylabel[5]="UFCD";
          char sw[3]="ws";
          printf("** A(%d,:) = ",i);
          for(LO j=0; j<(LO)A_indices_i.size(); j++){  
            printf("%6.4e(%d-%c%c) ",A_vals_i[j],A_indices_i[j],mylabel[myPointType[A_indices_i[j]]],sw[(int)edgeIsStrong[row_start+j]]);
          }
          printf("\n");
        }
#endif        
        
        
        SC a_ii            = SC_ZERO;
        SC pos_numerator   = SC_ZERO, neg_numerator   = SC_ZERO;
        SC pos_denominator = SC_ZERO, neg_denominator = SC_ZERO;
        // Find the diagonal and compute the sum ratio
        for(LO j=0; j<(LO)A_indices_i.size(); j++) {
          SC a_ik = A_vals_i[j]; 
          LO k = A_indices_i[j];
          
          // Diagonal
          if(i == k) { 
            a_ii = a_ik;
          }          
          // Only strong C-neighbors are in the denomintor
          // FIXME: myPointType needs to be ghosted
          if(myPointType[k] == C_PT && edgeIsStrong[row_start + j]) {
            if(STS::real(a_ik) > SC_ZERO) pos_denominator += a_ik;
            else neg_denominator += a_ik;
          }  
          
          // All neighbors are in the numerator
          // NOTE: As per PyAMG, this does not include the diagonal
          if(i != k) {
            if(STS::real(a_ik) > SC_ZERO) pos_numerator += a_ik;
            else neg_numerator += a_ik;
          }   
        }
        SC alpha = (neg_denominator == SC_ZERO) ? SC_ZERO : (neg_numerator / neg_denominator);
        SC beta  = (pos_denominator == SC_ZERO) ? SC_ZERO : (pos_numerator / pos_denominator);
        alpha /= -a_ii;        
        beta  /= -a_ii;

        // Loop over the entries
        for(LO p_j=0; p_j<(LO)P_indices_i.size(); p_j++){  
          LO P_col = pcol2cpoint[P_indices_i[p_j]];
          SC a_ij = SC_ZERO;
          
          // Find A_ij (if it is there)
          // FIXME: We can optimize this if we assume sorting
          for(LO a_j =0; a_j<(LO)A_indices_i.size(); a_j++) {
            if(A_indices_i[a_j] == P_col) {
              a_ij = A_vals_i[a_j];
              break;
            }
          }
          SC w_ij = (STS::real(a_ij) < 0 ) ? (alpha * a_ij) : (beta * a_ij);
#ifdef CMS_DEBUG
          SC alpha_or_beta = (STS::real(a_ij) < 0 ) ? alpha : beta;
          printf("P(%d,%d/%d) =  - %6.4e  * %6.4e  = %6.4e\n",i,P_indices_i[p_j],pcol2cpoint[P_indices_i[p_j]],alpha_or_beta,a_ij,w_ij);
#endif
          P_vals_i[p_j] = w_ij;          
        }//end for A_indices_i
      }//end else C_PT
    }//end for Numrows

    // Finish up
    PCrs->setAllValues(P_rowptr, P_colind, P_values);
    PCrs->fillComplete(/*domain*/coarseDomainMap, /*range*/A.getDomainMap());// Yes, we want A's domainMap here.

}


/* ************************************************************************* */
template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Coarsen_Ext_Plus_I(const Matrix & A,const GraphBase & graph,  RCP<const Map> & coarseColMap, RCP<const Map> & coarseDomainMap, LO num_c_points, LO num_f_points, Teuchos::Array<point_type> & myPointType, const Teuchos::Array<LO> & cpoint2pcol, const Teuchos::Array<LO> & pcol2cpoint, Teuchos::Array<size_t> & eis_rowptr, Teuchos::Array<bool> & edgeIsStrong, RCP<Matrix> & P) const {

    /* ============================================================= */
    /* Phase 3 : Extended+i Interpolation                            */
    /* De Sterck, Falgout, Nolting and Yang. "Distance-two           */
    /* interpolation for parallel algebraic multigrid", NLAA 2008    */
    /* 15:115-139                                                    */
    /* ============================================================= */    
    /* Definitions:                                                        */
    /* F = F-points                                                        */
    /* C = C-points                                                        */
    /* N_i = non-zero neighbors of node i                                  */
    /* S_i = {j\in N_i | j strongly influences i } [strong neighbors of i] */
    /* F_i^s = F \cap S_i [strong F-neighbors of i]                        */
    /* C_i^s = C \cap S_i [strong C-neighbors of i]                        */
    /* N_i^w = N_i\ (F_i^s \cup C_i^s) [weak neighbors of i]               */
    /*         This guy has a typo.  The paper had a \cap instead of \cup  */
    /*         I would note that this set can contain both F-points and    */
    /*         C-points.  They're just weak neighbors of this guy.         */
    /*         Note that N_i^w \cup F_i^s \cup C_i^s = N_i by construction */

    /* \hat{C}_i = C_i \cup (\bigcup_{j\inF_i^s} C_j)                      */
    /*         [C-neighbors and C-neighbors of strong F-neighbors of i]    */
    /*                                                                     */

    /* Rewritten Equation (19) on p. 123                                   */
    /* f_ik = \frac{\bar{a}_kj}{\sum{l\in \hat{C}_i\cup {i}} \bar{a}_kl    */
    /* w_ij = -\tilde{a}_ii^{-1} (a_ij + \sum_{k\inF_i^s} a_ik f_ik        */
    /*         for j in \hat{C}_i                                          */
    
    /* Rewritten Equation (20) on p. 124 [for the lumped diagonal]                                  */
    /* g_ik = \frac{\bar{a}_ki}{\sum{l\in \hat{C}_i\cup {i}} \bar{a}_kl                             */    
    /* \tilde{a}_ii = a_ii + \sum_{n\inN_i^w\setminus \hat{C}_i} a_in + \sum_{k\inF_i^s} a_ik g_ik  */

    using STS=Teuchos::ScalarTraits<SC>;
    LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
    SC SC_ZERO = STS::zero();

    // Initial (estimated) allocation
    // NOTE: If we only used Tpetra, then we could use these guys as is, but because Epetra, we can't, so there
    // needs to be a copy below.
    size_t Nrows = A.getNodeNumRows();
    double c_point_density = (double)num_c_points / (num_c_points+num_f_points);
    double mean_strong_neighbors_per_row = (double) graph.GetNodeNumEdges() / graph.GetNodeNumVertices();
    double mean_neighbors_per_row = (double)A.getNodeNumEntries() / Nrows;
    double nnz_per_row_est = c_point_density*mean_neighbors_per_row + (1-c_point_density)*mean_strong_neighbors_per_row*c_point_density*mean_neighbors_per_row;

    size_t nnz_est = std::max(Nrows,std::min((size_t)A.getNodeNumEntries(),(size_t)(nnz_per_row_est*Nrows)));
    Array<size_t> tmp_rowptr(Nrows+1);
    Array<LO> tmp_colind(nnz_est);

    // Algorithm (count+realloc)
    // For each row, i, 
    // - Count the number of elements in \hat{C}_j, aka [C-neighbors and C-neighbors of strong F-neighbors of i]   
    size_t ct=0;
    for(LO row=0; row < (LO) Nrows; row++) {
      ArrayView<const LO> indices;
      ArrayView<const SC> vals;
      std::set<LO> C_hat;

      if(myPointType[row] == C_PT) {
        // C-Points get a single 1 in their row
        C_hat.insert(cpoint2pcol[row]);
      }
      else {
        // C-neighbors of row (strong or not)
        A.getLocalRowView(row, indices, vals);
        for(LO j=0; j<indices.size(); j++)
          if(myPointType[indices[j]] == C_PT)
            C_hat.insert(cpoint2pcol[indices[j]]);
        
        // Strong neighbors of row      
        ArrayView<const LO> strong_neighbors = graph.getNeighborVertices(row);
        for(LO j=0; j<(LO)strong_neighbors.size(); j++) {
          if(myPointType[strong_neighbors[j]] == F_PT) {
            LO row2 = strong_neighbors[j];
            A.getLocalRowView(row2, indices, vals);
            // C-neighbors of strong F-neighbors
            for(LO k=0; k<indices.size(); k++) {
              if(myPointType[indices[k]] == C_PT)
                C_hat.insert(cpoint2pcol[indices[k]]);
            }         
          }
        }
      }// end else 

      // Realloc if needed
      if(ct + (size_t)C_hat.size() > (size_t)tmp_colind.size()) {
        tmp_colind.resize(std::max(ct+(size_t)C_hat.size(),(size_t)2*tmp_colind.size()));
      }
      
      // Copy
      std::copy(C_hat.begin(), C_hat.end(),tmp_colind.begin()+ct);
      ct+=C_hat.size();
      tmp_rowptr[row+1] = tmp_rowptr[row] + C_hat.size();
    }
    // Resize down
    tmp_colind.resize(tmp_rowptr[Nrows]);


    // Allocate memory & copy
    P = rcp(new CrsMatrixWrap(A.getRowMap(), coarseColMap, 0));
    RCP<CrsMatrix> PCrs   = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();
    ArrayRCP<size_t>  P_rowptr;
    ArrayRCP<LO>      P_colind;
    ArrayRCP<SC>      P_values;
    PCrs->allocateAllValues(tmp_rowptr[Nrows], P_rowptr, P_colind, P_values);
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_rowptr.size() !=P_rowptr.size(), Exceptions::RuntimeError,"ClassicalPFactory: Allocation size error (rowptr)");
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_colind.size() !=P_colind.size(), Exceptions::RuntimeError,"ClassicalPFactory: Allocation size error (colind)");

    // FIXME:  This can be short-circuited for Tpetra, if we decide we care
    std::copy(tmp_rowptr.begin(),tmp_rowptr.end(), P_rowptr.begin());
    std::copy(tmp_colind.begin(),tmp_colind.end(), P_colind.begin());

    // Algorithm (numeric)
    // For each row
    // - Pass though the row to compute the lumped diagonal.
    // - Pass through the row a second time to compute the entries
    
    /* \bar{a}_ij = {    0, if sign(a_ij) == sign(a_ii)                    */
    /*              { a_ij, otherwise                                      */


    // - Pass though the row to compute the lumped diagonal and the denominator of the last term
    //   in (19) and (20)
    Array<SC> denominator_sum_i(A.getNodeMaxNumRowEntries());    

    for(LO i=0; i < (LO)Nrows; i++) {
      /* Rewritten Equation (20) on p. 124 [for the lumped diagonal]                        */
      /* g_ik = \frac{\bar{a}_ki}{\sum{l\in \hat{C}_i\cup {i}} \bar{a}_kl                   */    
      /* \tilde{a}_ii = a_ii + \sum_{n\inN_i^w\setminus \hat{C}_i} a_in + \sum_{k\inF_i^s} a_ik g_ik  */

      if (myPointType[i] == C_PT) {
        // C Points get a single 1 in their row
        P_values[P_rowptr[i]] = Teuchos::ScalarTraits<SC>::one();
#ifdef CMS_DEBUG        
        // DEBUG
        printf("** A(%d,:) is a C-Point.\n",i);
#endif
      }
      else {
        SC atilde_ii = SC_ZERO;
        ArrayView<const LO> A_indices_i, A_indices_k;
        ArrayView<const SC> A_vals_i, A_vals_k;
        A.getLocalRowView(i, A_indices_i, A_vals_i);
        size_t row_start = eis_rowptr[i];
      
#ifdef CMS_DEBUG          
        // DEBUG
        {
          char mylabel[5]="UFCD";
          char sw[3]="ws";
          printf("** A(%d,:) = ",i);
          for(LO j=0; j<(LO)A_indices_i.size(); j++){  
            printf("%6.4e(%d-%c%c) ",A_vals_i[j],A_indices_i[j],mylabel[myPointType[A_indices_i[j]]],sw[(int)edgeIsStrong[row_start+j]]);
          }
          printf("\n");
        }
#endif

        for(LO j=0; j<(LO)A_indices_i.size(); j++){  
          // Init this to zero, for columns that aren't F-strong neighbors
          denominator_sum_i[j] = 0;
          
          // NOTE: Assumes row/col compatible local indexing on A 
          if(i == A_indices_i[j]) {
            // (a) a_ii
            atilde_ii += A_vals_i[j];
          }
          else if(!edgeIsStrong[row_start +j] && cpoint2pcol[A_indices_i[j]] != LO_INVALID) {
            // (c) Weak neighbor that is in \hat{C}_i (which means it it in P)
            atilde_ii += A_vals_i[j];          
          }
          else if(edgeIsStrong[row_start+j] && myPointType[A_indices_i[j]] == F_PT && A_vals_i[j] != SC_ZERO) {

            // (b) Strong F-neighbor
            // Added a short-circuit for identically zero a_ik here
            LO k     = A_indices_i[j];
            SC a_ik  = A_vals_i[j];
#ifdef CMS_DEBUG
            printf(" - A(%d,%d) is strong F-neighbor\n",i,k);
#endif
            // FIXME: This may try to get an off-proc row in parallel.
            A.getLocalRowView(k, A_indices_k, A_vals_k);
            
            // First, akk (because we need the sign)
            SC a_kk = SC_ZERO;
            for(LO l=0; l<(LO)A_indices_k.size(); l++) {  
              if(k == A_indices_k[l]) {
                a_kk = A_vals_k[l];
                break;
              }
            }
            int sign_akk = Sign(a_kk);
            
            SC a_ki_bar = SC_ZERO;
            SC a_kl_bar_sum = SC_ZERO;
            for(LO l=0; l<(LO)A_indices_k.size(); l++) {
              // Find Aki
              if(i == A_indices_k[l]) {
                SC a_ki = A_vals_k[l];              
                a_ki_bar = (Sign(a_ki) == sign_akk) ? SC_ZERO : a_ki;
                a_kl_bar_sum += a_ki_bar;
#ifdef CMS_DEBUG
                printf(" - - Found A(%d,%d) aki_bar = %6.4e\n",k,i,a_ki_bar);
#endif
              }
              else {
                // Is this in \hat{C}_i? (which means it it in P)
                if (cpoint2pcol[A_indices_k[l]] != LO_INVALID) {
                  SC a_kl = A_vals_k[l];
                  SC a_kl_bar = (Sign(a_kl) == sign_akk) ? SC_ZERO : a_kl;
#ifdef CMS_DEBUG
                  printf(" - - Found A(%d,%d) c-point akl_bar = %6.4e\n",k,A_indices_k[l],a_kl_bar);
#endif
                  a_kl_bar_sum += a_kl_bar;
                }
              }            
            }// end for A_indices_k 
            
            // Now, add it all together.
            // FIXME: This probably won't compile in complex
            denominator_sum_i[j] = a_kl_bar_sum;
            atilde_ii +=  a_ik * a_ki_bar / a_kl_bar_sum;
          }// end if Strong F-point neighbor of i
        }// end for A_indices_i 

#ifdef CMS_DEBUG        
        {
          // DEBUG
          printf("* denominator_sum_i =");
          for(LO jj=0; jj<(LO)A_indices_i.size(); jj++){  
            printf("%6.4e(%d-%d) ",denominator_sum_i[jj],jj,A_indices_i[jj]);
          }
          printf("\n");
        }
#endif
        
        
        /* Rewritten Equation (19) on p. 123                                   */
        /* f_ik = \frac{\bar{a}_kj}{\sum{l\in \hat{C}_i\cup {i}} \bar{a}_kl    */
        /* w_ij = -\tilde{a}_ii^{-1} (a_ij + \sum_{k\inF_i^s} a_ik f_ik        */
        /*         for j in \hat{C}_i                                          */     
        ArrayView<LO> P_indices_i = P_colind.view(P_rowptr[i],P_rowptr[i+1] - P_rowptr[i]);
        ArrayView<SC> P_vals_i     = P_values.view(P_rowptr[i],P_rowptr[i+1] - P_rowptr[i]);
        
        for(LO p_j=0; p_j <(LO)P_indices_i.size(); p_j++) {
          LO P_col = pcol2cpoint[P_indices_i[p_j]];
          SC a_ij = SC_ZERO;
          
          // Find A_ij (if it is there, which it might not be, because P is "extended")
          // FIXME: We can optimize this if we assume sorting
          for(LO a_j =0; a_j<(LO)A_indices_i.size(); a_j++) {
            if(A_indices_i[a_j] == P_col) {
              a_ij = A_vals_i[a_j];
              break;
            }
          }
          
          // Loop over all the F-strong neighbors in A to get that second term
          SC second_term=SC_ZERO;
          for(LO j=0; j<(LO)A_indices_i.size(); j++){  
            if(edgeIsStrong[row_start +j] && (myPointType[A_indices_i[j]] == F_PT) && (A_vals_i[j] != SC_ZERO)) {          
              // CMS Mod: Not counting self-connections here.
              if(i==A_indices_i[j]) continue;
   

              // FIXME: This is the off-rank deal again
              LO k    = A_indices_i[j];
              SC a_ik = A_vals_i[j];
              SC a_kj_bar = SC_ZERO;
              A.getLocalRowView(k, A_indices_k, A_vals_k);         
              
              // First, akk (because we need the sign)
              SC a_kk = SC_ZERO;
              for(LO l=0; l<(LO)A_indices_k.size(); l++) {  
                if(k == A_indices_k[l]) {
                  a_kk = A_vals_k[l];
                  break;
                }
              }
              int sign_akk = Sign(a_kk);
              
              for(LO l=0; l<(LO)A_indices_k.size(); l++) {
                if(A_indices_k[l] == P_col) {
                  SC a_kj = A_vals_k[l];              
                  a_kj_bar = (Sign(a_kj) == sign_akk) ? SC_ZERO : a_kj;
                  break;
                }
              }//end for A_indices_k
              second_term+= a_ik * a_kj_bar / denominator_sum_i[j];
            }//end if strong F-neighbor
          }//end A_indices_i
          
          SC w_ij = -(a_ij + second_term) / atilde_ii;
#ifdef CMS_DEBUG
          printf("P(%d,%d/%d) = (%6.4e + %6.4e) / %6.4e = %6.4e\n",i,P_indices_i[p_j],pcol2cpoint[P_indices_i[p_j]],a_ij,second_term,atilde_ii,w_ij);
#endif
          P_vals_i[p_j] = w_ij;
          
        }//end for P_indices
      }//end else is C_PT
    }//end for Nrows    



    // Finish up
    PCrs->setAllValues(P_rowptr, P_colind, P_values);
    PCrs->fillComplete(/*domain*/coarseDomainMap, /*range*/A.getDomainMap());// Yes, we want A's domainMap here.
}


/* ************************************************************************* */
template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  GenerateCoarseMap(const Map & fineMap, const Teuchos::Array<point_type> & fc_splitting, RCP<const Map> & coarseMap, Array<LO> & cpoint2pcol, Array<LO> & pcol2cpoint) const {
  LO num_c_points = 0;
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  // We need to make sure all non-C points are LO_INVALID
  cpoint2pcol.resize(0);
  cpoint2pcol.resize(fc_splitting.size(),LO_INVALID);
  
  for(LO i=0; i<(LO) fc_splitting.size(); i++)
    if(fc_splitting[i] == C_PT) {
      cpoint2pcol[i] = num_c_points;
      num_c_points++;
    }

  pcol2cpoint.resize(0);
  pcol2cpoint.resize(num_c_points);
    for(LO i=0; i<(LO)cpoint2pcol.size(); i++)
      if(cpoint2pcol[i] != LO_INVALID)
        pcol2cpoint[cpoint2pcol[i]] =i;

  // FIXME: Assumes scalar PDE
  std::vector<size_t> stridingInfo_(1);
  stridingInfo_[0]=1;
  GO domainGIDOffset = 0;

  coarseMap = StridedMapFactory::Build(fineMap.lib(),
                                       Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                       num_c_points,
                                       fineMap.getIndexBase(),
                                       stridingInfo_,
                                       fineMap.getComm(),
                                       domainGIDOffset); 
}

/* ************************************************************************* */
template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
GenerateStrengthFlags(const Matrix & A,const GraphBase & graph, Teuchos::Array<size_t> & eis_rowptr,Teuchos::Array<bool> & edgeIsStrong) const {
  // To make this easier, we'll create a bool array equal to the nnz in the matrix
  // so we know whether each edge is strong or not.  This will save us a bunch of
  // trying to match the graph and matrix later
  size_t Nrows = A.getNodeNumRows();
  eis_rowptr.resize(Nrows+1);

  if(edgeIsStrong.size() == 0) {
    // Preferred
    edgeIsStrong.resize(A.getNodeNumEntries(),false);
  }
  else {
    edgeIsStrong.resize(A.getNodeNumEntries(),false);
    edgeIsStrong.assign(A.getNodeNumEntries(),false);
  }
  
  eis_rowptr[0] = 0;
  for (LO i=0; i<(LO)Nrows; i++) {
    LO rowstart = eis_rowptr[i];
    ArrayView<const LO> A_indices;
    ArrayView<const SC> A_values;
    A.getLocalRowView(i, A_indices, A_values);
    LO A_size = (LO) A_indices.size();

    ArrayView<const LO> G_indices = graph.getNeighborVertices(i);
    LO G_size = (LO) G_indices.size();
    
    // Both of these guys should be in the same (sorted) order, but let's check
    bool is_ok=true;
    for(LO j=0; j<A_size-1; j++)
      if (A_indices[j] > A_indices[j+1]) { is_ok=false; break;}
    for(LO j=0; j<G_size-1; j++)
      if (G_indices[j] > G_indices[j+1]) { is_ok=false; break;}
    TEUCHOS_TEST_FOR_EXCEPTION(!is_ok, Exceptions::RuntimeError,"ClassicalPFactory: Exected A and Graph to be sorted");
    
    // Now cycle through and set the flags - if the edge is in G it is strong
    for(LO g_idx=0, a_idx=0; g_idx < G_size; g_idx++) {
      LO col = G_indices[g_idx];
      while (A_indices[a_idx] != col && a_idx < A_size) a_idx++;
      if(a_idx == A_size) {is_ok=false;break;}
      edgeIsStrong[rowstart+a_idx] = true;      
    }

    eis_rowptr[i+1] = eis_rowptr[i] + A_size;
  }
}



/* ************************************************************************* */
template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  DoGraphColoring(const GraphBase & graph, colors_view_type & myColors, LO & numColors) const {
  const ParameterList& pL = GetParameterList();
  
  using graph_t = typename LWGraph_kokkos::local_graph_type;
  using KernelHandle = KokkosKernels::Experimental::
    KokkosKernelsHandle<typename graph_t::row_map_type::value_type,
                        typename graph_t::entries_type::value_type,
                        typename graph_t::entries_type::value_type,
                        typename graph_t::device_type::execution_space,
                        typename graph_t::device_type::memory_space,
                        typename graph_t::device_type::memory_space>;
  KernelHandle kh;

  // Leave gc algorithm choice as the default
  kh.create_graph_coloring_handle();
  
  // Get the distance-1 graph coloring handle
  auto coloringHandle = kh.get_graph_coloring_handle();      
  
  // Set the distance-1 coloring algorithm to use
  if(pL.get<bool>("aggregation: deterministic") == true) {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_SERIAL );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: serial" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "serial") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_SERIAL );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: serial" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "vertex based") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_VB );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "vertex based bit array") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_VBBIT );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based bit array" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "vertex based color set") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_VBCS );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based color set" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "vertex based deterministic") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_VBD );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based deterministic" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "vertex based deterministic bit array") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_VBDBIT );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based deterministic bit array" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "edge based") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_EB );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: edge based" << std::endl;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Unrecognized distance 1 coloring algorithm");
  }
  
  // Create device views for graph rowptrs/colinds
  size_t numRows = graph.GetNodeNumVertices();
  auto graphLWK = dynamic_cast<const LWGraph_kokkos*>(&graph);
  auto graphLW  = dynamic_cast<const LWGraph*>(&graph);
  auto graphG   = dynamic_cast<const Graph*>(&graph);
  TEUCHOS_TEST_FOR_EXCEPTION(!graphLW && !graphLWK && !graphG,std::invalid_argument,"Graph is not a LWGraph or LWGraph_kokkos object");
    // Run d1 graph coloring
    // Assume that the graph is symmetric so row map/entries and col map/entries are the same

  if(graphLWK) {
    KokkosGraph::Experimental::graph_color(&kh, 
                                           numRows, 
                                           numRows, // FIXME: This should be the number of columns
                                           graphLWK->getRowPtrs(),
                                           graphLWK->getEntries(),
                                         true);
  }
  else if(graphLW) {
    auto rowptrs = graphLW->getRowPtrs();
    auto entries = graphLW->getEntries();
    // Copy rowptrs to a size_t, because kokkos-kernels doesn't like rowptrs as LO's
    Teuchos::Array<size_t> rowptrs_s(rowptrs.size());
    std::copy(rowptrs.begin(),rowptrs.end(),rowptrs_s.begin());
    Kokkos::View<const size_t*,Kokkos::LayoutLeft,Kokkos::HostSpace> rowptrs_v(rowptrs_s.data(),(size_t)rowptrs.size());
    Kokkos::View<const LO*,Kokkos::LayoutLeft,Kokkos::HostSpace> entries_v(entries.getRawPtr(),(size_t)entries.size());
    KokkosGraph::Experimental::graph_color(&kh, 
                                           numRows, 
                                           numRows, // FIXME: This should be the number of columns
                                           rowptrs_v,
                                           entries_v,
                                           true);    
  }
  else if(graphG) {  
    // FIXME:  This is a terrible, terrible hack, based on 0-based local indexing.
    RCP<const CrsGraph> graphC = graphG->GetGraph();
    size_t numEntries = graphC->getNodeNumEntries();
    ArrayView<const LO> indices;
    graphC->getLocalRowView(0,indices);
    Kokkos::View<size_t*,Kokkos::LayoutLeft,Kokkos::HostSpace> rowptrs_v("rowptrs_v",graphC->getNodeNumRows()+1);
    rowptrs_v[0]=0;
    for(LO i=0; i<(LO)graphC->getNodeNumRows()+1; i++) 
      rowptrs_v[i+1] = rowptrs_v[i] + graphC->getNumEntriesInLocalRow(i);
    Kokkos::View<const LO*,Kokkos::LayoutLeft,Kokkos::HostSpace> entries_v(&indices[0],numEntries);    
    KokkosGraph::Experimental::graph_color(&kh, 
                                           numRows, 
                                           numRows, // FIXME: This should be the number of columns
                                           rowptrs_v,
                                           entries_v,
                                           true);   
    
  }

  
  // Extract the colors and store them in the aggregates
  myColors = coloringHandle->get_vertex_colors();
  numColors = static_cast<LO>(coloringHandle->get_num_colors());
  
  //clean up coloring handle
  kh.destroy_graph_coloring_handle();

  
}// end DoGraphColoring
    
    






} //namespace MueLu



#define MUELU_CLASSICALPFACTORY_SHORT
#endif // MUELU_CLASSICALPFACTORY_DEF_HPP


