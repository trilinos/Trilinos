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
#ifndef MUELU_IPCFACTORY_DEF_HPP
#define MUELU_IPCFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <sstream>

#include "MueLu_IntrepidPCoarsenFactory_decl.hpp"

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Utilities.hpp"

#include "Teuchos_ScalarTraits.hpp"

// Intrepid Headers

//Intrepid_HGRAD_HEX_C1_FEM.hpp
//Intrepid_HGRAD_HEX_C2_FEM.hpp
//Intrepid_HGRAD_HEX_Cn_FEM.hpp
//Intrepid_HGRAD_HEX_I2_FEM.hpp
#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"
//Intrepid_HGRAD_LINE_Cn_FEM_JACOBI.hpp
//Intrepid_HGRAD_POLY_C1_FEM.hpp
//Intrepid_HGRAD_PYR_C1_FEM.hpp
//Intrepid_HGRAD_PYR_I2_FEM.hpp
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
//#include Intrepid_HGRAD_QUAD_C2_FEM.hpp
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
//Intrepid_HGRAD_TET_C1_FEM.hpp
//Intrepid_HGRAD_TET_C2_FEM.hpp
//Intrepid_HGRAD_TET_Cn_FEM.hpp
//Intrepid_HGRAD_TET_Cn_FEM_ORTH.hpp
//Intrepid_HGRAD_TET_COMP12_FEM.hpp
//Intrepid_HGRAD_TRI_C1_FEM.hpp
//Intrepid_HGRAD_TRI_C2_FEM.hpp
//Intrepid_HGRAD_TRI_Cn_FEM.hpp
//Intrepid_HGRAD_TRI_Cn_FEM_ORTH.hpp
//Intrepid_HGRAD_WEDGE_C1_FEM.hpp
//Intrepid_HGRAD_WEDGE_C2_FEM.hpp
//Intrepid_HGRAD_WEDGE_I2_FEM.hpp

namespace MueLu {


/*********************************************************************************************************/
namespace MueLuIntrepid {
inline std::string tolower(const std::string & str) {
  std::string data(str);
  std::transform(data.begin(), data.end(), data.begin(), [](unsigned char c) { return std::tolower(c); });
  return data;
}

template<class Scalar>
Teuchos::RCP<Intrepid2::Basis<Scalar,Intrepid2::FieldContainer<Scalar> > >  BasisFactory(const std::string & name) {
    using std::string;
    using Teuchos::rcp;
    string myerror("IntrepidBasisFactory: cannot parse string name '"+name+"'");
    

    // Syntax [HGRAD|HCURL|HDIV][_| ][HEX|LINE|POLY|PYR|QUAD|TET|TRI|WEDGE][_| ][C|I][1|2|n]

    // Get the derivative type
    size_t pos1 = name.find_first_of(" _");
    if(pos1==0) throw std::runtime_error(myerror);
    string deriv = tolower(name.substr(0,pos1));
    if(deriv!="hgrad" && deriv!="hcurl" && deriv!="hdiv") throw std::runtime_error(myerror);

    // Get the element type
    pos1++;
    size_t pos2 = name.find_first_of(" _",pos1);
    if(pos2==0) throw std::runtime_error(myerror);
    string el = tolower(name.substr(pos1,pos2-pos1));
    if(el!="hex" && el!="line" && el!="poly" && el!="pyr" && el!="quad" && el!="tet" && el!="tri" && el!="wedge") throw std::runtime_error(myerror);
    
    // Get the polynomial type    
    pos2++;
    string poly = tolower(name.substr(pos2,1));
    if(poly!="c" && poly!="i") throw std::runtime_error(myerror);

    // Get the degree
    pos2++;
    int degree=std::stoi(name.substr(pos2,1));
    if(degree<=0) throw std::runtime_error(myerror);

    // FIXME LATER: Allow for alternative point types for Kirby elements
    if(deriv=="hgrad" && el=="quad" && poly=="c"){
      if(degree==1) return rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<Scalar,Intrepid2::FieldContainer<Scalar> >());
      else          return rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<Scalar,Intrepid2::FieldContainer<Scalar> >(degree,Intrepid2::POINTTYPE_EQUISPACED));
    }
    else if(deriv=="hgrad" && el=="line" && poly=="c"){
      if(degree==1) return rcp(new Intrepid2::Basis_HGRAD_LINE_C1_FEM<Scalar,Intrepid2::FieldContainer<Scalar> >());
      else          return rcp(new Intrepid2::Basis_HGRAD_LINE_Cn_FEM<Scalar,Intrepid2::FieldContainer<Scalar> >(degree,Intrepid2::POINTTYPE_EQUISPACED));
    }

    // Error out
    throw std::runtime_error(myerror);
    return Teuchos::null;
}

/*********************************************************************************************************/
template <class Scalar, class ArrayScalar>
void IntrepidGetLoNodeInHi(const Teuchos::RCP<Intrepid2::Basis<Scalar,ArrayScalar> > &hi_basis,
			   const Teuchos::RCP<Intrepid2::Basis<Scalar,ArrayScalar> > &lo_basis,
			   std::vector<size_t> & lo_node_in_hi,
			   ArrayScalar & hi_DofCoords) {
  
  // Figure out which unknowns in hi_basis correspond to nodes on lo_basis. This varies by element type.
  size_t degree         = hi_basis->getDegree();
  lo_node_in_hi.resize(0);
  RCP<Intrepid2::DofCoordsInterface<ArrayScalar> > hi_dci;
  if(!rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<Scalar,ArrayScalar> >(hi_basis).is_null()) {
    // HGRAD QUAD Cn: Numbering as per the Kirby convention (straight across, bottom to top) 
    lo_node_in_hi.insert(lo_node_in_hi.end(),{0,degree, (degree+1)*(degree+1)-1, degree*(degree+1)});
    hi_dci = rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<Scalar,ArrayScalar> >(hi_basis);
  }
  else if(!rcp_dynamic_cast<Intrepid2::Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar> >(hi_basis).is_null()) {
    // HGRAD LINE Cn: Numbering as per the Kirby convention (straight across) 
    lo_node_in_hi.insert(lo_node_in_hi.end(),{0,degree});
    hi_dci = rcp_dynamic_cast<Intrepid2::Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar> >(hi_basis);
  } 
  else
    throw std::runtime_error("IntrepidPCoarsenFactory: Unknown element type");
  
  // Get coordinates of the hi_basis dof's
  hi_DofCoords.resize(hi_basis->getCardinality(),hi_basis->getBaseCellTopology().getDimension());
  hi_dci->getDofCoords(hi_DofCoords);
}


/*********************************************************************************************************/
template <class LocalOrdinal>
void BuildLoElemToNode(const Intrepid2::FieldContainer<LocalOrdinal> & hi_elemToNode,
		       const std::vector<bool> & hi_nodeIsOwned,
		       const std::vector<size_t> & lo_node_in_hi,
		       Intrepid2::FieldContainer<LocalOrdinal> & lo_elemToNode,
		       std::vector<bool> & lo_nodeIsOwned,
		       std::vector<LocalOrdinal> & hi_to_lo_map,
		       int & lo_numOwnedNodes) {
  typedef LocalOrdinal LO;
  using Teuchos::RCP;
  //  printf("CMS:BuildLoElemToNode: hi_elemToNode.rank() = %d hi_elemToNode.size() = %d\n",hi_elemToNode.rank(), hi_elemToNode.size());

  size_t numElem     = hi_elemToNode.dimension(0);
  size_t hi_numNodes = hi_nodeIsOwned.size();

  size_t lo_nperel = lo_node_in_hi.size();
  lo_elemToNode.resize(numElem, lo_nperel);

  // Build lo_elemToNode (in the hi local index ordering) and flag owned ones
  std::vector<bool> is_low_order(hi_numNodes,false);
  for(size_t i=0; i<numElem; i++)
    for(size_t j=0; j<lo_nperel; j++) {
      lo_elemToNode(i,j)  = hi_elemToNode(i,lo_node_in_hi[j]);
      is_low_order[hi_elemToNode(i,lo_node_in_hi[j])] = true; // This can overwrite and that is OK.
    }

  // Count the number of lo owned nodes, generating a local index for lo nodes
  lo_numOwnedNodes=0;
  size_t lo_numNodes=0;
  hi_to_lo_map.resize(hi_numNodes,Teuchos::OrdinalTraits<LO>::invalid());

  for(size_t i=0; i<hi_numNodes; i++)
    if(is_low_order[i]) {
      hi_to_lo_map[i] = lo_numNodes;
      lo_numNodes++;
      if(hi_nodeIsOwned[i]) lo_numOwnedNodes++;
    }

  // Flag the owned lo nodes
  lo_nodeIsOwned.resize(lo_numNodes,false);
  for(size_t i=0; i<hi_numNodes; i++) {
    if(is_low_order[i] && hi_nodeIsOwned[i])
      lo_nodeIsOwned[hi_to_lo_map[i]]=true;
  }  

  // Translate lo_elemToNode to a lo local index
  for(size_t i=0; i<numElem; i++)
    for(size_t j=0; j<lo_nperel; j++) 
      lo_elemToNode(i,j) = hi_to_lo_map[lo_elemToNode(i,j)];    

  
  // Check for the [E|T]petra column map ordering property, namely LIDs for owned nodes should all appear first.  
  // Since we're injecting from the higher-order mesh, it should be true, but we should add an error check & throw in case.
  bool map_ordering_test_passed=true;
  for(size_t i=0; i<lo_numNodes-1; i++)
    if(!lo_nodeIsOwned[i] && lo_nodeIsOwned[i+1]) 
      map_ordering_test_passed=false;
  
  if(!map_ordering_test_passed)
    throw std::runtime_error("MueLu::MueLuIntrepid::BuildLoElemToNode failed map ordering test");

}


/*********************************************************************************************************/
  template <class LocalOrdinal, class GlobalOrdinal, class Node> 
  void GenerateColMapFromImport(const Xpetra::Import<LocalOrdinal,GlobalOrdinal, Node> & hi_importer,const std::vector<LocalOrdinal> &hi_to_lo_map,const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> & lo_domainMap, const size_t & lo_columnMapLength, RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & lo_columnMap) {
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Node NO;
  typedef Xpetra::Map<LO,GO,NO> Map;
  typedef Xpetra::Vector<GO,LO,GO,NO> GOVector;

  GO go_invalid = Teuchos::OrdinalTraits<GO>::invalid();
  LO lo_invalid = Teuchos::OrdinalTraits<LO>::invalid();

  RCP<const Map> hi_domainMap = hi_importer.getSourceMap();
  RCP<const Map> hi_columnMap = hi_importer.getTargetMap();
  // Figure out the GIDs of my non-owned P1 nodes
  // HOW: We can build a GOVector(domainMap) and fill the values with either invalid() or the P1 domainMap.GID() for that guy.
  // Then we can use A's importer to get a GOVector(colMap) with that information.

  // NOTE: This assumes rowMap==colMap and [E|T]petra ordering of all the locals first in the colMap
  RCP<GOVector> dvec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(hi_domainMap);
  ArrayRCP<GO> dvec_data = dvec->getDataNonConst(0);
  for(size_t i=0; i<lo_domainMap.getNodeNumElements(); i++) {
    if(hi_to_lo_map[i]!=lo_invalid) dvec_data[i] = hi_domainMap->getGlobalElement(hi_to_lo_map[i]);
    else dvec_data[i] = go_invalid;
  }
  RCP<GOVector> cvec = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(hi_columnMap,true);
  cvec->doImport(*dvec,hi_importer,Xpetra::ADD);

  // Generate the lo_columnMap
  // HOW: We can use the local hi_to_lo_map from the GID's in cvec to generate the non-contiguous colmap ids.
  Array<GO> lo_col_data(lo_columnMapLength);
  ArrayRCP<GO> cvec_data = cvec->getDataNonConst(0);
  for(size_t i=0,idx=0; i<hi_columnMap->getNodeNumElements(); i++) {
    if(hi_to_lo_map[i]!=lo_invalid) {
      lo_col_data[idx] = cvec_data[i];
      idx++;
    }
  }  
  
  lo_columnMap = Xpetra::MapFactory<LO,GO,NO>::Build(lo_domainMap.lib(),Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),lo_col_data(),lo_domainMap.getIndexBase(),lo_domainMap.getComm());
}
		  
}//end MueLu::MueLuIntrepid namespace



/*********************************************************************************************************/

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GenerateLinearCoarsening_pn_kirby_to_p1(const Intrepid2::FieldContainer<LocalOrdinal> & hi_elemToNode, 
														 const std::vector<bool> & hi_nodeIsOwned,
														 const Intrepid2:: FieldContainer<double> hi_DofCoords,
														 const std::vector<size_t> &lo_node_in_hi,
														 const Intrepid2::Basis<double,Intrepid2::FieldContainer<double> > &lo_basis,
														 const std::vector<LocalOrdinal> & hi_to_lo_map,
														 const Teuchos::RCP<const Map> & lo_colMap, 
														 const Teuchos::RCP<const Map> & lo_domainMap, 
														 const Teuchos::RCP<const Map> & hi_map,
														 Teuchos::RCP<Matrix>& P) const{
  typedef Intrepid2::FieldContainer<double> FC;
  // Evaluate the linear basis functions at the Pn nodes
  size_t numFieldsHi = hi_elemToNode.dimension(1);
  size_t numFieldsLo = lo_basis.getCardinality();
  FC LoValues_at_HiDofs(numFieldsLo,numFieldsHi);
  lo_basis.getValues(LoValues_at_HiDofs, hi_DofCoords, Intrepid2::OPERATOR_VALUE);

  typedef typename Teuchos::ScalarTraits<SC>::halfPrecision SClo;
  typedef typename Teuchos::ScalarTraits<SClo>::magnitudeType MT;
  MT effective_zero = Teuchos::ScalarTraits<MT>::eps();

  // Allocate P
  P = rcp(new CrsMatrixWrap(hi_map,lo_colMap,0)); //FIXLATER: Need faster fill
  RCP<CrsMatrix> Pcrs   = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

  // Slow-ish fill
  size_t Nelem=hi_elemToNode.dimension(0);  
  std::vector<bool> touched(hi_map->getNodeNumElements(),false);
  Teuchos::Array<GO> col_gid(1); 
  Teuchos::Array<SC> val(1);
  for(size_t i=0; i<Nelem; i++) {
    for(size_t j=0; j<numFieldsHi; j++) {
      LO row_lid = hi_elemToNode(i,j);
      GO row_gid = hi_map->getGlobalElement(row_lid);
      if(hi_nodeIsOwned[row_lid] && !touched[row_lid]) {
	for(size_t k=0; k<numFieldsLo; k++) {
	  // Get the local id in P1's column map
	  LO col_lid = hi_to_lo_map[hi_elemToNode(i,lo_node_in_hi[k])];
	  col_gid[0] = {lo_colMap->getGlobalElement(col_lid)};
	  val[0]     = LoValues_at_HiDofs(k,j);
	  
	  // Skip near-zeros
	  if(Teuchos::ScalarTraits<SC>::magnitude(val[0]) >= effective_zero)
	    P->insertGlobalValues(row_gid,col_gid(),val());
	}
	touched[row_lid]=true;
      }
    }
  }
  P->fillComplete(lo_domainMap,hi_map);
}



  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("ipc: hi basis");
    SET_VALID_ENTRY("ipc: lo basis");
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
    
    validParamList->set< RCP<const FactoryBase> >("Nullspace",      Teuchos::null, "Generating factory of the nullspace");
    validParamList->set< RCP<const FactoryBase> >("ipc: element to node map",          Teuchos::null, "Generating factory of the element to node map");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "A");
    Input(fineLevel, "Nullspace");
    Input(fineLevel, "ipc: element to node map");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level &coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "P Coarsening", coarseLevel);

    std::string levelIDs = toString(coarseLevel.GetLevelID());

    const std::string prefix = "MueLu::IntrepidPCoarsenFactory(" + levelIDs + "): ";
    typedef Intrepid2::FieldContainer<LO> FCi;
    //    typedef Intrepid2::FieldContainer<SC> FC;
    //    typedef Intrepid2::Basis<SC,FC> Basis;
    // NOTE: This is hardwired to double on purpose.  See the note below.
    typedef Intrepid2::FieldContainer<double> FC;
    typedef Intrepid2::Basis<double,FC> Basis;

    GO go_invalid = Teuchos::OrdinalTraits<GO>::invalid();
    LO lo_invalid = Teuchos::OrdinalTraits<LO>::invalid();

    // Level Get
    RCP<Matrix> A     = Get< RCP<Matrix> >(fineLevel, "A");
    Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Acrs = dynamic_cast<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(*A);

    if (restrictionMode_) {
      SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);
      A = Utilities::Transpose(*A, true); // build transpose of A explicitly
    }

    // Build final prolongator
    RCP<Matrix> finalP;

    // Reuse pattern if available
    RCP<ParameterList> APparams = rcp(new ParameterList);
    if (coarseLevel.IsAvailable("AP reuse data", this)) {
      GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous AP data" << std::endl;

      APparams = coarseLevel.Get< RCP<ParameterList> >("AP reuse data", this);

      if (APparams->isParameter("graph"))
        finalP = APparams->get< RCP<Matrix> >("graph");
    }
    const ParameterList& pL = GetParameterList();

    /*******************/
    // FIXME LATER: Allow these to be manually specified instead of Intrepid
    // Get the Intrepid bases
    // NOTE: To make sure Stokhos works we only instantiate these guys with double.  There's a lot
    // of stuff in the guts of Intrepid2 that doesn't play well with Stokhos as of yet.
    RCP<Basis> hi_basis = MueLuIntrepid::BasisFactory<double>(pL.get<std::string>("ipc: hi basis"));
    RCP<Basis> lo_basis = MueLuIntrepid::BasisFactory<double>(pL.get<std::string>("ipc: lo basis"));

    // Get reference coordinates and the lo-to-hi injection list for the reference element
    std::vector<size_t> lo_node_in_hi;
    FC hi_DofCoords;
    MueLuIntrepid::IntrepidGetLoNodeInHi(hi_basis,lo_basis,lo_node_in_hi,hi_DofCoords);

    /*******************/    
    // Get the higher-order element-to-node map 
    const Teuchos::RCP<FCi> Pn_elemToNode = Get<Teuchos::RCP<FCi> >(fineLevel,"ipc: element to node map");

    // Calculate node ownership (the quick and dirty way)
    // NOTE: This exploits two things: 
    //  1) domainMap == rowMap
    //  2) Standard [e|t]petra ordering (namely the local unknowns are always numbered first).  
    // This routine does not work in general.    
    RCP<const Map> rowMap = A->getRowMap();
    RCP<const Map> colMap = Acrs.getColMap();
    RCP<const Map> domainMap = A->getDomainMap();
    int NumProc = rowMap->getComm()->getSize();
    assert(&*rowMap == &*domainMap);
    std::vector<bool> Pn_nodeIsOwned(colMap->getNodeNumElements(),false);
    LO num_owned_rows = 0;
    for(size_t i=0; i<rowMap->getNodeNumElements(); i++) {
      if(rowMap->getGlobalElement(i) == colMap->getGlobalElement(i)) {
	Pn_nodeIsOwned[i] = true;
	num_owned_rows++;
      }
    }

    // Generate lower-order element-to-node map
    FCi P1_elemToNode;
    std::vector<bool> P1_nodeIsOwned;
    std::vector<LO> hi_to_lo_map;
    int P1_numOwnedNodes;
    MueLuIntrepid::BuildLoElemToNode(*Pn_elemToNode,Pn_nodeIsOwned,lo_node_in_hi,P1_elemToNode,P1_nodeIsOwned,hi_to_lo_map,P1_numOwnedNodes);
    assert(hi_to_lo_map.size() == colMap->getNodeNumElements());

    // Generate the P1_domainMap
    // HOW: Since we know how many each proc has, we can use the non-uniform contiguous map constructor to do the work for us
    RCP<const Map> P1_domainMap = MapFactory::Build(rowMap->lib(),Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),P1_numOwnedNodes,rowMap->getIndexBase(),rowMap->getComm());
    Set(coarseLevel, "CoarseMap", P1_domainMap);       

    // Generate the P1_columnMap
    RCP<const Map> P1_colMap;
    if(NumProc==1) P1_colMap = P1_domainMap;
    else MueLuIntrepid::GenerateColMapFromImport<LO,GO,NO>(*Acrs.getCrsGraph()->getImporter(),hi_to_lo_map,*P1_domainMap,P1_nodeIsOwned.size(),P1_colMap);

    // Generate the coarsening
    GenerateLinearCoarsening_pn_kirby_to_p1(*Pn_elemToNode,Pn_nodeIsOwned,hi_DofCoords,lo_node_in_hi,*lo_basis,hi_to_lo_map,P1_colMap,P1_domainMap,A->getRowMap(),finalP);

    // Level Set
    if (!restrictionMode_) {
      // The factory is in prolongation mode
      Set(coarseLevel, "P",             finalP);

      APparams->set("graph", finalP);
      Set(coarseLevel, "AP reuse data", APparams);

      if (IsPrint(Statistics1)) {
        RCP<ParameterList> params = rcp(new ParameterList());
        params->set("printLoadBalancingInfo", true);
        params->set("printCommInfo",          true);
        GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*finalP, "P", params);
      }

    } else {
      // The factory is in restriction mode
      RCP<Matrix> R;
      {
        SubFactoryMonitor m2(*this, "Transpose P", coarseLevel);
        R = Utilities::Transpose(*finalP, true);
      }

      Set(coarseLevel, "R", R);

      if (IsPrint(Statistics1)) {
        RCP<ParameterList> params = rcp(new ParameterList());
        params->set("printLoadBalancingInfo", true);
        params->set("printCommInfo",          true);
        GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*R, "R", params);
      }
    }

  } //Build()

} //namespace MueLu

#endif // MUELU_IPCFACTORY_DEF_HPP

