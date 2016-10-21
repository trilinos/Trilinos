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
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_Utilities.hpp"

// Intrepid Headers

//Intrepid_HGRAD_HEX_C1_FEM.hpp
//Intrepid_HGRAD_HEX_C2_FEM.hpp
//Intrepid_HGRAD_HEX_Cn_FEM.hpp
//Intrepid_HGRAD_HEX_I2_FEM.hpp
//Intrepid_HGRAD_LINE_C1_FEM.hpp
//Intrepid_HGRAD_LINE_Cn_FEM.hpp
//Intrepid_HGRAD_LINE_Cn_FEM_JACOBI.hpp
//Intrepid_HGRAD_POLY_C1_FEM.hpp
//Intrepid_HGRAD_PYR_C1_FEM.hpp
//Intrepid_HGRAD_PYR_I2_FEM.hpp
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
//#include Intrepid_HGRAD_QUAD_C2_FEM.hpp
#include "Intrepid_HGRAD_QUAD_Cn_FEM.hpp"
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
Teuchos::RCP<Intrepid::Basis<Scalar,Intrepid::FieldContainer<Scalar> > >  BasisFactory(const std::string & name) {
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

    std::cout<<"String = "<<name<<" ";
    std::cout<<"Parsed as = ("<<deriv<<","<<el<<","<<poly<<","<<degree<<")"<<std::endl;

    
    // FIX: Allow for alternative point types for Kirby elements
    if(deriv=="hgrad" && el=="quad" && poly=="c"){
      if(degree==1) return rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Scalar,Intrepid::FieldContainer<Scalar> >());
      else          return rcp(new Intrepid::Basis_HGRAD_QUAD_Cn_FEM<Scalar,Intrepid::FieldContainer<Scalar> >(degree,Intrepid::POINTTYPE_EQUISPACED));
    }

    // Error out
    throw std::runtime_error(myerror);
    return Teuchos::null;
}

/*********************************************************************************************************/
template <class Scalar, class ArrayScalar>
void IntrepidGetLoNodeInHi(const Teuchos::RCP<Intrepid::Basis<Scalar,ArrayScalar> > &hi_basis,
			   const Teuchos::RCP<Intrepid::Basis<Scalar,ArrayScalar> > &lo_basis,
			   std::vector<size_t> & lo_node_in_hi,
			   ArrayScalar & hi_DofCoords) {
  
  // Figure out which unknowns in hi_basis correspond to nodes on lo_basis. This varies by element type.
  size_t degree         = hi_basis->getDegree();
  lo_node_in_hi.resize(0);
  RCP<Intrepid::DofCoordsInterface<ArrayScalar> > hi_dci;
  if(!rcp_dynamic_cast<Intrepid::Basis_HGRAD_QUAD_Cn_FEM<Scalar,ArrayScalar> >(hi_basis).is_null()) {
    // HGRAD QUAD Cn: Numbering as per the Kirby convention (straight across, bottom to top) 
    lo_node_in_hi.insert(lo_node_in_hi.end(),{0,degree, (degree+1)*(degree+1)-1, degree*(degree+1)});
    hi_dci = rcp_dynamic_cast<Intrepid::Basis_HGRAD_QUAD_Cn_FEM<Scalar,ArrayScalar> >(hi_basis);
  }
  else
    throw std::runtime_error("IntrepidPCoarsenFactory: Unknown element type");
  
  // Get coordinates of the hi_basis dof's
  hi_DofCoords.resize(hi_basis->getCardinality(),hi_basis->getBaseCellTopology().getDimension());
  hi_dci->getDofCoords(hi_DofCoords);
}


template <class LocalOrdinal>
void BuildLoElemToNode(const Intrepid::FieldContainer<LocalOrdinal> & hi_elemToNode,
		       const std::vector<bool> & hi_nodeIsOwned,
		       const std::vector<size_t> & lo_node_in_hi,
		       Intrepid::FieldContainer<LocalOrdinal> & lo_elemToNode,
		       std::vector<bool> & lo_nodeIsOwned,
		       std::vector<LocalOrdinal> & hi_to_lo_map,
		       int & lo_numOwnedNodes) {
  typedef LocalOrdinal LO;
  using Teuchos::RCP;

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
    if(hi_nodeIsOwned[i])
      lo_nodeIsOwned[hi_to_lo_map[i]]=true;
  }  

  // Translate lo_elemToNode to a lo local index
  for(size_t i=0; i<numElem; i++)
    for(size_t j=0; j<lo_nperel; j++) 
      lo_elemToNode(i,j) = hi_to_lo_map[lo_elemToNode(i,j)];    
}


		   



}//end MueLu::MueLuIntrepid namespace





/*********************************************************************************************************/

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GenerateLinearCoarsening_pn_kirby_to_p1(const Intrepid::FieldContainer<LocalOrdinal> & Pn_elemToNode, 
														 const std::vector<bool> & Pn_nodeIsOwned,
														 const Teuchos::RCP<Intrepid::Basis<Scalar,Intrepid::FieldContainer<Scalar> > > &PnBasis_rcp,
														 const Teuchos::RCP<Intrepid::Basis<Scalar,Intrepid::FieldContainer<Scalar> > > &P1Basis_rcp, 
														 const Teuchos::RCP<const Map> & P1_colMap, 
														 const Teuchos::RCP<const Map> & P1_domainMap, 
														 const Teuchos::RCP<const Map> & Pn_map,
														 Teuchos::RCP<Matrix>& P) const{
  typedef Intrepid::FieldContainer<Scalar> FC;

  // NOTE: For all of the basis functions we care about, the Pn Bais is a DofCoordsInterface.  It's just that "Basis" is not.
  Teuchos::RCP<Intrepid::DofCoordsInterface<FC> > PnDCI_rcp;
  // Sanity checks
  assert(Pn_elemToNode.dimension(1) == PnBasis_rcp->getCardinality());
  int degree = PnBasis_rcp->getDegree();

  // Figure out which unknowns in Pn correspond to nodes on P1.  This varies by element type 
  // NOTE: It would be nice if we could somehow integrate this into the Intrepid factory stuff.
  std::vector<int> p1_node_in_pn;
  if(!rcp_dynamic_cast<Intrepid::Basis_HGRAD_QUAD_Cn_FEM<SC,FC> >(PnBasis_rcp).is_null()) {
    // HGRAD QUAD Cn: Numbering as per the Kirby convention (straight across, bottom to top) 
    p1_node_in_pn.insert(p1_node_in_pn.end(),{0,degree, (degree+1)*(degree+1)-1, degree*(degree+1)});
    PnDCI_rcp = rcp_dynamic_cast<Intrepid::Basis_HGRAD_QUAD_Cn_FEM<SC,FC> >(PnBasis_rcp);
  }
  else
    throw std::runtime_error("IntrepidPCoarsenFactory: Unknown element type");


  // Get the reference coordinates for the Pn element  
  int numFieldsPn = PnBasis_rcp->getCardinality();
  int spaceDim    = PnBasis_rcp->getBaseCellTopology().getDimension();
  Intrepid::FieldContainer<SC> PnDofCoords(numFieldsPn,spaceDim);
  PnDCI_rcp->getDofCoords(PnDofCoords);

  // Evaluate the linear basis functions at the Pn nodes
  size_t numFieldsP1 = P1Basis_rcp->getCardinality();
  Intrepid::FieldContainer<SC> P1Values_at_PnDofs(numFieldsP1,numFieldsPn);
  P1Basis_rcp->getValues(P1Values_at_PnDofs, PnDofCoords, Intrepid::OPERATOR_VALUE);

  // Allocate P
  P = rcp(new CrsMatrixWrap(Pn_map, 0));
  RCP<CrsMatrix> Pcrs   = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

  // Extra slow fill

  // CMS: For this to work, P1_map has to be the *colmap* of P, we
  // need a separate domain map.
  size_t Nelem=Pn_elemToNode.dimension(0);  
  std::vector<bool> touched(Pn_map->getNodeNumElements(),false);
  Teuchos::Array<GO> col_gid(1); 
  Teuchos::Array<SC> val(1);
  for(size_t i=0; i<Nelem; i++) {
    for(int j=0; j<numFieldsPn; j++) {
      LO row_lid = Pn_elemToNode(i,j);
      GO row_gid = Pn_map->getGlobalElement(row_lid);
      for(size_t k=0; k<numFieldsP1; k++) {
        LO col_lid = Pn_elemToNode(i,p1_node_in_pn[k]);
	col_gid[0]= {P1_colMap->getGlobalElement(col_lid)};
	val[0] = P1Values_at_PnDofs(k,j);
	if(Pn_nodeIsOwned[row_lid] && !touched[row_lid]) {
	  P->insertGlobalValues(row_gid,col_gid(),val());
	  touched[row_lid]=true;
	}
      }
    }
  }
  P->fillComplete(P1_domainMap,Pn_map);
}



  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("ipc: hi basis");
    SET_VALID_ENTRY("ipc: lo basis");
    SET_VALID_ENTRY("inc: element to node map");
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void IntrepidPCoarsenFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "A");
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
    typedef Intrepid::FieldContainer<SC> FC;
    typedef Intrepid::FieldContainer<LO> FCi;
    typedef Intrepid::Basis<SC,FC> Basis;


    // Level Get
    RCP<Matrix> A     = Get< RCP<Matrix> >(fineLevel, "A");
    Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Acrs = dynamic_cast<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>&>(*A);

    if (restrictionMode_) {
      SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);
      A = Utilities::Transpose(*A, true); // build transpose of A explicitely
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
    // FIXME: Allow these to be manually specified instead of intrepid

    // Get the Intrepid bases
    RCP<Basis> hi_basis = MueLuIntrepid::BasisFactory<Scalar>(pL.get<std::string>("inc: hi basis"));
    RCP<Basis> lo_basis = MueLuIntrepid::BasisFactory<Scalar>(pL.get<std::string>("inc: lo basis"));

    // Get reference coordinates and the lo-to-hi injection list for the reference element
    std::vector<size_t> lo_node_in_hi;
    FC hi_DofCoords;
    MueLuIntrepid::IntrepidGetLoNodeInHi(hi_basis,lo_basis,lo_node_in_hi,hi_DofCoords);

    /*******************/    
    // Get the higher-order element-to-node map
    const Teuchos::RCP<FCi> &Pn_elemToNode = pL.get<Teuchos::RCP<FCi> >("inc: element to node map");

    // Calculate node ownership (the quick and dirty way)
    // NOTE: This exploits two things: 
    //  1) domainMap == rowMap
    //  2) Standard [e|t]petra ordering (namely the local unknowns are always numbered first).  
    // This routine does not work in general.
    RCP<const Map> rowMap = A->getRowMap();
    RCP<const Map> colMap = Acrs.getColMap();
    RCP<const Map> domainMap = A->getDomainMap();
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

#if 0
 BuildLoElemToNode(const Intrepid::FieldContainer<LocalOrdinal> & hi_elemToNode,
		       const std::vector<bool> & hi_nodeIsOwned,
		       const std::vector<size_t> & lo_node_in_hi,
		       Intrepid::FieldContainer<LocalOrdinal> & lo_elemToNode,
		       std::vector<bool> & lo_nodeIsOwned,
		       std::vector<LocalOrdinal> & hi_to_lo_map,
		   int & lo_numOwnedNodes);
#endif

    // Count the P1 owned nodes
    
        
    // Generate p1_map
    // NOTE: We need two maps here, a colmap and a domain map. 
    // The below algorithm assumes that A->getDomainMap() == A->getRowMap().
    // 1) The domain map can be generated by counting the number of owned LowOrder nodes in the RowMap and using the
    //    non-uniform contiguous map constructor.
    // 2) The LowOrder colmap can be used by building a Vector<GO> and filling that with either invalid() or the GID generated in (1).
    //    At this point we can use A's importer to get a Vector<GO> of the GID's corresponding to the LowOrder columns in A.
    // 3) With the Vector<GO> from (2) we can use the non-uniform non-contiguous constructor to generate the ColMap for P.

    // Step 1: Generate the P1 domain map
    // NTS: This is fixed for non-strided at the moment.  At some point, we need to expand this boy to handle multiple PDEs per node.

    // CMS: This is the wrong size, at current

    //    RCP<const Map> P1_domainMap = MapFactory::Build(rowMap->lib(),Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
    //						    num_owned_rows,rowMap->getIndexBase(),rowMap->getComm());


    // Step 2: LowOrder colmap can be used by building a Vector<GO> and filling that with either invalid() or the GID generated in (1).
    // At this point we can use A's importer to get a Vector<GO> of the GID's corresponding to the LowOrder columns in A.
    //    RCP<GOVector> pn_domainvec = VectorFactory::Build(domainMap,true);
    //    RCP<GOVector> pn_columnvec = VectorFactory::Build(colMap,true);
    //   for(size_t i=0; i<domainMap->getNodeNumElements(); i++) {
    //     if(Pn_nodeIsOwned[i]
    //   }

    Teuchos::RCP<Map> P1_colMap, P1_domainMap;
#if 0
 RCP<const Map > reducedMap = MapFactory::Build( A->getRowMap()->lib(),
                                                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                    gidsToImport, indexBase, A->getRowMap()->getComm()    );
#endif
    Set(coarseLevel, "CoarseMap", P1_domainMap);





    // Generate the coarsening
    GenerateLinearCoarsening_pn_kirby_to_p1(*Pn_elemToNode,Pn_nodeIsOwned,hi_basis,lo_basis,P1_colMap,P1_domainMap,A->getRowMap(),finalP);

#if 0
GenerateLinearCoarsening_pn_kirby_to_p1(const Intrepid::FieldContainer<int> & Pn_elemToNode, 
					const std::vector<bool> & Pn_nodeIsOwned,
					const Teuchos::RCP<Intrepid::Basis<Scalar,Intrepid::FieldContainer<Scalar> > > &PnBasis_rcp,
					const Teuchos::RCP<Intrepid::Basis<Scalar,Intrepid::FieldContainer<Scalar> > > &P1Basis_rcp, 
					const Teuchos::RCP<const Map> & P1_map, 
					const Teuchos::RCP<const Map> & Pn_map,
					Teuchos::RCP<Matrix>& P) const;
#endif




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

