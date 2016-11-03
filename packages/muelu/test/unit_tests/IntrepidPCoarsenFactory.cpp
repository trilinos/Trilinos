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
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_IntrepidPCoarsenFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid_FieldContainer.hpp"


namespace MueLuTests {

  /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GetLoNodeInHi, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef Intrepid::FieldContainer<SC> FC;

    out << "version: " << MueLu::Version() << std::endl;

    int max_degree=5;

    {
      // QUAD
      RCP<Intrepid::Basis_HGRAD_QUAD_C1_FEM<SC,FC> > lo = rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<SC,FC>());
      RCP<Intrepid::Basis_HGRAD_QUAD_Cn_FEM<SC,FC> > hi;
      for(int i=0;i<max_degree; i++) {
	hi = rcp(new Intrepid::Basis_HGRAD_QUAD_Cn_FEM<SC,FC>(i,Intrepid::POINTTYPE_EQUISPACED));

	std::vector<size_t> lo_node_in_hi;
	FC hi_dofCoords;
	MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<SC,FC>(hi,lo,lo_node_in_hi,hi_dofCoords);
	
	TEST_EQUALITY(hi_dofCoords.dimension(0),hi->getCardinality());	
	TEST_EQUALITY((size_t)hi_dofCoords.dimension(1),(size_t)hi->getBaseCellTopology().getDimension());
	TEST_EQUALITY(lo_node_in_hi.size(),(size_t)lo->getCardinality());
      }
    }
  }

  /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, BasisFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef Intrepid::FieldContainer<Scalar> FC;

    out << "version: " << MueLu::Version() << std::endl;
    
    // QUAD
    {bool test= rcp_dynamic_cast<Intrepid::Basis_HGRAD_QUAD_C1_FEM<SC,FC> >(MueLu::MueLuIntrepid::BasisFactory<SC>("hgrad_quad_c1")) !=Teuchos::null;TEST_EQUALITY(test,true);}
    {bool test= rcp_dynamic_cast<Intrepid::Basis_HGRAD_QUAD_Cn_FEM<SC,FC> >(MueLu::MueLuIntrepid::BasisFactory<SC>("hgrad_quad_c2")) !=Teuchos::null;TEST_EQUALITY(test,true);}
    {bool test= rcp_dynamic_cast<Intrepid::Basis_HGRAD_QUAD_Cn_FEM<SC,FC> >(MueLu::MueLuIntrepid::BasisFactory<SC>("hgrad_quad_c3")) !=Teuchos::null;TEST_EQUALITY(test,true);}

  }


  /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildLoElemToNode, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
  #   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef Intrepid::FieldContainer<SC> FC;
    typedef Intrepid::FieldContainer<LO> FCi;

    out << "version: " << MueLu::Version() << std::endl;
    int max_degree=5;

    {
      //QUAD
      // A one element test with Kirby-numbered nodes where the top edge is not owned      
      RCP<Intrepid::Basis_HGRAD_QUAD_C1_FEM<SC,FC> > lo = rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<SC,FC>());
      RCP<Intrepid::Basis_HGRAD_QUAD_Cn_FEM<SC,FC> > hi;
      for(int degree=2; degree < max_degree; degree++) {
	hi = rcp(new Intrepid::Basis_HGRAD_QUAD_Cn_FEM<SC,FC>(degree,Intrepid::POINTTYPE_EQUISPACED));
	int Nn = (degree+1)*(degree+1);
#ifdef DEBUG_BUILDLOE2N
	printf("** p=%d (nn=%d)**\n",degree,Nn);
#endif
	FCi hi_e2n(1,Nn), lo_e2n;
	std::vector<bool> hi_owned(Nn,false),lo_owned;
	std::vector<size_t> lo_node_in_hi;
	std::vector<LO> hi_to_lo_map;
	int lo_numOwnedNodes=0;
	FC hi_dofCoords;
	MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<SC,FC>(hi,lo,lo_node_in_hi,hi_dofCoords);

	for(int i=0; i<Nn; i++) {
	  hi_e2n(0,i)=i;
	  if(i < Nn-(degree+1)) hi_owned[i]=true;
	}

	MueLu::MueLuIntrepid::BuildLoElemToNode(hi_e2n,hi_owned,lo_node_in_hi,lo_e2n,lo_owned,hi_to_lo_map,lo_numOwnedNodes);

#ifdef DEBUG_BUILDLOE2N
	printf("hi_to_lo_map = ");
	for(int i=0;i<(int)hi_to_lo_map.size(); i++) {
	  if(hi_to_lo_map[i] == Teuchos::OrdinalTraits<LO>::invalid()) printf("- ");
	  else printf("%d ",hi_to_lo_map[i]);
	}

	printf("\nlo_e2n = ");
	for(int i=0;i<(int)lo_e2n.size(); i++)
	  printf("%d[%d] ",(int)lo_e2n[i],(int)lo_owned[i]);
	printf("\n");
#endif
	
	// Checks
	TEST_EQUALITY(lo_numOwnedNodes,2);

	size_t num_lo_nodes_located=0;
	for(size_t i=0;i<hi_to_lo_map.size(); i++) {
	  if(hi_to_lo_map[i] != Teuchos::OrdinalTraits<LO>::invalid())
	    num_lo_nodes_located++;
	}
	TEST_EQUALITY(lo_owned.size(),num_lo_nodes_located);
      }
    }//end QUAD

  }

  /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,GenerateColMapFromImport, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
  #   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO;

    out << "version: " << MueLu::Version() << std::endl;

    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO gst_invalid = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    GO lo_invalid = Teuchos::OrdinalTraits<LO>::invalid();
    int NumProcs = comm->getSize(); 
    //    int MyPID    = comm->getRank(); 
    // This function basically takes an existing domain->column importer (the "hi order" guy) and using the hi_to_lo_map, generates "lo order" version.  The domain map is already given to
    // us here, so we just need to make tha column one.

    // We'll test this by starting with some linearly distributed map for the "hi" domain map and a duplicated map as the "hi" column map.
    // We then do a "lo" domain map, which just grabs the first GID on each proc.
    GO numGlobalElements =100;
    Teuchos::Array<GO> hi_Cols(numGlobalElements);
    for(size_t i=0; i<(size_t)numGlobalElements; i++)
      hi_Cols[i] = i;


    RCP<Map> hi_domainMap   = MapFactory::Build(lib,numGlobalElements,0,comm);
    RCP<Map> hi_colMap      = MapFactory::Build(lib,gst_invalid,hi_Cols(),0,comm);
    RCP<Import> hi_importer = ImportFactory::Build(hi_domainMap,hi_colMap);
    //Xpetra::Import<LocalOrdinal,GlobalOrdinal, Node> hi_importer(hi_domainMap,hi_colMap);

    Teuchos::Array<GO> lo_Doms(1); lo_Doms[0]=hi_domainMap->getGlobalElement(0);
    RCP<Map> lo_domainMap= MapFactory::Build(lib,gst_invalid,lo_Doms(),0,comm);

    // Get the list of the first GID from each proc (aka the lo_colMap)
    std::vector<GO> send_colgids(1);
    std::vector<GO> recv_colgids(NumProcs); send_colgids[0] = hi_domainMap->getGlobalElement(0);
    comm->gatherAll(1*sizeof(GO),(char*)send_colgids.data(),NumProcs*sizeof(GO),(char*)recv_colgids.data());

    // Build the h2l map
    std::vector<LO> hi_to_lo_map(hi_colMap->getNodeNumElements(),lo_invalid); 
    for(size_t i=0; i<(size_t)NumProcs; i++) 
      hi_to_lo_map[recv_colgids[i]]=i;

    // Import
    size_t lo_columnMapLength = NumProcs; // One col per proc
    RCP<const Map> lo_colMap;
    MueLu::MueLuIntrepid::GenerateColMapFromImport(*hi_importer,hi_to_lo_map,*lo_domainMap,lo_columnMapLength,lo_colMap);
   
    // Final test
    for(size_t i=0; i<lo_colMap->getNodeNumElements(); i++)
      TEST_EQUALITY(recv_colgids[i],lo_colMap->getGlobalElement(i));
 
  }


 /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildP_PseudoPoisson, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
  #   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO; 
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;

    out << "version: " << MueLu::Version() << std::endl;

    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO gst_invalid = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    GO lo_invalid = Teuchos::OrdinalTraits<LO>::invalid();

    // Setup Levels
    Level fineLevel, coarseLevel;
    test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
    fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
    coarseLevel.SetFactoryManager(Teuchos::null);

    // Build a pseudo-poisson test matrix
    Intrepid::FieldContainer<LocalOrdinal> elem_to_node;
    RCP<Matrix> A = test_factory::Build1DPseudoPoissonHigherOrder(10,2,elem_to_node,lib);
    fineLevel.Set("A",A);
    fineLevel.Set("ipc: element to node map",rcp(&elem_to_node,false));

    // only one NS vector 
    LocalOrdinal NSdim = 1;
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),NSdim);
    nullSpace->randomize();
    fineLevel.Set("Nullspace",nullSpace);

    // ParameterList
    ParameterList Params;
    Params.set("ipc: hi basis","hgrad_line_c2");
    Params.set("ipc: lo basis","hgrad_line_c1");

    // Build P
    RCP<MueLu::IntrepidPCoarsenFactory<SC,LO,GO,NO> > IPCFact = rcp(new MueLu::IntrepidPCoarsenFactory<SC,LO,GO,NO>());
    IPCFact->SetParameterList(Params);
    coarseLevel.Request("P",IPCFact.get());  // request Ptent
    coarseLevel.Request("Nullspace",IPCFact.get());
    coarseLevel.Request("CoarseMap",IPCFact.get());
    coarseLevel.Request(*IPCFact);
    IPCFact->Build(fineLevel,coarseLevel);


    // Test P
    RCP<Matrix> P;
    coarseLevel.Get("P",P,IPCFact.get());
    RCP<CrsMatrix> Pcrs   = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();
    printf("CMS: Matrix P is %dx%d\n",(int)P->getRangeMap()->getGlobalNumElements(),(int)P->getDomainMap()->getGlobalNumElements());

    for(size_t i=0; i<P->getRowMap()->getNodeNumElements(); i++) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> values;
      Pcrs->getLocalRowView((LO)i,indices,values);

      //      printf("[%d] Prow[%d] = ",P->getRowMap()->getComm()->getRank(),(int)i);
      printf("[*] Prow[%d] = ",(int)P->getRowMap()->getGlobalElement(i));
      for(size_t j=0; j<(size_t)indices.size(); j++)
	printf("%d(%6.4e) ",(int)P->getColMap()->getGlobalElement(indices[j]),values[j]);
      printf("\n");
    }
    fflush(stdout);

  }


  /*********************************************************************************************************************/
#  define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,GetLoNodeInHi,Scalar,LO,GO,Node)  \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,BasisFactory,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,BuildLoElemToNode,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,GenerateColMapFromImport,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,BuildP_PseudoPoisson,Scalar,LO,GO,Node)


#include <MueLu_ETI_4arg.hpp>


} // namespace MueLuTests
