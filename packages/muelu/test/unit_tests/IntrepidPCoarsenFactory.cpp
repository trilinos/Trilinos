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
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_CreateXpetraPreconditioner.hpp"

#ifdef HAVE_MUELU_INTREPID2
#include "MueLu_IntrepidPCoarsenFactory.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_FieldContainer.hpp"



namespace MueLuTests {

  /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GetLoNodeInHi, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;

    typedef Intrepid2::FieldContainer<MT> FC;

    out << "version: " << MueLu::Version() << std::endl;

    int max_degree=5;

    {
      // QUAD
      RCP<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC> > lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC>());
      RCP<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC> > hi;
      for(int i=0;i<max_degree; i++) {
	hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC>(i,Intrepid2::POINTTYPE_EQUISPACED));

	std::vector<size_t> lo_node_in_hi;
	FC hi_dofCoords;
	MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<MT,FC>(hi,lo,lo_node_in_hi,hi_dofCoords);
	
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
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
    typedef Intrepid2::FieldContainer<MT> FC;

    out << "version: " << MueLu::Version() << std::endl;
    
    // QUAD
    {bool test= rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC> >(MueLu::MueLuIntrepid::BasisFactory<MT>("hgrad_quad_c1")) !=Teuchos::null;TEST_EQUALITY(test,true);}
    {bool test= rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC> >(MueLu::MueLuIntrepid::BasisFactory<MT>("hgrad_quad_c2")) !=Teuchos::null;TEST_EQUALITY(test,true);}
    {bool test= rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC> >(MueLu::MueLuIntrepid::BasisFactory<MT>("hgrad_quad_c3")) !=Teuchos::null;TEST_EQUALITY(test,true);}

  }


  /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildLoElemToNode, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
  #   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
    typedef Intrepid2::FieldContainer<MT> FC;
    typedef Intrepid2::FieldContainer<LO> FCi;

    out << "version: " << MueLu::Version() << std::endl;
    int max_degree=5;

    {
      //QUAD
      // A one element test with Kirby-numbered nodes where the top edge is not owned      
      RCP<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC> > lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC>());
      RCP<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC> > hi;
      for(int degree=2; degree < max_degree; degree++) {
	hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC>(degree,Intrepid2::POINTTYPE_EQUISPACED));
	int Nn = (degree+1)*(degree+1);

	FCi hi_e2n(1,Nn), lo_e2n;
	std::vector<bool> hi_owned(Nn,false),lo_owned;
	std::vector<size_t> lo_node_in_hi;
	std::vector<LO> hi_to_lo_map;
	int lo_numOwnedNodes=0;
	FC hi_dofCoords;
	MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<MT,FC>(hi,lo,lo_node_in_hi,hi_dofCoords);

	for(int i=0; i<Nn; i++) {
	  hi_e2n(0,i)=i;
	  if(i < Nn-(degree+1)) hi_owned[i]=true;
	}

	MueLu::MueLuIntrepid::BuildLoElemToNode(hi_e2n,hi_owned,lo_node_in_hi,lo_e2n,lo_owned,hi_to_lo_map,lo_numOwnedNodes);
	
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
    // This function basically takes an existing domain->column importer (the "hi order" guy) and using the hi_to_lo_map, 
    // generates "lo order" version.  The domain map is already given to us here, so we just need to make tha column one.

    // We'll test this by starting with some linearly distributed map for the "hi" domain map and a duplicated map as the "hi" column map.
    // We then do a "lo" domain map, which just grabs the first GID on each proc.
    GO numGlobalElements =100;
    Teuchos::Array<GO> hi_Cols(numGlobalElements);
    for(size_t i=0; i<(size_t)numGlobalElements; i++)
      hi_Cols[i] = i;

    RCP<Map> hi_domainMap   = MapFactory::Build(lib,numGlobalElements,0,comm);
    RCP<Map> hi_colMap      = MapFactory::Build(lib,gst_invalid,hi_Cols(),0,comm);
    RCP<Import> hi_importer = ImportFactory::Build(hi_domainMap,hi_colMap);

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
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TestPseudoPoisson(Teuchos::FancyOStream &out, int num_nodes, int degree, std::vector<Scalar> &pn_gold_in, std::vector<Scalar> &pn_gold_out,const std::string & hi_basis)
  {
  #   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO; 
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;

    out << "version: " << MueLu::Version() << std::endl;

    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO gst_invalid = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    GO lo_invalid = Teuchos::OrdinalTraits<LO>::invalid();
    int MyPID = comm->getRank();

    // Setup Levels
    Level fineLevel, coarseLevel;
    test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
    fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
    coarseLevel.SetFactoryManager(Teuchos::null);

    // Build a pseudo-poisson test matrix
    Intrepid2::FieldContainer<LocalOrdinal> elem_to_node;
    RCP<Matrix> A = test_factory::Build1DPseudoPoissonHigherOrder(num_nodes,degree,elem_to_node,lib);
    fineLevel.Set("A",A);
    fineLevel.Set("ipc: element to node map",rcp(&elem_to_node,false));

    // only one NS vector 
    LocalOrdinal NSdim = 1;
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),NSdim);
    nullSpace->setSeed(846930886);
    nullSpace->randomize();
    fineLevel.Set("Nullspace",nullSpace);

    // ParameterList
    ParameterList Params;
    Params.set("ipc: hi basis",hi_basis);
    Params.set("ipc: lo basis","hgrad_line_c1");

    // Build P
    RCP<MueLu::IntrepidPCoarsenFactory<SC,LO,GO,NO> > IPCFact = rcp(new MueLu::IntrepidPCoarsenFactory<SC,LO,GO,NO>());
    IPCFact->SetParameterList(Params);
    coarseLevel.Request("P",IPCFact.get());  // request Ptent
    coarseLevel.Request("Nullspace",IPCFact.get());
    coarseLevel.Request("CoarseMap",IPCFact.get());
    coarseLevel.Request(*IPCFact);
    IPCFact->Build(fineLevel,coarseLevel);

    // Get P
    RCP<Matrix> P;
    coarseLevel.Get("P",P,IPCFact.get());
    RCP<CrsMatrix> Pcrs   = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();
    if(!MyPID) printf("P size = %d x %d\n",(int)P->getRangeMap()->getGlobalNumElements(),(int)P->getDomainMap()->getGlobalNumElements());

    // Build serial comparison maps
    GO pn_num_global_dofs = A->getRowMap()->getGlobalNumElements();
    GO pn_num_serial_elements = !MyPID ? pn_num_global_dofs : 0;
    RCP<Map> pn_SerialMap = MapFactory::Build(lib,pn_num_global_dofs,pn_num_serial_elements,0,comm);

    GO p1_num_global_dofs = P->getDomainMap()->getGlobalNumElements();
    GO p1_num_serial_elements = !MyPID ? p1_num_global_dofs : 0;
    RCP<Map> p1_SerialMap = MapFactory::Build(lib, p1_num_global_dofs,p1_num_serial_elements,0,comm);

    RCP<Export> p1_importer = ExportFactory::Build(p1_SerialMap,P->getDomainMap());
    RCP<Export> pn_importer = ExportFactory::Build(A->getRowMap(),pn_SerialMap);

    // Allocate some vectors
    RCP<Vector> s_InVec = VectorFactory::Build(p1_SerialMap);
    RCP<Vector> p_InVec = VectorFactory::Build(P->getDomainMap());
    RCP<Vector> s_OutVec = VectorFactory::Build(pn_SerialMap);
    RCP<Vector> s_codeOutput = VectorFactory::Build(pn_SerialMap);
    RCP<Vector> p_codeOutput = VectorFactory::Build(A->getRowMap());


    // Fill serial GOLD vecs on Proc 0
    if(!MyPID) {
      for(size_t i=0; i<(size_t)pn_gold_in.size(); i++)
	s_InVec->replaceLocalValue(i,pn_gold_in[i]);

      for(size_t i=0; i<(size_t)pn_gold_out.size(); i++)
	s_OutVec->replaceLocalValue(i,pn_gold_out[i]);
    }

    // Migrate input data
    p_InVec->doExport(*s_InVec,*p1_importer,Xpetra::ADD);

    // Apply P
    P->apply(*p_InVec,*p_codeOutput);

    // Migrate Output data
    s_codeOutput->doExport(*p_codeOutput,*pn_importer,Xpetra::ADD);

    // Compare vs. GOLD
    s_codeOutput->update(-1.0,*s_OutVec,1.0);
    Teuchos::Array<MT> norm2(1);
    s_codeOutput->norm2(norm2());
    

    if(!MyPID) printf("Diff norm = %10.4e\n",norm2[0]);

  }


 /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildP_PseudoPoisson_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
    // GOLD vector collection
    std::vector<Scalar> p2_gold_in = {0,1,2,3,4,5,6,7,8,9};
    std::vector<Scalar> p2_gold_out= {0,1,2,3,4,5,6,7,8,9,
				  0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5};
    TestPseudoPoisson<Scalar,LocalOrdinal,GlobalOrdinal,Node>(out,p2_gold_in.size(),2,p2_gold_in,p2_gold_out,std::string("hgrad_line_c2"));
  }


/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildP_PseudoPoisson_p3, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
    // GOLD vector collection
    size_t total_num_points=10;
    int degree=3;
    std::vector<Scalar> p3_gold_in(total_num_points);
    std::vector<Scalar> p3_gold_out(total_num_points + (total_num_points-1) *(degree-1));
    for(size_t i=0; i<total_num_points; i++) {
      p3_gold_in[i] = i;
      p3_gold_out[i] = i;
    }

    size_t idx=total_num_points;
    for(size_t i=0; i<total_num_points-1; i++) {
      for(size_t j=0; j<(size_t)degree-1; j++) {
	p3_gold_out[idx] = i + ((double)j+1)/degree;
	idx++;
      }
    }

    TestPseudoPoisson<Scalar,LocalOrdinal,GlobalOrdinal,Node>(out,p3_gold_in.size(),3,p3_gold_in,p3_gold_out,std::string("hgrad_line_c3"));
  }

/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildP_PseudoPoisson_p4, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {

    // GOLD vector collection
    size_t total_num_points=10;
    int degree=4;
    std::vector<Scalar> gold_in(total_num_points);
    std::vector<Scalar> gold_out(total_num_points + (total_num_points-1) *(degree-1));
    for(size_t i=0; i<total_num_points; i++) {
      gold_in[i] = i;
      gold_out[i] = i;
    }

    size_t idx=total_num_points;
    for(size_t i=0; i<total_num_points-1; i++) {
      for(size_t j=0; j<(size_t)degree-1; j++) {
	gold_out[idx] = i + ((double)j+1)/degree;
	idx++;
      }
    }

    TestPseudoPoisson<Scalar,LocalOrdinal,GlobalOrdinal,Node>(out,gold_in.size(),4,gold_in,gold_out,std::string("hgrad_line_c4"));
  }


/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=2;
    std::string hi_basis("hgrad_line_c2");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    Intrepid2::FieldContainer<LocalOrdinal> elem_to_node;
    RCP<Matrix> A = test_factory::Build1DPseudoPoissonHigherOrder(num_nodes,degree,elem_to_node,lib);

    // Normalized RHS
    RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
    RHS1->setSeed(846930886);
    RHS1->randomize();
    Teuchos::Array<MT> norms(1);
    RHS1->norm2(norms);
    RHS1->scale(1/norms[0]);
    
    // Zero initial guess
    RCP<MultiVector> X1   = MultiVectorFactory::Build(A->getRowMap(), 1);
    X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

    // ParameterList
    ParameterList Params, level0;
    Params.set("multigrid algorithm","pcoarsen");
    Params.set("ipc: hi basis",hi_basis);
    Params.set("ipc: lo basis","hgrad_line_c1");
    Params.set("verbosity","high");
    Params.set("max levels",2);
    Params.set("coarse: max size",100);
    level0.set("ipc: element to node map",rcp(&elem_to_node,false));
    Params.set("level 0",level0);
      

    // Build hierarchy
    RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,Params);
  }

/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p3, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=3;
    std::string hi_basis("hgrad_line_c3");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    Intrepid2::FieldContainer<LocalOrdinal> elem_to_node;
    RCP<Matrix> A = test_factory::Build1DPseudoPoissonHigherOrder(num_nodes,degree,elem_to_node,lib);

    // Normalized RHS
    RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
    RHS1->setSeed(846930886);
    RHS1->randomize();
    Teuchos::Array<MT> norms(1);
    RHS1->norm2(norms);
    RHS1->scale(1/norms[0]);
    
    // Zero initial guess
    RCP<MultiVector> X1   = MultiVectorFactory::Build(A->getRowMap(), 1);
    X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

    // ParameterList
    ParameterList Params, level0;
    Params.set("multigrid algorithm","pcoarsen");
    Params.set("ipc: hi basis",hi_basis);
    Params.set("ipc: lo basis","hgrad_line_c1");
    Params.set("verbosity","high");
    Params.set("max levels",2);
    Params.set("coarse: max size",100);
    level0.set("ipc: element to node map",rcp(&elem_to_node,false));
    Params.set("level 0",level0);
      

    // Build hierarchy
    RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,Params);
  }

/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p4, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=4;
    std::string hi_basis("hgrad_line_c4");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    Intrepid2::FieldContainer<LocalOrdinal> elem_to_node;
    RCP<Matrix> A = test_factory::Build1DPseudoPoissonHigherOrder(num_nodes,degree,elem_to_node,lib);

    // Normalized RHS
    RCP<MultiVector> RHS1 = MultiVectorFactory::Build(A->getRowMap(), 1);
    RHS1->setSeed(846930886);
    RHS1->randomize();
    Teuchos::Array<MT> norms(1);
    RHS1->norm2(norms);
    RHS1->scale(1/norms[0]);
    
    // Zero initial guess
    RCP<MultiVector> X1   = MultiVectorFactory::Build(A->getRowMap(), 1);
    X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

    // ParameterList
    ParameterList Params, level0;
    Params.set("multigrid algorithm","pcoarsen");
    Params.set("ipc: hi basis",hi_basis);
    Params.set("ipc: lo basis","hgrad_line_c1");
    Params.set("verbosity","high");
    Params.set("max levels",2);
    Params.set("coarse: max size",100);
    level0.set("ipc: element to node map",rcp(&elem_to_node,false));
    Params.set("level 0",level0);
      

    // Build hierarchy
    RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,Params);
  }


  /*********************************************************************************************************************/
#  define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,GetLoNodeInHi,Scalar,LO,GO,Node)  \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,BasisFactory,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,BuildLoElemToNode,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,GenerateColMapFromImport,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,BuildP_PseudoPoisson_p2,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,BuildP_PseudoPoisson_p3,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,BuildP_PseudoPoisson_p4,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p2, Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p3, Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p4, Scalar, LO,GO,Node)


#include <MueLu_ETI_4arg.hpp>


} // namespace MueLuTests
#endif
