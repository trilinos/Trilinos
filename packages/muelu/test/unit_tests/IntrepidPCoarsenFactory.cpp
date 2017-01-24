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

#include "MueLu_TestHelpers_HO.hpp"
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
#include "MueLu_IntrepidPCoarsenFactory_def.hpp"   // Why does ETI suddenly decide to hate right here?
#include "Intrepid2_Types.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
//#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
//#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
#include "Kokkos_DynRankView.hpp"
#else
#include "Intrepid2_FieldContainer.hpp"
#endif

namespace MueLuTests {

  /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GetLoNodeInHi, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<MT,typename Node::device_type> FC;
    typedef typename Node::device_type::execution_space ES;
    typedef Intrepid2::Basis<ES,MT,MT> Basis;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
    typedef Intrepid2::Basis<MT,FC> Basis;
#endif

    out << "version: " << MueLu::Version() << std::endl;

    int max_degree=5;

    {
      // QUAD
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
      RCP<Basis> lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES,MT,MT>());
#else
      RCP<Basis> lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC>());
#endif

      for(int i=0;i<max_degree; i++) {
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
        RCP<Basis> hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES,MT,MT>(i,Intrepid2::POINTTYPE_EQUISPACED));
#else
        RCP<Basis> hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC>(i,Intrepid2::POINTTYPE_EQUISPACED));
#endif
        std::vector<size_t> lo_node_in_hi;
        FC hi_dofCoords;

#ifdef HAVE_MUELU_INTREPID2_REFACTOR 
        MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<MT,typename Node::device_type>(hi,lo,lo_node_in_hi,hi_dofCoords);
#else
        MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<MT,FC>(hi,lo,lo_node_in_hi,hi_dofCoords);
#endif  
        TEST_EQUALITY((size_t)hi_dofCoords.dimension(0),(size_t)hi->getCardinality());  
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
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef typename Node::device_type::execution_space ES;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
#endif
    out << "version: " << MueLu::Version() << std::endl;
    
    // QUAD
    int degree;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR 
    {bool test= rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES,MT,MT> >(MueLu::MueLuIntrepid::BasisFactory<MT,ES>("hgrad_quad_c1",degree)) !=Teuchos::null;TEST_EQUALITY(test,true);}
    {bool test= rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES,MT,MT> >(MueLu::MueLuIntrepid::BasisFactory<MT,ES>("hgrad_quad_c2",degree)) !=Teuchos::null;TEST_EQUALITY(test,true);}
    {bool test= rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES,MT,MT> >(MueLu::MueLuIntrepid::BasisFactory<MT,ES>("hgrad_quad_c3",degree)) !=Teuchos::null;TEST_EQUALITY(test,true);}
#else
    {bool test= rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC> >(MueLu::MueLuIntrepid::BasisFactory<MT>("hgrad_quad_c1",degree)) !=Teuchos::null;TEST_EQUALITY(test,true);}
    {bool test= rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC> >(MueLu::MueLuIntrepid::BasisFactory<MT>("hgrad_quad_c2",degree)) !=Teuchos::null;TEST_EQUALITY(test,true);}
    {bool test= rcp_dynamic_cast<Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC> >(MueLu::MueLuIntrepid::BasisFactory<MT>("hgrad_quad_c3",degree)) !=Teuchos::null;TEST_EQUALITY(test,true);}
#endif
  }


  /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildLoElemToNode, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
  #   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<MT,typename Node::device_type> FC;
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
    typedef typename Node::device_type::execution_space ES;
    typedef Intrepid2::Basis<ES,MT,MT> Basis;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
    typedef Intrepid2::FieldContainer<LO> FCi;
    typedef Intrepid2::Basis<MT,FC> Basis;
#endif

    out << "version: " << MueLu::Version() << std::endl;
    int max_degree=5;

    {
      //QUAD
      // A one element test with Kirby-numbered nodes where the top edge is not owned      
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
      RCP<Basis> lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES,MT,MT>());
#else
      RCP<Basis> lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC>());
#endif

      for(int degree=2; degree < max_degree; degree++) {
        int Nn = (degree+1)*(degree+1);
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
        RCP<Basis> hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES,MT,MT>(degree,Intrepid2::POINTTYPE_EQUISPACED));
        FCi hi_e2n("hi_e2n",1,Nn), lo_e2n;
#else
        RCP<Basis> hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC>(degree,Intrepid2::POINTTYPE_EQUISPACED));
        FCi hi_e2n(1,Nn), lo_e2n;
#endif
        std::vector<bool> hi_owned(Nn,false),lo_owned;
        std::vector<size_t> lo_node_in_hi;
        std::vector<LO> hi_to_lo_map;
        int lo_numOwnedNodes=0;
        FC hi_dofCoords;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR 
        MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<MT,typename Node::device_type>(hi,lo,lo_node_in_hi,hi_dofCoords);
#else
        MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<MT,FC>(hi,lo,lo_node_in_hi,hi_dofCoords);
#endif  

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
  /* How this guy works:
     num_p1_nodes - number of nodes in the p=1 mesh
     p1_gold_in   - input vector for the lo_basis 
     p2_gold_in   - output vector of linear interpolation from lo_basis to hi_basis
     
  */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TestPseudoPoisson(Teuchos::FancyOStream &out, int num_p1_nodes, int degree, std::vector<Scalar> &lo_gold_in, std::vector<Scalar> &hi_gold_out,const std::string & hi_basis, const std::string lo_basis = "hgrad_line_c1")
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
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
#else
    typedef Intrepid2::FieldContainer<LO> FCi;
#endif

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
    FCi elem_to_node;
    RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC,LO,GO,NO>(num_p1_nodes,degree,elem_to_node,lib);
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
    Params.set("ipc: lo basis",lo_basis);

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

    // Sanity
    if((int)P->getRangeMap()->getGlobalNumElements()!=(int)hi_gold_out.size())
      throw std::runtime_error("P range size does not match hi_gold_out");
    if((int)P->getDomainMap()->getGlobalNumElements()!=(int)lo_gold_in.size())
      throw std::runtime_error("P domain size does not match lo_gold_in");

    // Build serial comparison maps
    GO hi_num_global_dofs = A->getRowMap()->getGlobalNumElements();
    GO hi_num_serial_elements = !MyPID ? hi_num_global_dofs : 0;
    RCP<Map> hi_SerialMap = MapFactory::Build(lib,hi_num_global_dofs,hi_num_serial_elements,0,comm);

    GO lo_num_global_dofs = P->getDomainMap()->getGlobalNumElements();
    GO lo_num_serial_elements = !MyPID ? lo_num_global_dofs : 0;
    RCP<Map> lo_SerialMap = MapFactory::Build(lib, lo_num_global_dofs,lo_num_serial_elements,0,comm);

    RCP<Export> lo_importer = ExportFactory::Build(lo_SerialMap,P->getDomainMap());
    RCP<Export> hi_importer = ExportFactory::Build(A->getRowMap(),hi_SerialMap);

    // Allocate some vectors
    RCP<Vector> s_InVec = VectorFactory::Build(lo_SerialMap);
    RCP<Vector> p_InVec = VectorFactory::Build(P->getDomainMap());
    RCP<Vector> s_OutVec = VectorFactory::Build(hi_SerialMap);
    RCP<Vector> s_codeOutput = VectorFactory::Build(hi_SerialMap);
    RCP<Vector> p_codeOutput = VectorFactory::Build(A->getRowMap());


    // Fill serial GOLD vecs on Proc 0
    if(!MyPID) {
      for(size_t i=0; i<(size_t)lo_gold_in.size(); i++)
        s_InVec->replaceLocalValue(i,lo_gold_in[i]);

      for(size_t i=0; i<(size_t)hi_gold_out.size(); i++)
        s_OutVec->replaceLocalValue(i,hi_gold_out[i]);
    }

    // Migrate input data
    p_InVec->doExport(*s_InVec,*lo_importer,Xpetra::ADD);

    // Apply P
    P->apply(*p_InVec,*p_codeOutput);

    // Migrate Output data
    s_codeOutput->doExport(*p_codeOutput,*hi_importer,Xpetra::ADD);

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
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
#else
    typedef Intrepid2::FieldContainer<LO> FCi;
#endif


    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=2;
    std::string hi_basis("hgrad_line_c2");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    FCi elem_to_node;
    RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC,LO,GO,NO>(num_nodes,degree,elem_to_node,lib);

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
    //    Params.set("rap: fix zero diagonals",true);
    Params.set("ipc: hi basis",hi_basis);
    Params.set("ipc: lo basis","hgrad_line_c1");
    Params.set("verbosity","high");
    Params.set("max levels",2);
     if(lib==Xpetra::UseEpetra) Params.set("coarse: type","RELAXATION");// FIXME remove when we sort out the OAZ issue
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
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
#else
    typedef Intrepid2::FieldContainer<LO> FCi;
#endif

    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=3;
    std::string hi_basis("hgrad_line_c3");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    FCi elem_to_node;
    RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC,LO,GO,NO>(num_nodes,degree,elem_to_node,lib);

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
    //    Params.set("rap: fix zero diagonals",true);
    Params.set("ipc: hi basis",hi_basis);
    Params.set("ipc: lo basis","hgrad_line_c1");
    Params.set("verbosity","high");
    Params.set("max levels",2);
    if(lib==Xpetra::UseEpetra) Params.set("coarse: type","RELAXATION");// FIXME remove when we sort out the OAZ issue
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
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
#else
    typedef Intrepid2::FieldContainer<LO> FCi;
#endif

    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=4;
    std::string hi_basis("hgrad_line_c4");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    FCi elem_to_node;
    RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC,LO,GO,NO>(num_nodes,degree,elem_to_node,lib);

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
    //    Params.set("rap: fix zero diagonals",true);
    Params.set("ipc: hi basis",hi_basis);
    Params.set("ipc: lo basis","hgrad_line_c1");
    Params.set("verbosity","high");
    Params.set("max levels",2);
    if(lib==Xpetra::UseEpetra) Params.set("coarse: type","RELAXATION");// FIXME remove when we sort out the OAZ issue
    Params.set("coarse: max size",100);
    level0.set("ipc: element to node map",rcp(&elem_to_node,false));
    Params.set("level 0",level0);
      

    // Build hierarchy
    RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,Params);
  }


/*********************************************************************************************************************/
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class Basis>
bool test_representative_basis(Teuchos::FancyOStream &out, const std::string & name, Intrepid2::EPointType ptype, int max_degree)			       
  {
#   include "MueLu_UseShortNames.hpp"
    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;  
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<MT,typename Node::device_type> FC;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
#endif
    out << "version: " << MueLu::Version() << std::endl;
    
    for(int i=1; i<max_degree; i++) {
      FC hi_DofCoords;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
      RCP<Basis> hi = rcp(new Basis(i,ptype));
      Kokkos::Experimental::resize(hi_DofCoords,hi->getCardinality(),hi->getBaseCellTopology().getDimension());
      hi->getDofCoords(hi_DofCoords);
#else
      RCP<Basis> hi = rcp(new Basis(i,ptype));
      RCP<Intrepid2::DofCoordsInterface<FC> > hi_dci = rcp_dynamic_cast<Basis>(hi);
      hi_DofCoords.resize(hi->getCardinality(),hi->getBaseCellTopology().getDimension());
      hi_dci->getDofCoords(hi_DofCoords);
#endif
      for(int j=1; j<i; j++) {
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
	RCP<Basis> lo = rcp(new Basis(j,ptype));
#else
	RCP<Basis> lo = rcp(new Basis(j,ptype));
#endif

	// Get the candidates
	double threshold = 1e-10;
	std::vector<std::vector<size_t> > candidates;
	MueLu::MueLuIntrepid::GenerateRepresentativeBasisNodes<Basis,FC>(*lo,hi_DofCoords,threshold,candidates);

	// Correctness Test 1: Make sure that there are no duplicates in the representative lists / no low DOF has no candidates
	std::vector<bool> is_candidate(hi_DofCoords.dimension(0),false);
	bool no_doubles = true;
	for(int k=0; no_doubles && k<(int)candidates.size(); k++) {
	  if(candidates[k].size()==0) no_doubles=false;
	  for(int l=0; l<(int)candidates[k].size(); l++)	    
	    if(is_candidate[candidates[k][l]] == false) is_candidate[candidates[k][l]]=true;
	    else {no_doubles=false;break;}
	}
#if 1
	if(!no_doubles) {
	  printf("*** lo/hi = %d/%d ***\n",j,i);
	  for(int k=0; k<(int)candidates.size(); k++) {
	    printf("candidates[%d] = ",k);
	    for(int l=0; l<(int)candidates[k].size(); l++)
	      printf("%d ",(int)candidates[k][l]);
	    printf("\n");
	  }
	}
#endif	         
	if(!no_doubles) {
	  out<<"ERROR: "<<name<<" The 'no duplicates' test fails w/ lo/hi = "<<j<<"/"<<i<<std::endl;
	  return false;
	}		
      }
    }
    
    // Everything worked
    return true;
  }
 
/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_LINE_Equispaced, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;    
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef typename Node::device_type::execution_space ES;
    typedef Kokkos::DynRankView<MT,typename Node::device_type> FC;
    typedef Intrepid2::Basis_HGRAD_LINE_Cn_FEM<ES,MT,MT> Basis;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
    typedef Intrepid2::Basis_HGRAD_LINE_Cn_FEM<MT,FC> Basis;

#endif

    bool rv = test_representative_basis<Scalar,LocalOrdinal,GlobalOrdinal,Node,Basis>(out," GenerateRepresentativeBasisNodes_LINE_EQUISPACED",Intrepid2::POINTTYPE_EQUISPACED,10);
    TEST_EQUALITY(rv,true);
  }

/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_QUAD_Equispaced, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;    
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef typename Node::device_type::execution_space ES;
    typedef Kokkos::DynRankView<MT,typename Node::device_type> FC;
    typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES,MT,MT> Basis;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
    typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC> Basis;

#endif

    bool rv = test_representative_basis<Scalar,LocalOrdinal,GlobalOrdinal,Node,Basis>(out," GenerateRepresentativeBasisNodes_QUAD_EQUISPACED",Intrepid2::POINTTYPE_EQUISPACED,10);
    TEST_EQUALITY(rv,true);
  }

/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_QUAD_Spectral, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;    
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef typename Node::device_type::execution_space ES;
    typedef Kokkos::DynRankView<MT,typename Node::device_type> FC;
    typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES,MT,MT> Basis;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
    typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC> Basis;

#endif

    const Intrepid2::EPointType POINTTYPE_SPECTRAL = static_cast<Intrepid2::EPointType>(1);// Not sure why I have to do this...
    bool rv = test_representative_basis<Scalar,LocalOrdinal,GlobalOrdinal,Node,Basis>(out," GenerateRepresentativeBasisNodes_QUAD_SPECTRAL",POINTTYPE_SPECTRAL,10);
    TEST_EQUALITY(rv,true);
  }


/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_HEX_Equispaced, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;    
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef typename Node::device_type::execution_space ES;
    typedef Kokkos::DynRankView<MT,typename Node::device_type> FC;
    typedef Intrepid2::Basis_HGRAD_HEX_Cn_FEM<ES,MT,MT> Basis;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
    typedef Intrepid2::Basis_HGRAD_HEX_Cn_FEM<MT,FC> Basis;

#endif

    bool rv = test_representative_basis<Scalar,LocalOrdinal,GlobalOrdinal,Node,Basis>(out," GenerateRepresentativeBasisNodes_HEX_EQUISPACED",Intrepid2::POINTTYPE_EQUISPACED,8);
    TEST_EQUALITY(rv,true);
  }

/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_HEX_Spectral, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;    
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef typename Node::device_type::execution_space ES;
    typedef Kokkos::DynRankView<MT,typename Node::device_type> FC;
    typedef Intrepid2::Basis_HGRAD_HEX_Cn_FEM<ES,MT,MT> Basis;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
    typedef Intrepid2::Basis_HGRAD_HEX_Cn_FEM<MT,FC> Basis;
#endif

    const Intrepid2::EPointType POINTTYPE_SPECTRAL = static_cast<Intrepid2::EPointType>(1);// Not sure why I have to do this...
    bool rv = test_representative_basis<Scalar,LocalOrdinal,GlobalOrdinal,Node,Basis>(out," GenerateRepresentativeBasisNodes_HEX_SPECTRAL",POINTTYPE_SPECTRAL,8);
    TEST_EQUALITY(rv,true);
  }

/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, GenerateLoNodeInHighViaGIDs_QUAD_pn_to_p1, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {

    // NOTE: We need more tests for this that do pn to pm
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;    

#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef typename Node::device_type::execution_space ES;
    typedef Kokkos::DynRankView<MT,typename Node::device_type> FC;
    typedef Kokkos::DynRankView<LO,typename Node::device_type> FCi;
    typedef Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES,MT,MT> LoBasis;
    typedef Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES,MT,MT> HiBasis;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
    typedef Intrepid2::FieldContainer<LO> FCi;
    typedef Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC> LoBasis;
    typedef Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC> HiBasis;
#endif
    out << "version: " << MueLu::Version() << std::endl;
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    int max_degree = 10;
    double threshold = 1e-10;
    GO gst_invalid = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

    // Lo
    RCP<LoBasis> lo = rcp(new LoBasis());
    size_t numLo = lo->getCardinality();
    
    for(int i=0;i<max_degree; i++) {
      RCP<HiBasis> hi = rcp(new HiBasis(i,Intrepid2::POINTTYPE_EQUISPACED));
      size_t numHi    = hi->getCardinality();
      
      // The quad-only stuff 
      std::vector<size_t> lo_node_in_hi;
      FC hi_dofCoords;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR 
      MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<MT,typename Node::device_type>(hi,lo,lo_node_in_hi,hi_dofCoords);
      FCi hi_e2n("hi_e2n",1,numHi);
      FCi lo_e2n("lo_e2n",1,numLo);
#else
      MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<MT,FC>(hi,lo,lo_node_in_hi,hi_dofCoords);
      FCi hi_e2n(1,numHi);
      FCi lo_e2n(1,numLo);
#endif  

      // Dummy elem2node map
      Teuchos::Array<GO> hi_colids(numHi);
      for(size_t j=0; j<numHi; j++) {
	hi_e2n(j)    = j;
	hi_colids[j] = j;
      }

      // Dummy column map
      RCP<const Map> hi_colMap      = MapFactory::Build(lib,gst_invalid,hi_colids(),0,comm);

      // The dynamic stuff
      std::vector<std::vector<size_t> > candidates;
      MueLu::MueLuIntrepid::GenerateRepresentativeBasisNodes<LoBasis,FC>(*lo,hi_dofCoords,threshold,candidates);
      MueLu::MueLuIntrepid::GenerateLoNodeInHiViaGIDs<LO,GO,Node,FCi>(candidates,hi_e2n,hi_colMap,lo_e2n);

      // Compare and make sure we're cool
      bool node_diff = false;
      for(size_t j=0; j<numLo; j++)
	if(lo_node_in_hi[j]!=(size_t)lo_e2n(0,j)) node_diff=true;
#if 0
      printf("[%d] Comparison = ",i);
      for(size_t j=0; j<numLo; j++)
	printf("%d|%d ",(int)lo_node_in_hi[j],(int)lo_e2n(0,j));
      printf("\n");
#endif
      TEST_EQUALITY(node_diff,false);
      }
  }



  /*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildLoElemToNodeViaRepresentatives_QUAD_pn_to_p1, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
    // NOTE: We need more tests for this that do pn to pm
  #   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<MT,typename Node::device_type> FC;
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
    typedef typename Node::device_type::execution_space ES;
    typedef Intrepid2::Basis<ES,MT,MT> Basis;
#else
    typedef Intrepid2::FieldContainer<MT> FC;
    typedef Intrepid2::FieldContainer<LO> FCi;
    typedef Intrepid2::Basis<MT,FC> Basis;
#endif
    GO gst_invalid = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

    out << "version: " << MueLu::Version() << std::endl;
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    int max_degree=10;

    {
      //QUAD
      // A one element test with Kirby-numbered nodes where the top edge is not owned      
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
      RCP<Basis> lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<ES,MT,MT>());
#else
      RCP<Basis> lo = rcp(new Intrepid2::Basis_HGRAD_QUAD_C1_FEM<MT,FC>());
#endif

      for(int degree=2; degree < max_degree; degree++) {
        int Nn = (degree+1)*(degree+1);
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
        RCP<Basis> hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<ES,MT,MT>(degree,Intrepid2::POINTTYPE_EQUISPACED));
        FCi hi_e2n("hi_e2n",1,Nn), lo_e2n, lo_e2n_mk2;
#else
        RCP<Basis> hi = rcp(new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<MT,FC>(degree,Intrepid2::POINTTYPE_EQUISPACED));
        FCi hi_e2n(1,Nn), lo_e2n, lo_e2n_mk2;
#endif
        std::vector<bool> hi_owned(Nn,false),lo_owned, lo_owned_mk2;
        std::vector<size_t> lo_node_in_hi;
        std::vector<LO> hi_to_lo_map,  hi_to_lo_map_mk2;
        int lo_numOwnedNodes=0, lo_numOwnedNodes_mk2=0;
        FC hi_dofCoords;

	// El2node / ownership / colmap
	Teuchos::Array<GO> hi_colids(Nn);
        for(int i=0; i<Nn; i++) {
	  hi_colids[i] = i;
          hi_e2n(0,i)=i;
          if(i < Nn-(degree+1)) hi_owned[i]=true;
        }

	/*** Do stuff the injection way ***/
#ifdef HAVE_MUELU_INTREPID2_REFACTOR 
        MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<MT,typename Node::device_type>(hi,lo,lo_node_in_hi,hi_dofCoords);
#else
        MueLu::MueLuIntrepid::IntrepidGetLoNodeInHi<MT,FC>(hi,lo,lo_node_in_hi,hi_dofCoords);
#endif  
        MueLu::MueLuIntrepid::BuildLoElemToNode(hi_e2n,hi_owned,lo_node_in_hi,lo_e2n,lo_owned,hi_to_lo_map,lo_numOwnedNodes);

	/*** Do stuff the representative way ***/
	RCP<const Map> hi_colMap      = MapFactory::Build(lib,gst_invalid,hi_colids(),0,comm);
	FCi lo_elemToHiRepresentativeNode;
	double threshold = 1e-10;
	std::vector<std::vector<size_t> > candidates;
	MueLu::MueLuIntrepid::GenerateRepresentativeBasisNodes<Basis,FC>(*lo,hi_dofCoords,threshold,candidates);
	MueLu::MueLuIntrepid::GenerateLoNodeInHiViaGIDs(candidates,hi_e2n,hi_colMap,lo_elemToHiRepresentativeNode);
	MueLu::MueLuIntrepid::BuildLoElemToNodeViaRepresentatives(hi_e2n,hi_owned,lo_elemToHiRepresentativeNode,lo_e2n_mk2,lo_owned_mk2,hi_to_lo_map_mk2,lo_numOwnedNodes_mk2);

	// Compare stuff
        TEST_EQUALITY(lo_numOwnedNodes,2);
        TEST_EQUALITY(lo_numOwnedNodes_mk2,2);

	size_t num_lo_nodes_located=0;
        for(size_t i=0;i<hi_to_lo_map.size(); i++) {
          if(hi_to_lo_map[i] != Teuchos::OrdinalTraits<LO>::invalid())
            num_lo_nodes_located++;
        }
	TEST_EQUALITY(lo_owned.size(),num_lo_nodes_located);
	TEST_EQUALITY(lo_owned_mk2.size(),num_lo_nodes_located);

	for(size_t i=0; i<lo_e2n.dimension(0); i++) 
	  for(size_t j=0; j<lo_e2n.dimension(1); j++) 
	    TEST_EQUALITY(lo_e2n(i,j),lo_e2n_mk2(i,j));

	for(size_t i=0; i<(size_t) lo_owned.size(); i++) 
	  TEST_EQUALITY(lo_owned[i],lo_owned_mk2[i]);

	TEST_EQUALITY(hi_to_lo_map.size(),hi_to_lo_map_mk2.size());
	for(size_t i=0; i<(size_t) hi_to_lo_map.size(); i++) 
	  TEST_EQUALITY(hi_to_lo_map[i],hi_to_lo_map_mk2[i]);

      }
    }//end QUAD

  }


/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildP_PseudoPoisson_LINE_p3_to_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
    int hi_degree=3;
    int lo_degree=2;
    // GOLD vector collection
    // Note: Vectors are exodus-ordered, not Kirby-ordered
    size_t num_p1_points = 10;
    size_t num_hi_points = num_p1_points + (num_p1_points-1) *(hi_degree-1);
    size_t num_lo_points = num_p1_points + (num_p1_points-1) *(lo_degree-1);

    std::vector<Scalar> lo_gold_in(num_lo_points);
    std::vector<Scalar> hi_gold_out(num_hi_points);				    
    for(size_t i=0; i<num_p1_points; i++) {
      lo_gold_in[i] = i;
      hi_gold_out[i] = i;
    }

    size_t idx=num_p1_points;
    for(size_t i=0; i<num_p1_points-1; i++) {
      for(size_t j=0; j<(size_t)hi_degree-1; j++) {
        hi_gold_out[idx] = i + ((double)j+1)/hi_degree;
        idx++;
      }
    }

    idx=num_p1_points;
    for(size_t i=0; i<num_p1_points-1; i++) {
      for(size_t j=0; j<(size_t)lo_degree-1; j++) {
        lo_gold_in[idx] = i + ((double)j+1)/lo_degree;
        idx++;
      }
    }

    TestPseudoPoisson<Scalar,LocalOrdinal,GlobalOrdinal,Node>(out,num_p1_points,hi_degree,lo_gold_in,hi_gold_out,std::string("hgrad_line_c3"),std::string("hgrad_line_c2"));
  }

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildP_PseudoPoisson_LINE_p4_to_p3, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
    int hi_degree=4;
    int lo_degree=3;
    // GOLD vector collection
    // Note: Vectors are exodus-ordered, not Kirby-ordered
    size_t num_p1_points = 10;
    size_t num_hi_points = num_p1_points + (num_p1_points-1) *(hi_degree-1);
    size_t num_lo_points = num_p1_points + (num_p1_points-1) *(lo_degree-1);

    std::vector<Scalar> lo_gold_in(num_lo_points);
    std::vector<Scalar> hi_gold_out(num_hi_points);				    
    for(size_t i=0; i<num_p1_points; i++) {
      lo_gold_in[i] = i;
      hi_gold_out[i] = i;
    }

    size_t idx=num_p1_points;
    for(size_t i=0; i<num_p1_points-1; i++) {
      for(size_t j=0; j<(size_t)hi_degree-1; j++) {
        hi_gold_out[idx] = i + ((double)j+1)/hi_degree;
        idx++;
      }
    }

    idx=num_p1_points;
    for(size_t i=0; i<num_p1_points-1; i++) {
      for(size_t j=0; j<(size_t)lo_degree-1; j++) {
        lo_gold_in[idx] = i + ((double)j+1)/lo_degree;
        idx++;
      }
    }

    TestPseudoPoisson<Scalar,LocalOrdinal,GlobalOrdinal,Node>(out,num_p1_points,hi_degree,lo_gold_in,hi_gold_out,std::string("hgrad_line_c4"),std::string("hgrad_line_c3"));
  }

/*********************************************************************************************************************/
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory,BuildP_PseudoPoisson_LINE_p4_to_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
    int hi_degree=4;
    int lo_degree=2;
    // GOLD vector collection
    // Note: Vectors are exodus-ordered, not Kirby-ordered
    size_t num_p1_points = 10;
    size_t num_hi_points = num_p1_points + (num_p1_points-1) *(hi_degree-1);
    size_t num_lo_points = num_p1_points + (num_p1_points-1) *(lo_degree-1);

    std::vector<Scalar> lo_gold_in(num_lo_points);
    std::vector<Scalar> hi_gold_out(num_hi_points);				    
    for(size_t i=0; i<num_p1_points; i++) {
      lo_gold_in[i] = i;
      hi_gold_out[i] = i;
    }

    size_t idx=num_p1_points;
    for(size_t i=0; i<num_p1_points-1; i++) {
      for(size_t j=0; j<(size_t)hi_degree-1; j++) {
        hi_gold_out[idx] = i + ((double)j+1)/hi_degree;
        idx++;
      }
    }

    idx=num_p1_points;
    for(size_t i=0; i<num_p1_points-1; i++) {
      for(size_t j=0; j<(size_t)lo_degree-1; j++) {
        lo_gold_in[idx] = i + ((double)j+1)/lo_degree;
        idx++;
      }
    }

    TestPseudoPoisson<Scalar,LocalOrdinal,GlobalOrdinal,Node>(out,num_p1_points,hi_degree,lo_gold_in,hi_gold_out,std::string("hgrad_line_c4"),std::string("hgrad_line_c2"));
  }


/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p3_to_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
#else
    typedef Intrepid2::FieldContainer<LO> FCi;
#endif

    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=3;
    std::string hi_basis("hgrad_line_c3");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    FCi elem_to_node;
    RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC,LO,GO,NO>(num_nodes,degree,elem_to_node,lib);

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
    //    Params.set("rap: fix zero diagonals",true);
    Params.set("ipc: hi basis",hi_basis);
    Params.set("ipc: lo basis","hgrad_line_c2");
    Params.set("verbosity","high");
    Params.set("max levels",2);
    if(lib==Xpetra::UseEpetra) Params.set("coarse: type","RELAXATION");// FIXME remove when we sort out the OAZ issue
    Params.set("coarse: max size",100);
    level0.set("ipc: element to node map",rcp(&elem_to_node,false));
    Params.set("level 0",level0);
      

    // Build hierarchy
    RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,Params);
  }

/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p3, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
#else
    typedef Intrepid2::FieldContainer<LO> FCi;
#endif

    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=4;
    std::string hi_basis("hgrad_line_c4");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    FCi elem_to_node;
    RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC,LO,GO,NO>(num_nodes,degree,elem_to_node,lib);

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
    //    Params.set("rap: fix zero diagonals",true);
    Params.set("ipc: hi basis",hi_basis);
    Params.set("ipc: lo basis","hgrad_line_c3");
    Params.set("verbosity","high");
    Params.set("max levels",2);
    if(lib==Xpetra::UseEpetra) Params.set("coarse: type","RELAXATION");// FIXME remove when we sort out the OAZ issue
    Params.set("coarse: max size",100);
    level0.set("ipc: element to node map",rcp(&elem_to_node,false));
    Params.set("level 0",level0);
      

    // Build hierarchy
    RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,Params);
  }

/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
#else
    typedef Intrepid2::FieldContainer<LO> FCi;
#endif

    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=4;
    std::string hi_basis("hgrad_line_c4");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    FCi elem_to_node;
    RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC,LO,GO,NO>(num_nodes,degree,elem_to_node,lib);

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
    //    Params.set("rap: fix zero diagonals",true);
    Params.set("ipc: hi basis",hi_basis);
    Params.set("ipc: lo basis","hgrad_line_c2");
    Params.set("verbosity","high");
    Params.set("max levels",2);
    if(lib==Xpetra::UseEpetra) Params.set("coarse: type","RELAXATION");// FIXME remove when we sort out the OAZ issue
    Params.set("coarse: max size",100);
    level0.set("ipc: element to node map",rcp(&elem_to_node,false));
    Params.set("level 0",level0);
      

    // Build hierarchy
    RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,Params);
  }


/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p3_to_p2, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
#else
    typedef Intrepid2::FieldContainer<LO> FCi;
#endif

    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=4;
    std::string hi_basis("hgrad_line_c4");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    FCi elem_to_node;
    RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC,LO,GO,NO>(num_nodes,degree,elem_to_node,lib);

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
    ParameterList Params, level0, level1, level2;
    Params.set("multigrid algorithm","pcoarsen");
    //    Params.set("rap: fix zero diagonals",true);
    Params.set("verbosity","high");
    Params.set("max levels",3);
    if(lib==Xpetra::UseEpetra) Params.set("coarse: type","RELAXATION");// FIXME remove when we sort out the OAZ issue
    Params.set("coarse: max size",100);

    level0.set("ipc: element to node map",rcp(&elem_to_node,false));
    Params.set("level 0",level0);

    level1.set("ipc: hi basis",hi_basis);
    level1.set("ipc: lo basis","hgrad_line_c3");
    Params.set("level 1",level1);

    level2.set("ipc: hi basis","hgrad_line_c3");
    level2.set("ipc: lo basis","hgrad_line_c2");
    Params.set("level 2",level2);
      

    // Build hierarchy
    RCP<Hierarchy> tH = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,Params);
  }



/*********************************************************************************************************************/
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IntrepidPCoarsenFactory, CreatePreconditioner_p2_to_p1_sa, Scalar, LocalOrdinal, GlobalOrdinal, Node) 
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
#   if !defined(HAVE_MUELU_AMESOS) || !defined(HAVE_MUELU_IFPACK)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#   endif
#   if !defined(HAVE_MUELU_AMESOS2) || !defined(HAVE_MUELU_IFPACK2)
    MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#   endif

    typedef Scalar SC;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO; 
    typedef Node  NO;  
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
#ifdef HAVE_MUELU_INTREPID2_REFACTOR
    typedef Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> FCi;
#else
    typedef Intrepid2::FieldContainer<LO> FCi;
#endif

    out << "version: " << MueLu::Version() << std::endl;
    using Teuchos::RCP;
    int degree=2;
    std::string hi_basis("hgrad_line_c2");

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO num_nodes = 972;
    // Build a pseudo-poisson test matrix
    FCi elem_to_node;
    RCP<Matrix> A = TestHelpers::Build1DPseudoPoissonHigherOrder<SC,LO,GO,NO>(num_nodes,degree,elem_to_node,lib);

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
    ParameterList Params, level0, level1, level2;
    Params.set("multigrid algorithm","pcoarsen");
    //    Params.set("rap: fix zero diagonals",true);
    Params.set("verbosity","high");
    Params.set("max levels",3);
    if(lib==Xpetra::UseEpetra) Params.set("coarse: type","RELAXATION");// FIXME remove when we sort out the OAZ issue
    Params.set("coarse: max size",100);

    level0.set("ipc: element to node map",rcp(&elem_to_node,false));
    Params.set("level 0",level0);

    level1.set("multigrid algorithm","pcoarsen");
    level1.set("ipc: hi basis",hi_basis);
    level1.set("ipc: lo basis","hgrad_line_c1");
    Params.set("level 1",level1);

    level2.set("multigrid algorithm","sa");
    Params.set("level 2",level2);
      

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
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p4, Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_LINE_Equispaced,Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_QUAD_Equispaced,Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_QUAD_Spectral,  Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_HEX_Equispaced, Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateRepresentativeBasisNodes_HEX_Spectral,   Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, GenerateLoNodeInHighViaGIDs_QUAD_pn_to_p1, Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildLoElemToNodeViaRepresentatives_QUAD_pn_to_p1, Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_LINE_p3_to_p2,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_LINE_p4_to_p3,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, BuildP_PseudoPoisson_LINE_p4_to_p2,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p3_to_p2, Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p3, Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p2, Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p4_to_p3_to_p2, Scalar, LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory, CreatePreconditioner_p2_to_p1_sa, Scalar, LO,GO,Node)


#include <MueLu_ETI_4arg.hpp>


} // namespace MueLuTests
#endif
