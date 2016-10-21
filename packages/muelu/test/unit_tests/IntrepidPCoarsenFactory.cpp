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

    // This function basically takes an existing domain->column importer (the "hi order" guy) and using the hi_to_lo_map, generates "lo order" version.  The domain map is already given to
    // us here, so we just need to make tha column one.

    // We'll test this by starting with some linearly distributed map for the "hi" domain map and a duplicated map as the "hi" column map.
    // We then do a "lo" domain map, which just grabs the first GID on each proc.
    // Testing this should be easy.

#if 0
 GenerateColMapFromImport(const Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> & hi_importer,const std::vector<LocalOrdinal> &hi_to_lo_map,const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> & lo_domainMap, const size_t & lo_columnMapLength, RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & lo_columnMap);
#endif
  }

  /*********************************************************************************************************************/
#  define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,GetLoNodeInHi,Scalar,LO,GO,Node)  \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,BasisFactory,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IntrepidPCoarsenFactory,BuildLoElemToNode,Scalar,LO,GO,Node)


#include <MueLu_ETI_4arg.hpp>


} // namespace MueLuTests
