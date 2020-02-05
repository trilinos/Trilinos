/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_FECrsMatrix.hpp"
#include "Tpetra_FEMultiVector.hpp"
#include "Tpetra_Details_getNumDiags.hpp"
#include "KokkosCompat_View.hpp"
// TODO: add test where some nodes have zero rows
// TODO: add test where non-"zero" graph is used to build matrix; if no values are added to matrix, the operator effect should be zero. This tests that matrix values are initialized properly.
// TODO: add test where dynamic profile initially has no allocation, then entries are added. this will test new view functionality.

namespace { // (anonymous)

using Tpetra::TestingUtilities::getDefaultComm;
using Tpetra::createContigMapWithNode;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::rcp;
using Teuchos::outArg;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::tuple;
using Teuchos::NO_TRANS;
//using Teuchos::TRANS;
using Teuchos::CONJ_TRANS;
using std::endl;
typedef Tpetra::global_size_t GST;
typedef std::complex<double> scd;

template<class Scalar, class LO, class GO, class Node, class TOLERANCE >
bool compare_final_matrix_structure_impl(Teuchos::FancyOStream &out,Tpetra::CrsMatrix<Scalar,LO,GO,Node> & g1, Tpetra::CrsMatrix<Scalar,LO,GO,Node> & g2, TOLERANCE tol) {
  using std::endl;
 
  if (!g1.isFillComplete() || !g2.isFillComplete()) {out<<"Compare: FillComplete failed"<<endl;return false;}
  if (!g1.getRangeMap()->isSameAs(*g2.getRangeMap())) {out<<"Compare: RangeMap failed"<<endl;return false;}
  if (!g1.getRowMap()->isSameAs(*g2.getRowMap())) {out<<"Compare: RowMap failed"<<endl;return false;}
  if (!g1.getColMap()->isSameAs(*g2.getColMap())) {out<<"Compare: ColMap failed"<<endl;return false;}
  if (!g1.getDomainMap()->isSameAs(*g2.getDomainMap())) {out<<"Compare: DomainMap failed"<<endl;return false;}

  auto rowptr1 = g1.getLocalMatrix().graph.row_map;
  auto rowptr2 = g2.getLocalMatrix().graph.row_map;

  auto colind1 = g1.getLocalMatrix().graph.entries;
  auto colind2 = g2.getLocalMatrix().graph.entries;

  auto values1 = g1.getLocalMatrix().values;
  auto values2 = g2.getLocalMatrix().values;

  if (rowptr1.extent(0) != rowptr2.extent(0)) {out<<"Compare: rowptr extent failed"<<endl;return false;}      
  if (colind1.extent(0) != colind2.extent(0)) {out<<"Compare: colind extent failed"<<endl;return false;}      
  if (values1.extent(0) != values2.extent(0)) {out<<"Compare: values extent failed"<<endl;return false;}      

  bool success=true;
  TEST_COMPARE_ARRAYS(rowptr1,rowptr2);
  if (!success) {out<<"Compare: rowptr match failed"<<endl;return false;}

  TEST_COMPARE_ARRAYS(colind1,colind2);
  if (!success) {out<<"Compare: colind match failed"<<endl;return false;}

  // This is necessary to make sure that complex works (since Teuchos::ScalarTraits does not have a Kokkos::complex specialization)
  auto values1_h = Kokkos::create_mirror_view(values1);
  auto values2_h = Kokkos::create_mirror_view(values2);
  auto values1_av = Teuchos::av_reinterpret_cast<Scalar>(Kokkos::Compat::getArrayView(values1_h));
  auto values2_av = Teuchos::av_reinterpret_cast<Scalar>(Kokkos::Compat::getArrayView(values2_h));  
  TEST_COMPARE_FLOATING_ARRAYS(values1_av,values2_av,tol);
  if (!success) {out<<"Compare: values match failed"<<endl;return false;}

  return true;
}


template<class Scalar, class LO, class GO, class Node>
struct compare {
  static bool compare_final_matrix_structure(Teuchos::FancyOStream &out,Tpetra::CrsMatrix<Scalar,LO,GO,Node> & g1, Tpetra::CrsMatrix<Scalar,LO,GO,Node> & g2){ 
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    double errorTolSlack = 1.0e+2;
    const Mag tol = errorTolSlack * Teuchos::ScalarTraits<Scalar>::eps();
    return compare_final_matrix_structure_impl(out,g1,g2,tol);
  }
};

template<class LO, class GO, class Node>
struct compare<int,LO,GO,Node> {
  static bool compare_final_matrix_structure(Teuchos::FancyOStream &out,Tpetra::CrsMatrix<int,LO,GO,Node> & g1, Tpetra::CrsMatrix<int,LO,GO,Node> & g2) {
    return compare_final_matrix_structure_impl(out,g1,g2,0);
  }
};

template<class LO, class GO, class Node>
struct compare<long long,LO,GO,Node> {
  static bool compare_final_matrix_structure(Teuchos::FancyOStream &out,Tpetra::CrsMatrix<long long,LO,GO,Node> & g1, Tpetra::CrsMatrix<long long,LO,GO,Node> & g2) {
    return compare_final_matrix_structure_impl(out,g1,g2,0);
  }
};

template<class LO, class GO, class Node>
class GraphPack {
public:
  RCP<const Tpetra::Map<LO,GO,Node> > uniqueMap;
  RCP<const Tpetra::Map<LO,GO,Node> > overlapMap;
  std::vector<std::vector<GO> > element2node;

  // NOTE: This is hardwired for 1D bar elements
  Kokkos::View<GO*[2], typename Node::device_type > k_element2node;

  void print(int rank, std::ostream & out) {
    using std::endl;
    out << "["<<rank<<"] Unique Map  : ";
    for(size_t i=0; i<uniqueMap->getNodeNumElements(); i++)
      out << uniqueMap->getGlobalElement(i) << " ";
    out<<endl;      

    out << "["<<rank<<"] Overlap Map : ";
    for(size_t i=0; i<overlapMap->getNodeNumElements(); i++)
      out << overlapMap->getGlobalElement(i) << " ";
    out<<endl;

    out << "["<<rank<<"] element2node: ";
    for(size_t i=0; i<(size_t)element2node.size(); i++) {
      out <<"(";
      for(size_t j=0; j<(size_t)element2node[i].size(); j++)
        out<<element2node[i][j] << " ";
      out<<") ";
    }
    out<<endl;
  }
};


template<class LO, class GO, class Node>
void generate_fem1d_graph(size_t numLocalNodes, RCP<const Comm<int> > comm , GraphPack<LO,GO,Node> & pack) {
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();
  int rank    = comm->getRank();
  int numProc = comm->getSize();
  size_t numOverlapNodes = numLocalNodes; if(rank!=numProc-1) numOverlapNodes++;  if(rank!=0) numOverlapNodes++;
  size_t numLocalElements = (rank == numProc-1) ? numLocalNodes -1 : numLocalNodes;
  //  printf("CMS numOverlapNodes = %d numLocalElements = %d\n",numOverlapNodes,numLocalElements);

  pack.uniqueMap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocalNodes,comm);

  Teuchos::Array<GO> overlapIndices(numOverlapNodes);
  for(size_t i=0; i<numLocalNodes; i++) {
    overlapIndices[i] = pack.uniqueMap->getGlobalElement(i);
  }
  size_t last = numLocalNodes;
  if(rank != 0)           {overlapIndices[last] = overlapIndices[0] - 1; last++;}
  if(rank != numProc -1)  {overlapIndices[last] = overlapIndices[numLocalNodes-1] + 1; last++;}

  pack.overlapMap = rcp(new Tpetra::Map<LO,GO,Node>(INVALID,overlapIndices,0,comm));

  pack.element2node.resize(numLocalElements);
  for(size_t i=0; i<numLocalElements; i++) {
    pack.element2node[i].resize(2);
    pack.element2node[i][0] = pack.uniqueMap->getGlobalElement(i);
    pack.element2node[i][1] = pack.uniqueMap->getGlobalElement(i) + 1;
  }

  // Kokkos version of the element2node array
  Kokkos::resize(pack.k_element2node,numLocalElements);
  auto k_e2n = pack.k_element2node;
  auto l_umap = pack.uniqueMap->getLocalMap();
  Kokkos::parallel_for(Kokkos::RangePolicy<typename Node::execution_space>(0,numLocalElements), KOKKOS_LAMBDA(const size_t i) {
      GO gid = l_umap.getGlobalElement(i);
      k_e2n(i,0) = gid;
      k_e2n(i,1) = gid+1;
    });
  Kokkos::fence();

}

template<class Scalar>
std::vector<std::vector<Scalar> > generate_fem1d_element_values() {
  std::vector<std::vector<Scalar> > mat;
  mat.resize(2);
  mat[0].resize(2);
  mat[1].resize(2);
  mat[0][0] =  1;
  mat[0][1] = -1;
  mat[1][0] = -1;
  mat[1][1] =  1;

  return mat;
}

template<class ImplScalarType, class Node>
Kokkos::View<ImplScalarType[2][2], Kokkos::LayoutLeft, typename Node::device_type > generate_fem1d_element_values_kokkos() {  
  Kokkos::View<ImplScalarType[2][2], Kokkos::LayoutLeft, typename Node::device_type> mat ("fem1d_element_values");

  auto mat_h = Kokkos::create_mirror_view(mat);
  mat_h(0,0) =  1.0;
  mat_h(0,1) = -1.0;
  mat_h(1,0) = -1.0;
  mat_h(1,1) =  1.0;

  Kokkos::deep_copy(mat, mat_h);
  return mat;
}


//
// UNIT TESTS
//


////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FECrsMatrix, Assemble1D, LO, GO, Scalar, Node )
{
  using FEMAT = typename Tpetra::FECrsMatrix<Scalar,LO,GO,Node>;
  using CMAT = typename Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
  using FEG = typename Tpetra::FECrsGraph<LO,GO,Node>;

  // get a comm
  RCP<const Comm<int> > comm = getDefaultComm();
  
  // Generate a mesh
  size_t numLocal = 10;
  GraphPack<LO,GO,Node> pack;
  generate_fem1d_graph(numLocal,comm,pack);
  
  // Make the graph    
  // FIXME: We should be able to get away with 3 for StaticProfile here, but we need 4 since duplicates are
  // not being handled correctly. 
  RCP<FEG> graph = rcp(new FEG(pack.uniqueMap,pack.overlapMap,4));

  graph->beginFill();
  for(size_t i=0; i<(size_t)pack.element2node.size(); i++) {
    for(size_t j=0; j<pack.element2node[i].size(); j++) {
      GO gid_j = pack.element2node[i][j];
      for(size_t k=0; k<pack.element2node[i].size(); k++) {
        GO gid_k = pack.element2node[i][k];
        //        printf("Inserting (%d,%d)\n",gid_j,gid_k);
        graph->insertGlobalIndices(gid_j,1,&gid_k);
      }
    }
  }
  graph->endFill();


  // Generate the "local stiffness matrix"
  std::vector<std::vector<Scalar> > localValues = generate_fem1d_element_values<Scalar>();

  // Make the matrix two ways
  FEMAT mat1(graph); // Here we use graph as a FECrsGraph
  CMAT mat2(graph);  // Here we use graph as a CrsGraph in OWNED mode
  mat1.beginFill();
  for(size_t i=0; i<(size_t)pack.element2node.size(); i++) {
    for(size_t j=0; j<pack.element2node[i].size(); j++) {
      GO gid_j = pack.element2node[i][j];
      for(size_t k=0; k<pack.element2node[i].size(); k++) {
        GO gid_k = pack.element2node[i][k];
        mat1.sumIntoGlobalValues(gid_j,1,&localValues[j][k],&gid_k);
        mat2.sumIntoGlobalValues(gid_j,1,&localValues[j][k],&gid_k);
      }
    }
  }
  mat1.endFill();
  mat2.fillComplete();

  success = compare<Scalar,LO,GO,Node>::compare_final_matrix_structure(out,mat1,mat2);
  TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success)
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FECrsMatrix, Assemble1D_Kokkos, LO, GO, Scalar, Node )
{
  using exec_space = typename Node::execution_space;
  using range_type = Kokkos::RangePolicy<exec_space, LO>;
  using FEMAT = typename Tpetra::FECrsMatrix<Scalar,LO,GO,Node>;
  using CMAT = typename Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
  using FEG = typename Tpetra::FECrsGraph<LO,GO,Node>;
  using ImplScalarType = typename FEMAT::impl_scalar_type;

  // get a comm
  RCP<const Comm<int> > comm = getDefaultComm();
  
  // Generate a mesh
  size_t numLocal = 10;
  GraphPack<LO,GO,Node> pack;
  generate_fem1d_graph(numLocal,comm,pack);

  // Make the graph    
  // FIXME: We should be able to get away with 3 for StaticProfile here, but we need 4 since duplicates are
  // not being handled correctly. 
  RCP<FEG> graph = rcp(new FEG(pack.uniqueMap,pack.overlapMap,4));

  graph->beginFill();
  for(size_t i=0; i<(size_t)pack.element2node.size(); i++) {
    for(size_t j=0; j<pack.element2node[i].size(); j++) {
      GO gid_j = pack.element2node[i][j];
      for(size_t k=0; k<pack.element2node[i].size(); k++) {
        GO gid_k = pack.element2node[i][k];
        //printf("Inserting (%d,%d)\n",gid_j,gid_k);
        graph->insertGlobalIndices(gid_j,1,&gid_k);
      }
    }
  }
  graph->endFill();

  // Generate the "local stiffness matrix"
  std::vector<std::vector<Scalar> > localValues = generate_fem1d_element_values<Scalar>();
  auto kokkosValues = generate_fem1d_element_values_kokkos<ImplScalarType, Node>();

  // Make the matrix two ways
  FEMAT mat1(graph); // Here we use graph as a FECrsGraph
  CMAT mat2(graph);  // Here we use graph as a CrsGraph in OWNED mode
  mat1.beginFill();

  auto k_e2n = pack.k_element2node;
  auto localMat = mat1.getLocalMatrix();
  auto localMap = pack.overlapMap->getLocalMap();
  //get local map too

  Kokkos::parallel_for("assemble_1d",
		       range_type (0,k_e2n.extent(0)), 
		       KOKKOS_LAMBDA(const size_t i) {
    size_t extent = k_e2n.extent(1);
    for(size_t j=0; j < extent; j++) {
      LO lid_j = localMap.getLocalElement(k_e2n(i, j));
      for(size_t k=0; k < extent; k++) {
        LO lid_k = localMap.getLocalElement(k_e2n(i, k));
	ImplScalarType tmp = kokkosValues(j, k);
	localMat.sumIntoValues(lid_j, &lid_k, 1, &tmp, true, true);
      }
    }
  });
  mat1.endFill();
  
  for(size_t i=0; i<(size_t)pack.element2node.size(); i++) {
    for(size_t j=0; j<pack.element2node[i].size(); j++) {
      GO gid_j = pack.element2node[i][j];
      for(size_t k=0; k<pack.element2node[i].size(); k++) {
        GO gid_k = pack.element2node[i][k];
        mat2.sumIntoGlobalValues(gid_j,1,&localValues[j][k],&gid_k);
      }
    }
  }
  mat2.fillComplete();

  success = compare<Scalar,LO,GO,Node>::compare_final_matrix_structure(out,mat1,mat2);
  TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success)
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FECrsMatrix, Assemble1D_LocalIndex, LO, GO, Scalar, Node )
{
  using FEMAT = typename Tpetra::FECrsMatrix<Scalar,LO,GO,Node>;
  using CMAT = typename Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
  using FEG = typename Tpetra::FECrsGraph<LO,GO,Node>;

  // get a comm
  RCP<const Comm<int> > comm = getDefaultComm();
  
  // Generate a mesh
  size_t numLocal = 10;
  GraphPack<LO,GO,Node> pack;
  generate_fem1d_graph(numLocal,comm,pack);
  
  // Make the graph    
  RCP<FEG> graph = rcp(new FEG(pack.uniqueMap,pack.overlapMap,pack.overlapMap,3));

  graph->beginFill();
  for(size_t i=0; i<(size_t)pack.element2node.size(); i++) {
    for(size_t j=0; j<pack.element2node[i].size(); j++) {
      GO gid_j = pack.element2node[i][j];
      LO lid_j = pack.overlapMap->getLocalElement(gid_j);
      for(size_t k=0; k<pack.element2node[i].size(); k++) {
        GO gid_k = pack.element2node[i][k];
        LO lid_k = pack.overlapMap->getLocalElement(gid_k);
        //printf("[%d] Inserting gid (%d,%d) lid (%d,%d)\n",comm->getRank(),gid_j,gid_k,lid_j,lid_k);fflush(stdout);
        graph->insertLocalIndices(lid_j,1,&lid_k);
      }
    }
  }
  graph->endFill();


  // Generate the "local stiffness matrix"
  std::vector<std::vector<Scalar> > localValues = generate_fem1d_element_values<Scalar>();

  // Make the matrix two ways
  FEMAT mat1(graph); // Here we use graph as a FECrsGraph
  CMAT mat2(graph);  // Here we use graph as a CrsGraph in OWNED mode
  mat1.beginFill();
  for(size_t i=0; i<(size_t)pack.element2node.size(); i++) {
    for(size_t j=0; j<pack.element2node[i].size(); j++) {
      GO gid_j = pack.element2node[i][j];
      LO lid_j = pack.overlapMap->getLocalElement(gid_j);
      for(size_t k=0; k<pack.element2node[i].size(); k++) {
        GO gid_k = pack.element2node[i][k];
        LO lid_k = pack.overlapMap->getLocalElement(gid_k);
        mat1.sumIntoLocalValues(lid_j,1,&localValues[j][k],&lid_k);
        mat2.sumIntoGlobalValues(gid_j,1,&localValues[j][k],&gid_k);
      }
    }
  }
  mat1.endFill();
  mat2.fillComplete();

  success = compare<Scalar,LO,GO,Node>::compare_final_matrix_structure(out,mat1,mat2);
  TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success)
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FECrsMatrix, Assemble1D_LocalIndex_Kokkos, LO, GO, Scalar, Node )
{
  using exec_space = typename Node::execution_space;
  using range_type = Kokkos::RangePolicy<exec_space, LO>;
  using FEMAT = typename Tpetra::FECrsMatrix<Scalar,LO,GO,Node>;
  using CMAT = typename Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
  using FEG = typename Tpetra::FECrsGraph<LO,GO,Node>;
  using ImplScalarType = typename FEMAT::impl_scalar_type;

  // get a comm
  RCP<const Comm<int> > comm = getDefaultComm();
  
  // Generate a mesh
  size_t numLocal = 10;
  GraphPack<LO,GO,Node> pack;
  generate_fem1d_graph(numLocal,comm,pack);
  
  // Make the graph    
  RCP<FEG> graph = rcp(new FEG(pack.uniqueMap,pack.overlapMap,pack.overlapMap,3));

  graph->beginFill();
  for(size_t i=0; i<(size_t)pack.element2node.size(); i++) {
    for(size_t j=0; j<pack.element2node[i].size(); j++) {
      GO gid_j = pack.element2node[i][j];
      LO lid_j = pack.overlapMap->getLocalElement(gid_j);
      for(size_t k=0; k<pack.element2node[i].size(); k++) {
        GO gid_k = pack.element2node[i][k];
        LO lid_k = pack.overlapMap->getLocalElement(gid_k);
        //printf("[%d] Inserting gid (%d,%d) lid (%d,%d)\n",comm->getRank(),gid_j,gid_k,lid_j,lid_k);fflush(stdout);
        graph->insertLocalIndices(lid_j,1,&lid_k);
      }
    }
  }
  graph->endFill();

  // Generate the "local stiffness matrix"
  std::vector<std::vector<Scalar> > localValues = generate_fem1d_element_values<Scalar>();
  auto kokkosValues = generate_fem1d_element_values_kokkos<ImplScalarType, Node>();

  // Make the matrix two ways
  FEMAT mat1(graph); // Here we use graph as a FECrsGraph
  CMAT mat2(graph);  // Here we use graph as a CrsGraph in OWNED mode
  mat1.beginFill();
  auto k_e2n = pack.k_element2node;
  auto localMat = mat1.getLocalMatrix();
  auto localMap = pack.overlapMap->getLocalMap();

  Kokkos::parallel_for("assemble_1d_local_index", 
		       range_type (0, k_e2n.extent(0)),
		       KOKKOS_LAMBDA(const size_t i) {
    for(size_t j=0; j<k_e2n.extent(1); j++) {
      LO lid_j = localMap.getLocalElement(k_e2n(i, j));
      for(size_t k=0; k<k_e2n.extent(1); k++) {
        LO lid_k = localMap.getLocalElement(k_e2n(i, k));
	ImplScalarType tmp = kokkosValues(j, k);
	localMat.sumIntoValues(lid_j, &lid_k, 1, &tmp, true, true);
      }
    }
  });
  mat1.endFill();

  for(size_t i=0; i<(size_t)pack.element2node.size(); i++) {
    for(size_t j=0; j<pack.element2node[i].size(); j++) {
      GO gid_j = pack.element2node[i][j];
      for(size_t k=0; k<pack.element2node[i].size(); k++) {
        GO gid_k = pack.element2node[i][k];
        mat2.sumIntoGlobalValues(gid_j,1,&localValues[j][k],&gid_k);
      }
    }
  }

  mat2.fillComplete();

  success = compare<Scalar,LO,GO,Node>::compare_final_matrix_structure(out,mat1,mat2);
  TPETRA_GLOBAL_SUCCESS_CHECK(out,comm,success)
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FECrsMatrix, Assemble1D_Kokkos, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FECrsMatrix, Assemble1D, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FECrsMatrix, Assemble1D_LocalIndex, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FECrsMatrix, Assemble1D_LocalIndex_Kokkos, LO, GO, SCALAR, NODE ) 


TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

} // end namespace (anonymous)
