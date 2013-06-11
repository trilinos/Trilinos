// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_KokkosArrayKernelsUnitTestNew.hpp"

#include "Stokhos_Host_CrsMatrix.hpp"
#include "Stokhos_Host_BlockCrsMatrix.hpp"
#include "Stokhos_Host_StochasticProductTensor.hpp"
#include "Stokhos_Host_CrsProductTensor.hpp"
#include "Stokhos_Host_FlatSparse3Tensor.hpp"
#include "Stokhos_Host_FlatSparse3Tensor_kji.hpp"

#include "Stokhos_LexicographicBlockSparse3Tensor.hpp"
#include "Stokhos_Host_LexicographicBlockSparse3Tensor.hpp"

#include "KokkosArray_hwloc.hpp"
#include "KokkosArray_Cuda.hpp"

using namespace KokkosArrayKernelsUnitTest;

UnitTestSetup setup;

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsMatrixFree_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  typedef Stokhos::DefaultSparseMatOps SparseMatOps;
  bool test_block = true;

  success = test_crs_matrix_free<Scalar,Device,SparseMatOps>(
    setup, test_block, out);
}

#ifdef HAVE_STOKHOS_MKL

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsMatrixFree_HostMKL ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  typedef Stokhos::MKLSparseMatOps SparseMatOps;
  bool test_block = true;

  success = test_crs_matrix_free<Scalar,Device,SparseMatOps>(
    setup, test_block, out);
}

#endif

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsProductTensor_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;

  success = test_crs_product_tensor<Scalar,Stokhos::CrsProductTensor<Scalar,Device>,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, FlatSparse3Tensor_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;

  success = test_crs_product_tensor<Scalar,Stokhos::FlatSparse3Tensor<Scalar,Device>,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, FlatSparse3Tensor_kji_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;

  success = test_crs_product_tensor<Scalar,Stokhos::FlatSparse3Tensor_kji<Scalar,Device>,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, ProductTensorCijk ) {
  success = true;

  typedef double value_type;
  typedef KokkosArray::Host Device;
  typedef Stokhos::CrsProductTensor< value_type , Device > tensor_type ;

  tensor_type tensor =
   Stokhos::create_product_tensor<Device>( *setup.basis, *setup.Cijk );

  for (int i=0; i<setup.stoch_length; ++i) {
    const size_t iEntryBeg = tensor.entry_begin(i);
    const size_t iEntryEnd = tensor.entry_end(i);
    for (size_t iEntry = iEntryBeg ; iEntry < iEntryEnd ; ++iEntry ) {
      const size_t j = tensor.coord(iEntry,0);
      const size_t k = tensor.coord(iEntry,1);
      double c2 = tensor.value(iEntry);
      if (j == k) c2 *= 2.0;

      int ii = setup.inv_perm[i];
      int jj = setup.inv_perm[j];
      int kk = setup.inv_perm[k];
      double c = setup.Cijk->getValue(ii,jj,kk);

      if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
        out << "(" << ii << "," << jj << "," << kk << "):  " << c
            << " == " << c2 << " failed!" << std::endl;
        success = false;
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, FlatSparseCijk ) {
  success = true;

  typedef double value_type;
  typedef KokkosArray::Host Device;
  typedef Stokhos::FlatSparse3Tensor< value_type , Device > tensor_type ;
  typedef size_t size_type;

  tensor_type tensor =
   Stokhos::create_flat_sparse_3_tensor<Device>( *setup.basis, *setup.Cijk );

  for (int i=0; i<setup.stoch_length; ++i) {
    const size_type nk = tensor.num_k(i);
    const size_type kBeg = tensor.k_begin(i);
    const size_type kEnd = kBeg + nk;
    for (size_type kEntry = kBeg; kEntry < kEnd; ++kEntry) {
      const size_type k = tensor.k_coord(kEntry);
      const size_type nj = tensor.num_j(kEntry);
      const size_type jBeg = tensor.j_begin(kEntry);
      const size_type jEnd = jBeg + nj;
      for (size_type jEntry = jBeg; jEntry < jEnd; ++jEntry) {
        const size_type j = tensor.j_coord(jEntry);
        double c2 = tensor.value(jEntry);
        if (j == k) c2 *= 2.0;
        double c = setup.Cijk->getValue(i,j,k);
        if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
          out << "(" << i << "," << j << "," << k << "):  " << c
              << " == " << c2 << " failed!" << std::endl;
          success = false;
        }
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, FlatSparseCijk_kji ) {
  success = true;

  typedef double value_type;
  typedef KokkosArray::Host Device;
  typedef Stokhos::FlatSparse3Tensor_kji< value_type , Device > tensor_type ;
  typedef size_t size_type;

  tensor_type tensor =
   Stokhos::create_flat_sparse_3_tensor_kji<Device>(*setup.basis, *setup.Cijk);
  const size_type nk = tensor.num_k();

  for ( size_type k = 0; k < nk; ++k) {
    const size_type nj = tensor.num_j(k);
    const size_type jBeg = tensor.j_begin(k);
    const size_type jEnd = jBeg + nj;
    for (size_type jEntry = jBeg; jEntry < jEnd; ++jEntry) {
      const size_type j = tensor.j_coord(jEntry);
      const size_type ni = tensor.num_i(jEntry);
      const size_type iBeg = tensor.i_begin(jEntry);
      const size_type iEnd = iBeg + ni;
      for (size_type iEntry = iBeg; iEntry < iEnd; ++iEntry) {
        const size_type i = tensor.i_coord(iEntry);
        double c2 = tensor.value(iEntry);
        if (j == k) c2 *= 2.0;
        double c = setup.Cijk->getValue(i,j,k);
        if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
          out << "(" << i << "," << j << "," << k << "):  " << c
              << " == " << c2 << " failed!" << std::endl;
          success = false;
        }
      }
    }
  }
}

namespace KokkosArrayKernelsUnitTest {

  template< typename ScalarType , class Device >
  bool test_lexo_block_tensor(const UnitTestSetup& setup,
                              Teuchos::FancyOStream& out) {
    typedef ScalarType value_type ;
    typedef int ordinal_type;
    typedef KokkosArray::View< value_type** , KokkosArray::LayoutLeft ,
                               Device > block_vector_type ;
    typedef Stokhos::LexicographicBlockSparse3Tensor<value_type,Device> TensorType;
    typedef Stokhos::StochasticProductTensor< value_type , TensorType , Device > tensor_type ;

    typedef Stokhos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
    typedef typename matrix_type::graph_type graph_type ;

    //------------------------------
    // Generate input multivector:

    block_vector_type x =
      block_vector_type( "x" , setup.stoch_length , setup.fem_length );
    block_vector_type y =
      block_vector_type( "y" , setup.stoch_length , setup.fem_length );

    typename block_vector_type::HostMirror hx =
      KokkosArray::create_mirror( x );

    for ( int iColFEM = 0 ;   iColFEM < setup.fem_length ;   ++iColFEM ) {
      for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
        hx(iColStoch,iColFEM) =
          generate_vector_coefficient<ScalarType>(
            setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
      }
    }

    KokkosArray::deep_copy( x , hx );

    //------------------------------

    matrix_type matrix ;

    /*
    typedef UnitTestSetup::abstract_basis_type abstract_basis_type;
    typedef UnitTestSetup::basis_type basis_type;
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > less_type;
    typedef Stokhos::TotalOrderBasis<ordinal_type,value_type,less_type> product_basis_type;
    Teuchos::Array< Teuchos::RCP<const abstract_basis_type> > bases(setup.d);
    for (int i=0; i<setup.d; i++)
      bases[i] = rcp(new basis_type(setup.p,true));
    product_basis_type basis(bases, 1e-12);
    */
    const bool symmetric = true;
    Teuchos::RCP< Stokhos::LTBSparse3Tensor<ordinal_type, value_type> > Cijk =
      Stokhos::computeTripleProductTensorLTBBlockLeaf(*setup.basis, symmetric);

    matrix.block =
      Stokhos::create_stochastic_product_tensor< TensorType >( *setup.basis,
                                                               *Cijk );

    matrix.graph = KokkosArray::create_crsarray<graph_type>(
      std::string("test crs graph") , setup.fem_graph );

    matrix.values = block_vector_type(
      "matrix" , setup.stoch_length , setup.fem_graph_length );

    typename block_vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {
      for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
        const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM] ;

        for ( int k = 0 ; k < setup.stoch_length ; ++k ) {
          hM(k,iEntryFEM) =
            generate_matrix_coefficient<ScalarType>(
              setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , k );
          //hM(k,iEntryFEM) = 1.0;
        }
      }
    }

    KokkosArray::deep_copy( matrix.values , hM );

    //------------------------------

    Stokhos::multiply( matrix , x , y );

    typename block_vector_type::HostMirror hy = KokkosArray::create_mirror( y );
    KokkosArray::deep_copy( hy , y );

    bool success = setup.test_commuted_no_perm(hy, out);
    //bool success = true;
    return success;
  }

}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, LexoBlockTensor_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;

  success = test_lexo_block_tensor<Scalar,Device>(setup, out);
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  setup.setup();

  // Initialize host
  const std::pair<unsigned,unsigned> core_topo =
    KokkosArray::hwloc::get_core_topology();
  //const size_t core_capacity = KokkosArray::hwloc::get_core_capacity();

  const size_t gang_count = core_topo.first ;
  const size_t gang_worker_count = core_topo.second;
  KokkosArray::Host::initialize( gang_count , gang_worker_count );

#ifdef HAVE_KOKKOSARRAY_CUDA
  // Initialize Cuda
  KokkosArray::Cuda::initialize( KokkosArray::Cuda::SelectDevice(0) );
#endif

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  KokkosArray::Host::finalize();
#ifdef HAVE_KOKKOSARRAY_CUDA
  KokkosArray::Cuda::finalize();
#endif

  return ret;
}
