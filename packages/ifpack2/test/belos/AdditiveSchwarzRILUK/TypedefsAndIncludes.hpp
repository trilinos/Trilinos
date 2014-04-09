#ifndef __IFPACK2_TEST_ADDITIVESCHWARZ_RILUK_TYPEDEFS_AND_INCLUDES_HPP
#define __IFPACK2_TEST_ADDITIVESCHWARZ_RILUK_TYPEDEFS_AND_INCLUDES_HPP

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

typedef double scalar_type;
typedef scalar_type ST;
typedef int LO;
typedef int GO;

typedef Teuchos::ScalarTraits<scalar_type> STS;
typedef STS::magnitudeType magnitude_type;
typedef Teuchos::ScalarTraits<magnitude_type> STM;

typedef Kokkos::DefaultNode::DefaultNodeType node_type;
typedef Tpetra::Map<LO, GO, node_type> map_type;
typedef Tpetra::MultiVector<scalar_type, LO, GO, node_type> multivector_type;
typedef Tpetra::CrsMatrix<scalar_type, LO, GO, node_type> sparse_mat_type;

#endif // __IFPACK2_TEST_ADDITIVESCHWARZ_RILUK_TYPEDEFS_AND_INCLUDES_HPP
