// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_MEAN_BASED_PRECONDITIONER_HPP
#define STOKHOS_MEAN_BASED_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "SGPreconditioner.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Tpetra_MultiVector.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "Kokkos_CrsMatrix.hpp"
#include "mean_func.hpp"


namespace Sacado {
  namespace UQ {
    template <typename Storage> class PCE;
  }
}

namespace Kokkos {
namespace Example {    
  /*! 
   * \brief A stochastic preconditioner based on applying the inverse of the
   * mean.
   */
  template<class S, class LO, class GO, class N>
  class MeanBasedPreconditioner : 
    public SGPreconditioner<S, LO, GO, N> {
      
  public:

    //! Constructor 
    MeanBasedPreconditioner() {}
    
    //! Destructor
    virtual ~MeanBasedPreconditioner() {}

    //! Setup preconditioner - construct a MueLu hierarchy with the given matrix
    virtual 
    Teuchos::RCP<Tpetra::Operator<S,LO,GO,N> > 
    setupPreconditioner(const Teuchos::RCP<Tpetra::CrsMatrix<S,LO,GO,N> >& A, 
			const std::string& xmlFileName){

      typedef MueLu::TpetraOperator<S,LO,GO,N> PreconditionerType;
      RCP<PreconditionerType> mueluPreconditioner;
      mueluPreconditioner = MueLu::CreateTpetraPreconditioner(A,xmlFileName);
      return mueluPreconditioner;

    };
 

  private:
    
    //! Private to prohibit copying
    MeanBasedPreconditioner(const MeanBasedPreconditioner&);
    
    //! Private to prohibit copying
    MeanBasedPreconditioner& operator=(const MeanBasedPreconditioner&);
    

  }; // class MeanBasedPreconditioner

  template<class Storage, class LO, class GO, class N>
  class MeanBasedPreconditioner< Sacado::UQ::PCE<Storage>, LO, GO, N > :
    public SGPreconditioner<Sacado::UQ::PCE<Storage>, LO, GO, N> {
 
  public:
    typedef Sacado::UQ::PCE<Storage> Scalar;
 
    //! Constructor 
    MeanBasedPreconditioner() {}

    //! Destructor
    virtual ~MeanBasedPreconditioner() {}

    //! Setup preconditioner 
    virtual
    Teuchos::RCP<Tpetra::Operator<Scalar,LO,GO,N> >
    setupPreconditioner(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,N> >& A,
                        const std::string& xmlFileName){
      using Teuchos::ArrayView;
      using Teuchos::Array;
      typedef Tpetra::CrsMatrix<Scalar,LO,GO,N> MatrixType;
      typedef Tpetra::Map<LO,GO,N> Map;
      typedef MueLu::TpetraOperator<Scalar,LO,GO,N> PreconditionerType;

      typedef Kokkos::CrsMatrix<Scalar, LO, typename N::device_type, void, size_t> KokkosMatrixType;

      typedef typename KokkosMatrixType::StaticCrsGraphType KokkosGraphType;
      typedef typename KokkosMatrixType::values_type KokkosMatrixValuesType; 
      typedef typename KokkosMatrixValuesType::cijk_type cijk_type;

      RCP< const Map > rmap = A->getRowMap();
      RCP< const Map > cmap = A->getColMap();

      KokkosMatrixType  kokkos_matrix = A->getKokkosMatrix();
      size_t nnz = kokkos_matrix.nnz();
      size_t ncols = kokkos_matrix.numCols();
      KokkosGraphType kokkos_graph = kokkos_matrix.graph;

      cijk_type cijk = kokkos_matrix.values.cijk();

      KokkosMatrixValuesType matrix_values = kokkos_matrix.values;
      KokkosMatrixValuesType mean_matrix_values =
        KokkosMatrixValuesType("values", cijk, nnz, 1); 

      Kokkos::Example::GetMeanValsFunc <KokkosMatrixValuesType > getmeanfunc(mean_matrix_values, matrix_values); 
      Kokkos::parallel_for( nnz , getmeanfunc );


      KokkosMatrixType mean_kokkos_matrix("mean-matrix", ncols, mean_matrix_values, kokkos_graph);
      RCP < MatrixType > M = rcp( new MatrixType(rmap, cmap, mean_kokkos_matrix));

      RCP< PreconditionerType > mueluPreconditioner;
      mueluPreconditioner = MueLu::CreateTpetraPreconditioner(M,xmlFileName);
      return mueluPreconditioner;
    };

  private:

    //! Private to prohibit copying
    MeanBasedPreconditioner(const MeanBasedPreconditioner&);

    //! Private to prohibit copying
    MeanBasedPreconditioner& operator=(const MeanBasedPreconditioner&);


  }; // class MeanBasedPreconditioner



}
}

#endif // STOKHOS_MEAN_BASED_PRECONDITIONER_HPP
