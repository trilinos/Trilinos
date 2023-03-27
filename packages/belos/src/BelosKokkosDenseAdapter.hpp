// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Jennifer A. Loe (jloe@sandia.gov)
//                or Heidi K. Thornquist (hkthorn@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
#ifndef BELOS_KOKKOS_DENSE_MAT_TRAITS_HPP
#define BELOS_KOKKOS_DENSE_MAT_TRAITS_HPP

/*! \file BelosKokkosDenseAdapter.hpp
  \brief Full specialization of Belos::DenseMatTraits for Kokkos::DualView
  with arbitrary scalarType, executionSpace, etc. All views are expected to
  have 2 dimensions. 
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

//Kokkos BLAS files:
#include "KokkosBlas1_scal.hpp"
#include "KokkosBlas1_axpby.hpp"

namespace Belos {

  //! Full specialization of Belos::DenseMatTraits for Teuchos::SerialDenseMatrix<int,ST>.
  template<class MatrixType>
  class DenseMatTraits<Kokkos::DualView<MatrixType::DataType,MatrixType::Layout,MatrixType::Device,MatrixType::MemoryTraits>>{
  public:
    
    //@{ \name Creation methods

    /*! \brief Creates a new empty \c Kokkos::DualView<Scalar**,vargs...> with no dimension.

    \return Reference-counted pointer to a new dense matrix of type \c Kokkos::DualView<Scalar**,vargs...>.
    */
    //static Teuchos::RCP<Kokkos::DualView<Scalar**,vargs...>> Create() { 
    static Teuchos::RCP<MatrixType> Create() { 
      return Teuchos::rcp(new MatrixType("BelosDenseView"));
    }     

    /*! \brief Creates a new empty \c Kokkos::DualView<Scalar**,vargs...> containing 
     *     and \c numRows rows and \c numCols columns.
     *         Will be initialized to zeros if last parameter is true.

    \return Reference-counted pointer to a new dense matrix of type \c Kokkos::DualView<Scalar**,vargs...>.
    */
    static Teuchos::RCP<Kokkos::DualView<Scalar**,vargs...>> Create( const int numRows, const int numCols, bool initZero = true) { 
      if(initZero){
        return Teuchos::rcp(new Kokkos::DualView<Scalar**,vargs...>("BelosDenseView",numRows,numCols));
      }
      else {
        //return Teuchos::rcp(new Kokkos::DualView<Scalar**,vargs...>(Kokkos::WithoutInitializing,"BelosDenseView",numRows,numCols));
        return Teuchos::rcp(new Kokkos::DualView<Scalar**,vargs...>(Kokkos::WithoutInitializing,numRows,numCols));
      }
    }     

    //! \brief Returns a raw pointer to the data on the host.
    static Scalar* GetRawHostPtr(const Kokkos::DualView<Scalar**,vargs...> & dm ) { 
    //LAPACK could use the host ptr to modify entries, so mark as modified.
    //TODO: Is there a better way to handle this?
        dm.modify_host();
        return dm.h_view.data();
    }

    //! \brief Returns an RCP to a Kokkos::DualView<Scalar**,vargs...> which has a subview of the given Kokkos::DualView<Scalar**,vargs...>.
    //        Row and column indexing is zero-based.
    static Teuchos::RCP<Kokkos::DualView<Scalar**,vargs...>> 
    Subview( Kokkos::DualView<Scalar**,vargs...> & source, int numRows, int numCols, int startRow=0, int startCol=0){
      return Kokkos::subview(source, 
          Kokkos::pair<int,int>(startRow,startRow+numRows), Kokkos::pair<int,int>(startCol,startCol+numCols)); 
      //Does subview work on dual views? Brian says yes. 
    }     

    static Teuchos::RCP<const Kokkos::DualView<Scalar**,vargs...>> 
    SubviewConst( const Kokkos::DualView<Scalar**,vargs...>& source, int numRows, int numCols, int startRow=0, int startCol=0){
      return Kokkos::subview(source, 
          Kokkos::pair<int,int>(startRow,startRow+numRows), Kokkos::pair<int,int>(startCol,startCol+numCols)); 
          //TODO: Check const-ness semantics here?
    }     

    //! \brief Returns a deep copy of the requested subview.
    static Teuchos::RCP<Kokkos::DualView<Scalar**,vargs...>> 
    SubviewCopy( const Kokkos::DualView<Scalar**,vargs...>& source, int numRows, int numCols, int startRow=0, int startCol=0){ 
      //Maaybe we could get away with just a host copy here??
      //Hmmm... but it says we return a dual view. 
      //Maybe it should return a dual view with only host data copied in. Require sync to work on device. 
      //This is related to the functionality of the Assign function. 
      auto tmpViewRCP = Teuchos::rcp(new Kokkos::DualView<Scalar**,vargs...>
                                      (Kokkos::WithoutInitializing,"BelosDenseView",numRows,numCols));
      Kokkos::deep_copy(*tmpViewRCP, Kokkos::subview(source, 
          Kokkos::pair<int,int>(startRow,startRow+numRows), Kokkos::pair<int,int>(startCol,startCol+numCols))); 
      return tmpViewRCP;
    }     
    //@}

    //@{ \name Attribute methods

    //! \brief Obtain the number of rows of \c dm.
    static int GetNumRows( const Kokkos::DualView<Scalar**,vargs...>& dm ) { 
      return dm.extent(0); 
    }     

    //! \brief Obtain the number of columns of \c dm.
    static int GetNumCols( const Kokkos::DualView<Scalar**,vargs...>& dm ) { 
      return dm.extent(1); 
    }     

    //! \brief Obtain the stride between the columns of \c dm.
    static int GetStride( const Kokkos::DualView<Scalar**,vargs...>& dm ) { 
      int strides[2];
      dm.stride(strides);
      return strides[0]; //TODO check this is the right dimension. 
    }

    //@}

    //@{ \name Shaping methods

    /* \brief Reshaping method for changing the size of \c dm to have \c numRows rows and \c numCols columns.
     *        All values will be initialized to zero if the final argument is true. 
     *        If the final argument is fale, the previous entries in 
     *        the matrix will be maintained. For new entries that did not exist in the previous matrix, values will
     *        contain noise from memory. 
    */
    static void Reshape( Kokkos::DualView<Scalar**,vargs...>& dm, const int numRows, const int numCols, bool initZero = false) {
      if(initZero){
        dm.realloc(numRows,numCols); //changes size of both host and device view.
        Kokkos::deep_copy(dm.h_view, 0.0);
        dm.modify_host();
      }
      else{
        dm.resize(numRows,numCols); //keeps values in old array.
      }
    }     

    //@}

    //@{ \name Data access methods

    //! \brief Access a reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static Scalar & Value( Kokkos::DualView<Scalar**,vargs...>& dm, const int i, const int j )
    { 
    //Mark as modified on host, since we don't know if it will be. 
      dm.modify_host();
      return dm.h_view(i,j);
    }

    //! \brief Access a const reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static const Scalar & Value( const Kokkos::DualView<Scalar**,vargs...>& dm, const int i, const int j ) { 
      return dm.h_view(i,j);
      //TODO check const semantics here?
    }

    static void SyncHostToDevice(Kokkos::DualView<Scalar**,vargs...> &)
    { }

    static void SyncDeviceToHost(Kokkos::DualView<Scalar**,vargs...> &)
    { }
    //@}
    //@{ \name Operator methods
    
    //!  \brief Adds sourceDM to thisDM and returns answer in thisDM.
    static void Add( Kokkos::DualView<Scalar**,vargs...>& thisDM, const Kokkos::DualView<Scalar**,vargs...>& sourceDM) {
      KokkosBlas::axpy(1.0,sourceDM.h_view, thisDM.h_view); //axpy(alpha,x,y), y = y + alpha*x
      thisDM.modify_host();
    }

    //!  \brief Fill all entries with \c value. Value is zero if not specified.
    static void PutScalar( Kokkos::DualView<Scalar**,vargs...>& dm, Scalar value = Teuchos::ScalarTraits<Scalar>::zero()){ 
      Kokkos::deep_copy( dm.h_view, value);
      dm.modify_host();
    }

    //!  \brief Multiply all entries by a scalar. DM = value.*DM
    static void Scale( Kokkos::DualView<Scalar**,vargs...>& dm, Scalar value) { 
      KokkosBlas::scal( dm.h_view, value, dm.h_view);
      dm.modify_host();
    }

    //!  \brief Fill the Kokkos::DualView with random entries.
    //!   Entries are assumed to be the same on each MPI rank (each matrix copy). 
    static void Randomize( Kokkos::DualView<Scalar**,vargs...>& dm)
    { }

    //TODO really need random synced. 

    //!  \brief Copies entries of source to dest (deep copy). 
    static void Assign( Kokkos::DualView<Scalar**,vargs...>& dest, const Kokkos::DualView<Scalar**,vargs...>& source) { 
      //Kokkos::deep_copy(dest,source); //Brian Kelley says this works on dual views.
      //TODO: But do we really want to do this? What if we only need the host pieces copied right now?
      Kokkos::deep_copy(dest.h_view,source.h_view); 
      dest.modify_host();
    }

    //!  \brief Returns the Frobenius norm of the dense matrix.
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormFrobenius( Kokkos::DualView<Scalar**,vargs...>& dm) { 
      using KAT = Kokkos::ArithTraits<Scalar>;
      using mag_t = typename KAT::mag_type;
      mag_t frobNorm;
      Kokkos::parallel_reduce(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {dm.extent(0), dm.extent(1)}),
          KOKKOS_LAMBDA(size_t i, size_t j, mag_t& lfrobNorm)
          {
          mag_t absVal = KAT::abs(dm(i, j));
          lfrobNorm += absVal * absVal;
          }, frobNorm);
      return Kokkos::sqrt(frobNorm);
    }
    //@}
  };
  
} // namespace Belos

#endif // end file BELOS_TEUCHOS_DENSE_MAT_TRAITS_HPP
