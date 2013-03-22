// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2011 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

/**
  \file   Amesos2_MultiVecAdapter.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Wed Feb 10 14:53:29 2010

  \brief  A templated adapter/wrapper class for Trilinos Multivector
          type classes.  Provides the functions necessary for Amesos2
          to function.  Other wrapper methods may of course be added
          for special cases.
*/

#ifndef AMESOS2_MULTIVEC_ADAPTER_DECL_HPP
#define AMESOS2_MULTIVEC_ADAPTER_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Tpetra_Map.hpp>

#include "Amesos2_TypeDecl.hpp"

namespace Amesos2 {


  /**
   * \brief A templated MultiVector class adapter for Amesos2
   *
   * Specializations of this templated class provide a unified interface
   * to MultiVector types for Amesos2.  Any specializations are expected to
   * implement the following methods:
   *
   * <br>
   * <b>Implementation Requirements:</b>
   * <ul>
   * <li>Wrapper constructor
   * \code MultiVecAdapter<MultiVecType>(const Teuchos::RCP<MultiVecType>& mat); \endcode
   * </li>
   * 
   * <li> Method to get locality of multivec, either globally or locally indexed.
   *
   * \code
   * bool isLocallyIndexed() const;
   * bool isGloballyIndexed() const;
   *
   * Teuchos::RCP<Tpetra::Map<LO,GO,Node> > getMap() const;
   * \endcode
   * </li>
   *
   * <li> Method to get multi-vector communicator.
   *
   * \code
   * const Teuchos::RCP<const Teuchos::Comm<int> >& getComm() const;
   * \endcode
   * </li>
   *
   * <li> Methods to get the local and global length of vectors and number of
   * vectors.
   *
   * \code
   * size_t getLocalLength();
   *
   * size_t  getLocalNumVectors();
   *
   * global_size_type getGlobalLength();
   *
   * global_size_type getGlobalNumVectors();
   * \endcode
   * </li>
   *
   * <li> Get access to multi-vector stride
   *
   * \code
   * size_t getStride() const;
   * \endcode
   * </li>
   *
   * <li> Ask whether this multi-vector has constant stride between vectors on
   * this node.
   *
   * \code
   * bool isConstantStride() const;
   * \endcode
   * </li>
   *
   * <li> Vector access methods.
   *
   * \code
   * Teuchos::RCP<const Tpetra::Vector<Scalar,LO,GO,Node> >
   * getVector( size_t j ) const;
   *
   * Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO,Node> >
   * getVectorNonConst( size_t j );
   * \endcode
   * </li>
   *
   * <li> Data access methods
   *
   * \code
   * void get1dCopy( const Teuchos::ArrayView<scalar_type>& A, size_t lda) const;
   *
   * void get2dCopy( Teuchos::ArrayView<const Teuchos::ArrayView<scalar_type> > A ) const;
   * \endcode
   * </li>
   *
   * <li> Method to export an array of new values into the global multi-vector.
   *
   * \code
   * template<typename Value_t>
   * void globalize( const Teuchos::ArrayView<Value_t>& newVals )
   * \endcode
   * </li>
   * 
   * <li> Get a description of this adapter.
   *
   * \code
   * std::string description() const;
   * \endcode
   * </li>
   *
   * <li> Print the multivec to the \c os output stream.
   *
   * \code
   * void describe(
   *   Teuchos::FancyOStream& os,
   *   const Teuchos::EVerbosityLevel verblevel) const;
   * \endcode
   * </li>
   *
   * \ingroup amesos2_multivec_adapters 
   */
  template <class MV>
  struct MultiVecAdapter {};

  /** \brief Factory creation method for MultiVecAdapters
   *
   * Developers should favor this method for creating Amesos2
   * MultiVector adapters over using the constructors.
   *
   * \relates MultiVecAdapter
   */
  template <class MV>
  Teuchos::RCP<MultiVecAdapter<MV> >
  createMultiVecAdapter(Teuchos::RCP<MV> mv){
    using Teuchos::rcp;
    
    if(mv.is_null()) return Teuchos::null;
    return( rcp(new MultiVecAdapter<MV>(mv)) );
  }

  template <class MV>
  Teuchos::RCP<const MultiVecAdapter<MV> >
  createConstMultiVecAdapter(Teuchos::RCP<const MV> mv){
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    
    if(mv.is_null()) return Teuchos::null;
    return( rcp(new MultiVecAdapter<MV>(Teuchos::rcp_const_cast<MV,const MV>(mv))).getConst() );
  }


  ///////////////////////////////////////////////////////////
  // Utilities for getting and putting data from MultiVecs //
  ///////////////////////////////////////////////////////////
  
  namespace Util {
    
    /*
     * If the multivector scalar type and the desired scalar tpye are
     * the same, then we can do a simple straight copy.
     */
    template <typename MV>
    struct same_type_get_copy {
      static void apply(const Teuchos::Ptr<const MV>& mv,
			const Teuchos::ArrayView<typename MV::scalar_t>& v,
			const size_t ldx,
			Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map );
    };

    /*
     * In the case where the scalar type of the multi-vector and the
     * corresponding S type are different, then we need to first get a
     * copy of the scalar values, then convert each one into the S
     * type before inserting into the vals array.
     */
    template <typename MV, typename S>
    struct diff_type_get_copy {
      static void apply(const Teuchos::Ptr<const MV>& mv,
			const Teuchos::ArrayView<S>& v,
			const size_t& ldx,
			Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map );
    };

    /** \internal
     * 
     * \brief Helper class for getting 1-D copies of multivectors
     *
     * Handles datatype conversion when appropriate.
     */
    template <class MV, typename S>
    struct get_1d_copy_helper {
      static void do_get(const Teuchos::Ptr<const MV>& mv,
			 const Teuchos::ArrayView<S>& vals,
			 const size_t ldx,
			 Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map );
      
      static void do_get(const Teuchos::Ptr<const MV>& mv,
			 const Teuchos::ArrayView<S>& vals,
			 const size_t ldx,
			 EDistribution distribution, typename MV::global_ordinal_t indexBase = 0);

      static void do_get(const Teuchos::Ptr<const MV>& mv,
			 const Teuchos::ArrayView<S>& vals,
			 const size_t ldx);
    };

    /*
     * If the multivector scalar type and the desired scalar tpye are
     * the same, then we can do a simple straight copy.
     */
    template <typename MV>
    struct same_type_data_put {
      static void apply(const Teuchos::Ptr<MV>& mv,
			const Teuchos::ArrayView<typename MV::scalar_t>& data,
			const size_t ldx,
			Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map );
    };

    /*
     * In the case where the scalar type of the multi-vector and the
     * corresponding S type are different, then we need to first get a
     * copy of the scalar values, then convert each one into the S
     * type before inserting into the vals array.
     */
    template <typename MV, typename S>
    struct diff_type_data_put {
      static void apply(const Teuchos::Ptr<MV>& mv,
			const Teuchos::ArrayView<S>& data,
			const size_t& ldx,
			Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map );
    };
  
    /** \internal
     * 
     * \brief Helper class for putting 1-D data arrays into multivectors
     *
     * Handles dataype conversion when necessary before putting the data.
     */
    template <class MV, typename S>
    struct put_1d_data_helper {
      static void do_put(const Teuchos::Ptr<MV>& mv,
			 const Teuchos::ArrayView<S>& data,
			 const size_t ldx,
			 Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map );
      
      static void do_put(const Teuchos::Ptr<MV>& mv,
			 const Teuchos::ArrayView<S>& data,
			 const size_t ldx,
			 EDistribution distribution, typename MV::global_ordinal_t indexBase = 0);

      static void do_put(const Teuchos::Ptr<MV>& mv,
			 const Teuchos::ArrayView<S>& data,
			 const size_t ldx);
    };
  }
} // end namespace Amesos2

#include "Amesos2_TpetraMultiVecAdapter_decl.hpp"
#ifdef HAVE_AMESOS2_EPETRA
#  include "Amesos2_EpetraMultiVecAdapter_decl.hpp"
#endif

#endif  // AMESOS2_MULTIVEC_ADAPTER_DECL_HPP
