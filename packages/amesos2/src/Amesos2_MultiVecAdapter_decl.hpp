// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Amesos2_VectorTraits.hpp"

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

    /** \internal
     *
     * \brief Helper struct for getting pointers to the MV
     * data - only used when number of vectors = 1
     * and single MPI process
     */
    template <typename MV, typename V>
    struct vector_pointer_helper {

      typedef typename VectorTraits<V>::ptr_scalar_type ptr_return_type ;

      static ptr_return_type * get_pointer_to_vector ( const Teuchos::Ptr< MV> &mv ) ;

      static ptr_return_type * get_pointer_to_vector ( Teuchos::Ptr< MV> &mv ) ;

      static ptr_return_type * get_pointer_to_vector ( const Teuchos::Ptr< const MV > &mv ) ;

      static ptr_return_type * get_pointer_to_vector ( Teuchos::Ptr< const MV > &mv ) ;
    };

    /*
     * If the multivector scalar type and the desired scalar tpye are
     * the same, then we can do a simple straight copy.
     */
    template <typename MV>
    struct same_type_get_copy {
      static void apply(const Teuchos::Ptr<const MV> mv,
                        const Teuchos::ArrayView<typename MV::scalar_t>& v,
                        const size_t ldx,
                        Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                        EDistribution distribution );
    };

    /*
     * In the case where the scalar type of the multi-vector and the
     * corresponding S type are different, then we need to first get a
     * copy of the scalar values, then convert each one into the S
     * type before inserting into the vals array.
     */
    template <typename MV, typename S>
    struct diff_type_get_copy {
      static void apply(const Teuchos::Ptr<const MV> mv,
                        const Teuchos::ArrayView<S>& v,
                        const size_t ldx,
                        Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                        EDistribution distribution );
    };

    /** \internal
     *
     * \brief Helper class for getting 1-D copies of multivectors
     *
     * Handles datatype conversion when appropriate.
     */
    template <class MV, typename S>
    struct get_1d_copy_helper {
      static void
      do_get (const Teuchos::Ptr<const MV>& mv,
              const Teuchos::ArrayView<S>& vals,
              const size_t ldx,
              Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
              EDistribution distribution = ROOTED);

      static void
      do_get (const Teuchos::Ptr<const MV>& mv,
              const Teuchos::ArrayView<S>& vals,
              const size_t ldx,
              EDistribution distribution,
              typename MV::global_ordinal_t indexBase = 0);

      static void
      do_get (const Teuchos::Ptr<const MV>& mv,
              const Teuchos::ArrayView<S>& vals,
              const size_t ldx);
    };

    /*
      do_get

      Return type (bool):
        true:  The input kokkos_vals view is now pointing directly to the adapter's data (same memory and type).
               If this is x for an Ax=b solve, you don't need 'do_put x' after the solve since you modified the adapter directly.
        false: The input kokkos_vals view is now resized to match the adapter's size.
               kokkos_vals will only have the adapter values deep_copied if bInitialize is true (see below).
               If this is x for an Ax=b solve, you must call 'do_put x' after the solve to deep copy back to the adapter.

      Inputs
        bInitialize (bool): tells the adapter whether kokkos_vals needs to have the values of the adapter.
          true:  We require kokkos_vals to have the same size and values of the adapter.
                 For b in Ax=b solves, set bInitialize true because you need the size and values of the adapter.
          false: We require kokkos_vals to have the same size as the adapter but we don't need the values.
                 For x in Ax=b solves, set bInitialize false because you just need the size, not the values.

          Note: When this method returns true, meaning direct assignment of the view occurred,
                bInitialize is not used because you already have the values whether you need them or not.

        kokkos_vals (View<scalar_t**>): The view which will contain the x or b data.
          Do not allocate the size of kokkos_vals, let the do_get method do it for you.
          This is because kokkos_vals may be set to point directly to the adapter memory
          and then any pre-allocation of size will have been wasted.
    */
    template <class MV, typename KV>
    struct get_1d_copy_helper_kokkos_view {
      static bool
      do_get (bool bInitialize,
              const Teuchos::Ptr<const MV>& mv,
              KV& kokkos_vals,
              const size_t ldx,
              Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
              EDistribution distribution = ROOTED);

      static bool
      do_get (bool bInitialize,
              const Teuchos::Ptr<const MV>& mv,
              KV& kokkos_vals,
              const size_t ldx,
              EDistribution distribution,
              typename MV::global_ordinal_t indexBase = 0);

      static bool
      do_get (bool bInitialize,
              const Teuchos::Ptr<const MV>& mv,
              KV& kokkos_vals,
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
                        Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                        EDistribution distribution );
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
                        const size_t ldx,
                        Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                        EDistribution distribution );
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
                         Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                         EDistribution distribution = ROOTED);

      static void do_put(const Teuchos::Ptr<MV>& mv,
                         const Teuchos::ArrayView<S>& data,
                         const size_t ldx,
                         EDistribution distribution, typename MV::global_ordinal_t indexBase = 0);

      static void do_put(const Teuchos::Ptr<MV>& mv,
                         const Teuchos::ArrayView<S>& data,
                         const size_t ldx);
    };

    template <class MV, typename KV>
    struct put_1d_data_helper_kokkos_view {
      static void do_put(const Teuchos::Ptr<MV>& mv,
                         KV& kokkos_data,
                         const size_t ldx,
                         Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                         EDistribution distribution = ROOTED);

      static void do_put(const Teuchos::Ptr<MV>& mv,
                         KV& kokkos_data,
                         const size_t ldx,
                         EDistribution distribution, typename MV::global_ordinal_t indexBase = 0);

      static void do_put(const Teuchos::Ptr<MV>& mv,
                         KV& kokkos_data,
                         const size_t ldx);
    };
  }
} // end namespace Amesos2

#include "Amesos2_TpetraMultiVecAdapter_decl.hpp"
#include "Amesos2_KokkosMultiVecAdapter_decl.hpp"
#ifdef HAVE_AMESOS2_EPETRA
#  include "Amesos2_EpetraMultiVecAdapter_decl.hpp"
#endif

#endif  // AMESOS2_MULTIVEC_ADAPTER_DECL_HPP
