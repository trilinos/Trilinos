// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_Util.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu May 27 13:11:13 CDT 2010

  \brief  Utility functions for Amesos2
*/

#ifndef AMESOS2_UTIL_HPP
#define AMESOS2_UTIL_HPP

#include <cstdio>
#include <fstream>
#include <iostream>

#include "Amesos2_config.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_FancyOStream.hpp"

#include <Tpetra_Map.hpp>
#include <Tpetra_DistObject_decl.hpp>
#include <Tpetra_ComputeGatherMap.hpp> // added for gather map... where is the best place??

#include "Amesos2_TypeDecl.hpp"
#include "Amesos2_Meta.hpp"
#include "Amesos2_Kokkos_View_Copy_Assign.hpp"

#ifdef HAVE_AMESOS2_METIS
#include "metis.h" // to discuss, remove from header?
#endif

namespace Amesos2 {

  namespace Util {

    /**
     * \internal
     * \defgroup amesos2_util Amesos2 Utilities
     * @{
     */

    using Teuchos::RCP;
    using Teuchos::ArrayView;

    /**
     * \brief Gets a Tpetra::Map described by the EDistribution.
     *
     * \param distribution The distribution that the returned map will conform to
     * \param num_global_elements A global_size_t value that gives the number of
     *                     global elements in the map.
     * \param comm         The communicator to create the map on.
     *
     * \tparam LO          The local ordinal type
     * \tparam GO          The global ordinal type
     * \tparam GS          The global size type
     * \tparam Node        The Kokkos node type
     *
     * \ingroup amesos2_utils
     */

    template <typename LO, typename GO, typename GS, typename Node>
    const Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >
    getGatherMap( const Teuchos::RCP< const Tpetra::Map<LO,GO,Node> > &map );


    template <typename LO, typename GO, typename GS, typename Node>
    const Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >
    getDistributionMap(EDistribution distribution,
                       GS num_global_elements,
                       const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                       GO indexBase = 0,
                       const Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >& map = Teuchos::null);

    /**
     * Transposes the compressed sparse matrix representation.
     *
     * \ingroup amesos2_utils
     */
    template <typename Scalar,
              typename GlobalOrdinal,
              typename GlobalSizeT>
    void transpose(ArrayView<Scalar> vals,
                   ArrayView<GlobalOrdinal> indices,
                   ArrayView<GlobalSizeT> ptr,
                   ArrayView<Scalar> trans_vals,
                   ArrayView<GlobalOrdinal> trans_indices,
                   ArrayView<GlobalSizeT> trans_ptr);

    /**
     * \brief Scales a 1-D representation of a multivector.
     *
     * \param [in/out] vals The values of the multi-vector.  On exit will contain the scaled values.
     * \param [in] l The length of each vector in the multivector
     * \param [in] ld The leading dimension of the multivector
     * \param [in] s Contains the scaling factors of the diagonal scaling matrix
     *
     * The first vector will be scaled by \c s[0] , the second vector
     * by \c s[1] , etc.
     *
     * \ingroup amesos2_utils
     */
    template <typename Scalar1, typename Scalar2>
    void scale(ArrayView<Scalar1> vals, size_t l,
               size_t ld, ArrayView<Scalar2> s);

    /**
     * \brief Scales a 1-D representation of a multivector.
     *
     * \param [in/out] vals The values of the multi-vector.  On exit will contain the scaled values.
     * \param [in] l The length of each vector in the multivector
     * \param [in] ld The leading dimension of the multivector
     * \param [in] s Contains the scaling factors of the diagonal scaling matrix
     *
     * Scales each vector by diag(s), with the scaling multiplication
     * being performed by the `binary_op' parameter.  BinaryOp is some
     * class that defines a \c operator() method as
     *
     * \code
     * Scalar1 operator()(Scalar1 x, Scalar2 y){  }
     * \endcode
     *
     * \ingroup amesos2_utils
     */
    template <typename Scalar1, typename Scalar2, class BinaryOp>
    void scale(ArrayView<Scalar1> vals, size_t l,
               size_t ld, ArrayView<Scalar2> s, BinaryOp binary_op);


    /// Prints a line of 70 "-"s on std::cout.
    void printLine( Teuchos::FancyOStream &out );

    // Helper function used to convert Kokkos::complex pointer
    // to std::complex pointer; needed for optimized code path
    // when retrieving the CRS raw pointers
    template < class T0, class T1 >
    struct getStdCplxType
    {
      using common_type = typename std::common_type<T0,T1>::type;
      using type = common_type;
    };

    template < class T0, class T1 >
    struct getStdCplxType< T0, T1* >
    {
      using common_type = typename std::common_type<T0,T1>::type;
      using type = common_type;
    };

#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_AMESOS2_KOKKOS)
    template < class T0 >
    struct getStdCplxType< T0, Kokkos::complex<T0>* >
    {
      using type = std::complex<T0>;
    };

    template < class T0 , class T1 >
    struct getStdCplxType< T0, Kokkos::complex<T1>* >
    {
      using common_type = typename std::common_type<T0,T1>::type;
      using type = std::complex<common_type>;
    };
#endif

    //////////////////////////////////
    // Matrix/MultiVector Utilities //
    //////////////////////////////////



    /**
     * \brief A generic base class for the CRS and CCS helpers.
     *
     * S, GO, and GS are the desired types.  They are also the types
     * of the respective input parameters.  Matrix is expected to be
     * an Amesos2 MatrixAdapter type.
     *
     * The \c Op template parameter is a function-like class that
     * provides a static \c apply() function.
     *
     * \ingroup amesos2_util
     */

    template<class M, typename KV_S, typename KV_GO, typename KV_GS, class Op>
    struct same_gs_helper_kokkos_view
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
                         KV_S& nzvals,
                         KV_GO& indices,
                         KV_GS& pointers,
                         typename KV_GS::value_type& nnz,
                         const Teuchos::Ptr<
                           const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
                                             typename M::node_t> > map,
                         EDistribution distribution,
                         EStorage_Ordering ordering)
      {
        Op::template apply_kokkos_view<KV_S, KV_GO, KV_GS>(mat, nzvals,
          indices, pointers, nnz, map, distribution, ordering);
      }
    };

    template<class M, typename KV_S, typename KV_GO, typename KV_GS, class Op>
    struct diff_gs_helper_kokkos_view
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
                         KV_S& nzvals,
                         KV_GO& indices,
                         KV_GS& pointers,
                         typename KV_GS::value_type& nnz,
                         const Teuchos::Ptr<
                           const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
                                             typename M::node_t> > map,
                         EDistribution distribution,
                         EStorage_Ordering ordering)
      {
        typedef typename M::global_size_t mat_gs_t;
        typedef typename Kokkos::View<mat_gs_t*, Kokkos::HostSpace> KV_TMP;
        size_t i, size = (pointers.extent(0) > 0 ? pointers.extent(0) : 1); // making sure it is at least 1, even for empty local matrix
        KV_TMP pointers_tmp(Kokkos::ViewAllocateWithoutInitializing("pointers_tmp"), size);

        mat_gs_t nnz_tmp = 0;
        Op::template apply_kokkos_view<KV_S, KV_GO, KV_TMP>(mat, nzvals,
          indices, pointers_tmp, nnz_tmp, Teuchos::ptrInArg(*map), distribution, ordering);
        nnz = Teuchos::as<typename KV_GS::value_type>(nnz_tmp);

        typedef typename KV_GS::value_type view_gs_t;
        if (pointers.extent(0) == 1) {
          Kokkos::deep_copy(pointers, 0);
        } else {
          auto host_pointers = Kokkos::create_mirror_view(pointers);
          for (i = 0; i < pointers.extent(0); ++i){
            host_pointers(i) = Teuchos::as<view_gs_t>(pointers_tmp(i));
          }
          Kokkos::deep_copy(pointers, host_pointers);
        }
        nnz = Teuchos::as<view_gs_t>(nnz_tmp);
      }
    };

    template<class M, typename KV_S, typename KV_GO, typename KV_GS, class Op>
    struct same_go_helper_kokkos_view
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
                         KV_S& nzvals,
                         KV_GO& indices,
                         KV_GS& pointers,
                         typename KV_GS::value_type& nnz,
                         const Teuchos::Ptr<
                           const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
                                             typename M::node_t> > map,
                         EDistribution distribution,
                         EStorage_Ordering ordering)
      {
        typedef typename M::global_size_t mat_gs_t;
        typedef typename KV_GS::value_type view_gs_t;
        std::conditional_t<std::is_same_v<view_gs_t,mat_gs_t>,
          same_gs_helper_kokkos_view<M, KV_S, KV_GO, KV_GS, Op>,
          diff_gs_helper_kokkos_view<M, KV_S, KV_GO, KV_GS, Op> >::do_get(mat, nzvals, indices,
                                                                                pointers, nnz, map,
                                                                                distribution, ordering);
      }
    };

    template<class M, typename KV_S, typename KV_GO, typename KV_GS, class Op>
    struct diff_go_helper_kokkos_view
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
                         KV_S& nzvals,
                         KV_GO& indices,
                         KV_GS& pointers,
                         typename KV_GS::value_type& nnz,
                         const Teuchos::Ptr<
                           const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
                                             typename M::node_t> > map,
                         EDistribution distribution,
                         EStorage_Ordering ordering)
      {
        typedef typename M::global_ordinal_t mat_go_t;
        typedef typename M::global_size_t mat_gs_t;
        typedef typename Kokkos::View<mat_go_t*, Kokkos::HostSpace> KV_TMP;
        size_t i, size = indices.extent(0);
        KV_TMP indices_tmp(Kokkos::ViewAllocateWithoutInitializing("indices_tmp"), size);

        typedef typename KV_GO::value_type view_go_t;
        typedef typename KV_GS::value_type view_gs_t;
        std::conditional_t<std::is_same_v<view_gs_t,mat_gs_t>,
          same_gs_helper_kokkos_view<M, KV_S, KV_TMP, KV_GS, Op>,
          diff_gs_helper_kokkos_view<M, KV_S, KV_TMP, KV_GS, Op> >::do_get(mat, nzvals, indices_tmp,
                                                                                 pointers, nnz, map,
                                                                                 distribution, ordering);
        auto host_indices = Kokkos::create_mirror_view(indices);
        for (i = 0; i < size; ++i){
          host_indices(i) = Teuchos::as<view_go_t>(indices_tmp(i));
        }
        Kokkos::deep_copy(indices, host_indices);
      }
    };

    template<class M, typename KV_S, typename KV_GO, typename KV_GS, class Op>
    struct same_scalar_helper_kokkos_view
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
                         KV_S& nzvals,
                         KV_GO& indices,
                         KV_GS& pointers,
                         typename KV_GS::value_type& nnz,
                         const Teuchos::Ptr<
                           const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
                                             typename M::node_t> > map,
                         EDistribution distribution,
                         EStorage_Ordering ordering)
      {
        typedef typename M::global_ordinal_t mat_go_t;
        typedef typename KV_GO::value_type view_go_t;
        std::conditional_t<std::is_same_v<view_go_t, mat_go_t>,
          same_go_helper_kokkos_view<M, KV_S, KV_GO, KV_GS, Op>,
          diff_go_helper_kokkos_view<M, KV_S, KV_GO, KV_GS, Op> >::do_get(mat, nzvals, indices,
                                                                                pointers, nnz, map,
                                                                                distribution, ordering);
      }
    };

    template<class M, typename KV_S, typename KV_GO, typename KV_GS, class Op>
    struct diff_scalar_helper_kokkos_view
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
                         KV_S& nzvals,
                         KV_GO& indices,
                         KV_GS& pointers,
                         typename KV_GS::value_type& nnz,
                         const Teuchos::Ptr<
                           const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
                                             typename M::node_t> > map,
                         EDistribution distribution,
                         EStorage_Ordering ordering)
      {
        typedef typename M::global_ordinal_t mat_go_t;
#if KOKKOS_VERSION >= 40799
        typedef typename KokkosKernels::ArithTraits<typename M::scalar_t>::val_type mat_scalar_t;
#else
        typedef typename Kokkos::ArithTraits<typename M::scalar_t>::val_type mat_scalar_t;
#endif
        typedef typename Kokkos::View<mat_scalar_t*, Kokkos::HostSpace> KV_TMP;
        size_t i, size = nzvals.extent(0);
        KV_TMP nzvals_tmp(Kokkos::ViewAllocateWithoutInitializing("nzvals_tmp"), size);

        typedef typename KV_S::value_type view_scalar_t;
        typedef typename KV_GO::value_type view_go_t;
        std::conditional_t<std::is_same_v<view_go_t, mat_go_t>,
          same_go_helper_kokkos_view<M, KV_TMP, KV_GO, KV_GS, Op>,
          diff_go_helper_kokkos_view<M, KV_TMP, KV_GO, KV_GS, Op> >::do_get(mat, nzvals_tmp, indices,
                                                                                  pointers, nnz, map,
                                                                                  distribution, ordering);

        auto host_nzvals = Kokkos::create_mirror_view(nzvals);
        for (i = 0; i < size; ++i){
          host_nzvals(i) = Teuchos::as<view_scalar_t>(nzvals_tmp(i));
        }
        Kokkos::deep_copy(nzvals, host_nzvals);
      }
    };


    template<class Matrix, typename KV_S, typename KV_GO, typename KV_GS, class Op>
    struct get_cxs_helper_kokkos_view
    {
      static void do_get(const Teuchos::Ptr<const Matrix> mat,
                         KV_S& nzvals,
                         KV_GO& indices,
                         KV_GS& pointers,
                         typename KV_GS::value_type& nnz,
                         EDistribution distribution,
                         EStorage_Ordering ordering=ARBITRARY,
                         typename KV_GO::value_type indexBase = 0)
      {
        typedef typename Matrix::local_ordinal_t lo_t;
        typedef typename Matrix::global_ordinal_t go_t;
        typedef typename Matrix::global_size_t gs_t;
        typedef typename Matrix::node_t node_t;

        const Teuchos::RCP<const Tpetra::Map<lo_t,go_t,node_t> > map
          = getDistributionMap<lo_t,go_t,gs_t,node_t>(distribution,
                                                      Op::get_dimension(mat),
                                                      mat->getComm(),
                                                      indexBase,
                                                      Op::getMapFromMatrix(mat) //getMap must be the map returned, NOT rowmap or colmap
                                                      );
        do_get(mat, nzvals, indices, pointers, nnz, Teuchos::ptrInArg(*map), distribution, ordering);
      }

      /**
       * Basic function overload that uses the matrix's row/col map as
       * returned by Op::getMap().
       */
      static void do_get(const Teuchos::Ptr<const Matrix> mat,
                         KV_S& nzvals,
                         KV_GO& indices,
                         KV_GS& pointers,
                         typename KV_GS::value_type& nnz,
                         EDistribution distribution, // Does this one need a distribution argument??
                         EStorage_Ordering ordering=ARBITRARY)
      {
        const Teuchos::RCP<const Tpetra::Map<typename Matrix::local_ordinal_t,
                                             typename Matrix::global_ordinal_t,
                                             typename Matrix::node_t> > map
          = Op::getMap(mat);
        do_get(mat, nzvals, indices, pointers, nnz, Teuchos::ptrInArg(*map), distribution, ordering);
      }

      /**
       * Function overload that takes an explicit map to use for the
       * representation's distribution.
       */
      static void do_get(const Teuchos::Ptr<const Matrix> mat,
                         KV_S& nzvals,
                         KV_GO& indices,
                         KV_GS& pointers,
                         typename KV_GS::value_type& nnz,
                         const Teuchos::Ptr<
                           const Tpetra::Map<typename Matrix::local_ordinal_t,
                                             typename Matrix::global_ordinal_t,
                                             typename Matrix::node_t> > map,
                         EDistribution distribution,
                         EStorage_Ordering ordering=ARBITRARY)
      {
        typedef typename Matrix::scalar_t mat_scalar;
        typedef typename KV_S::value_type view_scalar_t;

        std::conditional_t<std::is_same_v<mat_scalar,view_scalar_t>,
          same_scalar_helper_kokkos_view<Matrix,KV_S,KV_GO,KV_GS,Op>,
          diff_scalar_helper_kokkos_view<Matrix,KV_S,KV_GO,KV_GS,Op> >::do_get(mat,
                                                                                     nzvals, indices,
                                                                                     pointers, nnz,
                                                                                     map,
                                                                                     distribution, ordering);
      }
    };

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    /*
     * These two function-like classes are meant to be used as the \c
     * Op template parameter for the \c get_cxs_helper template class.
     */
    template<class Matrix>
    struct get_ccs_func
    {      
      template<typename KV_S, typename KV_GO, typename KV_GS>
      static void apply_kokkos_view(const Teuchos::Ptr<const Matrix> mat,
                        KV_S& nzvals,
                        KV_GO& rowind,
                        KV_GS& colptr,
                        typename Matrix::global_size_t& nnz,
                        const Teuchos::Ptr<
                          const Tpetra::Map<typename Matrix::local_ordinal_t,
                                            typename Matrix::global_ordinal_t,
                                            typename Matrix::node_t> > map,
                        EDistribution distribution,
                        EStorage_Ordering ordering)
      {
        mat->getCcs_kokkos_view(nzvals, rowind, colptr, nnz, map, ordering, distribution);
      }

      static
      const Teuchos::RCP<const Tpetra::Map<typename Matrix::local_ordinal_t,
                                           typename Matrix::global_ordinal_t,
                                           typename Matrix::node_t> >
      getMapFromMatrix(const Teuchos::Ptr<const Matrix> mat)
      {
        return mat->getMap(); // returns Teuchos::null if mat is Epetra_CrsMatrix
      }

      static
      const Teuchos::RCP<const Tpetra::Map<typename Matrix::local_ordinal_t,
                                           typename Matrix::global_ordinal_t,
                                           typename Matrix::node_t> >
      getMap(const Teuchos::Ptr<const Matrix> mat)
      {
        return mat->getColMap();
      }

      static
      typename Matrix::global_size_t
      get_dimension(const Teuchos::Ptr<const Matrix> mat)
      {
        return mat->getGlobalNumCols();
      }
    };

    template<class Matrix>
    struct get_crs_func
    {
     template<typename KV_S, typename KV_GO, typename KV_GS>
      static void apply_kokkos_view(const Teuchos::Ptr<const Matrix> mat,
                        KV_S& nzvals,
                        KV_GO& colind,
                        KV_GS& rowptr,
                        typename Matrix::global_size_t& nnz,
                        const Teuchos::Ptr<
                          const Tpetra::Map<typename Matrix::local_ordinal_t,
                                            typename Matrix::global_ordinal_t,
                                            typename Matrix::node_t> > map,
                        EDistribution distribution,
                        EStorage_Ordering ordering)
      {
        mat->getCrs_kokkos_view(nzvals, colind, rowptr, nnz, map, ordering, distribution);
      }

      static
      const Teuchos::RCP<const Tpetra::Map<typename Matrix::local_ordinal_t,
                                           typename Matrix::global_ordinal_t,
                                           typename Matrix::node_t> >
      getMapFromMatrix(const Teuchos::Ptr<const Matrix> mat)
      {
        return mat->getMap(); // returns Teuchos::null if mat is Epetra_CrsMatrix
      }

      static
      const Teuchos::RCP<const Tpetra::Map<typename Matrix::local_ordinal_t,
                                           typename Matrix::global_ordinal_t,
                                           typename Matrix::node_t> >
      getMap(const Teuchos::Ptr<const Matrix> mat)
      {
        return mat->getRowMap();
      }

      static
      typename Matrix::global_size_t
      get_dimension(const Teuchos::Ptr<const Matrix> mat)
      {
        return mat->getGlobalNumRows();
      }
    };
#endif  // DOXYGEN_SHOULD_SKIP_THIS

    /**
     * \brief A generic helper class for getting a CCS representation
     * of a Matrix.
     *
     * The template types \c S , \c GO , and \c GS (scalar, global
     * ordinal, and global size type, respectively) are the types that
     * you would like to get from the Matrix, regardless of what types
     * are actually housed in the matrix.  Type conversions will be
     * performed when necessary.
     *
     * \subsection get_ccs_helper_example Example:
     *
     * Say for example that you have a matrix that has \c
     * complex<double> scalar values, \c int global ordinals, and
     * unsigned long size type, but your solver has a special complex
     * data type that it defines and has size type \c int.  As long as
     * the Teuchos::ValueTypeConversionTraits class is specialized for
     * conversion between the \c complex<double> and the solver's
     * complex double type, then you can use this helper to help with
     * this conversion.  We assume that we want the global matrix
     * representation at the root processor (\c Rooted), and the row
     * indices can be in an arbitrary order (\c Arbitrary):
     *
     * \code
     * // with unsigned long size_t
     * typedef Tpetra::CrsMatrix<std::complex<double>, int, int> mat_t;
     * mat_t my_mat;
     * // initialize mt_mat
     * Array<solver_complex> nzvals(nnz);
     * Array<int> rowind(nnz);
     * Array<int> colptr(numcols+1);
     * get_ccs_helper<mat_t,solver_complex,int,int>::do_get(mat,nzvals,rowind,rowptr,nnz,Rooted,Arbitrary);
     * \endcode
     *
     * \sa \ref get_crs_helper
     * \ingroup amesos2_util
     */
    template<class Matrix, typename KV_S, typename KV_GO, typename KV_GS>
    struct get_ccs_helper_kokkos_view : get_cxs_helper_kokkos_view<Matrix,KV_S,KV_GO,KV_GS,get_ccs_func<Matrix> >
    {};

    /**
     * \brief Similar to get_ccs_helper , but used to get a CRS
     * representation of the given matrix.
     *
     * \sa \ref get_ccs_helper
     * \ingroup amesos2_util
     */
    template<class Matrix, typename KV_S, typename KV_GO, typename KV_GS>
    struct get_crs_helper_kokkos_view : get_cxs_helper_kokkos_view<Matrix,KV_S,KV_GO,KV_GS,get_crs_func<Matrix> >
    {};
    /* End Matrix/MultiVector Utilities */


    ////////////////////////////////////////
    //           Definitions              //
    ////////////////////////////////////////


    template <typename LO, typename GO, typename GS, typename Node>
    const Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >
    getGatherMap( const Teuchos::RCP< const Tpetra::Map<LO,GO,Node> > &map )
    {
      //RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream( Teuchos::null ); // may need to pass an osstream to computeGatherMap for debugging cases...
      Teuchos::RCP< const Tpetra::Map<LO,GO,Node> > gather_map = Tpetra::Details::computeGatherMap(map, Teuchos::null);
      return gather_map;
    }


    template <typename LO, typename GO, typename GS, typename Node>
    const Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >
    getDistributionMap(EDistribution distribution,
                       GS num_global_elements,
                       const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                       GO indexBase,
                       const Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >& map)
    {
        // TODO: Need to add indexBase to cases other than ROOTED
        //  We do not support these maps in any solver now.
      switch( distribution ){
      case DISTRIBUTED:
      case DISTRIBUTED_NO_OVERLAP:
        return Tpetra::createUniformContigMapWithNode<LO,GO, Node>(num_global_elements, comm);
      case GLOBALLY_REPLICATED:
        return Tpetra::createLocalMapWithNode<LO,GO, Node>(num_global_elements, comm);
      case ROOTED:
        {
          int rank = Teuchos::rank(*comm);
          size_t my_num_elems = Teuchos::OrdinalTraits<size_t>::zero();
          if( rank == 0 ) { my_num_elems = num_global_elements; }

          return rcp(new Tpetra::Map<LO,GO, Node>(num_global_elements,
                                                  my_num_elems, indexBase, comm));
        }
      case CONTIGUOUS_AND_ROOTED:
        {
          const Teuchos::RCP<const Tpetra::Map<LO,GO,Node> > gathermap
          = getGatherMap<LO,GO,GS,Node>( map );  //getMap must be the map returned, NOT rowmap or colmap
          return gathermap;
        }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true,
                            std::logic_error,
                            "Control should never reach this point.  "
                            "Please contact the Amesos2 developers." );
      }
    }

    template <typename Scalar,
              typename GlobalOrdinal,
              typename GlobalSizeT>
    void transpose(Teuchos::ArrayView<Scalar> vals,
                   Teuchos::ArrayView<GlobalOrdinal> indices,
                   Teuchos::ArrayView<GlobalSizeT> ptr,
                   Teuchos::ArrayView<Scalar> trans_vals,
                   Teuchos::ArrayView<GlobalOrdinal> trans_indices,
                   Teuchos::ArrayView<GlobalSizeT> trans_ptr)
    {
      /* We have a compressed-row storage format of this matrix.  We
       * transform this into a compressed-column format using a
       * distribution-counting sort algorithm, which is described by
       * D.E. Knuth in TAOCP Vol 3, 2nd ed pg 78.
       */

#ifdef HAVE_AMESOS2_DEBUG
      typename Teuchos::ArrayView<GlobalOrdinal>::iterator ind_it, ind_begin, ind_end;
      ind_begin = indices.begin();
      ind_end = indices.end();
      size_t min_trans_ptr_size = *std::max_element(ind_begin, ind_end) + 1;
      TEUCHOS_TEST_FOR_EXCEPTION( Teuchos::as<size_t>(trans_ptr.size()) < min_trans_ptr_size,
                          std::invalid_argument,
                          "Transpose pointer size not large enough." );
      TEUCHOS_TEST_FOR_EXCEPTION( trans_vals.size() < vals.size(),
                          std::invalid_argument,
                          "Transpose values array not large enough." );
      TEUCHOS_TEST_FOR_EXCEPTION( trans_indices.size() < indices.size(),
                          std::invalid_argument,
                          "Transpose indices array not large enough." );
#else
      typename Teuchos::ArrayView<GlobalOrdinal>::iterator ind_it, ind_end;
#endif
      // Count the number of entries in each column
      Teuchos::Array<GlobalSizeT> count(trans_ptr.size(), 0);
      ind_end = indices.end();
      for( ind_it = indices.begin(); ind_it != ind_end; ++ind_it ){
        ++(count[(*ind_it) + 1]);
      }
      // Accumulate
      typename Teuchos::Array<GlobalSizeT>::iterator cnt_it, cnt_end;
      cnt_end = count.end();
      for( cnt_it = count.begin() + 1; cnt_it != cnt_end; ++cnt_it ){
        *cnt_it = *cnt_it + *(cnt_it - 1);
      }
      // This becomes the array of column pointers
      trans_ptr.assign(count);

      /* Move the nonzero values into their final place in nzval, based on the
       * counts found previously.
       *
       * This sequence deviates from Knuth's algorithm a bit, following more
       * closely the description presented in Gustavson, Fred G. "Two Fast
       * Algorithms for Sparse Matrices: Multiplication and Permuted
       * Transposition" ACM Trans. Math. Softw. volume 4, number 3, 1978, pages
       * 250--269, http://doi.acm.org/10.1145/355791.355796.
       *
       * The output indices end up in sorted order
       */

      GlobalSizeT size = ptr.size();
      for( GlobalSizeT i = 0; i < size - 1; ++i ){
        GlobalOrdinal u = ptr[i];
        GlobalOrdinal v = ptr[i + 1];
        for( GlobalOrdinal j = u; j < v; ++j ){
          GlobalOrdinal k = count[indices[j]];
          trans_vals[k] = vals[j];
          trans_indices[k] = i;
          ++(count[indices[j]]);
        }
      }
    }


    template <typename Scalar1, typename Scalar2>
    void
    scale(Teuchos::ArrayView<Scalar1> vals, size_t l,
          size_t ld, Teuchos::ArrayView<Scalar2> s)
    {
      size_t vals_size = vals.size();
#ifdef HAVE_AMESOS2_DEBUG
      size_t s_size = s.size();
      TEUCHOS_TEST_FOR_EXCEPTION( s_size < l,
                          std::invalid_argument,
                          "Scale vector must have length at least that of the vector" );
#endif
      size_t i, s_i;
      for( i = 0, s_i = 0; i < vals_size; ++i, ++s_i ){
        if( s_i == l ){
          // bring i to the next multiple of ld
          i += ld - s_i;
          s_i = 0;
        }
        vals[i] *= s[s_i];
      }
    }

    template <typename Scalar1, typename Scalar2, class BinaryOp>
    void
    scale(Teuchos::ArrayView<Scalar1> vals, size_t l,
          size_t ld, Teuchos::ArrayView<Scalar2> s,
          BinaryOp binary_op)
    {
      size_t vals_size = vals.size();
#ifdef HAVE_AMESOS2_DEBUG
      size_t s_size = s.size();
      TEUCHOS_TEST_FOR_EXCEPTION( s_size < l,
                          std::invalid_argument,
                          "Scale vector must have length at least that of the vector" );
#endif
      size_t i, s_i;
      for( i = 0, s_i = 0; i < vals_size; ++i, ++s_i ){
        if( s_i == l ){
          // bring i to the next multiple of ld
          i += ld - s_i;
          s_i = 0;
        }
        vals[i] = binary_op(vals[i], s[s_i]);
      }
    }

    template<class row_ptr_view_t, class cols_view_t, class per_view_t>
    void
    reorder(row_ptr_view_t & row_ptr, cols_view_t & cols,
            per_view_t & perm, per_view_t & peri, size_t & nnz,
            bool permute_matrix)
    {
      #ifndef HAVE_AMESOS2_METIS
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "Cannot reorder for cuSolver because no METIS is available.");
      #else
        typedef typename cols_view_t::value_type ordinal_type;
        typedef typename row_ptr_view_t::value_type size_type;

        // begin on host where we'll run metis reorder
        auto host_row_ptr = Kokkos::create_mirror_view(row_ptr);
        auto host_cols = Kokkos::create_mirror_view(cols);
        Kokkos::deep_copy(host_row_ptr, row_ptr);
        Kokkos::deep_copy(host_cols, cols);

        // strip out the diagonals - metis will just crash with them included.
        // make space for the stripped version
        typedef Kokkos::View<idx_t*, Kokkos::HostSpace> host_metis_array;
        const ordinal_type size = row_ptr.size() - 1;
        size_type max_nnz = host_row_ptr(size);
        host_metis_array host_strip_diag_row_ptr(
          Kokkos::ViewAllocateWithoutInitializing("host_strip_diag_row_ptr"),
          size+1);
        host_metis_array host_strip_diag_cols(
          Kokkos::ViewAllocateWithoutInitializing("host_strip_diag_cols"),
          max_nnz);

        size_type new_nnz = 0;
        for(ordinal_type i = 0; i < size; ++i) {
          host_strip_diag_row_ptr(i) = new_nnz;
          for(size_type j = host_row_ptr(i); j < host_row_ptr(i+1); ++j) {
            if (i != host_cols(j)) {
              host_strip_diag_cols(new_nnz++) = host_cols(j);
            }
          }
        }
        host_strip_diag_row_ptr(size) = new_nnz;

        // we'll get original permutations on host
        host_metis_array host_perm(
          Kokkos::ViewAllocateWithoutInitializing("host_perm"), size);
        host_metis_array host_peri(
          Kokkos::ViewAllocateWithoutInitializing("host_peri"), size);

        // If we want to remove metis.h included in this header we can move this
        // to the cpp, but we need to decide how to handle the idx_t declaration.
        idx_t metis_size = size;
        int err = METIS_NodeND(&metis_size, host_strip_diag_row_ptr.data(), host_strip_diag_cols.data(),
          NULL, NULL, host_perm.data(), host_peri.data());

        TEUCHOS_TEST_FOR_EXCEPTION(err != METIS_OK, std::runtime_error,
          "METIS_NodeND failed to sort matrix.");

        // put the permutations on our saved device ptrs
        // these will be used to permute x and b when we solve
        typedef typename cols_view_t::execution_space exec_space_t;
        auto device_perm = Kokkos::create_mirror_view(exec_space_t(), host_perm);
        auto device_peri = Kokkos::create_mirror_view(exec_space_t(), host_peri);
        deep_copy(device_perm, host_perm);
        deep_copy(device_peri, host_peri);

        // also set the permutation which may need to convert the type from
        // metis to the native ordinal_type
        deep_copy_or_assign_view(perm, device_perm);
        deep_copy_or_assign_view(peri, device_peri);

        if (permute_matrix) {
          // we'll permute matrix on device to a new set of arrays
          row_ptr_view_t new_row_ptr(
            Kokkos::ViewAllocateWithoutInitializing("new_row_ptr"), row_ptr.size());
          cols_view_t new_cols(
            Kokkos::ViewAllocateWithoutInitializing("new_cols"), cols.size() - new_nnz/2);

          // permute row indices
          Kokkos::RangePolicy<exec_space_t> policy_row(0, row_ptr.size());
          Kokkos::parallel_scan(policy_row, KOKKOS_LAMBDA(
          ordinal_type i, size_type & update, const bool &final) {
          if(final) {
            new_row_ptr(i) = update;
          }
          if(i < size) {
              ordinal_type count = 0;
              const ordinal_type row = device_perm(i);
              for(ordinal_type k = row_ptr(row); k < row_ptr(row + 1); ++k) {
                const ordinal_type j = device_peri(cols(k)); /// col in A
                count += (i >= j); /// lower triangular
              }
              update += count;
            }
          });

          // permute col indices
          Kokkos::RangePolicy<exec_space_t> policy_col(0, size);
          Kokkos::parallel_for(policy_col, KOKKOS_LAMBDA(ordinal_type i) {
            const ordinal_type kbeg = new_row_ptr(i);
            const ordinal_type row = device_perm(i);
            const ordinal_type col_beg = row_ptr(row);
            const ordinal_type col_end = row_ptr(row + 1);
            const ordinal_type nk = col_end - col_beg;
            for(ordinal_type k = 0, t = 0; k < nk; ++k) {
              const ordinal_type tk = kbeg + t;
              const ordinal_type sk = col_beg + k;
              const ordinal_type j = device_peri(cols(sk));
              if(i >= j) {
                new_cols(tk) = j;
                ++t;
              }
            }
          });

          // finally set the inputs to the new sorted arrays
          row_ptr = new_row_ptr;
          cols = new_cols;
        }

        nnz = new_nnz;
      #endif // HAVE_AMESOS2_METIS
    }

    template<class values_view_t, class row_ptr_view_t,
      class cols_view_t, class per_view_t>
    void
    reorder_values(values_view_t & values, const row_ptr_view_t & orig_row_ptr,
      const row_ptr_view_t & new_row_ptr,
      const cols_view_t & orig_cols, const per_view_t & perm, const per_view_t & peri,
      size_t nnz)
    {
        typedef typename cols_view_t::value_type ordinal_type;
        typedef typename cols_view_t::execution_space exec_space_t;

        auto device_perm = Kokkos::create_mirror_view(exec_space_t(), perm);
        auto device_peri = Kokkos::create_mirror_view(exec_space_t(), peri);
        deep_copy(device_perm, perm);
        deep_copy(device_peri, peri);

        const ordinal_type size = orig_row_ptr.size() - 1;

        auto host_orig_row_ptr = Kokkos::create_mirror_view(orig_row_ptr);
        auto new_nnz = host_orig_row_ptr(size); // TODO: Maybe optimize this by caching

        values_view_t new_values(
          Kokkos::ViewAllocateWithoutInitializing("new_values"), values.size() - new_nnz/2);

        // permute col indices
        Kokkos::RangePolicy<exec_space_t> policy_col(0, size);
        Kokkos::parallel_for(policy_col, KOKKOS_LAMBDA(ordinal_type i) {
          const ordinal_type kbeg = new_row_ptr(i);
          const ordinal_type row = device_perm(i);
          const ordinal_type col_beg = orig_row_ptr(row);
          const ordinal_type col_end = orig_row_ptr(row + 1);
          const ordinal_type nk = col_end - col_beg;
          for(ordinal_type k = 0, t = 0; k < nk; ++k) {
            const ordinal_type tk = kbeg + t;
            const ordinal_type sk = col_beg + k;
            const ordinal_type j = device_peri(orig_cols(sk));
            if(i >= j) {
              new_values(tk) = values(sk);
              ++t;
            }
          }
        });

        values = new_values;
    }

    template<class array_view_t, class per_view_t>
    void
    apply_reorder_permutation(const array_view_t & array,
      array_view_t & permuted_array, const per_view_t & permutation) {
        if(permuted_array.extent(0) != array.extent(0) || permuted_array.extent(1) != array.extent(1)) {
          permuted_array = array_view_t(
            Kokkos::ViewAllocateWithoutInitializing("permuted_array"),
            array.extent(0), array.extent(1));
        }
        typedef typename array_view_t::execution_space exec_space_t;
        Kokkos::RangePolicy<exec_space_t> policy(0, array.extent(0));
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(size_t i) {
          for(size_t j = 0; j < array.extent(1); ++j) {
            permuted_array(i, j) = array(permutation(i), j);
          }
        });
    }


    // used to read matrix with gapped GIDs for test/example
    template <typename GO, typename Scalar>
    bool
    readEntryFromFile (GO& gblRowInd, GO& gblColInd, Scalar& val, const std::string& s)
    {
      if (s.size () == 0 || s.find ("%") != std::string::npos) {
        return false; // empty line or comment line
      }
      std::istringstream in (s);

      if (! in) {
        return false;
      }
      in >> gblRowInd;
      if (! in) {
        return false;
      }
      in >> gblColInd;
      if (! in) {
        return false;
      }
      in >> val;
      return true;
    }

    template<class map_type, class MAT>
    Teuchos::RCP<MAT>
    readCrsMatrixFromFile (const std::string& matrixFilename,
                           Teuchos::RCP<Teuchos::FancyOStream> & fos,
                           const Teuchos::RCP<const map_type>& rowMap,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const bool convert_to_zero_base,
                           const int header_size)
    {
      using Scalar = typename MAT::scalar_type;
      using GO = typename MAT::global_ordinal_type;

      using counter_type = std::map<GO, size_t>;
      using pair_type = std::pair<const GO, size_t>;
      using Teuchos::RCP;

      auto comm = rowMap->getComm ();
      const int myRank = comm->getRank ();

      std::ifstream inFile;
      int opened = 0;
      if (myRank == 0)
      {
        try {
          inFile.open (matrixFilename);
          if (inFile) {
            opened = 1;
          }
        }
        catch (...) {
          opened = 0;
        }
      }
      Teuchos::broadcast<int, int> (*comm, 0, Teuchos::outArg (opened));
      TEUCHOS_TEST_FOR_EXCEPTION
        (opened == 0, std::runtime_error, "readCrsMatrixFromFile: "
         "Failed to open file \"" << matrixFilename << "\" on Process 0.");

      RCP<MAT> A;
      if (myRank == 0)
      {
        std::string line;

        // Skip the first N lines.  This is a hack, specific to the input file in question.
        //*fos << "  Reading matrix market file. Skip " << header_size << "  header lines" << std::endl;
        for ( int i = 0; i < header_size; ++i ) {
          std::getline (inFile, line);
        }

        counter_type counts;
        Teuchos::Array<Scalar> vals;
        Teuchos::Array<GO> gblRowInds;
        Teuchos::Array<GO> gblColInds;
        while (inFile) {
          std::getline (inFile, line);
          GO gblRowInd {};
          GO gblColInd {};
          Scalar val {};
          const bool gotLine = readEntryFromFile (gblRowInd, gblColInd, val, line);
          if (gotLine) {
            //*fos << "  read mtx rank: " << myRank << " |  gblRowInd = " << gblRowInd << "  gblColInd = " << gblColInd << std::endl;
            if ( convert_to_zero_base ) {
              gblRowInd -= 1 ;
              gblColInd -= 1 ;
            }
            counts[gblRowInd]++;
            vals.push_back(val);
            gblRowInds.push_back(gblRowInd);
            gblColInds.push_back(gblColInd);
          }
        }

        // Max number of entries in any row
        auto pr = std::max_element(
          std::begin(counts),
          std::end(counts),
          [] (pair_type const& p1, pair_type const& p2){ return p1.second < p2.second; }
        );
        size_t maxCount = (counts.empty()) ? size_t(0) : pr->second;
        A = Teuchos::rcp(new MAT(rowMap, maxCount));
        for (typename Teuchos::Array<GO>::size_type i=0; i<gblRowInds.size(); i++) {
          A->insertGlobalValues (gblRowInds[i], gblColInds(i,1), vals(i,1));
        }
      } else {
        A = Teuchos::rcp(new MAT(rowMap, 0));
      }

      A->fillComplete (domainMap, rangeMap);
      return A;
    }
    /** @} */

  } // end namespace Util

} // end namespace Amesos2

#endif  // #ifndef AMESOS2_UTIL_HPP
