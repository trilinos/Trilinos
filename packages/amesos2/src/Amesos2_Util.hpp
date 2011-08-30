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
  \file   Amesos2_Util.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu May 27 13:11:13 CDT 2010

  \brief  Utility functions for Amesos2
*/

#ifndef AMESOS2_UTIL_HPP
#define AMESOS2_UTIL_HPP

#include "Amesos2_config.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_FancyOStream.hpp>

#ifdef HAVE_AMESOS2_EPETRA
#  include <Epetra_Map.h>
#endif

#include <Tpetra_Map.hpp>

#include "Amesos2_TypeDecl.hpp"
#include "Amesos2_Meta.hpp"


namespace Amesos2 {

  namespace Util {

    /**
     * \internal
     * \defgroup amesos2_util Amesos2 Utilities
     * @{
     */

    using Teuchos::RCP;
    using Teuchos::ArrayView;

    using Meta::is_same;
    using Meta::if_then_else;

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
    getDistributionMap(EDistribution distribution,
		       GS num_global_elements,
		       const Teuchos::RCP<const Teuchos::Comm<int> >& comm);
    

#ifdef HAVE_AMESOS2_EPETRA
    /**
     * \brief Transform an Epetra_Map object into a Tpetra::Map
     *
     * \ingroup amesos2_utils
     */
    template <typename LO, typename GO, typename GS, typename Node>
    RCP<Tpetra::Map<LO,GO,Node> >
    epetra_map_to_tpetra_map(const Epetra_BlockMap& map);

    /**
     * \brief Transform a Tpetra::Map object into an Epetra_Map
     *
     * \ingroup amesos2_utils
     */
    template <typename LO, typename GO, typename GS, typename Node>
    RCP<Epetra_Map>
    tpetra_map_to_epetra_map(const Tpetra::Map<LO,GO,Node>& map);

    /**
     * \brief Transform an Epetra_Comm object into a Teuchos::Comm object
     *
     * \ingroup amesos2_utils
     */
    const RCP<const Teuchos::Comm<int> > to_teuchos_comm(RCP<const Epetra_Comm> c);

    /**
     * \brief Transfrom a Teuchos::Comm object into an Epetra_Comm object
     *
     * \ingroup amesos2_utils
     */
    const RCP<const Epetra_Comm> to_epetra_comm(RCP<const Teuchos::Comm<int> > c);
#endif	// HAVE_AMESOS2_EPETRA

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


    //////////////////////////////////
    // Matrix/MultiVector Utilities //
    //////////////////////////////////
    
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    /*
     * The following represents a general way of getting a CRS or CCS
     * representation of a matrix with implicit type conversions.  The
     * \c get_crs_helper and \c get_ccs_helper classes are templated
     * on 4 types:
     *
     * - A Matrix type (conforming to the Amesos2 MatrixAdapter interface)
     * - A scalar type
     * - A global ordinal type
     * - A global size type
     *
     * The last three template types correspond to the input argument
     * types.  For example, if the scalar type is \c double , then we
     * require that the \c nzvals argument is a \c
     * Teuchos::ArrayView<double> type.
     *
     * These helpers perform any type conversions that must be
     * performed to go between the Matrix's types and the input types.
     * If no conversions are necessary the static functions can be
     * effectively inlined.
     */

    template <class M, typename S, typename GO, typename GS, class Op>
    struct same_gs_helper
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
			 const ArrayView<typename M::scalar_t> nzvals,
			 const ArrayView<typename M::global_ordinal_t> indices,
			 const ArrayView<GS> pointers,
			 GS& nnz,
			 const Teuchos::Ptr<
			   const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
			                     typename M::node_t> > map,
			 EStorage_Ordering ordering)
      {
	Op::apply(mat, nzvals, indices, pointers, nnz, map, ordering);
      }
    };

    template <class M, typename S, typename GO, typename GS, class Op>
    struct diff_gs_helper
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
			 const ArrayView<typename M::scalar_t> nzvals,
			 const ArrayView<typename M::global_ordinal_t> indices,
			 const ArrayView<GS> pointers,
			 GS& nnz,
			 const Teuchos::Ptr<
			   const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
			                     typename M::node_t> > map,
			 EStorage_Ordering ordering)
      {
	typedef typename M::global_size_t mat_gs_t;
	typename ArrayView<GS>::size_type i, size = pointers.size();
	Teuchos::Array<mat_gs_t> pointers_tmp(size);
	mat_gs_t nnz_tmp = 0;

	Op::apply(mat, nzvals, indices, pointers_tmp, nnz_tmp, map, ordering);

	for (i = 0; i < size; ++i){
	  pointers[i] = Teuchos::as<GS>(pointers_tmp[i]);
	}
	nnz = Teuchos::as<GS>(nnz_tmp);
      }
    };

    template <class M, typename S, typename GO, typename GS, class Op>
    struct same_go_helper
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
			 const ArrayView<typename M::scalar_t> nzvals,
			 const ArrayView<GO> indices,
			 const ArrayView<GS> pointers,
			 GS& nnz,
			 const Teuchos::Ptr<
			   const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
			                     typename M::node_t> > map,
			 EStorage_Ordering ordering)
      {
	typedef typename M::global_size_t mat_gs_t;
	if_then_else<is_same<GS,mat_gs_t>::value,
	  same_gs_helper<M,S,GO,GS,Op>,
	  diff_gs_helper<M,S,GO,GS,Op> >::type::do_get(mat, nzvals, indices,
						       pointers, nnz, map,
						       ordering);
      }
    };

    template <class M, typename S, typename GO, typename GS, class Op>
    struct diff_go_helper
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
			 const ArrayView<typename M::scalar_t> nzvals,
			 const ArrayView<GO> indices,
			 const ArrayView<GS> pointers,
			 GS& nnz,
			 const Teuchos::Ptr<
			   const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
			                     typename M::node_t> > map,
			 EStorage_Ordering ordering)
      {
	typedef typename M::global_ordinal_t mat_go_t;
	typedef typename M::global_size_t mat_gs_t;
	typename ArrayView<GO>::size_type i, size = indices.size();
	Teuchos::Array<mat_go_t> indices_tmp(size);

	if_then_else<is_same<GS,mat_gs_t>::value,
	  same_gs_helper<M,S,GO,GS,Op>,
	  diff_gs_helper<M,S,GO,GS,Op> >::type::do_get(mat, nzvals, indices_tmp,
						       pointers, nnz, map,
						       ordering);

	for (i = 0; i < size; ++i){
	  indices[i] = Teuchos::as<GO>(indices_tmp[i]);
	}
      }
    };

    template <class M, typename S, typename GO, typename GS, class Op>
    struct same_scalar_helper
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
			 const ArrayView<S> nzvals,
			 const ArrayView<GO> indices,
			 const ArrayView<GS> pointers,
			 GS& nnz,
			 const Teuchos::Ptr<
			   const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
			                     typename M::node_t> > map,
			 EStorage_Ordering ordering)
      {
	typedef typename M::global_ordinal_t mat_go_t;
	if_then_else<is_same<GO,mat_go_t>::value,
	  same_go_helper<M,S,GO,GS,Op>,
	  diff_go_helper<M,S,GO,GS,Op> >::type::do_get(mat, nzvals, indices,
						       pointers, nnz, map,
						       ordering);
      }
    };

    template <class M, typename S, typename GO, typename GS, class Op>
    struct diff_scalar_helper
    {
      static void do_get(const Teuchos::Ptr<const M> mat,
			 const ArrayView<S> nzvals,
			 const ArrayView<GO> indices,
			 const ArrayView<GS> pointers,
			 GS& nnz,
			 const Teuchos::Ptr<
			   const Tpetra::Map<typename M::local_ordinal_t,
                                             typename M::global_ordinal_t,
			                     typename M::node_t> > map,
			 EStorage_Ordering ordering)
      {
	typedef typename M::scalar_t mat_scalar_t;
	typedef typename M::global_ordinal_t mat_go_t;
	typename ArrayView<S>::size_type i, size = nzvals.size();
	Teuchos::Array<mat_scalar_t> nzvals_tmp(size);

	if_then_else<is_same<GO,mat_go_t>::value,
	  same_go_helper<M,S,GO,GS,Op>,
	  diff_go_helper<M,S,GO,GS,Op> >::type::do_get(mat, nzvals_tmp, indices,
						       pointers, nnz, map,
						       ordering);

	for (i = 0; i < size; ++i){
	  nzvals[i] = Teuchos::as<S>(nzvals_tmp[i]);
	}
      }
    };
#endif	// DOXYGEN_SHOULD_SKIP_THIS

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
    template<class Matrix, typename S, typename GO, typename GS, class Op>
    struct get_cxs_helper
    {
      static void do_get(const Teuchos::Ptr<const Matrix> mat,
			 const ArrayView<S> nzvals,
			 const ArrayView<GO> indices,
			 const ArrayView<GS> pointers,
			 GS& nnz, EDistribution distribution,
			 EStorage_Ordering ordering=ARBITRARY)
      {
	typedef typename Matrix::local_ordinal_t lo_t;
	typedef typename Matrix::global_ordinal_t go_t;
	typedef typename Matrix::global_size_t gs_t;
	typedef typename Matrix::node_t node_t;
	
	const Teuchos::RCP<const Tpetra::Map<lo_t,go_t,node_t> > map
	  = getDistributionMap<lo_t,go_t,gs_t,node_t>(distribution,
						      Op::get_dimension(mat),
						      mat->getComm());
	do_get(mat, nzvals, indices, pointers, nnz, Teuchos::ptrInArg(*map), ordering);
      }

      /**
       * Basic function overload that uses the matrix's row/col map as
       * returned by Op::getMap().
       */
      static void do_get(const Teuchos::Ptr<const Matrix> mat,
			 const ArrayView<S> nzvals,
			 const ArrayView<GO> indices,
			 const ArrayView<GS> pointers,
			 GS& nnz, EStorage_Ordering ordering=ARBITRARY)
      {
	const Teuchos::RCP<const Tpetra::Map<typename Matrix::local_ordinal_t,
	                                     typename Matrix::global_ordinal_t,
	                                     typename Matrix::node_t> > map
	  = Op::getMap(mat);
	do_get(mat, nzvals, indices, pointers, nnz, Teuchos::ptrInArg(*map), ordering);
      }

      /**
       * Function overload that takes an explicit map to use for the
       * representation's distribution.
       */
      static void do_get(const Teuchos::Ptr<const Matrix> mat,
			 const ArrayView<S> nzvals,
			 const ArrayView<GO> indices,
			 const ArrayView<GS> pointers,
			 GS& nnz,
			 const Teuchos::Ptr<
			   const Tpetra::Map<typename Matrix::local_ordinal_t,
                                             typename Matrix::global_ordinal_t,
			                     typename Matrix::node_t> > map,
			 EStorage_Ordering ordering=ARBITRARY)
      {
	typedef typename Matrix::scalar_t mat_scalar;
	if_then_else<is_same<mat_scalar,S>::value,
	  same_scalar_helper<Matrix,S,GO,GS,Op>,
	  diff_scalar_helper<Matrix,S,GO,GS,Op> >::type::do_get(mat,
								nzvals, indices,
								pointers, nnz,
								map, ordering);
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
      static void apply(const Teuchos::Ptr<const Matrix> mat,
			const ArrayView<typename Matrix::scalar_t> nzvals,
			const ArrayView<typename Matrix::global_ordinal_t> rowind,
			const ArrayView<typename Matrix::global_size_t> colptr,
			typename Matrix::global_size_t& nnz,
			const Teuchos::Ptr<
			  const Tpetra::Map<typename Matrix::local_ordinal_t,
                                            typename Matrix::global_ordinal_t,
			                    typename Matrix::node_t> > map,
			EStorage_Ordering ordering)
      {
	mat->getCcs(nzvals, rowind, colptr, nnz, map, ordering);
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
      static void apply(const Teuchos::Ptr<const Matrix> mat,
			const ArrayView<typename Matrix::scalar_t> nzvals,
			const ArrayView<typename Matrix::global_ordinal_t> colind,
			const ArrayView<typename Matrix::global_size_t> rowptr,
			typename Matrix::global_size_t& nnz,
			const Teuchos::Ptr<
			  const Tpetra::Map<typename Matrix::local_ordinal_t,
			                    typename Matrix::global_ordinal_t,
			                    typename Matrix::node_t> > map,
			EStorage_Ordering ordering)
      {
	mat->getCrs(nzvals, colind, rowptr, nnz, map, ordering);
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
#endif	// DOXYGEN_SHOULD_SKIP_THIS

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
    template<class Matrix, typename S, typename GO, typename GS>
    struct get_ccs_helper : get_cxs_helper<Matrix,S,GO,GS,get_ccs_func<Matrix> >
    {};

    /**
     * \brief Similar to get_ccs_helper , but used to get a CRS
     * representation of the given matrix.
     *
     * \sa \ref get_ccs_helper 
     * \ingroup amesos2_util
     */
    template<class Matrix, typename S, typename GO, typename GS>
    struct get_crs_helper : get_cxs_helper<Matrix,S,GO,GS,get_crs_func<Matrix> >
    {};

    /* End Matrix/MultiVector Utilities */


    ////////////////////////////////////////
    //           Definitions              //
    ////////////////////////////////////////

    template <typename LO, typename GO, typename GS, typename Node>
    const Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >
    getDistributionMap(EDistribution distribution,
		       GS num_global_elements,
		       const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
    {
      switch( distribution ){
      case DISTRIBUTED:
      case DISTRIBUTED_NO_OVERLAP:
	return Tpetra::createUniformContigMap<LO,GO>(num_global_elements, comm);
	break;
      case GLOBALLY_REPLICATED:
	return Tpetra::createLocalMap<LO,GO>(num_global_elements, comm);
	break;
      case ROOTED:
	{
	  int rank = Teuchos::rank(*comm);
	  size_t my_num_elems = Teuchos::OrdinalTraits<size_t>::zero();
	  if( rank == 0 ) my_num_elems = num_global_elements;
	  return Tpetra::createContigMap<LO,GO>(num_global_elements,
						my_num_elems, comm);
	  break;
	}
      default:
	TEST_FOR_EXCEPTION( true,
			    std::logic_error,
			    "Control should never reach this point.  "
			    "Please contact the Amesos2 developers." );
	break;
      }
    }

#ifdef HAVE_AMESOS2_EPETRA
    template <typename LO, typename GO, typename GS, typename Node>
    Teuchos::RCP<Tpetra::Map<LO,GO,Node> >
    epetra_map_to_tpetra_map(const Epetra_BlockMap& map)
    {
      using Teuchos::as;
      using Teuchos::rcp;

      int num_my_elements = map.NumMyElements();
      Teuchos::Array<int> my_global_elements(num_my_elements);
      map.MyGlobalElements(my_global_elements.getRawPtr());

      typedef Tpetra::Map<LO,GO,Node> map_t;
      RCP<map_t> tmap = rcp(new map_t(Teuchos::OrdinalTraits<GS>::invalid(),
				      my_global_elements(),
				      as<GO>(map.IndexBase()),
				      to_teuchos_comm(Teuchos::rcpFromRef(map.Comm()))));
      return tmap;
    }

    template <typename LO, typename GO, typename GS, typename Node>
    Teuchos::RCP<Epetra_Map>
    tpetra_map_to_epetra_map(const Tpetra::Map<LO,GO,Node>& map)
    {
      using Teuchos::as;
      
      Teuchos::Array<GO> elements_tmp;
      elements_tmp = map.getNodeElementList();
      int num_my_elements = elements_tmp.size();
      Teuchos::Array<int> my_global_elements(num_my_elements);
      for (int i = 0; i < num_my_elements; ++i){
	my_global_elements[i] = as<int>(elements_tmp[i]);
      }

      using Teuchos::rcp;
      RCP<Epetra_Map> emap = rcp(new Epetra_Map(-1,
						num_my_elements,
						my_global_elements.getRawPtr(),
						as<GO>(map.getIndexBase()),
						*to_epetra_comm(map.getComm())));
      return emap;
    }
#endif	// HAVE_AMESOS2_EPETRA

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
      TEST_FOR_EXCEPTION( Teuchos::as<size_t>(trans_ptr.size()) < min_trans_ptr_size,
			  std::invalid_argument,
			  "Transpose pointer size not large enough." );
      TEST_FOR_EXCEPTION( trans_vals.size() < vals.size(),
			  std::invalid_argument,
			  "Transpose values array not large enough." );
      TEST_FOR_EXCEPTION( trans_indices.size() < indices.size(),
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
      TEST_FOR_EXCEPTION( s_size < l,
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
      TEST_FOR_EXCEPTION( s_size < l,
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

    /** @} */

  } // end namespace Util

} // end namespace Amesos2

#endif	// #ifndef AMESOS2_UTIL_HPP
