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


#ifndef AMESOS2_MULTIVECADAPTER_DEF_HPP
#define AMESOS2_MULTIVECADAPTER_DEF_HPP

#include "Amesos2_TpetraMultiVecAdapter_def.hpp"
// EpetraMultiVecAdapter_def.hpp not included because the specialization is not a template

#include "Amesos2_Util.hpp"     // for getDistributionMap

namespace Amesos2{

  namespace Util {

    ///////////////////////////////
    // Pointer-getting utilities //
    ///////////////////////////////

    template <typename MV, typename V>
    typename vector_pointer_helper<MV, V>::ptr_return_type *
    vector_pointer_helper<MV, V>::get_pointer_to_vector ( const Teuchos::Ptr< MV > &mv ) {
      return mv->getMVPointer_impl();
    }

    template <typename MV, typename V>
    typename vector_pointer_helper<MV, V>::ptr_return_type *
    vector_pointer_helper<MV, V>::get_pointer_to_vector ( Teuchos::Ptr< MV > &mv ) {
      return mv->getMVPointer_impl();
    }

    template <typename MV, typename V>
    typename vector_pointer_helper<MV, V>::ptr_return_type *
    vector_pointer_helper<MV, V>::get_pointer_to_vector ( const Teuchos::Ptr< const MV > &mv ) {
      return mv->getMVPointer_impl();
    }

    template <typename MV, typename V>
    typename vector_pointer_helper<MV, V>::ptr_return_type *
    vector_pointer_helper<MV, V>::get_pointer_to_vector ( Teuchos::Ptr< const MV > &mv ) {
      return mv->getMVPointer_impl();
    }


    ////////////////////////////
    // Copy-getting utilities //
    ////////////////////////////

    /*
     * If the multivector scalar type and the desired scalar tpye are
     * the same, then we can do a simple straight copy.
     */
    template <typename MV>
    void same_type_get_copy<MV>::apply(const Teuchos::Ptr<const MV>& mv,
                                       const Teuchos::ArrayView<typename MV::scalar_t>& v,
                                       const size_t ldx,
                                       Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                       EDistribution distribution )
    {
      mv->get1dCopy (v, ldx, distribution_map, distribution);
    }

    template <typename MV, typename KV>
    void same_type_get_copy_kokkos_view<MV, KV>::apply(const Teuchos::Ptr<const MV>& mv,
                                       KV& kokkos_view,
                                       const size_t ldx,
                                       Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                       EDistribution distribution )
    {
      mv->get1dCopy_kokkos_view (kokkos_view, ldx, distribution_map, distribution);
    }

    /*
     * In the case where the scalar type of the multi-vector and the
     * corresponding S type are different, then we need to first get a
     * copy of the scalar values, then convert each one into the S
     * type before inserting into the vals array.
     */
    template <typename MV, typename S>
    void diff_type_get_copy<MV,S>::
    apply (const Teuchos::Ptr<const MV>& mv,
           const Teuchos::ArrayView<S>& v,
           const size_t& ldx,
           Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
           EDistribution distribution )
    {
      typedef typename MV::scalar_t mv_scalar_t;
      typedef typename Teuchos::Array<mv_scalar_t>::size_type size_type;

      TEUCHOS_TEST_FOR_EXCEPTION(
        mv.getRawPtr () == NULL, std::invalid_argument,
        "Amesos2::diff_type_get_copy::apply: mv is null.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        distribution_map.getRawPtr () == NULL, std::invalid_argument,
        "Amesos2::diff_type_get_copy::apply: distribution_map is null.");

      const size_type vals_length = v.size ();
      Teuchos::Array<mv_scalar_t> vals_tmp (vals_length);

      mv->get1dCopy (vals_tmp (), ldx, distribution_map, distribution);
      for (size_type i = 0; i < vals_length; ++i) {
        v[i] = Teuchos::as<S> (vals_tmp[i]);
      }
    }

    template <typename MV, typename KV>
    void diff_type_get_copy_kokkos_view<MV,KV>::
    apply (const Teuchos::Ptr<const MV>& mv,
           KV& kokkos_view,
           const size_t& ldx,
           Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
           EDistribution distribution )
    {
      typedef typename MV::scalar_t mv_scalar_t;

      TEUCHOS_TEST_FOR_EXCEPTION(
        mv.getRawPtr () == NULL, std::invalid_argument,
        "Amesos2::diff_type_get_copy_kokkos_view::apply: mv is null.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        distribution_map.getRawPtr () == NULL, std::invalid_argument,
        "Amesos2::diff_type_get_copy_kokkos_view::apply: distribution_map is null.");

      // Create a View to hold the mismatched type data
      // This should have same layout and memory space as kokkos_view and only
      // differs in the value_type which matches the MV.
      Kokkos::View<mv_scalar_t**, typename KV::array_layout, typename KV::execution_space>
        mv_data(Kokkos::ViewAllocateWithoutInitializing("mv_data"), kokkos_view.extent(0), kokkos_view.extent(1));

      mv->get1dCopy_kokkos_view(mv_data, ldx, distribution_map, distribution);

      // Now copy element by element in KV::execution_space
      // Note this code is not currently tested in the default Tacho Amesos2_Solver_Test
      // setup so is not tested anywhere. Manually tested by forcing diff_type_get_copy_kokkos_view
      // to be always on. MDM-TODO make sure it's tested

      // MDM-TODO inquire if deep copy is ok instead of loop (would bypass Teuchos::as)
      //  Kokkos::deep_copy(kokkos_view, mv_data);
      Kokkos::parallel_for(
        Kokkos::RangePolicy<typename KV::execution_space, int> (0, kokkos_view.size()),
        KOKKOS_LAMBDA (size_t n) {
        kokkos_view.data()[n] = Teuchos::as<typename KV::value_type>(mv_data.data()[n]);
      });
    }

    /** \internal
     *
     * \brief Helper class for getting 1-D copies of multivectors
     *
     * Handles datatype conversion when appropriate.
     */
    template <class MV, typename S>
    void get_1d_copy_helper<MV,S>::
    do_get (const Teuchos::Ptr<const MV>& mv,
            const Teuchos::ArrayView<S>& vals,
            const size_t ldx,
            Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
            EDistribution distribution)
    {
      // Dispatch to the copy function appropriate for the type
      if_then_else<is_same<typename MV::scalar_t,S>::value,
        same_type_get_copy<MV>,
        diff_type_get_copy<MV,S> >::type::apply (mv, vals, ldx, distribution_map, distribution);
    }

    template <class MV, typename S>
    void get_1d_copy_helper<MV,S>::
    do_get (const Teuchos::Ptr<const MV>& mv,
            const Teuchos::ArrayView<S>& vals,
            const size_t ldx,
            EDistribution distribution,
            typename MV::global_ordinal_t indexBase)
    {
      typedef typename MV::local_ordinal_t lo_t;
      typedef typename MV::global_ordinal_t go_t;
      typedef typename MV::global_size_t gs_t;
      typedef typename MV::node_t node_t;

      TEUCHOS_TEST_FOR_EXCEPTION(
        mv.getRawPtr () == NULL, std::invalid_argument,
        "Amesos2::get_1d_copy_helper::do_get(5 args): mv is null.");

      Teuchos::RCP<const Tpetra::Map<lo_t,go_t,node_t> > map
        = Amesos2::Util::getDistributionMap<lo_t,go_t,gs_t,node_t> (distribution,
                                                                    mv->getGlobalLength (),
                                                                    mv->getComm (),
                                                                    indexBase,
                                                                    mv->getMap());

      do_get (mv, vals, ldx, Teuchos::ptrInArg (*map), distribution);
    }

    template <class MV, typename S>
    void get_1d_copy_helper<MV,S>::do_get(const Teuchos::Ptr<const MV>& mv,
                                          const Teuchos::ArrayView<S>& vals,
                                          const size_t ldx)
    {
      typedef Tpetra::Map<typename MV::local_ordinal_t,
                          typename MV::global_ordinal_t,
                          typename MV::node_t> map_type;
      TEUCHOS_TEST_FOR_EXCEPTION(
        mv.getRawPtr () == NULL, std::invalid_argument,
        "Amesos2::get_1d_copy_helper::do_get(3 args): mv is null.");

      Teuchos::RCP<const map_type> map = mv->getMap ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        map.is_null (), std::invalid_argument,
        "Amesos2::get_1d_copy_helper::do_get(3 args): mv->getMap() is null.");


      do_get (mv, vals, ldx, Teuchos::ptrInArg (*map), ROOTED); // ROOTED the default here for now
    }

    template <class MV, typename KV>
    void get_1d_copy_helper_kokkos_view<MV,KV>::
    do_get (const Teuchos::Ptr<const MV>& mv,
            KV& kokkos_vals,
            const size_t ldx,
            Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
            EDistribution distribution)
    {
      // Dispatch to the copy function appropriate for the type
      if_then_else<is_same<typename MV::scalar_t,typename KV::value_type>::value,
        same_type_get_copy_kokkos_view<MV,KV>, // MDM-TODO Temp change this back to same_type
        diff_type_get_copy_kokkos_view<MV,KV> >::type::apply (mv, kokkos_vals, ldx, distribution_map, distribution);
    }

    template <class MV, typename KV>
    void get_1d_copy_helper_kokkos_view<MV,KV>::
    do_get (const Teuchos::Ptr<const MV>& mv,
            KV& kokkos_vals,
            const size_t ldx,
            EDistribution distribution,
            typename MV::global_ordinal_t indexBase)
    {
      typedef typename MV::local_ordinal_t lo_t;
      typedef typename MV::global_ordinal_t go_t;
      typedef typename MV::global_size_t gs_t;
      typedef typename MV::node_t node_t;

      TEUCHOS_TEST_FOR_EXCEPTION(
        mv.getRawPtr () == NULL, std::invalid_argument,
        "Amesos2::get_1d_copy_helper_kokkos_view::do_get(5 args): mv is null.");

      Teuchos::RCP<const Tpetra::Map<lo_t,go_t,node_t> > map
        = Amesos2::Util::getDistributionMap<lo_t,go_t,gs_t,node_t> (distribution,
                                                                    mv->getGlobalLength (),
                                                                    mv->getComm (),
                                                                    indexBase,
                                                                    mv->getMap());

      do_get (mv, kokkos_vals, ldx, Teuchos::ptrInArg (*map), distribution);
    }

    template <class MV, typename KV>
    void get_1d_copy_helper_kokkos_view<MV,KV>::do_get(const Teuchos::Ptr<const MV>& mv,
                                          KV& kokkos_vals,
                                          const size_t ldx)
    {
      typedef Tpetra::Map<typename MV::local_ordinal_t,
                          typename MV::global_ordinal_t,
                          typename MV::node_t> map_type;
      TEUCHOS_TEST_FOR_EXCEPTION(
        mv.getRawPtr () == NULL, std::invalid_argument,
        "Amesos2::get_1d_copy_helper::do_get(3 args): mv is null.");

      Teuchos::RCP<const map_type> map = mv->getMap ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        map.is_null (), std::invalid_argument,
        "Amesos2::get_1d_copy_helper_kokkos_view::do_get(3 args): mv->getMap() is null.");

      do_get (mv, kokkos_vals, ldx, Teuchos::ptrInArg (*map), ROOTED); // ROOTED the default here for now
    }


    ///////////////////////////
    // Copy-puting utilities //
    ///////////////////////////

    template <typename MV>
    void same_type_data_put<MV>::apply(const Teuchos::Ptr<MV>& mv,
                                       const Teuchos::ArrayView<typename MV::scalar_t>& data,
                                       const size_t ldx,
                                       Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                       EDistribution distribution )
    {
      mv->put1dData (data, ldx, distribution_map, distribution);
    }

    template <typename MV, typename KV>
    void same_type_data_put_kokkos_view<MV,KV>::apply(const Teuchos::Ptr<MV>& mv,
                                       KV& kokkos_data,
                                       const size_t ldx,
                                       Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                       EDistribution distribution )
    {
      mv->template put1dData_kokkos_view<KV>(kokkos_data, ldx, distribution_map, distribution);
    }

    /*
     * In the case where the scalar type of the multi-vector and the
     * corresponding S type are different, then we need to first get a
     * copy of the scalar values, then convert each one into the S
     * type before inserting into the vals array.
     */
    template <typename MV, typename S>
    void diff_type_data_put<MV,S>::apply(const Teuchos::Ptr<MV>& mv,
                                         const Teuchos::ArrayView<S>& data,
                                         const size_t ldx,
                                         Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                         EDistribution distribution )
    {
      typedef typename MV::scalar_t mv_scalar_t;
      typedef typename Teuchos::Array<mv_scalar_t>::size_type size_type;

      TEUCHOS_TEST_FOR_EXCEPTION(
        mv.getRawPtr () == NULL, std::invalid_argument,
        "Amesos2::diff_type_data_put(4 args): mv is null.");

      const size_type vals_length = data.size ();
      Teuchos::Array<mv_scalar_t> data_tmp (vals_length);

      for (size_type i = 0; i < vals_length; ++i) {
        data_tmp[i] = Teuchos::as<mv_scalar_t> (data[i]);
      }

      mv->put1dData (data_tmp (), ldx, distribution_map, distribution);
    }

    template <typename MV, typename KV>
    void diff_type_data_put_kokkos_view<MV,KV>::apply(const Teuchos::Ptr<MV>& mv,
                                         KV& kokkos_data,
                                         const size_t ldx,
                                         Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                         EDistribution distribution )
    {
      typedef typename MV::scalar_t mv_scalar_t;

      TEUCHOS_TEST_FOR_EXCEPTION(
        mv.getRawPtr () == NULL, std::invalid_argument,
        "Amesos2::diff_type_data_put_kokkos_view(4 args): mv is null.");

      typedef Kokkos::View<mv_scalar_t**, typename KV::array_layout, typename KV::execution_space> matrix_kokkos_view_t;
      matrix_kokkos_view_t matrix_kokkos_data (Kokkos::ViewAllocateWithoutInitializing("data_tmp"),
        kokkos_data.extent(0), kokkos_data.extent(1));

      // MDM-TODO inquire if deep copy is ok instead of loop (would bypass Teuchos::as)
      // Kokkos::deep_copy(matrix_kokkos_data, kokkos_data);
      Kokkos::parallel_for(
        Kokkos::RangePolicy<typename KV::execution_space, int> (0, kokkos_data.size()),
        KOKKOS_LAMBDA (size_t n) {
        // MDM-TODO can we remove the data() access without looping in 2D?
        matrix_kokkos_data.data()[n] = Teuchos::as<mv_scalar_t>(kokkos_data.data()[n]);
      });

      mv->template put1dData_kokkos_view<matrix_kokkos_view_t> (matrix_kokkos_data, ldx, distribution_map, distribution);
    }

    /** \internal
     *
     * \brief Helper class for putting 1-D data arrays into multivectors
     *
     * Handles dataype conversion when necessary before putting the data.
     */
    template <class MV, typename S>
    void put_1d_data_helper<MV,S>::do_put(const Teuchos::Ptr<MV>& mv,
                                          const Teuchos::ArrayView<S>& data,
                                          const size_t ldx,
                                          Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                          EDistribution distribution )
    {
      // Dispatch to the copy function appropriate for the type
      if_then_else<is_same<typename MV::scalar_t,S>::value,
        same_type_data_put<MV>,
        diff_type_data_put<MV,S> >::type::apply(mv, data, ldx, distribution_map, distribution);
    }

    template <class MV, typename S>
    void put_1d_data_helper<MV,S>::do_put(const Teuchos::Ptr<MV>& mv,
                                          const Teuchos::ArrayView<S>& data,
                                          const size_t ldx,
                                          EDistribution distribution,  typename MV::global_ordinal_t indexBase)
    {
      typedef typename MV::local_ordinal_t lo_t;
      typedef typename MV::global_ordinal_t go_t;
      typedef typename MV::global_size_t gs_t;
      typedef typename MV::node_t node_t;

      const Teuchos::RCP<const Tpetra::Map<lo_t,go_t,node_t> > map
        = Amesos2::Util::getDistributionMap<lo_t,go_t,gs_t,node_t>(distribution,
                                                                   mv->getGlobalLength(),
                                                                   mv->getComm(), 
                                                                   indexBase,
                                                                   mv->getMap());

      do_put(mv, data, ldx, Teuchos::ptrInArg(*map), distribution);
    }

    template <class MV, typename S>
    void put_1d_data_helper<MV,S>::do_put (const Teuchos::Ptr<MV>& mv,
                                           const Teuchos::ArrayView<S>& data,
                                           const size_t ldx)
    {
      const Teuchos::RCP<const Tpetra::Map<typename MV::local_ordinal_t,
        typename MV::global_ordinal_t,
        typename MV::node_t> > map
        = mv->getMap();
      do_put (mv, data, ldx, Teuchos::ptrInArg (*map), ROOTED); // Default as ROOTED for now
    }

    template <class MV, typename KV>
    void put_1d_data_helper_kokkos_view<MV,KV>::do_put(const Teuchos::Ptr<MV>& mv,
                                          KV& kokkos_data,
                                          const size_t ldx,
                                          Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                          EDistribution distribution )
    {
      // Dispatch to the copy function appropriate for the type
      if_then_else<is_same<typename MV::scalar_t,typename KV::value_type>::value,
        same_type_data_put_kokkos_view<MV,KV>,
        diff_type_data_put_kokkos_view<MV,KV> >::type::apply(mv, kokkos_data, ldx, distribution_map, distribution);
    }

    template <class MV, typename KV>
    void put_1d_data_helper_kokkos_view<MV,KV>::do_put(const Teuchos::Ptr<MV>& mv,
                                          KV& kokkos_data,
                                          const size_t ldx,
                                          EDistribution distribution,  typename MV::global_ordinal_t indexBase)
    {
      typedef typename MV::local_ordinal_t lo_t;
      typedef typename MV::global_ordinal_t go_t;
      typedef typename MV::global_size_t gs_t;
      typedef typename MV::node_t node_t;

      const Teuchos::RCP<const Tpetra::Map<lo_t,go_t,node_t> > map
        = Amesos2::Util::getDistributionMap<lo_t,go_t,gs_t,node_t>(distribution,
                                                                   mv->getGlobalLength(),
                                                                   mv->getComm(),
                                                                   indexBase,
                                                                   mv->getMap());

      do_put(mv, kokkos_data, ldx, Teuchos::ptrInArg(*map), distribution);
    }

    template <class MV, typename KV>
    void put_1d_data_helper_kokkos_view<MV,KV>::do_put (const Teuchos::Ptr<MV>& mv,
                                           KV& kokkos_data,
                                           const size_t ldx)
    {
      const Teuchos::RCP<const Tpetra::Map<typename MV::local_ordinal_t,
        typename MV::global_ordinal_t,
        typename MV::node_t> > map
        = mv->getMap();
      do_put (mv, kokkos_data, ldx, Teuchos::ptrInArg (*map), ROOTED); // Default as ROOTED for now
    }

  } // end namespace Util

} // end namespace Amesos2

#endif  // AMESOS2_EPETRAMULTIVECADAPTER_DEF
