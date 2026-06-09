// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_MULTIVECADAPTER_DEF_HPP
#define AMESOS2_MULTIVECADAPTER_DEF_HPP

#include "Amesos2_TpetraMultiVecAdapter_def.hpp"
#include "Amesos2_KokkosMultiVecAdapter_def.hpp"

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
    void same_type_get_copy<MV>::apply(const Teuchos::Ptr<const MV> mv,
                                       const Teuchos::ArrayView<typename MV::scalar_t>& v,
                                       const size_t ldx,
                                       Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                       EDistribution distribution )
    {
      mv->get1dCopy (v, ldx, distribution_map, distribution);
    }

    /*
     * In the case where the scalar type of the multi-vector and the
     * corresponding S type are different, then we need to first get a
     * copy of the scalar values, then convert each one into the S
     * type before inserting into the vals array.
     */
    template <typename MV, typename S>
    void diff_type_get_copy<MV,S>::
    apply (const Teuchos::Ptr<const MV> mv,
           const Teuchos::ArrayView<S>& v,
           const size_t ldx,
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
      std::conditional_t<std::is_same_v<typename MV::scalar_t,S>,
        same_type_get_copy<MV>,
        diff_type_get_copy<MV,S> >::apply (mv, vals, ldx, distribution_map, distribution);
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

    // KV
    template <class MV, typename KV>
    bool diff_type_get_view<MV, KV>::
    apply (bool bInitialize,
           const Teuchos::Ptr<const MV>& mv,
           KV& kokkos_vals,
           const size_t ldx,
           Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
           EDistribution distribution)
    {
      using input_scalar_type  = typename MV::scalar_t;
      using output_scalar_type = typename KV::non_const_value_type;
      using execution_space = typename KV::execution_space;
      using memory_space    = typename execution_space::memory_space;

      // copy to workspace (of same type as input scalar type)
      Kokkos::View<input_scalar_type**, Kokkos::LayoutLeft, memory_space> output_view;
      bool bAssigned = mv->get1dCopy_kokkos_view(bInitialize, output_view, ldx, distribution_map, distribution);

      // copy back (of the output scalar type)
      const size_t nrows = output_view.extent(0);
      const size_t ncols = output_view.extent(1);
      Kokkos::resize(kokkos_vals, ldx, ncols);

      #if 1
      auto h_output_view = Kokkos::create_mirror_view(output_view);
      auto h_kokkos_vals = Kokkos::create_mirror_view(kokkos_vals);
      Kokkos::deep_copy(output_view, h_output_view);
      for (size_t j=0; j<ncols; j++) {
        for (size_t i=0; i<nrows; i++) h_kokkos_vals(i,j) = Teuchos::as<output_scalar_type> (h_output_view(i,j));
      }
      Kokkos::deep_copy(kokkos_vals, h_kokkos_vals);
      #else
      //Kokkos::parallel_for("Amesos2::MultiVecAdapter::diff_type_get_view::apply", Kokkos::RangePolicy<execution_space>(0, nrows),
      //  KOKKOS_LAMBDA(const int i) { for (size_t j=0; j<ncols; j++) kokkos_vals(i,j) = Teuchos::as<output_scalar_type> (output_view(i,j)); });
      #endif

      // made a copy, instead of assigning input pointer (so safe to shallow-copy, which won't over-write the original memory space)
      bAssigned = false;
      return bAssigned;
    }

    template <class MV, class KV>
    bool same_type_get_view<MV, KV>::
    apply (bool bInitialize,
           const Teuchos::Ptr<const MV>& mv,
           KV& kokkos_vals,
           const size_t ldx,
           Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
           EDistribution distribution)
    {
      return mv->get1dCopy_kokkos_view(bInitialize, kokkos_vals, ldx, distribution_map, distribution);
    }

    // > entry point
    template <class MV, typename KV>
    bool get_1d_copy_helper_kokkos_view<MV,KV>::
    do_get (bool bInitialize,
            const Teuchos::Ptr<const MV>& mv,
            KV& kokkos_vals,
            const size_t ldx,
            Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
            EDistribution distribution)
    {
      // TODO: check input-scalar (host vs non-host)
      //using input_scalar_type  = typename MV::scalar_t;
      using input_scalar_type  = typename MV::host_value_t;
      using output_scalar_type = typename KV::non_const_value_type;

      // call appropriate apply based on if the input and output scalar types are the same or different
      bool bAssigned =  std::conditional_t<std::is_same_v<input_scalar_type, output_scalar_type>,
               same_type_get_view<MV, KV>,
               diff_type_get_view<MV, KV>>::apply (bInitialize, mv, kokkos_vals, ldx, distribution_map, distribution);
      return bAssigned;
    }

    template <class MV, typename KV>
    bool get_1d_copy_helper_kokkos_view<MV,KV>::
    do_get (bool bInitialize,
            const Teuchos::Ptr<const MV>& mv,
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

      return do_get (bInitialize, mv, kokkos_vals, ldx, Teuchos::ptrInArg (*map), distribution);
    }

    template <class MV, typename KV>
    bool get_1d_copy_helper_kokkos_view<MV,KV>::
    do_get (bool bInitialize,
            const Teuchos::Ptr<const MV>& mv,
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
      //NOTE: KokkosMultiVecAdapter returns null map
      //TEUCHOS_TEST_FOR_EXCEPTION(
      //  map.is_null (), std::invalid_argument,
      //  "Amesos2::get_1d_copy_helper_kokkos_view::do_get(3 args): mv->getMap() is null.");

      return do_get (bInitialize, mv, kokkos_vals, ldx, Teuchos::ptrInArg (*map), ROOTED); // ROOTED the default here for now
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
      std::conditional_t<std::is_same_v<typename MV::scalar_t,S>,
        same_type_data_put<MV>,
        diff_type_data_put<MV,S> >::apply(mv, data, ldx, distribution_map, distribution);
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

    // KV
    template <class MV, class KV>
    void same_type_put_view<MV,KV>::apply(const Teuchos::Ptr<MV>& mv,
                                          KV& kokkos_data,
                                          const size_t ldx,
                                          Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                          EDistribution distribution )
    {
      mv->put1dData_kokkos_view(kokkos_data, ldx, distribution_map, distribution);
    }

    template <class MV, class KV>
    void diff_type_put_view<MV,KV>::apply(const Teuchos::Ptr<MV>& mv,
                                          KV& kokkos_data,
                                          const size_t ldx,
                                          Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                          EDistribution distribution )
    {
      // TODO: check input-scalar (host vs non-host)
      //using output_scalar_type  = typename MV::scalar_t;
      using output_scalar_type  = typename MV::host_value_t;
      using execution_space = typename KV::execution_space;
      using memory_space    = typename execution_space::memory_space;

      // copy to input workspace (of same type as output scalar type)
      const size_t nrows = kokkos_data.extent(0);
      const size_t ncols = kokkos_data.extent(1);
      Kokkos::View<output_scalar_type**, Kokkos::LayoutLeft, memory_space> output_view ("output_view", ldx, ncols);

      #if 1
      auto h_output_view = Kokkos::create_mirror_view(output_view);
      auto h_kokkos_data = Kokkos::create_mirror_view(kokkos_data);
      Kokkos::deep_copy(h_kokkos_data, kokkos_data);
      for (size_t j=0; j<ncols; j++) {
        for (size_t i=0; i<nrows; i++) h_output_view(i,j) = Teuchos::as<output_scalar_type> (h_kokkos_data(i,j));
      }
      Kokkos::deep_copy(output_view, h_output_view);
      #else
      // TODO: how to custom-cast within parallel-for
      Kokkos::parallel_for("Amesos2::MultiVecAdapter::diff_type_get_view::apply", Kokkos::RangePolicy<execution_space>(0, nrows),
        KOKKOS_LAMBDA(const int i) { for (size_t j=0; j<ncols; j++) output_view(i,j) = Teuchos::as<output_scalar_type> (kokkos_data(i,j)); });
      #endif

      // put to output (of output scalar type)
      mv->put1dData_kokkos_view(output_view, ldx, distribution_map, distribution);
    }

    template <class MV, typename KV>
    void put_1d_data_helper_kokkos_view<MV,KV>::do_put(const Teuchos::Ptr<MV>& mv,
                                          KV& kokkos_data,
                                          const size_t ldx,
                                          Teuchos::Ptr<const Tpetra::Map<typename MV::local_ordinal_t, typename MV::global_ordinal_t, typename MV::node_t> > distribution_map,
                                          EDistribution distribution )
    {
      // call appropriate apply based on if the input and output scalar types are the same or different
      using input_scalar_type   = typename KV::non_const_value_type;
      using output_scalar_type  = typename MV::scalar_t;
      std::conditional_t<std::is_same_v<input_scalar_type, output_scalar_type>,
               same_type_put_view<MV, KV>,
               diff_type_put_view<MV, KV>>::apply (mv, kokkos_data, ldx, distribution_map, distribution);
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

#endif  // AMESOS2_MULTIVECADAPTER_DEF_HPP
