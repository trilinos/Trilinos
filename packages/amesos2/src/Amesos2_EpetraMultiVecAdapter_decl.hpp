// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_EpetraMultiVecAdapter_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Tue Jul 20 23:34:52 CDT 2010

  \brief  Amesos2::MultiVecAdapter specialization for the
          Epetra_MultiVector class.
*/

#ifndef AMESOS2_EPETRA_MULTIVEC_ADAPTER_DECL_HPP
#define AMESOS2_EPETRA_MULTIVEC_ADAPTER_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>

#include <Epetra_MultiVector.h>

#include "Amesos2_MultiVecAdapter_decl.hpp"
#include "Amesos2_Kokkos_View_Copy_Assign.hpp"

namespace Amesos2 {

  /**
   * \brief Amesos2 adapter for the Epetra_MultiVector class.
   *
   * \ingroup amesos2_multivec_adapters
   */
  template <>
  class MultiVecAdapter<Epetra_MultiVector>
  {
  public:

    // public type definitions
    typedef double                                                scalar_t;
    typedef int                                            local_ordinal_t;
    typedef Tpetra::Map<>::global_ordinal_type            global_ordinal_t;
    typedef size_t                                           global_size_t;
    typedef Tpetra::Map<>::node_type                                node_t;
    typedef Epetra_MultiVector                                  multivec_t;

    friend Teuchos::RCP<MultiVecAdapter<multivec_t> > createMultiVecAdapter<>(Teuchos::RCP<multivec_t>);
    friend Teuchos::RCP<const MultiVecAdapter<multivec_t> > createConstMultiVecAdapter<>(Teuchos::RCP<const multivec_t>);


    static const char* name;


  protected:
    /// Copy constructor
    MultiVecAdapter( const MultiVecAdapter<multivec_t>& adapter );

    /**
     * \brief Initialize an adapter from a multi-vector RCP.
     *
     * \param m An RCP pointing to the multi-vector which is to be wrapped.
     */
    MultiVecAdapter( const Teuchos::RCP<multivec_t>& m );


  public:

    ~MultiVecAdapter() = default;


    /// Checks whether this multi-vector is local to the calling node.
    bool isLocallyIndexed() const;

    bool isGloballyIndexed() const;


    Teuchos::RCP<Epetra_MultiVector> clone() const;

    /**
     * \brief Get a Tpetra::Map that describes this MultiVector
     *
     * Not part of the MultiVecAdapter interface, but useful for other
     * adaptations.
     */
    Teuchos::RCP<const Tpetra::Map<local_ordinal_t, global_ordinal_t, node_t> >
    getMap() const;

    /// Returns the Teuchos::Comm object associated with this multi-vector
    const Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;


    /// Get the length of vectors local to the calling node
    size_t getLocalLength() const;


    /// Get the number of vectors on this node
    size_t getLocalNumVectors() const;


    /// Get the length of vectors in the global space
    global_size_t getGlobalLength() const;


    /// Get the number of global vectors
    size_t getGlobalNumVectors() const;


    /// Return the stride between vectors on this node
    size_t getStride() const;


    /// Return \c true if this MV has constant stride between vectors on this node
    bool isConstantStride() const;


    /// Const vector access
    Teuchos::RCP<const Tpetra::Vector<scalar_t,local_ordinal_t,global_ordinal_t,node_t> >
    getVector( size_t j ) const;


    /**
     * \brief Nonconst vector access
     *
     * \note Vectors returned hold a copy of the data in the multi-vector.  So
     * any changes to the returned vector will not be represented in the
     * underlying multi-vector.
     */
    Teuchos::RCP<Tpetra::Vector<scalar_t,local_ordinal_t,global_ordinal_t,node_t> >
    getVectorNonConst( size_t j );


    /// Return pointer to vector when number of vectors == 1 and num procs == 1
    double * getMVPointer_impl() const;

    /**
     * \brief Copies the multi-vector's data into the user-provided vector.
     *
     *  Each multi-vector is \c lda apart in memory.
     */
    void get1dCopy( const Teuchos::ArrayView<scalar_t>& A,
                    size_t lda,
                    Teuchos::Ptr<
                    const Tpetra::Map<local_ordinal_t,
                    global_ordinal_t,
                    node_t> > distribution_map,
        EDistribution distribution) const;

    template<typename KV>
    bool get1dCopy_kokkos_view(bool bInitialize, KV & A,
                    size_t lda,
                    Teuchos::Ptr<
                      const Tpetra::Map<local_ordinal_t,
                      global_ordinal_t,
                      node_t> > distribution_map,
                    EDistribution distribution) const {
      Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace> host_new_data;
      get1dCopy_kokkos_view_host(host_new_data, lda, distribution_map, distribution);
      bool bAssigned;
      deep_copy_or_assign_view(bInitialize, A, host_new_data, bAssigned);
      return false; // currently Epetra and prior get1dCopy_kokkos_view_host call cannot get direct assignment so always return false
    }

    void get1dCopy_kokkos_view_host(
                    Kokkos::View<scalar_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> & new_data,
                    size_t lda,
                    Teuchos::Ptr<
                      const Tpetra::Map<local_ordinal_t,
                      global_ordinal_t,
                      node_t> > distribution_map,
                    EDistribution) const;

    /**
     * \brief Extracts a 1 dimensional view of this multi-vector's data
     *
     * Guarantees that the view returned will reside in contiguous storage.
     *
     * \warning
     * It is recommended to use the \c get1dCopy function, from a
     * data-hiding perspective. Use if you know what you are doing.
     *
     * \param local if \c true , each node will get a view of the vectors it is
     * in possession of.  The default, \c false , will give each calling node a
     * view of the global multi-vector.
     *
     * \note This function is not declared \c const as it normally would be,
     * since it must modify local copies of the vector data before returning the
     * result.
     */
    Teuchos::ArrayRCP<scalar_t> get1dViewNonConst( bool local = false );


    /**
     * \brief Export \c newVals into the global MultiVector space.
     *
     * \note we assume the vectors in newVals have the same leading dimension as
     * those in \c this
     *
     * \tparam Value_t The type of the data values that are being put into \c mv_
     *
     * \param newVals The values to be exported into the global space.
     */
    void put1dData( const Teuchos::ArrayView<const scalar_t>& new_data,
                    size_t lda,
                    Teuchos::Ptr<
                    const Tpetra::Map<local_ordinal_t,
                    global_ordinal_t,
                    node_t> > source_map,
        EDistribution distribution );

    template<typename KV>
    void put1dData_kokkos_view(
                    KV & new_data,
                    size_t lda,
                      Teuchos::Ptr<
                      const Tpetra::Map<local_ordinal_t,
                      global_ordinal_t,
                      node_t> > source_map,
                    EDistribution distribution ) {
      Kokkos::View<scalar_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> host_new_data(
        Kokkos::ViewAllocateWithoutInitializing("host_new_data"),
        new_data.extent(0), new_data.extent(1));
      Kokkos::deep_copy(host_new_data, new_data);
      put1dData_kokkos_view_host(host_new_data, lda, source_map, distribution);
    }

    void put1dData_kokkos_view_host(
                    Kokkos::View<scalar_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> & new_data,
                    size_t lda,
                    Teuchos::Ptr<
                    const Tpetra::Map<local_ordinal_t,
                    global_ordinal_t,
                    node_t> > source_map,
        EDistribution distribution );

    /// Get a short description of this adapter class
    std::string description() const;


    /// Print a description of this adapter to the Fancy Output Stream.
    void describe( Teuchos::FancyOStream& os,
                   const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;


  private:

    /// The multi-vector this adapter wraps
    Teuchos::RCP<multivec_t> mv_;

    mutable Teuchos::RCP<Epetra_Import> importer_;
    mutable Teuchos::RCP<Epetra_Export> exporter_;

    mutable Teuchos::RCP<const Epetra_BlockMap> mv_map_;

  };                              // end class MultiVecAdapter<NewMultiVec>

} // end namespace Amesos2


#endif // AMESOS2_EPETRA_MULTIVEC_ADAPTER_DECL_HPP
