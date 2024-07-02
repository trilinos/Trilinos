// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_KokkosMultiVecAdapter_decl.hpp
  \author
  \date

  \brief  Amesos2::MultiVecAdapter specialization for the
          Kokkos::View class.
*/

#ifndef AMESOS2_KOKKOS_MULTIVEC_ADAPTER_DECL_HPP
#define AMESOS2_KOKKOS_MULTIVEC_ADAPTER_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
#include <Tpetra_Core.hpp>

#include "Amesos2_MultiVecAdapter_decl.hpp"
#include "Amesos2_Kokkos_View_Copy_Assign.hpp"

namespace Amesos2 {

  /**
   * \brief Amesos2 adapter for the Tpetra::MultiVector class.
   *
   * \ingroup amesos2_multivec_adapters
   */
  template< typename Scalar,
            typename ExecutionSpace >
  class MultiVecAdapter<Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >
  {
  public:
    // public type definitions
    typedef Tpetra::Map<>::node_type                                      node_t;
    typedef Kokkos::View<Scalar**,Kokkos::LayoutLeft, ExecutionSpace> multivec_t;
    typedef int                                                  local_ordinal_t;
    typedef Tpetra::Map<>::global_ordinal_type                  global_ordinal_t;
    typedef size_t                                                 global_size_t;
    typedef Scalar                                                      scalar_t;

    typedef Kokkos::View<scalar_t**, Kokkos::LayoutLeft, ExecutionSpace> kokkos_view_t;

    friend Teuchos::RCP<MultiVecAdapter<multivec_t> > createMultiVecAdapter<> (Teuchos::RCP<multivec_t>);
    friend Teuchos::RCP<const MultiVecAdapter<multivec_t> > createConstMultiVecAdapter<> (Teuchos::RCP<const multivec_t>);

    static const char* name;


  protected:
    // Do not allow direct construction of MultiVecAdapter's.  Only
    // allow construction through the non-member friend functions.

    /// Copy constructor
    MultiVecAdapter( const MultiVecAdapter<multivec_t>& adapter );

    /**
     * \brief Initialize an adapter from a multi-vector RCP.
     *
     * \param m An RCP pointing to the multi-vector which is to be wrapped.
     */
    MultiVecAdapter( const Teuchos::RCP<multivec_t>& m );


  public:

    ~MultiVecAdapter()
    { }

    /// Checks whether this multivector is local to the calling node.
    bool isLocallyIndexed() const
    {
      if(getComm()->getSize() == 1){
        return true;
      } // There may be other conditions to check
      return false;
    }

    bool isGloballyIndexed() const;

    Teuchos::RCP<const Tpetra::Map<
                         local_ordinal_t,
                         global_ordinal_t,
                         node_t > >
    getMap() const
    {
      return Teuchos::null; // serial only for Kokkos adapter right now
    }

    /// Returns the Teuchos::Comm object associated with this multi-vector
    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const
    {
      return Tpetra::getDefaultComm(); // serial only for Kokkos adapter right now
    }

    /// Get the length of vectors local to the calling node
    size_t getLocalLength() const
    {
      return mv_->extent(0);
    }


    /// Get the number of vectors on this node
    size_t getLocalNumVectors() const
    {
      return mv_->extent(1);
    }


    /// Get the length of vectors in the global space
    global_size_t getGlobalLength() const
    {
      return mv_->extent(0);
    }


    /// Get the number of global vectors
    global_size_t getGlobalNumVectors() const
    {
      return mv_->extent(1);
    }


    /// Return the stride between vectors on this node
    size_t getStride() const
    {
      return mv_->getStride();
    }


    /// Return \c true if this MV has constant stride between vectors on this node
    bool isConstantStride() const
    {
      return mv_->isConstantStride();
    }

    /// Return pointer to vector when number of vectors == 1 and single MPI process
    Scalar * getMVPointer_impl() const;


    Teuchos::RCP<multivec_t> clone() const;

    /**
     * \brief Copies the multivector's data into the user-provided vector.
     *
     *  Each vector of the multivector is placed \c lda apart in the
     *  given ArrayView.  Giving a distribution map is useful in the
     *  case where the data needs to end up on different processors
     *  than it currently resides.  For example, the SuperLU_DIST
     *  interface may receive a B multivector that is distributed
     *  across 13 processors, but only 12 of those 13 processors are
     *  in SuperLU_DIST's processor grid.  The rows of the multivector
     *  need to then be distributed amongst the 12 that are in the
     *  grid.
     *
     *  \param [in/out] A  user-supplied storage for multi-vector data
     *  \param [in] lda    user-supplied spacing for consecutive vectors
     *                     in \c A
     *  \param [in] distribution_map is a Tpetra::Map that describes the
     *                     desired distribution of the multivector's
     *                     data accross the calling processors.  The map
     *                     describes where the 'rows' of the multivector
     *                     will end up.
     *
     *  \param [in] distribution
     *
     *  \throw std::runtime_error Thrown if the space available in \c A
     *  is not large enough given \c lda , the value of \c global_copy ,
     *  and the number of vectors in \c this.
     */
    void
    get1dCopy ( const Teuchos::ArrayView<scalar_t>& av,
                size_t lda,
                Teuchos::Ptr<
                  const Tpetra::Map<local_ordinal_t,
                  global_ordinal_t,
                  node_t> > distribution_map,
                EDistribution distribution) const;

    template<typename KV>
    bool
    get1dCopy_kokkos_view (
      bool bInitialize, KV& kokkos_view,
      [[maybe_unused]] size_t lda,
      [[maybe_unused]] Teuchos::Ptr<
        const Tpetra::Map<local_ordinal_t,
        global_ordinal_t,
        node_t> 
      > distribution_map,
      [[maybe_unused]] EDistribution distribution
    ) const {
      bool bAssigned; // deep_copy_or_assign_view sets true if assigned (no deep copy)
      deep_copy_or_assign_view(bInitialize, kokkos_view, *mv_, bAssigned);
      return bAssigned;
    }

    /**
     * \brief Extracts a 1 dimensional view of this MultiVector's data
     *
     * Guarantees that the view returned will reside in contiguous storage.
     *
     * \warning
     * It is recommended to use the \c get1dCopy function, from a
     * data-hiding perspective. Use if you know what you are doing.
     *
     * \param local if \c true , each node will get a view of the vectors it is
     * in possession of.  The default, \c false , will give each calling node a
     * view of the global multivector.
     */
    Teuchos::ArrayRCP<scalar_t> get1dViewNonConst (bool local = false);

    /**
     * \brief Export data into the global MultiVector space.
     *
     * \param new_data The data to be exported into \c this.
     * \param source_map describes how the input array data is distributed
     *                 accross processors.  This data will be redistributed
     *                 to match the map of the adapted multivector.
     *  \param [in] distribution
     */
    void
    put1dData (const Teuchos::ArrayView<const scalar_t>& new_data,
              size_t lda,
              Teuchos::Ptr<
                const Tpetra::Map<local_ordinal_t,
                global_ordinal_t,
                node_t> > distribution_map,
              EDistribution) const;

    template<typename KV>
    void
    put1dData_kokkos_view (
      KV& kokkos_new_data,
      [[maybe_unused]] size_t lda,
      [[maybe_unused]] Teuchos::Ptr<
        const Tpetra::Map<local_ordinal_t,
        global_ordinal_t,
        node_t>
      > distribution_map,
      EDistribution
    ) const {
      deep_copy_or_assign_view(*mv_, kokkos_new_data);
    }


    //! Get a short description of this adapter class
    std::string description () const;

    //! Print a description of this adapter to the given output stream
    void
    describe (Teuchos::FancyOStream& os,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;

  private:

    //! The multivector which this adapter wraps
    Teuchos::RCP<multivec_t> mv_;

    //! The Tpetra::Export specialization used by this class.
    typedef Tpetra::Export<local_ordinal_t, global_ordinal_t, node_t> export_type;

    //! The Tpetra::Import specialization used by this class.
    typedef Tpetra::Import<local_ordinal_t, global_ordinal_t, node_t> import_type;

    /// \brief Used for data redistribution from the user's input
    ///   MultiVector to the solver's input MultiVector.
    ///
    /// This is an Export, because the user's input data need not
    /// necessarily be one-to-one, but the solver's input data must
    /// (presumably) always be one-to-one.
    mutable Teuchos::RCP<export_type> exporter_;

    /// \brief Used for data redistribution from the solver's output
    ///   MultiVector to the user's output MultiVector.
    ///
    /// This is an Import, because the user's output data need not
    /// necessarily be one-to-one, but the solver's output data must
    /// (presumably) always be one-to-one.
    mutable Teuchos::RCP<import_type> importer_;

  }; // end class MultiVecAdapter<Tpetra::MultiVector>

} // end namespace Amesos2


#endif // AMESOS2_KOKKOS_MULTIVEC_ADAPTER_DECL_HPP
