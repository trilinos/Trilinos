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
  \file   Amesos2_TpetraMultiVecAdapter_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Wed May 26 19:49:10 CDT 2010

  \brief  Amesos2::MultiVecAdapter specialization for the
	  Tpetra::MultiVector class.
*/

#ifndef AMESOS2_TPETRA_MULTIVEC_ADAPTER_DECL_HPP
#define AMESOS2_TPETRA_MULTIVEC_ADAPTER_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
#include <Tpetra_MultiVector.hpp>

#include "Amesos2_MultiVecAdapter_decl.hpp"

namespace Amesos2 {

  /**
   * \brief Amesos2 adapter for the Tpetra::MultiVector class.
   *
   * \ingroup amesos2_multivec_adapters
   */
  template< typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    class    Node >
  class MultiVecAdapter<Tpetra::MultiVector<Scalar,
					    LocalOrdinal,
					    GlobalOrdinal,
					    Node> >
  {
  public:

    // public type definitions
    typedef Tpetra::MultiVector<Scalar,
				LocalOrdinal,
				GlobalOrdinal,
				Node>       multivec_t;
    typedef Scalar                          scalar_t;
    typedef LocalOrdinal                    local_ordinal_t;
    typedef GlobalOrdinal                   global_ordinal_t;
    typedef Node                            node_t;
    typedef Tpetra::global_size_t           global_size_t;

    friend Teuchos::RCP<MultiVecAdapter<multivec_t> > createMultiVecAdapter<>(Teuchos::RCP<multivec_t>);
    friend Teuchos::RCP<const MultiVecAdapter<multivec_t> > createConstMultiVecAdapter<>(Teuchos::RCP<const multivec_t>);

    static const char* name;


  protected:
    // Do not allow direct construction of MultiVecAdapter's.  Only
    // allow construction throw the non-member friend functions.
    
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
    }

    // TODO
    bool isGloballyIndexed() const;


    const Teuchos::RCP<const Tpetra::Map<
			 local_ordinal_t,
			 global_ordinal_t,
			 node_t > >&
    getMap() const
    {
      return mv_->getMap();
    }

    /// Returns the Teuchos::Comm object associated with this multi-vector
    const Teuchos::RCP<const Teuchos::Comm<int> >& getComm() const
    {
      return getMap()->getComm();
    }

    /// Get the length of vectors local to the calling node
    size_t getLocalLength() const
    {
      return Teuchos::as<size_t>(mv_->getLocalLength());
    }


    /// Get the number of vectors on this node
    size_t getLocalNumVectors() const
    {
      return mv_->getNumVectors();
    }


    /// Get the length of vectors in the global space
    global_size_t getGlobalLength() const
    {
      return mv_->getGlobalLength();
      //return getMap()->getMaxAllGlobalIndex() + 1;
    }


    /// Get the number of global vectors
    global_size_t getGlobalNumVectors() const
    {
      return Teuchos::as<global_size_t>(mv_->getNumVectors());
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


    /// Const vector access
    Teuchos::RCP<const Tpetra::Vector<scalar_t,local_ordinal_t,global_ordinal_t,node_t> >
    getVector( size_t j ) const
    {
      return mv_->getVector(j);
    }


    /// Nonconst vector access
    Teuchos::RCP<Tpetra::Vector<scalar_t,local_ordinal_t,global_ordinal_t,node_t> >
    getVectorNonConst( size_t j )
    {
      return mv_->getVectorNonConst(j);
    }


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
     *  \throw std::runtime_error Thrown if the space available in \c A
     *  is not large enough given \c lda , the value of \c global_copy ,
     *  and the number of vectors in \c this.
     */
    void get1dCopy( const Teuchos::ArrayView<scalar_t>& A,
		    size_t lda,
		    Teuchos::Ptr<
		      const Tpetra::Map<local_ordinal_t,
		                        global_ordinal_t,
		                        node_t> > distribution_map ) const;

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
    Teuchos::ArrayRCP<scalar_t> get1dViewNonConst( bool local = false );

    /**
     * \brief Export data into the global MultiVector space.
     *
     * \param new_data The data to be exported into \c this.
     * \param source_map describes how the input array data is distributed
     *                 accross processors.  This data will be redistributed
     *                 to match the map of the adapted multivector.
     */
    void put1dData( const Teuchos::ArrayView<const scalar_t>& new_data,
		    size_t lda,
		    Teuchos::Ptr<
		      const Tpetra::Map<local_ordinal_t,
		                        global_ordinal_t,
		                        node_t> > source_map );


    /// Get a short description of this adapter class
    std::string description() const;


    /// Print a description of this adapter to the Fancy Output Stream.
    void describe( Teuchos::FancyOStream& os,
		   const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default ) const;


  private:

    /// The multivector this adapter wraps
    Teuchos::RCP<multivec_t> mv_;


    /// Used for transferring between local and global multivectors
    mutable Teuchos::RCP<Tpetra::Import<local_ordinal_t,
					global_ordinal_t,
					node_t> > importer_;
    mutable Teuchos::RCP<Tpetra::Export<local_ordinal_t,
					global_ordinal_t,
					node_t> > exporter_;

    /// The adapted multivector's map
    mutable Teuchos::RCP<const Tpetra::Map<local_ordinal_t,
					   global_ordinal_t,
					   node_t> > mv_map_;

  };                              // end class MultiVecAdapter<Tpetra::MultiVector>

} // end namespace Amesos2


#endif // AMESOS2_TPETRA_MULTIVEC_ADAPTER_DECL_HPP
