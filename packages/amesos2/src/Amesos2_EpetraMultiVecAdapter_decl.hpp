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

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>

#include <Epetra_MultiVector.h>

#include "Amesos2_MultiVecAdapter_decl.hpp"

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
    typedef int                                           global_ordinal_t;
    typedef size_t                                           global_size_t;
    typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  node_t;
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

    ~MultiVecAdapter()
    { }


    /// Checks whether this multi-vector is local to the calling node.
    bool isLocallyIndexed() const;

    bool isGloballyIndexed() const;


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
     * \note Vectors returned hold a copy of the date in the multi-vector.  So
     * any changes to the returned vector will not be represented in the
     * underlying multi-vector.
     */
    Teuchos::RCP<Tpetra::Vector<scalar_t,local_ordinal_t,global_ordinal_t,node_t> >
    getVectorNonConst( size_t j );


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
		    node_t> > distribution_map ) const;


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
		    node_t> > source_map );



    /// Get a short description of this adapter class
    std::string description() const;


    /// Print a description of this adapter to the Fancy Output Stream.
    void describe( Teuchos::FancyOStream& os,
		   const Teuchos::EVerbosityLevel verbLevel ) const;


  private:

    /// The multi-vector this adapter wraps
    Teuchos::RCP<multivec_t> mv_;

    mutable Teuchos::RCP<Epetra_Import> importer_;
    mutable Teuchos::RCP<Epetra_Export> exporter_;

    mutable Teuchos::RCP<const Epetra_BlockMap> mv_map_;

  };                              // end class MultiVecAdapter<NewMultiVec>

} // end namespace Amesos2


#endif // AMESOS2_EPETRA_MULTIVEC_ADAPTER_DECL_HPP
