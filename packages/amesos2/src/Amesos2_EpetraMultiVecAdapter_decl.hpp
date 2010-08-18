// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2010 Sandia Corporation
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

#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>

#include <Epetra_MultiVector.h>

#include "Amesos2_MultiVecAdapter_decl.hpp"

namespace Amesos {


/**
 * \brief Amesos2 adapter for the Epetra_MultiVector class.
 */
template <>
class MultiVecAdapter<Epetra_MultiVector>
{
public:

  // public type definitions
  typedef double                                                scalar_type;
  typedef int                                            local_ordinal_type;
  typedef int                                           global_ordinal_type;
  typedef size_t                                           global_size_type;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  node_type;
  typedef Epetra_MultiVector                                  multivec_type;

  static const char* name;


  /// Default Constructor
  MultiVecAdapter()
    : mv_(Teuchos::null)
    , l_mv_(Teuchos::null)
    , l_l_mv_(Teuchos::null)
    { }


  /// Copy constructor
  MultiVecAdapter( const MultiVecAdapter<multivec_type>& adapter )
    : mv_(adapter.mv_)
    , l_mv_(adapter.l_mv_)
    , l_l_mv_(adapter.l_l_mv_)
    { }

  /**
   * \brief Initialize an adapter from a multi-vector RCP.
   *
   * \param m An RCP pointing to the multi-vector which is to be wrapped.
   */
  MultiVecAdapter( const Teuchos::RCP<multivec_type>& m )
    : mv_(m)
    , l_mv_(Teuchos::null)
    , l_l_mv_(Teuchos::null)
    { }


  ~MultiVecAdapter()
    { }


  /**
   * \brief Scales values of \c this by a factor of \c alpha
   *
   * \f$ this = alpha * this\f$
   *
   * \param [in] alpha scalar factor
   *
   * \return A reference to \c this
   */
  MultiVecAdapter<multivec_type>& scale( const scalar_type alpha );


  /**
   * \brief Updates the values of \c this to \f$this = alpha*this + beta*B\f$
   *
   * \param [in] beta scalar coefficient of \c B
   * \param [in] B additive MultiVector
   * \param [in] alpha scalar coefficient of \c this
   *
   * \return A reference to \c this
   */
  MultiVecAdapter<multivec_type>& update(
    const scalar_type beta,
    const MultiVecAdapter<multivec_type>& B,
    const scalar_type alpha );


  /// Checks whether this multi-vector is local to the calling node.
  bool isLocal() const;


  /// Returns the Teuchos::Comm object associated with this multi-vector
  const Teuchos::RCP<const Teuchos::Comm<int> >& getComm() const;


  /// Get the length of vectors local to the calling node
  size_t getLocalLength() const;


  /// Get the number of vectors on this node
  size_t getLocalNumVectors() const;


  /// Get the length of vectors in the global space
  global_size_type getGlobalLength() const;


  /// Get the number of global vectors
  size_t getGlobalNumVectors() const;


  /// Return the stride between vectors on this node
  size_t getStride() const;


  /// Return \c true if this MV has constant stride between vectors on this node
  bool isConstantStride() const;


  /// Const vector access
  Teuchos::RCP<const Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >
  getVector( size_t j ) const;


  /**
   * \brief Nonconst vector access
   *
   * \note Vectors returned hold a copy of the date in the multi-vector.  So
   * any changes to the returned vector will not be represented in the
   * underlying multi-vector.
   */
  Teuchos::RCP<Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >
  getVectorNonConst( size_t j );


  /**
   * \brief Copies the multi-vector's data into the user-provided vector.
   *
   *  Each multi-vector is \c lda apart in memory.
   */
  void get1dCopy( const Teuchos::ArrayView<scalar_type>& A, size_t lda );


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
  Teuchos::ArrayRCP<scalar_type> get1dViewNonConst( bool local = false );

  /**
   * \brief Get a 2-D copy of the multi-vector.
   *
   * Copies a 2-D representation of the multi-vector into the user-supplied
   * array-of-arrays.
   *
   * \param A user-supplied storage for the 2-D copy.
   *
   * \throw std::length_error Thrown if the size of the given \c A is not
   * equal to the number of vectors in \c this .
   *
   * \throw std::length_error Thrown if the size of the arrays in \c A is not
   * equal to the global vector length.  Assumes all arrays in \c A are equal
   * in length.
   *
   * \throw std::runtime_error Thrown if there is an error copying vector
   * values into the user storage.
   */
  void get2dCopy( Teuchos::ArrayView<const Teuchos::ArrayView<scalar_type> > A );

  /**
   * \brief Extracts a 2 dimensional view of this multi-vector's data.
   *
   * Guarantees that the view returned will reside in contiguous storage.  The
   * data is not \c const , so it may be altered if desired.
   *
   * \warning
   * It is recommended to use the \c get2dCopy function, from a
   * data-hiding perspective. Use if you know what you are doing.
   *
   * \param local if \c true , each node will get a view of the vectors it is
   * in possession of.  The default, \c false , will give each calling node a
   * view of the global multi-vector.
   *
   * \return An array-of-arrays view of this multi-vector's data
   *
   * \note This function is not declared \c const as it normally would be,
   * since it must modify local copies of the vector data before returning the
   * result.
   */
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> >
  get2dViewNonConst( bool local = false );


  /**
   * \brief Exports any local changes to this MultiVector to the global space.
   *
   * \post The values in \c l_mv_ will be equal to the values in \c mv_
   */
  void globalize();


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
  template<typename Value_t>
  void globalize( const Teuchos::ArrayView<Value_t>& newVals );


  /// Get a short description of this adapter class
  std::string description() const;


  /// Print a description of this adapter to the Fancy Output Stream.
  void describe( Teuchos::FancyOStream& os,
    const Teuchos::EVerbosityLevel verbLevel ) const;


private:

  /**
   * \brief "Localizes" the wrapped \c mv_ .
   *
   * It defines the private maps \c o_map_ and \c l_map_ and imports global
   * data into the local node.  If \c mv_ is not distributed, this method does
   * nothing.
   *
   * It is intended to set things up properly for calls to \c get1dCopy() and
   * \c get1dView().
   *
   * \sa get1dCopy(), get1dView()
   */
  void localize();


  /**
   * \brief Get a Tpetra::Map that describes this MultiVector
   *
   * Not part of the MultiVecAdapter interface, but useful for other
   * adaptations.
   */
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> >
  getMap() const;


  /// The multi-vector this adapter wraps
  Teuchos::RCP<multivec_type> mv_;

  /**
   * \brief local multi-vector.
   *
   * Contains a local view of the entire multi-vector.
   */
  Teuchos::RCP<multivec_type> l_mv_;

  /**
   * \brief local-local multi-vector.
   *
   * Holds only a representation of the vectors local to the calling processor.
   */
  Teuchos::RCP<multivec_type> l_l_mv_;

  Teuchos::RCP<Epetra_BlockMap> o_map_, l_map_;

  Teuchos::RCP<Epetra_Import> importer_;
  Teuchos::RCP<Epetra_Export> exporter_;

};                              // end class MultiVecAdapter<NewMultiVec>

} // end namespace Amesos


#endif // AMESOS2_EPETRA_MULTIVEC_ADAPTER_DECL_HPP
