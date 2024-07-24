// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKVECTOR_DECL_HPP
#define TPETRA_BLOCKVECTOR_DECL_HPP

#include "Tpetra_BlockVector_fwd.hpp"
#include "Tpetra_BlockMultiVector.hpp"
#include "Tpetra_Vector.hpp"

namespace Tpetra {

/// \class BlockVector
/// \brief Vector for multiple degrees of freedom per mesh point
/// \author Mark Hoemmen
///
/// \tparam Scalar The type of each entry of the block vector.  (You
///   can use real-valued or complex-valued types here, unlike in
///   Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the documentation of
///   the first template parameter of Map for requirements.
/// \tparam GO The type of global indices.  See the documentation of
///   the second template parameter of Map for requirements.
/// \tparam Node The Kokkos Node type.  See the documentation of the
///   third template parameter of Map for requirements.
///
/// BlockVector is like ::Tpetra::MultiVector, but its interface
/// supports multiple degrees of freedom per mesh point.  You can
/// specify a mesh point by its local or global index, and read or
/// write the values at that point.  Every mesh point must have the
/// same number of degrees of freedom.  We call the number of degrees
/// of freedom per mesh point the <i>block size</i>.
///
/// BlockVector is a special case of BlockMultiVector, for
/// "multivectors" that are not "multi."  That is, a BlockVector has a
/// single vector (column).  Please refer to the documentation of
/// BlockMultiVector for details.
template<class Scalar,
         class LO,
         class GO,
         class Node>
class BlockVector : public BlockMultiVector<Scalar, LO, GO, Node> {
private:
  typedef BlockMultiVector<Scalar, LO, GO, Node> base_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

public:
  //! \name Typedefs to facilitate template metaprogramming.
  //@{

  //! The type of entries in the vector.
  typedef typename base_type::scalar_type scalar_type;
  //! The implementation type of entries in the vector.
  typedef typename base_type::impl_scalar_type impl_scalar_type;
  //! The type of local indices.
  typedef typename base_type::local_ordinal_type local_ordinal_type;
  //! The type of global indices.
  typedef typename base_type::global_ordinal_type global_ordinal_type;
  //! The Kokkos Node type.
  typedef typename base_type::node_type node_type;
  //! The Kokkos Device type.
  typedef typename Node::device_type device_type;

  //! The specialization of Tpetra::Map that this class uses.
  typedef Tpetra::Map<LO, GO, Node> map_type;
  //! The specialization of Tpetra::MultiVector that this class uses.
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
  //! The specialization of Tpetra::Vector that this class uses.
  typedef Tpetra::Vector<Scalar, LO, GO, Node> vec_type;

  /// \brief "Block view" of all degrees of freedom at a mesh point.
  ///
  /// A "block view" lets you address all degrees of freedom at a mesh
  /// point.  You don't have to use this class to access the degrees
  /// of freedom.  If you do choose to use this class, it implements
  /// operator()(LO i), so you can access and modify its entries.
  ///
  /// The preferred way to refer to the little_vec_type and
  /// const_little_vec_type types, is to get them from the typedefs
  /// below.  This is because different specializations of BlockVector
  /// reserve the right to use different types to implement
  /// little_vec_type or const_little_vec_type.  This gives us a
  /// porting strategy to move from "classic" Tpetra to the Kokkos
  /// refactor version.
  typedef typename base_type::little_vec_type 
                              little_vec_type;
  typedef typename base_type::little_host_vec_type 
                              little_host_vec_type;

  /// \brief "Const block view" of all degrees of freedom at a mesh point.
  ///
  /// This is just like little_vec_type, except that you can't modify
  /// its entries.
  typedef typename base_type::const_little_vec_type 
                              const_little_vec_type;
  typedef typename base_type::const_little_host_vec_type 
                              const_little_host_vec_type;

  //@}
  //! \name Constructors
  //@{

  /// \brief Default constructor.
  ///
  /// Creates an empty BlockVector.  An empty BlockVector has zero
  /// rows, and block size zero.
  BlockVector ();

  //! Copy constructor (shallow copy).
  BlockVector (const BlockVector<Scalar, LO, GO, Node>&) = default;

  //! Move constructor (shallow move).
  BlockVector (BlockVector<Scalar, LO, GO, Node>&&) = default;

  //! Copy assigment (shallow copy).
  BlockVector<Scalar, LO, GO, Node>&
  operator= (const BlockVector<Scalar, LO, GO, Node>&) = default;

  //! Move assigment (shallow move).
  BlockVector<Scalar, LO, GO, Node>&
  operator= (BlockVector<Scalar, LO, GO, Node>&&) = default;

  //! "Copy constructor" with option to deep copy.
  BlockVector (const BlockVector<Scalar, LO, GO, Node>& in,
               const Teuchos::DataAccess copyOrView);

  /// \brief Constructor that takes a mesh Map and a block size.
  ///
  /// \param meshMap [in] Map that describes the distribution of mesh
  ///   points (rather than the distribution of unknowns for those
  ///   mesh points).
  ///
  /// \param blockSize [in] The number of degrees of freedom per mesh
  ///   point.  We assume that this is the same for all mesh points in
  ///   the above Map.
  ///
  /// The <i>mesh Map</i> describes the distribution of mesh points.
  /// Its corresponding <i>point Map</i> describes the distribution of
  /// degrees of freedom corresponding to those mesh points.  If you
  /// have already computed the point Map corresponding to the above
  /// mesh Map, then it is more efficient to call the three-argument
  /// constructor below, that takes both the mesh Map and the point
  /// Map.
  ///
  /// There are two ways to get the point Map corresponding to a given
  /// mesh Map and block size.  You may either call the class method
  /// makePointMap() (inherited from the parent class
  /// BlockMultiVector), or you may call this two-argument
  /// constructor, and then call getPointMap().
  ///
  /// The point Map enables reinterpretation of a BlockVector as a
  /// standard Tpetra::Vector, or as a Tpetra::MultiVector with one
  /// column.  This lets users solve linear systems with Trilinos'
  /// solvers and preconditioners, that expect vectors as
  /// Tpetra::MultiVector or Tpetra::Vector instances.
  BlockVector (const map_type& meshMap, const LO blockSize);

  /// \brief Constructor that takes a mesh Map, a point Map, and a
  ///   block size.
  ///
  /// See the documentation of the two-argument constructor above.
  BlockVector (const map_type& meshMap,
               const map_type& pointMap,
               const LO blockSize);

  /// \brief View an existing Vector or MultiVector.
  ///
  /// \param X_mv [in/out] The Vector or MultiVector to view.  It MUST
  ///   have view semantics; otherwise this constructor throws.  Its
  ///   Map must be the same (in the sense of isSameAs) as the point
  ///   Map corresponding to the given mesh Map and block size.  If
  ///   this is a MultiVector, it must have only one column.
  ///
  /// \param meshMap [in] The mesh Map to use for interpreting the
  ///   given MultiVector or Vector (in place) as a BlockVector.
  ///
  /// \param blockSize [in] The number of degrees of freedom per mesh
  ///   point.  We assume that this is the same for all mesh points.
  BlockVector (const mv_type& X_mv,
               const map_type& meshMap,
               const LO blockSize);

  /// \brief View an existing Vector.
  ///
  /// \param X_vec [in/out] The Vector view.  It MUST have view
  ///   semantics; otherwise this constructor throws.  Its Map must be
  ///   the same (in the sense of isSameAs) as the point Map
  ///   corresponding to the given mesh Map and block size.
  ///
  /// \param meshMap [in] The mesh Map to use for interpreting the
  ///   given Vector (in place) as a BlockVector.
  ///
  /// \param blockSize [in] The number of degrees of freedom per mesh
  ///   point.  We assume that this is the same for all mesh points.
  BlockVector (const vec_type& X_vec,
               const map_type& meshMap,
               const LO blockSize);

  /// \brief View an existing BlockVector using a different mesh
  ///   Map, supplying the corresponding point Map.
  ///
  /// This method corresponds to MultiVector's "offset view" constructor.
  BlockVector (const BlockVector<Scalar, LO, GO, Node>& X,
               const map_type& newMeshMap,
               const map_type& newPointMap,
               const size_t offset = 0);

  /// \brief View an existing BlockVector using a different mesh
  ///   Map; compute the new point Map.
  ///
  /// This method corresponds to MultiVector's "offset view" constructor.
  BlockVector (const BlockVector<Scalar, LO, GO, Node>& X,
               const map_type& newMeshMap,
               const size_t offset = 0);

  //@}
  //! \name Access to Maps, the block size, and a Vector view.
  //@{

  /// \brief Get a Tpetra::Vector that views this BlockVector's data.
  ///
  /// This is how you can give a BlockVector to Trilinos' solvers and
  /// preconditioners.
  vec_type getVectorView ();

  //@}
  //! \name Fine-grained data access
  //@{

  /// \brief Replace all values at the given mesh point, using a local index.
  ///
  /// \param localRowIndex [in] Local index of the mesh point.
  /// \param vals [in] Input values with which to replace whatever
  ///   existing values are at the mesh point.
  ///
  /// \return true if successful, else false.  This method will
  ///   <i>not</i> succeed if the given local index of the mesh point
  ///   is invalid on the calling process.
  ///
  /// \note This method, the other "replace" and "sumInto" methods,
  ///   and the view methods, are marked const.  This is because they
  ///   do not change pointers.  They do, of course, change the values
  ///   in the BlockVector, but that does not require marking the
  ///   methods as nonconst.
  bool replaceLocalValues (const LO localRowIndex, const Scalar vals[]);

  /// \brief Replace all values at the given mesh point, using a global index.
  ///
  /// \param globalRowIndex [in] Global index of the mesh point.
  /// \param vals [in] Input values with which to sum into whatever
  ///   existing values are at the mesh point.
  ///
  /// \return true if successful, else false.  This method will
  ///   <i>not</i> succeed if the given global index of the mesh point
  ///   is invalid on the calling process.
  bool replaceGlobalValues (const GO globalRowIndex, const Scalar vals[]);

  /// \brief Sum into all values at the given mesh point, using a local index.
  ///
  /// \param localRowIndex [in] Local index of the mesh point.
  /// \param vals [in] Input values with which to replace whatever
  ///   existing values are at the mesh point.
  ///
  /// \return true if successful, else false.  This method will
  ///   <i>not</i> succeed if the given local index of the mesh point
  ///   is invalid on the calling process.
  bool sumIntoLocalValues (const LO localRowIndex, const Scalar vals[]);

  /// \brief Sum into all values at the given mesh point, using a global index.
  ///
  /// \param globalRowIndex [in] Global index of the mesh point.
  /// \param vals [in] Input values with which to replace whatever
  ///   existing values are at the mesh point.
  ///
  /// \return true if successful, else false.  This method will
  ///   <i>not</i> succeed if the given global index of the mesh point
  ///   is invalid on the calling process.
  bool sumIntoGlobalValues (const GO globalRowIndex, const Scalar vals[]);

  /// \brief Get a view of the degrees of freedom at the given mesh point,
  ///   using a local index.
  ///
  /// The preferred way to refer to little_vec_type is to get it from
  /// BlockVector's typedef.  This is because different
  /// specializations of BlockVector reserve the right to use
  /// different types to implement little_vec_type.  This gives us a
  /// porting strategy to move from "classic" Tpetra to the Kokkos
  /// refactor version.
  const_little_host_vec_type getLocalBlockHost (const LO localRowIndex,
                                                Access::ReadOnlyStruct) const;
  little_host_vec_type getLocalBlockHost (const LO localRowIndex,
                                          Access::OverwriteAllStruct);
  little_host_vec_type getLocalBlockHost (const LO localRowIndex,
                                          Access::ReadWriteStruct);
  //@}
};

} // namespace Tpetra

#endif // TPETRA_BLOCKMULTIVECTOR_DECL_HPP
