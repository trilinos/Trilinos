// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_OVERLAPPINGPARTITIONER_DECL_HPP
#define IFPACK2_OVERLAPPINGPARTITIONER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Partitioner.hpp"
#include "Tpetra_RowGraph.hpp"

namespace Ifpack2 {

/// \class OverlappingPartitioner
/// \brief Create overlapping partitions of a local graph.
/// \tparam GraphType Specialization of Tpetra::CrsGraph or
///   Tpetra::RowGraph.
///
/// This class enables the extension of nonoverlapping partitions to
/// an arbitrary amount of constant overlap.  "Overlap" here refers to the
/// overlap among the local parts, and not the overlap among the
/// processes.
///
/// For partitions with varying overlap, use UserPartitioner, which inherits
/// from OverlappingPartitioner.
///
/// Supported parameters are:
/// - "partitioner: local parts" (<tt>local_ordinal_type</tt>): the
///   number of local partitions.  Default is one (one local
///   partition, which means "do not partition locally").
/// - "partitioner: overlap" (int): the amount of overlap between
///   partitions.  Default is zero (no overlap).
///
/// This class implements common functionality for derived classes
/// like LinearPartitioner.  Graphs given as input to any derived
/// class must contain no singletons.  Usually, this means that the
/// graph is the graph of a Tpetra::RowMatrix that has been filtered
/// using SingletonFilter.
template<class GraphType>
class OverlappingPartitioner : public Partitioner<GraphType> {
public:
  typedef typename GraphType::local_ordinal_type local_ordinal_type;
  typedef typename GraphType::global_ordinal_type global_ordinal_type;
  typedef typename GraphType::node_type node_type;
  typedef typename GraphType::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  typedef typename GraphType::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
 
  typedef Tpetra::RowGraph<local_ordinal_type, global_ordinal_type, node_type> row_graph_type;

  //! Constructor.
  OverlappingPartitioner (const Teuchos::RCP<const row_graph_type>& graph);

  //! Destructor.
  virtual ~OverlappingPartitioner();

  /// \brief Number of computed local partitions.
  ///
  /// This is a (signed) \c int because negative values are
  /// significant.  We can't use \c local_ordinal_type here, because
  /// that type may be either signed or unsigned.
  int numLocalParts () const;

  //! The number of levels of overlap.
  int overlappingLevel() const;

  /// \brief Local index of the nonoverlapping partition of the given row.
  ///
  /// \param MyRow [in] Local index of the row.
  local_ordinal_type operator () (const local_ordinal_type MyRow) const;

  //! Local index of the overlapping partition of the j-th vertex in partition i.
  local_ordinal_type
  operator () (const local_ordinal_type i, const local_ordinal_type j) const;

  //! the number of rows contained in the given partition.
  size_t numRowsInPart (const local_ordinal_type Part) const;

  //! Fill \c List with the local indices of the rows in the (overlapping) partition \c Part.
  void rowsInPart (const local_ordinal_type Part, Teuchos::ArrayRCP<local_ordinal_type>& List) const;

  //! A view of the local indices of the nonoverlapping partitions of each local row.
  virtual Teuchos::ArrayView<const local_ordinal_type> 
  nonOverlappingPartition () const;

  //! Set all the parameters for the partitioner.
  /*! The supported parameters are:
   * - \c "partitioner: overlap" (int, default = 0).
   * - \c "partitioner: local parts" (int, default = 1).
   * - \c "partitioner: print level" (bool, default = false).
   */
  virtual void setParameters (Teuchos::ParameterList& List);

  //! Set all the parameters for the partitioner.
  /*! This function is used by derived classes to set their own
   * parameters. These classes should not derive SetParameters(),
   * so that common parameters can be set just once.
   */
  virtual void setPartitionParameters (Teuchos::ParameterList& List) = 0;

  //! Computes the partitions. Returns 0 if successful.
  virtual void compute();

  //! Computes the partitions. Returns 0 if successful.
  virtual void computePartitions() = 0;

  //! Computes the partitions. Returns 0 if successful.
  virtual void computeOverlappingPartitions();
  
  //! Returns true if partitions have been computed successfully.
  virtual bool isComputed() const;

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual std::ostream& print(std::ostream& os) const;

  //! @name Implementation of Teuchos::Describable 
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object. 
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}
 
protected:

  /// \brief Number of local subgraphs.
  ///
  /// This is a (signed) \c int because negative values are
  /// significant.  We can't use \c local_ordinal_type here, because
  /// that type may be either signed or unsigned.
  int NumLocalParts_;

  /// \brief Mapping from local row to partition number.
  ///
  /// \c Partition_[i] contains the local index of the nonoverlapping
  /// partition to which local row i belongs.  If the application has
  /// explicitly defined \c Parts_, then \c Partition_ is unused.
  Teuchos::Array<local_ordinal_type> Partition_;

  /// \brief Mapping from partition to all rows it contains.
  ///
  /// \c Used with 'partitioner: parts' (or 'partitioner: global ID parts')
  /// \c Parts_[i][j] is the local (or global) index of the j-th row contained in the
  /// (overlapping) partition i.
  Teuchos::Array<Teuchos::ArrayRCP<local_ordinal_type> > Parts_;

  //! The graph to be partitioned
  Teuchos::RCP<const row_graph_type> Graph_;

  //! Level of overlap.
  int OverlappingLevel_;

  //! If \c true,  the graph has been successfully partitioned.
  bool IsComputed_;

  //! If \c true, information are reported to stdout.
  bool verbose_;

  /// \brief If \c true, only add row to partition (block) if doing so won't add new columns to the column map.
  ///
  /// Default is false.  This option should decrease computational cost of the apply, but convergence
  /// will likely degrade.
  bool maintainSparsity_;
}; // class Ifpack2::OverlappingPartitioner

}// namespace Ifpack2

#endif // IFPACK2_OVERLAPPINGPARTITIONER_DECL_HPP
