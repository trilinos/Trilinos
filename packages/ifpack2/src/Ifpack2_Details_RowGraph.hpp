// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_ROWGRAPH_HPP
#define IFPACK2_DETAILS_ROWGRAPH_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_RowGraph.hpp"

namespace Ifpack2 {
namespace Details {

/// \class RowGraph
/// \brief All Ifpack2 implementations of Tpetra::RowGraph must
///   inherit from this class.
/// \tparam GraphType Tpetra::RowGraph specialization.
///
/// \warning This class is an implementation detail of Ifpack2.  Users
///   should not rely on its interface.
///
/// This class exists to facilitate Tpetra interface changes.  See
/// e.g., GitHub Issue #2630.
template<class GraphType>
class RowGraph :
    virtual public Tpetra::RowGraph<typename GraphType::local_ordinal_type,
                                    typename GraphType::global_ordinal_type,
                                    typename GraphType::node_type> {
public:
  //! \name Typedefs
  //@{
  typedef typename GraphType::local_ordinal_type local_ordinal_type;
  typedef typename GraphType::global_ordinal_type global_ordinal_type;
  typedef typename GraphType::node_type node_type;
  typedef typename GraphType::local_inds_host_view_type local_inds_host_view_type;
  typedef typename GraphType::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
  typedef typename GraphType::global_inds_host_view_type global_inds_host_view_type;
  typedef typename GraphType::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  //@}
  //! \name Destructor
  //@{

  //! Destructor (virtual for memory safety of derived classes)
  virtual ~RowGraph () = default;

  //@}
  /// \name Work-around implementations of deprecated virtual methods
  ///
  /// These methods exist to smooth the path for fixing GitHub Issue
  /// #2630.  This is why their existence depends on a Tpetra macro.
  //@{

};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_ROWGRAPH_HPP
