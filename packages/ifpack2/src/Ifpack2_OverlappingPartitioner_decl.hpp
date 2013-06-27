/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
//@HEADER
*/

#ifndef IFPACK2_OVERLAPPINGPARTITIONER_DECL_HPP
#define IFPACK2_OVERLAPPINGPARTITIONER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Partitioner.hpp"
#include "Tpetra_RowGraph.hpp"

/* \brief Ifpack2_OverlappingPartitioner: A class to create overlapping
    partitions of a local graph.

Class Ifpack2_OverlappingPartitioner enables the extension of 
non-overlapping partitions to an arbitrary value of overlap.
Note that overlap refers to the overlap among \e local parts,
and not the overlap among the processes.

Supported parameters are:
- \c "partitioner: local parts": the required number of parts;
- \c "partitioner: overlap": the required amount of overlap is set in 
  parameter. Default = 0 (integer).
- \c "partitioner: verbose": if \c true, information are reported on
  cout. Nothing is reported otherwise.

This class is a semi-virtual class, that contains the basic utilities
for derived classes Ifpack2_LinearPartitioner, Ifpack2_GreedyPartitioner,
Ifpack2_METISPartitioner, and Ifpack2_EquationPartitioner. Graphs in
input to one of these classes are supposed to contain no singletons.
Usually, this means that the graph is derived from an Tpetra_RowMatrix,
that has been filtered using Ifpack2_SingletonFilter.

\date Last update: Aug-12.
*/  

namespace Ifpack2 {

template<class GraphType>
class OverlappingPartitioner : public Partitioner<GraphType> {

public:
  typedef typename GraphType::local_ordinal_type LocalOrdinal;
  typedef typename GraphType::global_ordinal_type GlobalOrdinal;
  typedef typename GraphType::node_type Node;

  //! Constructor.
  OverlappingPartitioner(const Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> >& Graph);

  //! Destructor.
  virtual ~OverlappingPartitioner();

  //! Returns the number of computed local partitions.
  LocalOrdinal numLocalParts() const;

  //! Returns the overlapping level.
  size_t overlappingLevel() const;

  //! Returns the local non-overlapping partition ID of the specified row.
  /*! Returns the non-overlapping partition ID of the specified row.
   \param 
   MyRow - (In) local row numbe

   \return
   Local ID of non-overlapping partition for \t MyRow.
   */
  LocalOrdinal operator() (LocalOrdinal MyRow) const;

  //! Returns the local overlapping partition ID of the j-th node in partition i.
  LocalOrdinal operator() (LocalOrdinal i, LocalOrdinal j) const;

  //! Returns the number of rows contained in specified partition.
  size_t numRowsInPart(const LocalOrdinal Part) const;

  //! Copies into List the rows in the (overlapping) partition Part.
  void rowsInPart(const LocalOrdinal Part,  Teuchos::ArrayRCP<LocalOrdinal> & List) const;

  //! Returns an ArrayRCP to the integer vector containing the non-overlapping partition ID of each local row.
  virtual Teuchos::ArrayView<const LocalOrdinal>  nonOverlappingPartition() const;

  //! Sets all the parameters for the partitioner.
  /*! The supported parameters are:
   * - \c "partitioner: overlap" (int, default = 0).
   * - \c "partitioner: local parts" (int, default = 1).
   * - \c "partitioner: print level" (int, default = 0).
   */
  virtual void setParameters(Teuchos::ParameterList& List);

  //! Sets all the parameters for the partitioner.
  /*! This function is used by derived classes to set their own
   * parameters. These classes should not derive SetParameters(),
   * so that common parameters can be set just once.
   */
  virtual void setPartitionParameters(Teuchos::ParameterList& List) = 0;

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


  //! @name Overridden from Teuchos::Describable 
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object. 
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}
 
protected:

  //! Number of local subgraphs.  This is a localordinal because negatives mean something.
  LocalOrdinal NumLocalParts_;
  //! Partition_[i] contains the ID of non-overlapping part it belongs to
  Teuchos::Array<LocalOrdinal> Partition_; 
  //! Parts_[i][j] is the ID of the j-th row contained in the (overlapping) 
  // partition i
  Teuchos::Array<Teuchos::ArrayRCP<LocalOrdinal> > Parts_;
  //! Reference to the graph to be partitioned
  Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> > Graph_;
  //! Overlapping level.
  size_t OverlappingLevel_;
  //! If \c true,  the graph has been successfully partitioned.
  bool IsComputed_;
  //! If \c true, information are reported on cout.
  bool verbose_;

}; // class Ifpack2::OverlappingPartitioner

}// namespace Ifpack2

#endif // IFPACK2_OVERLAPPINGPARTITIONER_DECL_HPP
