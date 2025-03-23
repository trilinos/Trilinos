//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

#ifndef _Isorropia_TpetraRedistributor_hpp_
#define _Isorropia_TpetraRedistributor_hpp_

#include <Isorropia_Redistributor.hpp>
//#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <EEP_Isorropia_Tpetra.hpp>
#include <EEP_Isorropia_TpetraPartitioner.hpp>

#ifdef USE_UTILS
#include <../../utils/ispatest_epetra_utils.hpp>
#include <../../utils/ispatest_lbeval_utils.hpp>
using namespace ispatest;
#endif

#include <Teuchos_RCP.hpp>

#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Import_decl.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_CrsGraph_decl.hpp>

namespace Isorropia {

namespace Tpetra {

/** @ingroup partitioning_grp partitioning_rcp_grp partitioning_ptr_grp
     Class which is constructed with a Partitioner instance, and
     provides several methods for redistributing Epetra objects
     given the partitioning computed by the Partitioner object.
*/

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Redistributor : public Isorropia::Redistributor {
public:

  /** @ingroup partitioning_rcp_grp
      This constructor calls the Isorropia::Epetra::Partitioner::partition
      method on the @c partitioner if it has not already been called.
 
      \param partitioner (in) this input partitioner determines the new partitioning
            to be created when Isorropia::Epetra::Redistributor::redistribute is called
   */
  Redistributor(Teuchos::RCP< Isorropia::Tpetra::Partitioner<LocalOrdinal, GlobalOrdinal, Node> > partitioner);

  /** 
       Destructor
   */
  virtual ~Redistributor();

  /** @ingroup partitioning_rcp_grp
      Method to accept a Epetra_CrsGraph object, and
      return a redistributed Epetra_CrsGraph object.

      \param input_graph (in) the graph for which we want a new graph that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param callFillComplete (in) The new graph is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_graph has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed graph 
  */
  Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>
    redistribute(const Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> input_graph, bool callFillComplete= true);

  /** @ingroup partitioning_ptr_grp
      Method to accept a Epetra_CrsGraph object, and
      return a redistributed Epetra_CrsGraph object.

      \param input_graph (in) the graph for which we want a new graph that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputGraphPtr (out) pointer to the new redistributed graph 

      \param callFillComplete (in) The new graph is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_graph has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

  */
  void redistribute(const Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> input_graph, Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> outputGraphPtr, bool callFillComplete= true);

private:
  /** @ingroup partitioning_grp
      Create an importer object to be used in the redistribution
      \param src_map (in) the map describing the pattern of the import operation
   */
  void create_importer(const Teuchos::RCP<const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> src_map); // Epetra_BlockMap EEP

  Teuchos::RCP< Isorropia::Tpetra::Partitioner<LocalOrdinal, GlobalOrdinal, Node> > partitioner_;
  Teuchos::RCP< ::Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > importer_;
  Teuchos::RCP< ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > target_map_;

}; //class Redistributor

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Redistributor<LocalOrdinal, GlobalOrdinal, Node>::Redistributor(Teuchos::RCP<Isorropia::Tpetra::Partitioner<LocalOrdinal, GlobalOrdinal, Node>> partitioner)
  : partitioner_(partitioner),
  importer_(),
  target_map_()
{
  if (!partitioner_->alreadyComputed()) {
    partitioner_->partition();
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Redistributor<LocalOrdinal, GlobalOrdinal, Node>::~Redistributor()
{
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>
Redistributor<LocalOrdinal, GlobalOrdinal, Node>::redistribute(const Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> input_graph, bool callFillComplete)
{
  Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>outputGraphPtr(nullptr);
  redistribute(input_graph, outputGraphPtr, callFillComplete);

  return Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>(outputGraphPtr);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Redistributor<LocalOrdinal, GlobalOrdinal, Node>::redistribute(const Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> input_graph, Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> outputGraphPtr, bool callFillComplete)
{
  create_importer( input_graph->getRowMap() ); // EEP___
#if 0 // EEP___
  // First obtain the length of each of my new rows

  int myOldRows = input_graph.getLocalNumRows(); // NumMyRows();
  int myNewRows = target_map_->getLocalNumEntries(); // NumMyElements();

  double *nnz = new double [myOldRows];
  for (int i=0; i < myOldRows; i++){
    nnz[i] = 0; // input_graph.NumMyIndices(i); // EEP___
  }

  ::Tpetra::Vector<LocalOrdinal, GlobalOrdinal, Node, double> oldRowSizes(Copy, input_graph.RowMap(), nnz);

  if (myOldRows)
    delete [] nnz;

  ::Tpetra::Vector<LocalOrdinal, GlobalOrdinal, Node, double> newRowSizes(*target_map_);

  newRowSizes.Import(oldRowSizes, *importer_, Insert);

  int *rowSize = new int [myNewRows];
  for (int i=0; i< myNewRows; i++){
    rowSize[i] = static_cast<int>(newRowSizes[i]);
  }

  // Receive new rows, send old rows

  outputGraphPtr = new ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(Copy, *target_map_, rowSize, true);

  if (myNewRows)
    delete [] rowSize;

  outputGraphPtr->Import(input_graph, *importer_, Insert);

  // Set the new domain map such that
  // (a) if old DomainMap == old RangeMap, preserve this property,
  // (b) otherwise, let the new DomainMap be the old DomainMap 
  const Epetra_BlockMap *newDomainMap;
  if (input_graph.DomainMap().SameAs(input_graph.RangeMap()))
     newDomainMap = &(outputGraphPtr->RangeMap());
  else
     newDomainMap = &(input_graph.DomainMap());

  if (callFillComplete && (!outputGraphPtr->Filled()))
    outputGraphPtr->FillComplete(*newDomainMap, *target_map_);
#endif // EEP

  return;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Redistributor<LocalOrdinal, GlobalOrdinal, Node>::create_importer(const Teuchos::RCP<const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> src_map)
{

  if (!Teuchos::is_null(partitioner_) && partitioner_->numProperties() >
                                src_map->getComm()->getSize()) {
    throw std::runtime_error/*Isorropia::Exception*/("Cannot redistribute: Too many parts for too few processors.");
  }

  if (Teuchos::is_null(target_map_) && !Teuchos::is_null(partitioner_))
      target_map_ = partitioner_->createNewMap();

  importer_ = Teuchos::rcp(new ::Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>(target_map_, src_map));

}
  
}//namespace Tpetra
}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

