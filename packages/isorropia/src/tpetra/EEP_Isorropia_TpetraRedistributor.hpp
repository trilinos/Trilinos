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
  std::cout << "EEP Entering Redistributor<>::constructor()" << std::endl;
  if (!partitioner_->alreadyComputed()) {
    std::cout << "In Entering Redistributor<>::constructor(), pos 001" << std::endl;
    partitioner_->partition();
  }
  std::cout << "EEP Leaving Redistributor<>::constructor()" << std::endl;
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
  std::cout << "EEP Entering Redistributor<>::redistribute(2)" << std::endl;
  Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>outputGraphPtr(nullptr);
  redistribute(input_graph, outputGraphPtr, callFillComplete);

  std::cout << "EEP Leaving Redistributor<>::redistribute(2)" << std::endl;
  return Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>(outputGraphPtr);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Redistributor<LocalOrdinal, GlobalOrdinal, Node>::redistribute(const Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> input_graph, Teuchos::RCP<::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> outputGraphPtr, bool callFillComplete) // EEP____ check lots of changes
{
  std::cout << "EEP Entering Redistributor<>::redistribute(3)" << std::endl;
  create_importer( input_graph->getRowMap() ); // EEP__

  // First obtain the length of each of my new rows

  int myOldRows = input_graph->getLocalNumRows(); // NumMyRows();
  int myNewRows = target_map_->getLocalNumElements(); // NumMyElements();

  std::cout << "EEP In Redistributor<>::redistribute(3), pos 000"
            << ": myOldRows = " << myOldRows
            << ", myNewRows = " << myNewRows
            << std::endl;
  
  double *nnz = new double [myOldRows];
  for (int i=0; i < myOldRows; i++){
    nnz[i] = input_graph->getNumEntriesInLocalRow(i); // input_graph.NumMyIndices(i); // EEP__
  }

  std::cout << "EEP In Redistributor<>::redistribute(3), pos 001"
            << std::endl;

  ::Teuchos::ArrayView<double> tmpArray1(nnz, myOldRows);
  ::Tpetra::Vector<double, LocalOrdinal, GlobalOrdinal, Node> oldRowSizes(input_graph->getRowMap(), tmpArray1);

  std::cout << "EEP In Redistributor<>::redistribute(3), pos 002"
            << std::endl;

  //if (myOldRows)
  //  delete [] nnz;

  ::Tpetra::Vector<double, LocalOrdinal, GlobalOrdinal, Node> newRowSizes(target_map_);

  newRowSizes.doImport(oldRowSizes, *importer_, ::Tpetra::INSERT);
  auto tmpView = newRowSizes.getLocalViewHost(::Tpetra::Access::ReadOnly);

  std::cout << "EEP In Redistributor<>::redistribute(3), pos 003"
            << std::endl;

  size_t *rowSize = new size_t [myNewRows];
  for (int i=0; i< myNewRows; i++){
    rowSize[i] = static_cast<size_t>(tmpView.data()[i]);
  }

  std::cout << "EEP In Redistributor<>::redistribute(3), pos 004"
            << ": rowSize =";
  for (int i(0); i < myNewRows; ++i) {
    std::cout << " " << rowSize[i];
  }
  std::cout << std::endl;
  
  // Receive new rows, send old rows

  ::Teuchos::ArrayView<size_t> tmpArray2(rowSize, myNewRows);
  std::cout << "EEP In Redistributor<>::redistribute(3), pos 004.2"
            << ": tmpArray2 =";
  for (int i(0); i < myNewRows; ++i) {
    std::cout << " " << tmpArray2[i];
  }
  std::cout << std::endl;
  std::cout << "EEP In Redistributor<>::redistribute(3), pos 004.3"
            << ": *target_map_ = " << *target_map_
            << std::endl;
  outputGraphPtr = Teuchos::rcp( new ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>(target_map_, tmpArray2) ); // , true);

  //if (myNewRows)
  //  delete [] rowSize;

  std::cout << "EEP In Redistributor<>::redistribute(3), pos 005"
	    << ": *outputGraphPtr = " << *outputGraphPtr
            << std::endl;

  std::cout << "EEP In Redistributor<>::redistribute(3), pos 005.2"
	    << ": outputGraphPtr->getRangeMap() = " << outputGraphPtr->getRangeMap()
            << std::endl;

  outputGraphPtr->doImport(*input_graph, *importer_, ::Tpetra::INSERT);
  std::cout << "EEP In Redistributor<>::redistribute(3), pos 005.3"
	    << ": outputGraphPtr->getRangeMap() = " << outputGraphPtr->getRangeMap()
            << std::endl;
  outputGraphPtr->fillComplete(); // EEP____ check
  std::cout << "EEP In Redistributor<>::redistribute(3), pos 005.4"
	    << ": outputGraphPtr->getRangeMap() = " << outputGraphPtr->getRangeMap()
            << std::endl;
  outputGraphPtr->resumeFill(); // EEP____ check
  
  std::cout << "EEP In Redistributor<>::redistribute(3), pos 006"
	    << ": outputGraphPtr->getRangeMap() = " << outputGraphPtr->getRangeMap()
	    << std::endl;

  // Set the new domain map such that
  // (a) if old DomainMap == old RangeMap, preserve this property,
  // (b) otherwise, let the new DomainMap be the old DomainMap 
  Teuchos::RCP<const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> newDomainMap;
  if ( input_graph->getDomainMap()->isSameAs( *(input_graph->getRangeMap()) ) ) {
    std::cout << "EEP In Redistributor<>::redistribute(3), pos 006.a" << std::endl;
    newDomainMap = outputGraphPtr->getRangeMap();
  }
  else {
    std::cout << "EEP In Redistributor<>::redistribute(3), pos 006.b" << std::endl;
    newDomainMap = input_graph->getDomainMap();
  }

  std::cout << "EEP In Redistributor<>::redistribute(3), pos 007"
            << std::endl;

  if (callFillComplete) {
    std::cout << "EEP In Redistributor<>::redistribute(3), pos 008" << std::endl;
    if (outputGraphPtr->isFillComplete() == false) {
      std::cout << "EEP In Redistributor<>::redistribute(3), pos 009"
	        << ": newDomainMap = " << newDomainMap
	        << ", target_map_ = " << target_map_
		<< std::endl;
      std::cout << "EEP In Redistributor<>::redistribute(3), information on newDomainMap"
                << ": isOneToOne() = " << newDomainMap->isOneToOne()
                << ", GlobalNumElements() = " << newDomainMap->getGlobalNumElements()
                << ", LocalNumElements() = " << newDomainMap->getLocalNumElements()
                << ", IndexBase() = " << newDomainMap->getIndexBase()
                << ", MinLocalIndex() = " << newDomainMap->getMinLocalIndex()
                << ", MaxLocalIndex() = " << newDomainMap->getMaxLocalIndex()
                << ", MinGlobalIndex() = " << newDomainMap->getMinGlobalIndex()
                << ", MaxGlobalIndex() = " << newDomainMap->getMaxGlobalIndex()
                << ", MinAllGlobalIndex() = " << newDomainMap->getMinAllGlobalIndex()
                << ", MaxAllGlobalIndex() = " << newDomainMap->getMaxAllGlobalIndex()
                << std::endl;
      std::cout << "EEP In Redistributor<>::redistribute(3), information on target_map_"
                << ": isOneToOne() = " << target_map_->isOneToOne()
                << ", GlobalNumElements() = " << target_map_->getGlobalNumElements()
                << ", LocalNumElements() = " << target_map_->getLocalNumElements()
                << ", IndexBase() = " << target_map_->getIndexBase()
                << ", MinLocalIndex() = " << target_map_->getMinLocalIndex()
                << ", MaxLocalIndex() = " << target_map_->getMaxLocalIndex()
                << ", MinGlobalIndex() = " << target_map_->getMinGlobalIndex()
                << ", MaxGlobalIndex() = " << target_map_->getMaxGlobalIndex()
                << ", MinAllGlobalIndex() = " << target_map_->getMinAllGlobalIndex()
                << ", MaxAllGlobalIndex() = " << target_map_->getMaxAllGlobalIndex()
                << std::endl;
      outputGraphPtr->fillComplete(newDomainMap, target_map_); // EEP____ check
    }
  }

  std::cout << "EEP Leaving Redistributor<>::redistribute(3)" << std::endl;
  return;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Redistributor<LocalOrdinal, GlobalOrdinal, Node>::create_importer(const Teuchos::RCP<const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> src_map)
{
  std::cout << "EEP Entering Redistributor<>::create_importer()" << std::endl;

  if (!Teuchos::is_null(partitioner_) && partitioner_->numProperties() >
                                src_map->getComm()->getSize()) {
    throw std::runtime_error/*Isorropia::Exception*/("Cannot redistribute: Too many parts for too few processors.");
  }

  if (Teuchos::is_null(target_map_) && !Teuchos::is_null(partitioner_)) {
      std::cout << "EEP In Redistributor<>::create_importer(), pos 000" << std::endl;
      target_map_ = partitioner_->createNewMap();
  }

  importer_ = Teuchos::rcp(new ::Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>(target_map_, src_map));

  std::cout << "EEP Leaving Redistributor<>::create_importer()" << std::endl;
}
  
}//namespace Tpetra
}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

