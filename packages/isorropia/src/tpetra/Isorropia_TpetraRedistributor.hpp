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
#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_ISORROPIA_TPETRA


#include <Kokkos_DefaultNode.hpp>
#include <Isorropia_TpetraPartitioner.hpp>

//#include <Tpetra_CrsGraph_decl.hpp>
//#include <Tpetra_RowMatrix.hpp>
//#include <Tpetra_MultiVector_decl.hpp>


namespace Isorropia {

namespace Tpetra {


/** @ingroup partitioning_grp partitioning_rcp_grp partitioning_ptr_grp
     Class which is constructed with a Partitioner instance, and
     provides several methods for redistributing Tpetra objects
     given the partitioning computed by the Partitioner object.
*/

template <typename Node = ::Tpetra::Map<int, int>::node_type >
class Redistributor : public Isorropia::Redistributor {
public:

  /** @ingroup partitioning_rcp_grp
      This constructor calls the Isorropia::Tpetra::Partitioner::partition
      method on the @c partitioner if it has not already been called.

      \param[in] partitioner this input partitioner determines the new partitioning
            to be created when Isorropia::Tpetra::Redistributor::redistribute is called
   */
  Redistributor(Teuchos::RCP<Isorropia::Tpetra::Partitioner<Node> > partitioner);

  /** @ingroup partitioning_ptr_grp
      This constructor calls the Isorropia::Tpetra::Partitioner::partition
      method on the @c partitioner if it has not already been called.

      \param[in] partitioner this input partitioner determines the new partitioning
            to be created when Isorropia::Tpetra::Redistributor::redistribute is called
   */
  Redistributor(Isorropia::Tpetra::Partitioner<Node> *partitioner);

  /**
       Destructor
   */
  virtual ~Redistributor();

  /** @ingroup partitioning_rcp_grp
      Method to accept a Tpetra::CrsGraph object, and
      return a redistributed Tpetra::CrsGraph object.

      \param[in] input_graph the graph for which we want a new graph that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param[in] callFillComplete The new graph is FillComplete'd if callFillComplete is @c true.
      In that case, the range map is set to equal the row map.
      The domain map will equal the range map, unless the
      input_graph has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed graph
  */

  Teuchos::RCP< ::Tpetra::CrsGraph<int,int,Node> >
  redistribute(const ::Tpetra::CrsGraph<int,int,Node>& input_graph, bool callFillComplete= true);

  /** @ingroup partitioning_ptr_grp
      Method to accept a Tpetra::CrsGraph object, and
      return a redistributed Tpetra::CrsGraph object.

      \param[in] input_graph the graph for which we want a new graph that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param[out] outputGraphPtr pointer to the new redistributed graph

      \param[in] callFillComplete The new graph is FillComplete'd if callFillComplete is @c true.
      In that case, the range map is set to equal the row map.
      The domain map will equal the range map, unless the
      input_graph has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

  */
  void redistribute(const ::Tpetra::CrsGraph<int,int,Node>& input_graph,
                    ::Tpetra::CrsGraph<int,int,Node> * &outputGraphPtr, bool callFillComplete= true);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Tpetra::CrsMatrix object, and
      return a redistributed Tpetra::CrsMatrix object.

      \param[in] input_matrix the matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param[in] callFillComplete The new matrix is FillComplete'd if callFillComplete is @c true.
      In that case, the range map is set to equal the row map.
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed matrix
  */
  Teuchos::RCP< ::Tpetra::CrsMatrix<double,int,int,Node> >
  redistribute(const ::Tpetra::CrsMatrix<double,int,int,Node>& input_matrix, bool callFillComplete= true);

  /** @ingroup partitioning_ptr_grp
      Method to accept a Tpetra::CrsMatrix object, and
      return a redistributed Tpetra::CrsMatrix object.

      \param[in] inputMatrix the matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param[out] outputMatrix pointer to the new redistributed matrix

      \param[in] callFillComplete The new matrix is FillComplete'd if callFillComplete is @c true.
      In that case, the range map is set to equal the row map.
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.
  */

  void redistribute(const ::Tpetra::CrsMatrix<double,int,int,Node>& inputMatrix,
                    ::Tpetra::CrsMatrix<double,int,int,Node> * &outputMatrix, bool callFillComplete= true);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Tpetra::Vector object, and
      return a redistributed Tpetra::Vector object.

      \param[in] input_vector the vector for which we want a new vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \return a reference counted pointer to the new redistributed vector
  */
  Teuchos::RCP< ::Tpetra::Vector<double,int,int,Node> >
  redistribute(const ::Tpetra::Vector<double,int,int,Node>& input_vector);

  /** @ingroup partitioning_ptr_grp
      Method to accept a Tpetra::Vector object, and
      return a redistributed Tpetra::Vector object.

      \param[in] inputVector the vector for which we want a new vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param[out] outputVector pointer to the new redistributed vector
  */
  void
  redistribute(const ::Tpetra::Vector<double,int,int,Node>& inputVector,
               ::Tpetra::Vector<double,int,int,Node> * &outputVector);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Tpetra::MultiVector object, and
      return a redistributed Tpetra::MultiVector object.

      \param[in] input_vector the multi vector for which we want a new multi vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \return a reference counted pointer to the new redistributed multi vector

  */
  Teuchos::RCP< ::Tpetra::MultiVector<double,int,int,Node> >
  redistribute(const ::Tpetra::MultiVector<double,int,int,Node>& input_vector);


  /** @ingroup partitioning_ptr_grp
      Method to accept a Tpetra::MultiVector object, and
      return a redistributed Tpetra::MultiVector object.

      \param[in] inputVector the multi vector for which we want a new multi vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param[out] outputVector a reference counted pointer to the new redistributed multi vector

  */
  void redistribute(const ::Tpetra::MultiVector<double,int,int,Node>& inputVector,
                          ::Tpetra::MultiVector<double,int,int,Node> * &outputVector);

  /** @ingroup partitioning_grp
      Reverse redistribute an Tpetra::Vector.

      \param[in] input_vector a vector that is distributed according to the partitioner that was used to create this Redistributor

      \param[out] output_vector a copy of the @c input_vector which has been redistributed according
                    to the reverse of the partitioner that was used to create this Redistributor

  */
  void
  redistribute_reverse(const ::Tpetra::Vector<double,int,int,Node> & input_vector, ::Tpetra::Vector<double,int,int,Node>& output_vector);

  /** @ingroup partitioning_grp
      Reverse redistribute an Tpetra::MultiVector.

      \param[in] input_vector a multi vector that is distributed according to the partitioner that was used to create this Redistributor

      \param[out] output_vector a copy of the @c input_vector which has been redistributed according
                    to the reverse of the partitioner that was used to create this Redistributor
  */
  void
  redistribute_reverse(const ::Tpetra::MultiVector<double,int,int,Node>& input_vector, ::Tpetra::MultiVector<double,int,int,Node>& output_vector);
private:
  /** @ingroup partitioning_grp
      Create an importer object to be used in the redistribution
      \param[in] src_map the map describing the pattern of the import operation
   */
  void create_importer(const ::Tpetra::Map<int,int,Node>& src_map);

  Teuchos::RCP<Isorropia::Tpetra::Partitioner<Node> > partitioner_;
  Teuchos::RCP< ::Tpetra::Import<int,int,Node> > importer_;
  Teuchos::RCP< ::Tpetra::Map<int,int,Node> > target_map_;

  bool created_importer_;

}; //class Redistributor

}//namespace Tpetra

}//namespace Isorropia

#endif //HAVE_ISORROPIA_TPETRA

#endif

