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

#ifndef _Isorropia_TpetraZoltanLib_hpp_
#define _Isorropia_TpetraZoltanLib_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>


#include <zoltan_cpp.h>

#ifdef HAVE_ISORROPIA_TPETRA

#include <Isorropia_TpetraCostDescriber.hpp>
#include <Isorropia_TpetraLibrary.hpp>
#include <Isorropia_TpetraQueryObject.hpp>

#include <Tpetra_CrsGraph_decl.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace Isorropia {

namespace Tpetra {

template <typename Node = ::Tpetra::Map<int, int>::node_type >
class ZoltanLibClass : public Library<Node> {
public:

  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph, int inputType=Library<Node>::unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph,
                 Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords, int inputType=Library<Node>::unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph,
                 Teuchos::RCP<CostDescriber<Node> > costs, int inputType=Library<Node>::unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph, Teuchos::RCP<CostDescriber<Node> > costs,
                 Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords, Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > weights,
                 int inputType=Library<Node>::unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix, int inputType=Library<Node>::unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix,
                 Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords, int inputType=Library<Node>::unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix,
                 Teuchos::RCP<CostDescriber<Node> > costs, int inputType=Library<Node>::unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix, Teuchos::RCP<CostDescriber<Node> > costs,
                 Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords, Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > weights,
                 int inputType=Library<Node>::unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords, int inputType=Library<Node>::unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords,
                 Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > weights, int inputType=Library<Node>::unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const ::Tpetra::Map<int,int,Node> > input_map, int inputType=Library<Node>::unspecified_input_);




  /** Method to partition the object that the ZoltanLibClass was contructed with.

      \param[in] paramlist  Parameters to govern partitioning.

      \param[out]  newPartitions The new partition for each of my objects, in
                                   local ID order.  The objects may be rows or
                               non-zeroes (for
                               CrsGraph and RowMatrix input) or coordinates (for
                               MultiVector input).  Partition numbers can range from
                               zero to numProcs-1.
      \param[out]  exportsSize  The number of my objects that will be exported to
                              another process under the new partitioning.  This is
                             also the number of elements in newPartitions that are
                             not equal to my process rank.
      \param[out]  imports   A list of the global IDs of the objects that will be
                            imported to my process under the new partitioning
   */

  virtual int
  repartition(Teuchos::ParameterList& paramlist,
              std::vector<int>& newPartitions,
              int& exportsSize,
              std::vector<int>& imports);

  /** Method to color the object that the ZoltanLibClass was contructed with.

      \param[in] paramlist  Parameters to govern coloring.

      \param[out]  colorAssignment A list of integers indicating the coloring of
                              the object, in local ID order.
  */
  virtual int
  color(Teuchos::ParameterList& paramlist,
        std::vector<int>& colorAssignment);

  /** Method to order the object that the ZoltanLibClass was contructed with.

      \param[in] paramlist  Parameters to govern ordering .

      \param[out]  orderAssignment A list of integers indicating the ordering of
                              the object, in local ID order.
  */
  virtual int
  order(Teuchos::ParameterList& paramlist,
        std::vector<int>& orderAssignment);

protected:
  virtual int precompute();
  virtual int postcompute();
  void computeCost();
  void preCheckPartition();

  void setParameterList(Teuchos::ParameterList& zoltanParamList);

private:
  Teuchos::ParameterList zoltanParamList_;
  std::string partMethod_; // stores partitioning method used, perhaps should be in TpetraLibrary?
  Zoltan *zz_;
  Teuchos::RCP<ZoltanLib::QueryObject<Node> > queryObject_;
  int num_obj_;

};//class ZoltanLibClass

}//namespace Tpetra
}//namespace Isorropia

#endif //HAVE_ISORROPIA_TPETRA

#endif

