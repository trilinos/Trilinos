//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
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


template <class Node=Kokkos::DefaultNode::DefaultNodeType>
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

