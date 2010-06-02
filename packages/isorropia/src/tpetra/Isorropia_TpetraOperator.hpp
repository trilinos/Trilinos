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

#ifndef _Isorropia_TpetraOperator_hpp_
#define _Isorropia_TpetraOperator_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_Operator.hpp>

#ifdef HAVE_ISORROPIA_TPETRA
#include <Isorropia_TpetraLibrary.hpp>
#include <Isorropia_TpetraCostDescriber.hpp>
#include <Tpetra_CrsGraph_decl.hpp>

namespace Isorropia {

namespace Tpetra {


/** An implementation of the Partitioner interface that operates on
    Epetra matrices and linear systems.

*/

template <class Node=Kokkos::DefaultNode::DefaultNodeType>
class Operator : virtual public Isorropia::Operator 
{
public:

  /** Constructor that accepts an ::Tpetra::CrsGraph object

     \param input_graph Matrix-graph object for which a new operation
        is to be computed. 

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 
    
     \param[in] base index base

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */

  Operator(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph,
              const Teuchos::ParameterList& paramlist, int base);



  /** Destructor */
  virtual ~Operator();

  /** setParameters() is an internal method which handles
      the parameters from a Teuchos::ParameterList object. 
   */
  void setParameters(const Teuchos::ParameterList& paramlist);

  virtual void compute(bool force_compute) = 0 ;

  /** Query whether compute_operation() has already been called.
   */
  bool alreadyComputed() const 
  {
    return operation_already_computed_;
  }

  int numProperties() const {
    return (numberOfProperties_);
  }

  int numLocalProperties() const {
    return (localNumberOfProperties_);
  }

  /** Return the new partition ID for a given element that
     resided locally in the old operation.
  */
  virtual const int& operator[](int myElem) const;

  /** Return the number of elements in a given partition.
  */
  virtual int numElemsWithProperty(int property) const;

  /** Fill user-allocated list (of length len) with the
      global element ids to be located in the given partition.
  */
  virtual void elemsWithProperty(int property,
			 int* elementList,
			 int len) const;

  virtual int extractPropertiesCopy(int len,
				    int& size,
				    int* array) const ;

  virtual int extractPropertiesView(int& size,
				    const int*& array) const;

private:

  void paramsToUpper(Teuchos::ParameterList &, int &changed, bool rmUnderscore=true);
  void stringToUpper(std::string &s, int &changed, bool rmUnderscore=false);
  int numberOfProperties_;
  int localNumberOfProperties_;
  std::vector<int> numberElemsByProperties_;

protected:
  Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph_;
  Teuchos::RCP<const ::Tpetra::Map<int,int,Node> > input_map_;

  Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix_;
  Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords_;
  Teuchos::RCP<Isorropia::Tpetra::CostDescriber<Node> > costs_;
  Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > weights_;


  Teuchos::ParameterList paramlist_;

  int exportsSize_;
  std::vector<int> imports_;
  std::vector<int> properties_;

  bool operation_already_computed_;

  int global_num_vertex_weights_;
  int global_num_graph_edge_weights_;
  int global_num_hg_edge_weights_;

  Teuchos::RCP<Library<Node> > lib_;

  int base_;
  
  void computeNumberOfProperties();
};//class Operator

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

