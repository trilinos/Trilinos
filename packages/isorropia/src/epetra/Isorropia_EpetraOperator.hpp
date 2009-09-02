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

#ifndef _Isorropia_EpetraOperator_hpp_
#define _Isorropia_EpetraOperator_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_Operator.hpp>
#include <Isorropia_EpetraLibrary.hpp>

#ifdef HAVE_EPETRA
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;

namespace Isorropia {

namespace Epetra {
  class CostDescriber;

/** An implementation of the Partitioner interface that operates on
    Epetra matrices and linear systems.

*/

class Operator : virtual public Isorropia::Operator {
public:

  /** Constructor that accepts an Epetra_CrsGraph object

     \param input_graph Matrix-graph object for which a new operation
        is to be computed. 
    */

  Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph, int base);

  /** Constructor that accepts an Epetra_CrsGraph object

     \param input_graph Matrix-graph object for which a new operation
        is to be computed. 

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */

  Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
              const Teuchos::ParameterList& paramlist, int base);

  /** Constructor that accepts an Epetra_CrsGraph object and a CostDescriber

     \param input_graph Matrix-graph object for which a new operation
        is to be computed. 

     \param costs CostDescriber object which allows for user-specified
       weights 

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator (Teuchos::RCP<const Epetra_CrsGraph> input_graph,
              Teuchos::RCP<CostDescriber> costs,
              const Teuchos::ParameterList& paramlist, int base);

  /**
     Constructor that accepts an Epetra_RowMatrix object.

     \param input_matrix Matrix object for which a new operation is
        to be computed. 

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	   const Teuchos::ParameterList& paramlist, int base);

  /**
     Constructor that accepts an Epetra_RowMatrix object and a
     CostDescriber.

     \param input_matrix Matrix object for which a new operation is
        to be computed. 

     \param costs CostDescriber object which allows for user-specified
       weights 

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	   Teuchos::RCP<CostDescriber> costs,
	   const Teuchos::ParameterList& paramlist, int base);

  /**
     Constructor that accepts an Epetra_MultiVector object and a
     ParameterList

     \param coords The input geometric coordinates (represented in a
                    multivector)

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator(Teuchos::RCP<const Epetra_MultiVector> coords,
	   const Teuchos::ParameterList& paramlist, int base);

  /**
     Constructor that accepts an Epetra_MultiVector object and a
     ParameterList

     \param coords The input geometric coordinates (represented in a
                    multivector)

     \param weights A one or more dimensional weight for each of the
                    The input geometric coordinates 

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator(Teuchos::RCP<const Epetra_MultiVector> coords,
           Teuchos::RCP<const Epetra_MultiVector> weights,
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
  bool alreadyComputed() const {
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

  void paramsToUpper(Teuchos::ParameterList &, int &changed);
  void stringToUpper(std::string &s, int &changed);
  int numberOfProperties_;
  int localNumberOfProperties_;
  std::vector<int> numberElemsByProperties_;

protected:
  Teuchos::RCP<const Epetra_BlockMap> input_map_;
  Teuchos::RCP<const Epetra_CrsGraph> input_graph_;
  Teuchos::RCP<const Epetra_RowMatrix> input_matrix_;
  Teuchos::RCP<const Epetra_MultiVector> input_coords_;
  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs_;
  Teuchos::RCP<const Epetra_MultiVector> weights_;

  Teuchos::ParameterList paramlist_;

  int exportsSize_;
  std::vector<int> imports_;
  std::vector<int> properties_;

  bool operation_already_computed_;

  int global_num_vertex_weights_;
  int global_num_graph_edge_weights_;
  int global_num_hg_edge_weights_;

  Teuchos::RCP<Library> lib_;

  int base_;

  void computeNumberOfProperties();
};//class Operator

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

