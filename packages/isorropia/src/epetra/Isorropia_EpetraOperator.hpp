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

//   Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph, int base);

//   Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, int base);

//   Operator(Teuchos::RCP<const Epetra_MultiVector> input_coords, int base);

//   Operator(Teuchos::RCP<const Epetra_BlockMap> input_map, int base);

//   Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph, 
// 	   Teuchos::RCP<const Epetra_MultiVector> input_coords, int base);

//   Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, 
// 	   Teuchos::RCP<const Epetra_MultiVector> input_coords, int base);


  /** Constructor that accepts an Epetra_CrsGraph object

     \param input_graph Matrix-graph object for which a new operation
        is to be computed. 

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 
    
     \param[in] base index base

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */

  Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
              const Teuchos::ParameterList& paramlist, int base);


  /** Constructor that accepts an Epetra_BlockMap object

     \param[in] input_map  BlockMap object for which a new operation
        is to be computed. 

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

     \param[in] base
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */

  Operator(Teuchos::RCP<const Epetra_BlockMap> input_map,
              const Teuchos::ParameterList& paramlist, int base);



  /** Constructor that accepts an Epetra_CrsGraph object
      and an Epetra_MultiVector object of coordinates

     \param input_graph Matrix-graph object for which a new operation
        is to be computed. 

     \param[in] input_coords The input geometric coordinates (represented in a
                    multivector)

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

     \param[in] base

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */

  Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	   Teuchos::RCP<const Epetra_MultiVector> input_coords,
           const Teuchos::ParameterList& paramlist, int base);


  /** Constructor that accepts an Epetra_CrsGraph object, a CostDescriber,
      an Epetra_MultiVector object of coordinates, and an 
      Epetra_MultiVector object of coordinate weights.

     \param input_graph Matrix-graph object for which a new operation
        is to be computed. 

     \param costs CostDescriber object which allows for user-specified
       weights 

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

     \param[in] base

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator (Teuchos::RCP<const Epetra_CrsGraph> input_graph,
            Teuchos::RCP<CostDescriber> costs,
              const Teuchos::ParameterList& paramlist, int base);


  /** Constructor that accepts an Epetra_CrsGraph object, a CostDescriber,
      an Epetra_MultiVector object of coordinates, and an 
      Epetra_MultiVector object of coordinate weights.

     \param input_graph Matrix-graph object for which a new operation
        is to be computed. 

     \param costs CostDescriber object which allows for user-specified
       weights 

     \param coords The input geometric coordinates (represented in a
                    multivector)

     \param weights A one or more dimensional weight for each of the
                    The input geometric coordinates

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

     \param[in] base

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator (Teuchos::RCP<const Epetra_CrsGraph> input_graph,
            Teuchos::RCP<CostDescriber> costs,
	    Teuchos::RCP<const Epetra_MultiVector> coords,
	    Teuchos::RCP<const Epetra_MultiVector> weights,
              const Teuchos::ParameterList& paramlist, int base);


  /**
     Constructor that accepts an Epetra_RowMatrix object.

     \param input_matrix Matrix object for which a new operation is
        to be computed. 

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

     \param[in] base Index base

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	   const Teuchos::ParameterList& paramlist, int base);


  /**
     Constructor that accepts an Epetra_RowMatrix object and an 
     Epetra_MultiVector object of coordinates.

     \param input_matrix Matrix object for which a new operation is
        to be computed. 

     \param coords The input geometric coordinates (represented in a
                    multivector)

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

     \param[in] base

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	   Teuchos::RCP<const Epetra_MultiVector> coords,
	   const Teuchos::ParameterList& paramlist, int base);


  /**
     Constructor that accepts an Epetra_RowMatrix object and a
     CostDescriber.

     \param input_matrix Matrix object for which a new operation is
        to be computed. 

     \param costs CostDescriber object which allows for user-specified
       weights 

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

     \param[in] base

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	   Teuchos::RCP<CostDescriber> costs,
	   const Teuchos::ParameterList& paramlist, int base);


  /**
     Constructor that accepts an Epetra_RowMatrix object, a
     CostDescriber, an Epetra_MultiVector object of coordinates, 
     and an Epetra_MultiVector object of coordinate weights.

     \param input_matrix Matrix object for which a new operation is
        to be computed. 

     \param costs CostDescriber object which allows for user-specified
       weights 

     \param coords The input geometric coordinates (represented in a                
                    multivector)

     \param weights A one or more dimensional weight for each of the                
                    The input geometric coordinates

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

     \param[in] base

  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the operation. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  */
  Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	   Teuchos::RCP<CostDescriber> costs,
	   Teuchos::RCP<const Epetra_MultiVector> coords,
	   Teuchos::RCP<const Epetra_MultiVector> weights,
	   const Teuchos::ParameterList& paramlist, int base);


  /**
     Constructor that accepts an Epetra_MultiVector object and a
     ParameterList

     \param coords The input geometric coordinates (represented in a
                    multivector)

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 

     \param[in] base index base

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

     \param[in] paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. 
     
     \param[in] base index base
 
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

  /** Get the cost object
   */

  Teuchos::RCP<Isorropia::Epetra::CostDescriber> & getCosts() { return costs_; }

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

  void paramsToUpper(Teuchos::ParameterList &, int &changed, bool rmUnderscore=true);
  void stringToUpper(std::string &s, int &changed, bool rmUnderscore=false);
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


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

