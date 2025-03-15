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

#ifndef _Isorropia_TpetraOperator_hpp_
#define _Isorropia_TpetraOperator_hpp_

//#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

//#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_Operator.hpp>
//#include <Isorropia_EpetraLibrary.hpp>

//#include <Isorropia_EpetraOperator.hpp>
//#include <Isorropia_Exception.hpp>
#include <EEP_Isorropia_Tpetra.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

//#include <Epetra_Comm.h>
//#include <Epetra_Map.h>
//#include <Epetra_Import.h>
//#include <Epetra_Vector.h>
//#include <Epetra_MultiVector.h>
#include <Tpetra_CrsGraph_decl.hpp>
//#include <Epetra_CrsMatrix.h>
//#include <Epetra_LinearProblem.h>

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

//#ifdef MIN
//#undef MIN
//#endif

//#define MIN(a,b) ((a)<(b)?(a):(b))

namespace Isorropia {

namespace Tpetra {
  //class CostDescriber;

/** An implementation of the Partitioner interface that operates on
    Epetra matrices and linear systems.

*/
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
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

  Operator(Teuchos::RCP< const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > input_graph, // EEP__
           const Teuchos::ParameterList& paramlist, int base)
  : input_graph_(input_graph), /*input_matrix_(0),
    input_coords_(0),
    costs_(0), weights_(0),*/
    operation_already_computed_(false),
    /*lib_(0),*/ base_(base) // EEP__
  {
    input_map_ = Teuchos::rcp(&(input_graph->getRowMap()), false);

    if (paramlist.name() != "EmptyParameterList") {
      setParameters(paramlist);
    }
  }

#if 0 // EEP
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
#endif // EEP
  /** Query whether compute_operation() has already been called.
   */
  bool alreadyComputed() const { // EEP__
    return operation_already_computed_;
  }
#if 0 // EEP
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
#endif // EEP
private:

  //void paramsToUpper(Teuchos::ParameterList &, int &changed, bool rmUnderscore=true);
  //void stringToUpper(std::string &s, int &changed, bool rmUnderscore=false);
  int numberOfProperties_;
  int localNumberOfProperties_;
  std::vector<int> numberElemsByProperties_;

protected:
  Teuchos::RCP< const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > input_map_;
  Teuchos::RCP< const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > input_graph_;
  //Teuchos::RCP< const Tpetra::RowMatrix<> > input_matrix_;
  //Teuchos::RCP< const Tpetra::MultiVector<> > input_coords_;
  //Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs_;
  //Teuchos::RCP<const Tpetra::MultiVector> weights_;

  Teuchos::ParameterList paramlist_;

  int exportsSize_;
  std::vector<int> imports_;
  std::vector<int> properties_;

  bool operation_already_computed_;

  int global_num_vertex_weights_;
  int global_num_graph_edge_weights_;
  int global_num_hg_edge_weights_;

  //Teuchos::RCP<Library> lib_; // EEP__

  int base_;

  void computeNumberOfProperties();
};//class Operator

#if 0 // EEP
Operator::
Operator(Teuchos::RCP<const Epetra_BlockMap> input_map,
	 const Teuchos::ParameterList& paramlist, int base)
  : input_map_(input_map),input_graph_(0), input_matrix_(0),
    input_coords_(0),
    costs_(0), weights_(0),
    operation_already_computed_(false),
    lib_(0), base_(base)
{
  if(paramlist.name() != "EmptyParameterList")
  {
    setParameters(paramlist);
  }
}

Operator::
Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	 Teuchos::RCP<const Epetra_MultiVector> input_coords,
	 const Teuchos::ParameterList& paramlist, int base)
    : input_graph_(input_graph), input_matrix_(0),
      input_coords_(input_coords),
      costs_(0), weights_(0),
      operation_already_computed_(false),
      lib_(0), base_(base)
{
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);

  if(paramlist.name() != "EmptyParameterList")
  {
     setParameters(paramlist);
  }

}



Operator::
Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	 Teuchos::RCP<CostDescriber> costs,
	 const Teuchos::ParameterList& paramlist, int base)
  : input_graph_(input_graph), input_matrix_(0),
    input_coords_(0),
    costs_(costs), weights_(0),
    operation_already_computed_(false),
    lib_(0), base_(base)
{

  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);

  if(paramlist.name() != "EmptyParameterList")
  {
    setParameters(paramlist);
  }
}

Operator::
Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	 Teuchos::RCP<CostDescriber> costs,
	 Teuchos::RCP<const Epetra_MultiVector> input_coords,
	 Teuchos::RCP<const Epetra_MultiVector> weights,
	 const Teuchos::ParameterList& paramlist, int base)
  : input_graph_(input_graph), input_matrix_(0),
    input_coords_(input_coords),
    costs_(costs), weights_(weights),
    operation_already_computed_(false),
    lib_(0), base_(base)
{

  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);

  if(paramlist.name() != "EmptyParameterList")
  {
     setParameters(paramlist);
  }
}

// Operator::
// Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, int base)
//   : input_graph_(0),
//     input_matrix_(input_matrix),
//     input_coords_(0),
//     costs_(0),
//     weights_(0),
//     operation_already_computed_(false),
//     lib_(0),
//     base_(base)
// {
//   input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()), false);
// }

// Operator::
// Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, 
// 	 Teuchos::RCP<const Epetra_MultiVector> input_coords, int base)
//   : input_graph_(0),
//     input_matrix_(input_matrix),
//     input_coords_(input_coords),
//     costs_(0),
//     weights_(0),
//     operation_already_computed_(false),
//     lib_(0),
//     base_(base)
// {
//   input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()), false);
// }

Operator::
Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	 const Teuchos::ParameterList& paramlist, int base)
  : input_graph_(0), input_matrix_(input_matrix),
    input_coords_(0),
    costs_(0), weights_(0),
    operation_already_computed_(false),
    lib_(0),
    base_(base)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);

  if(paramlist.name() != "EmptyParameterList")
  {
     setParameters(paramlist);
  }
}

Operator::
Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	 Teuchos::RCP<const Epetra_MultiVector> input_coords,
	 const Teuchos::ParameterList& paramlist, int base)
  : input_graph_(0), input_matrix_(input_matrix),
    input_coords_(input_coords),
    costs_(0), weights_(0),
    operation_already_computed_(false),
    lib_(0),
    base_(base)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);

  if(paramlist.name() != "EmptyParameterList")
  {
     setParameters(paramlist);
  }
}

Operator::
Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	 Teuchos::RCP<CostDescriber> costs,
	 const Teuchos::ParameterList& paramlist, int base)
  : input_graph_(0),
    input_matrix_(input_matrix),
    input_coords_(0),
    costs_(costs),
    weights_(0),
    operation_already_computed_(false),
    lib_(0),
    base_(base)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);

  if(paramlist.name() != "EmptyParameterList")
  {
     setParameters(paramlist);
  }
}

Operator::
Operator(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	 Teuchos::RCP<CostDescriber> costs,
         Teuchos::RCP<const Epetra_MultiVector> input_coords,
         Teuchos::RCP<const Epetra_MultiVector> weights,
	 const Teuchos::ParameterList& paramlist, int base)
  : input_graph_(0),
    input_matrix_(input_matrix),
    input_coords_(input_coords),
    costs_(costs),
    weights_(weights),
    operation_already_computed_(false),
    lib_(0),
    base_(base)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);

  if(paramlist.name() != "EmptyParameterList")
  {
     setParameters(paramlist);
  }
}

// Operator::
// Operator(Teuchos::RCP<const Epetra_MultiVector> input_coords, int base)
//   : input_graph_(0),
//     input_matrix_(0),
//     input_coords_(input_coords),
//     costs_(0),
//     weights_(0),
//     operation_already_computed_(false),
//     lib_(0),
//     base_(base)
// {
//   input_map_ = Teuchos::rcp(&(input_coords->Map()), false);
// }

Operator::
Operator(Teuchos::RCP<const Epetra_MultiVector> input_coords,
	 const Teuchos::ParameterList& paramlist, int base)
  : input_graph_(0),
    input_matrix_(0),
    input_coords_(input_coords),
    costs_(0),
    weights_(0),
    operation_already_computed_(false),
    lib_(0),
    base_(base)
{
  input_map_ = Teuchos::rcp(&(input_coords->Map()),false);

  if(paramlist.name() != "EmptyParameterList")
  {
     setParameters(paramlist);
  }
}

Operator::
Operator(Teuchos::RCP<const Epetra_MultiVector> input_coords,
         Teuchos::RCP<const Epetra_MultiVector> weights,
	 const Teuchos::ParameterList& paramlist, int base)
  : input_graph_(0),
    input_matrix_(0),
    input_coords_(input_coords),
    costs_(0),
    weights_(weights),
    operation_already_computed_(false),
    lib_(0),
    base_(base)
{
  input_map_ = Teuchos::rcp(&(input_coords->Map()),false);

  if(paramlist.name() != "EmptyParameterList")
  {
     setParameters(paramlist);
  }
}

Operator::~Operator()
{
}

void Operator::setParameters(const Teuchos::ParameterList& paramlist)
{
  int changed;
  paramlist_ = paramlist;
  paramsToUpper(paramlist_, changed);

  Teuchos::ParameterList &sublist = paramlist_.sublist("ZOLTAN");
  std::string str_sym="STRUCTURALLY SYMMETRIC";
  sublist.set("GRAPH_SYMMETRIZE", "TRANSPOSE");
  if (paramlist_.isParameter(str_sym)
      && paramlist_.get<std::string>(str_sym) == "YES") {
    sublist.set("GRAPH_SYMMETRIZE", "NONE");
  }

}

const int& Operator::operator[](int myElem) const
{
  return (properties_[myElem]);
}

int Operator::numElemsWithProperty(int property) const
{
  if ((unsigned) property < numberElemsByProperties_.size())
    return numberElemsByProperties_[property];
  return (0);
}

void
Operator::elemsWithProperty(int property, int* elementList, int len) const
{
  int length = 0;
  std::vector<int>::const_iterator elemsIter;

  for (elemsIter = properties_.begin() ; (length < len) && (elemsIter != properties_.end()) ;
       elemsIter ++) {
    if (*elemsIter == property)
      elementList[length++] = elemsIter - properties_.begin();
  }
}
#endif // EEP

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Operator<LocalOrdinal, GlobalOrdinal, Node>::computeNumberOfProperties()
{
  std::vector<int>::const_iterator elemsIter;
  std::vector<int>::iterator numberIter;

  int max = 0;
  numberElemsByProperties_.assign(properties_.size(), 0);

  numberIter = numberElemsByProperties_.begin();
  for(elemsIter = properties_.begin() ; elemsIter != properties_.end() ; elemsIter ++) {
    int property;
    property = *elemsIter;
    if (max < property) {
      max = property;
      int toAdd = max - numberElemsByProperties_.size() + 1;
      if (toAdd > 0) {
	numberElemsByProperties_.insert(numberElemsByProperties_.end(), toAdd, 0);
	numberIter = numberElemsByProperties_.begin();
      }
    }
    (*(numberIter + property)) ++;
  }

  input_map_->Comm().MaxAll(&max, &numberOfProperties_, 1);

  numberOfProperties_ = numberOfProperties_ - base_ + 1;

  localNumberOfProperties_ = max - base_ + 1;
}

#if 0 // EEP
void Operator::stringToUpper(std::string &s, int &changed, bool rmUnderscore)
{
  std::string::iterator siter;
  changed = 0;

  for (siter = s.begin(); siter != s.end() ; siter++)
  {
    if (islower(*siter)){
      *siter = toupper(*siter);
      changed++;
    }
    if (rmUnderscore && (*siter == '_')) {
      *siter = ' ';
      changed ++;
    }
  }
}

void Operator::paramsToUpper(Teuchos::ParameterList &plist, int &changed, bool rmUnderscore)
{
  changed = 0;

  // get a list of all parameter names in the list

  std::vector<std::string> paramNames ;
  Teuchos::ParameterList::ConstIterator pIter;

  pIter = plist.begin();

  while (1){
    //////////////////////////////////////////////////////////////////////
    // Compiler considered this while statement an error
    // for ( pIter = plist.begin() ; pIter != plist.end() ; pIter++ ){
    // }
    //////////////////////////////////////////////////////////////////////
    if (pIter == plist.end()) break;
    const std::string & nm = plist.name(pIter);
    paramNames.push_back(nm);
    pIter++;
  }

  // Change parameter names and values to upper case

  for (unsigned int i=0; i < paramNames.size(); i++){

    std::string origName(paramNames[i]);
    int paramNameChanged = 0;
    stringToUpper(paramNames[i], paramNameChanged, rmUnderscore);

    if (plist.isSublist(origName)){
      Teuchos::ParameterList &sublist = plist.sublist(origName);

      int sublistChanged=0;
      paramsToUpper(sublist, sublistChanged, false);

      if (paramNameChanged){

        // this didn't work, so I need to remove the old sublist
        // and create a new one
        //
        //sublist.setName(paramNames[i]);

        Teuchos::ParameterList newlist(sublist);
        plist.remove(origName);
        plist.set(paramNames[i], newlist);
      }
    }
    else if (plist.isType<std::string>(origName)){

      std::string paramVal(plist.get<std::string>(origName));

      int paramValChanged=0;
      stringToUpper(paramVal, paramValChanged);

      if (paramNameChanged || paramValChanged){
        if (paramNameChanged){
          plist.remove(origName);
        }
        plist.set(paramNames[i], paramVal);
        changed++;
      }
    }
  } // next parameter or sublist
}

int Operator::extractPropertiesCopy(int len,
				    int& size,
				    int* array) const
{
  int ierr;
  const int *ptr;

  ierr = extractPropertiesView(size, ptr);
  if (ierr)
    return (ierr);

  size = MIN(size, len);
  memcpy (array, ptr, size * sizeof(int));

  return (0);
}

int Operator::extractPropertiesView(int& size,
				    const int*& array) const
{
  size = properties_.size();
  if (size)
    array = &(properties_[0]);
  else
    array = NULL;
  return (0);
}
#endif // EEP

}//namespace Epetra
}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

