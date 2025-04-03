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

#include <EEP_Isorropia_TpetraCostDescriber.hpp>
#include <Isorropia_Operator.hpp>
#include <EEP_Isorropia_TpetraLibrary.hpp>

#include <EEP_Isorropia_Tpetra.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_CrsGraph_decl.hpp>

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

#ifdef MIN
#undef MIN
#endif

#define MIN(a,b) ((a)<(b)?(a):(b))

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
  : input_graph_(input_graph),
    costs_(0),
    operation_already_computed_(false),
    lib_(0), base_(base) // EEP__
  {
    input_map_ = input_graph->getRowMap(); // Teuchos::rcp(&(input_graph->getRowMap()), false); // EEP___

    if (paramlist.name() != "EmptyParameterList") {
      setParameters(paramlist);
    }
  }

  /** Destructor */
  virtual ~Operator();

  /** setParameters() is an internal method which handles
      the parameters from a Teuchos::ParameterList object. 
   */
  void setParameters(const Teuchos::ParameterList& paramlist);

  /** Get the cost object
   */

  Teuchos::RCP< Isorropia::Tpetra::CostDescriber< LocalOrdinal, GlobalOrdinal, Node> > & getCosts() { return costs_; }

  virtual void compute(bool force_compute) = 0 ;

  /** Query whether compute_operation() has already been called.
   */
  bool alreadyComputed() const { // EEP__
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
  virtual const int& operator[](int myElem) const; // EEP__ forced by virtual = 0 in Isorropia_Operator

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
  Teuchos::RCP< const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > input_map_;
  Teuchos::RCP< const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > input_graph_;
  Teuchos::RCP< Isorropia::Tpetra::CostDescriber<LocalOrdinal, GlobalOrdinal, Node> > costs_;

  Teuchos::ParameterList paramlist_;

  int exportsSize_;
  std::vector<int> imports_;
  std::vector<int> properties_;

  bool operation_already_computed_;

  int global_num_vertex_weights_;
  int global_num_graph_edge_weights_;
  int global_num_hg_edge_weights_;

  Teuchos::RCP< Library<LocalOrdinal, GlobalOrdinal, Node> > lib_; // EEP__

  int base_;

  void computeNumberOfProperties();
}; // class Operator

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Operator<LocalOrdinal, GlobalOrdinal, Node>::~Operator()
{
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Operator<LocalOrdinal, GlobalOrdinal, Node>::setParameters(const Teuchos::ParameterList& paramlist)
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

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
const int& Operator<LocalOrdinal, GlobalOrdinal, Node>::operator[](int myElem) const // EEP__ forced by virtual = 0 in Isorropia_Operator
{
  return (properties_[myElem]);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int Operator<LocalOrdinal, GlobalOrdinal, Node>::numElemsWithProperty(int property) const
{
  if ((unsigned) property < numberElemsByProperties_.size())
    return numberElemsByProperties_[property];
  return (0);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Operator<LocalOrdinal, GlobalOrdinal, Node>::elemsWithProperty(int property, int* elementList, int len) const
{
  int length = 0;
  std::vector<int>::const_iterator elemsIter;

  for (elemsIter = properties_.begin() ; (length < len) && (elemsIter != properties_.end()) ;
       elemsIter ++) {
    if (*elemsIter == property)
      elementList[length++] = elemsIter - properties_.begin();
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Operator<LocalOrdinal, GlobalOrdinal, Node>::computeNumberOfProperties()
{
  std::cout << "EEP Entering Operator<>::computeNumberOfProperties()" << std::endl;

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

  max = numberOfProperties_; //input_map_->getComm()->maxAll(&max, &numberOfProperties_, 1); // EEP____: generalize to MPI np > 1
  std::cout << "EEP In Operator<>::computeNumberOfProperties()"
            << ": max = " << max
            << std::endl;
  
  numberOfProperties_ = numberOfProperties_ - base_ + 1;

  localNumberOfProperties_ = max - base_ + 1;

  std::cout << "EEP Leaving Operator<>::computeNumberOfProperties()"
            << ", base_ = " << base_
            << ", numberOfProperties_ = " << numberOfProperties_
            << ", localNumberOfProperties_ = " << localNumberOfProperties_
	    << std::endl;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Operator<LocalOrdinal, GlobalOrdinal, Node>::stringToUpper(std::string &s, int &changed, bool rmUnderscore)
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

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
Operator<LocalOrdinal, GlobalOrdinal, Node>::paramsToUpper(Teuchos::ParameterList &plist, int &changed, bool rmUnderscore)
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

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int
Operator<LocalOrdinal, GlobalOrdinal, Node>::extractPropertiesCopy(int len,
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

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int
Operator<LocalOrdinal, GlobalOrdinal, Node>::extractPropertiesView(int& size,
				    const int*& array) const
{
  size = properties_.size();
  if (size)
    array = &(properties_[0]);
  else
    array = NULL;
  return (0);
}

}//namespace Epetra
}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

