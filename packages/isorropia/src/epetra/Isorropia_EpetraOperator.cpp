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

#include <Isorropia_EpetraOperator.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#else /* HAVE_EPETRA */
#error "This module needs Epetra"
#endif /* HAVE_EPETRA */



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

#ifdef HAVE_EPETRA

namespace Epetra {

// Operator::
// Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph, int base)
//   : input_graph_(input_graph),
//     input_matrix_(0),
//     input_coords_(0),
//     costs_(0),
//     weights_(0),
//     operation_already_computed_(false),
//     lib_(0),
//     base_(base)
// {
//   input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
// }

// Operator::
// Operator(Teuchos::RCP<const Epetra_BlockMap> input_map, int base)
//   : input_map_(input_map),
//     input_graph_(0),
//     input_matrix_(0),
//     input_coords_(0),
//     costs_(0),
//     weights_(0),
//     operation_already_computed_(false),
//     lib_(0),
//     base_(base)
// {

// }

// Operator::
// Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph, 
// 	 Teuchos::RCP<const Epetra_MultiVector> input_coords,int base)
//   : input_graph_(input_graph),
//     input_matrix_(0),
//     input_coords_(input_coords),
//     costs_(0),
//     weights_(0),
//     operation_already_computed_(false),
//     lib_(0),
//     base_(base)
// {
//   input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
// }


Operator::
Operator(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	 const Teuchos::ParameterList& paramlist, int base)
  : input_graph_(input_graph), input_matrix_(0),
    input_coords_(0),
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


void
Operator::computeNumberOfProperties()
{
  std::vector<int>::const_iterator elemsIter;
  std::vector<int>::iterator numberIter;
  const Epetra_Comm& input_comm = input_map_->Comm();

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

  input_comm.MaxAll(&max, &numberOfProperties_, 1);

  numberOfProperties_ = numberOfProperties_ - base_ + 1;

  localNumberOfProperties_ = max - base_ + 1;
}

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



} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

