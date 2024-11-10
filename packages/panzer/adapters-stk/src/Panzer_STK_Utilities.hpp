// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef PANZER_HAVE_EPETRA_STACK

#ifndef __Panzer_STK_Utilities_hpp__
#define __Panzer_STK_Utilities_hpp__

#include "Panzer_STK_Interface.hpp"

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

namespace panzer {
  class GlobalIndexer;
}

namespace panzer_stk {

/** Write a vector to the cell data of a STK mesh. This will look up
  * the cell field <code>prefix+fieldName+postfix</code>, which is assumed
  * to be in the STK mesh object. If not an assertion exeption will be thrown.
  *
  * \param[in] mesh STK mesh object
  * \param[in] data Vector of doubles equatl to the total number of elements on this processor
  * \param[in] fieldName Name of field to be written (must be a STK field)
  */
void write_cell_data(panzer_stk::STK_Interface & mesh,const std::vector<double> & data,const std::string & fieldName);

void write_solution_data(const panzer::GlobalIndexer& dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_MultiVector & x,const std::string & prefx="",const std::string & postfix="");
void write_solution_data(const panzer::GlobalIndexer& dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_Vector & x,const std::string & prefix="",const std::string & postfix="");

/** Using a container, compute the sorted permutation vector
  * do not modifiy the original container.
  *
  * Motivated by this board on StackOverflow:
  * http://stackoverflow.com/questions/4523220/sorting-a-vector-of-double-precision-reals-and-obtain-their-order
  */
template <typename RAContainer,class Compare>
void sorted_permutation(const RAContainer & cont,std::vector<std::size_t> & permutation,const Compare & comp);

/** Using a container, compute the sorted permutation vector
  * do not modifiy the original container.
  *
  * Motivated by this board on StackOverflow:
  * http://stackoverflow.com/questions/4523220/sorting-a-vector-of-double-precision-reals-and-obtain-their-order
  */
template <typename RAContainer>
void sorted_permutation(const RAContainer & cont,std::vector<std::size_t> & permutation);

}

namespace panzer_stk {
// utility class used by the sorted permutation objects
template <typename RAContainer,typename Compare>
struct PermFunctor {
   PermFunctor(const RAContainer & cont,const Compare & comp)
      : compare(comp), values(cont) {}
   PermFunctor(const PermFunctor & p)
      : compare(p.compare), values(p.values) {}

   bool operator()(std::size_t a,std::size_t b) const
   { return compare(values[a],values[b]); }

private:
   const Compare & compare;
   const RAContainer & values;

   PermFunctor();
};

template <typename RAContainer>
void sorted_permutation(const RAContainer & cont,std::vector<std::size_t> & permutation)
{
   std::less<typename RAContainer::value_type> comp;
   sorted_permutation(cont,permutation,comp);
}

template <typename RAContainer,class Compare>
void sorted_permutation(const RAContainer & cont,std::vector<std::size_t> & permutation,const Compare & comp)
{
   PermFunctor<RAContainer,Compare> pf(cont,comp);

   permutation.resize(cont.size());
   for(std::size_t i=0;i<cont.size();i++)
      permutation[i] = i;

   std::sort(permutation.begin(),permutation.end(),pf);
}

}

#endif

#endif // PANZER_HAVE_EPETRA_STACK