// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_STK_Utilities_hpp__
#define __Panzer_STK_Utilities_hpp__

#include "Panzer_STK_Interface.hpp"

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

namespace panzer {
   template <typename LO,typename GO> class DOFManagerFEI;
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

void write_solution_data(const panzer::DOFManagerFEI<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_MultiVector & x,const std::string & prefx="",const std::string & postfix="");
void write_solution_data(const panzer::DOFManagerFEI<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_Vector & x,const std::string & prefix="",const std::string & postfix="");

void read_solution_data(const panzer::DOFManagerFEI<int,int> & dofMngr,const panzer_stk::STK_Interface & mesh,Epetra_MultiVector & x);
void read_solution_data(const panzer::DOFManagerFEI<int,int> & dofMngr,const panzer_stk::STK_Interface & mesh,Epetra_Vector & x);

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
