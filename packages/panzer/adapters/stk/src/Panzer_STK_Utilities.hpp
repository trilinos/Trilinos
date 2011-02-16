#ifndef __Panzer_STK_Utilities_hpp__
#define __Panzer_STK_Utilities_hpp__

#include "Panzer_STK_Interface.hpp"

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

namespace panzer {
   template <typename LO,typename GO> class DOFManager;
}

namespace panzer_stk { 

void write_solution_data(const panzer::DOFManager<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_MultiVector & x);
void write_solution_data(const panzer::DOFManager<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_Vector & x);

void read_solution_data(const panzer::DOFManager<int,int> & dofMngr,const panzer_stk::STK_Interface & mesh,Epetra_MultiVector & x);
void read_solution_data(const panzer::DOFManager<int,int> & dofMngr,const panzer_stk::STK_Interface & mesh,Epetra_Vector & x);

/** Using a container, computed the sorted permutation vector
  * do not modifiy the original container.
  *
  * Motivated by this board on StackOverflow: 
  * http://stackoverflow.com/questions/4523220/sorting-a-vector-of-double-precision-reals-and-obtain-their-order
  */ 
template <typename RAContainer,class Compare>
void sorted_permutation(const RAContainer & cont,std::vector<std::size_t> & permutation,const Compare & comp);

/** Using a container, computed the sorted permutation vector
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
