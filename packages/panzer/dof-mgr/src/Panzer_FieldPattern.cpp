// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_FieldPattern.hpp"

namespace panzer {

FieldPattern::~FieldPattern() {}

int FieldPattern::numberIds() const
{
   int count = 0;
   int dim = getDimension();

   // compute number of IDs
   for(int i=0;i<dim+1;i++) {
      for(int sc=0;sc<getSubcellCount(i);sc++)
         count += getSubcellIndices(i,sc).size();
   }

   return count;
}

void FieldPattern::print(std::ostream & os) const
{
   int dim = getDimension()+1;
   os << "FieldPattern: " << dim << " Subcell types" << std::endl;
   for(int i=0;i<dim;i++) {
      int subcells = getSubcellCount(i);
      os << "FieldPattern:    " << subcells << " subcells of type " << i << std::endl;
       
      for(int j=0;j<subcells;j++) {
         const std::vector<int> & indices = getSubcellIndices(i,j);
         os << "FieldPattern:       subcell " << j << " = [ ";
         for(std::size_t k=0;k<indices.size();k++)
            os << indices[k] << " "; 
         os << "]" << std::endl;
      }
   }
}

bool FieldPattern::sameGeometry(const FieldPattern & fp) const
{
   bool equal = true;

   // test same dimension
   std::size_t dim = getDimension();
   equal &= (dim==(std::size_t) fp.getDimension());
   
   // check sub cells
   for(std::size_t d=0;d<dim;d++)
      equal &= getSubcellCount(d)==fp.getSubcellCount(d);

   return equal;
}

bool FieldPattern::consistentSubcells() const
{
   bool consistent = true;

   std::size_t dim = getDimension();
   for(std::size_t d=0;d<dim+1;d++) {
      int numSC = getSubcellCount(d);
      std::size_t sz = getSubcellIndices(d,0).size();
      for(int i=1;i<numSC;i++) {
         consistent &= (sz==getSubcellIndices(d,i).size());
      }
   }

   return consistent;
}

bool FieldPattern::equals(const FieldPattern & fp) const
{
   // same geometry is required
   if(not this->sameGeometry(fp))
      return false;
   
   // check to make sure subcell indices are equal
   int dimension = this->getDimension(); 
   for(int d=0;d<dimension+1;d++) {
      for(int sc=0;sc<this->getSubcellCount(d);sc++) {
         const std::vector<int> & myVector = this->getSubcellIndices(d,sc);
         const std::vector<int> & argVector = fp.getSubcellIndices(d,sc);

         // check size of vectors
         if(myVector.size()!=argVector.size()) 
            return false;
         
         // check content of vectors
         bool eq = std::equal(myVector.begin(),myVector.end(),argVector.begin());
         if(not eq) 
            return false;
      }
   }

  return true;
}

std::ostream & operator<<(std::ostream & os,const FieldPattern & fp)
{
   fp.print(os);
   return os;
}

}
