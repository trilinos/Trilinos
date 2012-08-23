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
