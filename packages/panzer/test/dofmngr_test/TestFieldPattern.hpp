#ifndef __TestFieldPattern_hpp__
#define __TestFieldPattern_hpp__

#include "dofmngr_v2/Panzer_FieldPattern.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

class TestFieldPattern : public FieldPattern {
public:
   TestFieldPattern() {}

   /* This function has no functionality in this case.
    * If called it will throw an assertion failure
    */
   virtual void getSubcellClosureIndices(int dim,int cellIndex,std::vector<int> & indices) const
   { TEUCHOS_ASSERT(false); }

   virtual int getSubcellCount(int dim) const
   {  return subcellIndices[dim].size(); }

   virtual const std::vector<int> & getSubcellIndices(int dim,int cellIndex) const
   {  return subcellIndices[dim][cellIndex]; }

   virtual int getDimension() const
   { return subcellIndices.size()-1; }

   std::vector<std::vector<int> > & operator[](int v)
   { return subcellIndices[v]; } 

public:
   std::vector<std::vector<std::vector<int> > > subcellIndices;
};

}

#endif
