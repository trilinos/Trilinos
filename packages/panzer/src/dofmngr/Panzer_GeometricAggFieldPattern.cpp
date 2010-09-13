#include "Panzer_GeometricAggFieldPattern.hpp"

using Teuchos::RCP;

namespace panzer {

GeometricAggFieldPattern::GeometricAggFieldPattern()
   : patternBuilt_(false), dimension_(0)
{}

GeometricAggFieldPattern::GeometricAggFieldPattern(std::vector<Teuchos::RCP<const FieldPattern> > & patterns) 
   : patternBuilt_(false), dimension_(0)
{ 
   buildPattern(patterns); 
}

GeometricAggFieldPattern::GeometricAggFieldPattern(const Teuchos::RCP<const FieldPattern> & pattern) 
   : patternBuilt_(false), dimension_(0)
{ 
   buildPattern(pattern); 
}

void GeometricAggFieldPattern::buildPattern(const std::vector<Teuchos::RCP<const FieldPattern> > & patterns)
{
   std::size_t numPat = patterns.size();

   // must be at least one field to do something
   if(numPat<1) {
      std::cout << "No patterns!" << std::endl;
      return;
   }

   bool sameGeometry=true;
   for(std::size_t i=1;i<patterns.size();i++)
      sameGeometry &= patterns[0]->sameGeometry(*patterns[i]);       
   TEST_FOR_EXCEPTION(not sameGeometry,std::logic_error,
             "GeometricAggFieldPattern::buildPattern(): Patterns must "
             "have the same geometry!");

   // grab the dimension
   dimension_ = patterns[0]->getDimension();
   patternData_.resize(dimension_+1);

   // build space for subcells
   std::vector<int> subcellCount(dimension_+1);
   for(std::size_t d=0;d<dimension_+1;d++) {
      subcellCount[d] = patterns[0]->getSubcellCount(d);
      patternData_[d].resize(subcellCount[d]);
   }

   // build geometric pattern: increment it logically
   // over all the subcells.
   int counter = 0;
   for(std::size_t d=0;d<dimension_+1;d++) {
      for(int s=0;s<subcellCount[d];s++) {
         std::vector<int> & current = patternData_[d][s];
         for(std::size_t p=0;p<patterns.size();p++) {
            RCP<const FieldPattern> field = patterns[p];
            std::size_t num = field->getSubcellIndices(d,s).size();
 
            if(current.size()<num) { 
               for(int i=num-current.size();i>0;i--,counter++) 
                  current.push_back(counter);
            }
         } 
      }
   }

   // record that the pattern has been built
   patternBuilt_ = true;
}

void GeometricAggFieldPattern::buildPattern(const Teuchos::RCP<const FieldPattern> & pattern)
{
   std::vector<Teuchos::RCP<const FieldPattern> > patterns;
   patterns.push_back(pattern);
   buildPattern(patterns);
}

int GeometricAggFieldPattern::getSubcellCount(int dim) const
{
   if(patternBuilt_) return patternData_[dim].size();

   TEST_FOR_EXCEPTION(true,std::logic_error,
            "GeometricAggFieldPattern::getSubcellCount() cannot be called before "
            "GeometricAggFieldPattern::buildPattern()");
}

const std::vector<int> & GeometricAggFieldPattern::getSubcellIndices(int dim,int cellIndex) const
{
   if(patternBuilt_) return patternData_[dim][cellIndex];

   TEST_FOR_EXCEPTION(true,std::logic_error,
            "GeometricAggFieldPattern::getSubcellIndices() cannot be called before "
            "GeometricAggFieldPattern::buildPattern()");
}

int GeometricAggFieldPattern::getDimension() const
{
   if(patternBuilt_) return dimension_;

   TEST_FOR_EXCEPTION(true,std::logic_error,
            "GeometricAggFieldPattern::getDimension() cannot be called before "
            "GeometricAggFieldPattern::buildPattern()");
}

}
