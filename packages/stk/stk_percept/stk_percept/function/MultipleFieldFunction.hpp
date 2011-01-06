#ifndef stk_percept_MultipleFieldFunction_hpp
#define stk_percept_MultipleFieldFunction_hpp

#include <cmath>
#include <math.h>
#include <map>

#include <typeinfo>

#include <stk_percept/function/MDArray.hpp>

#include <stk_percept/function/FunctionOperator.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/Function.hpp>
#include <stk_percept/function/internal/HasValue.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/ElementOp.hpp>
#include <stk_percept/function/BucketOp.hpp>

#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/norm/IntrepidManager.hpp>

namespace stk
{
  namespace percept
  {

    /** This class compbines several Fields into a single function, e.g., density_field*temperature_field/pressure_field
     *  It uses Intrepid to compute basis functions then evaluates all the required fields and stores them in MDArray's 
     *  ready for the operator() to compute the actual function.
     */
    class MultipleFieldFunction : public FieldFunction
    {
    public:
//       MultipleFieldFunction(const char *name, mesh::FieldBase *field, mesh::BulkData *bulk, SearchType searchType = SIMPLE_SEARCH,
//                              Dimensions domain_dimensions = Dimensions(),
//                              Dimensions codomain_dimensions = Dimensions(),
//                              unsigned integration_order = 0) : FieldFunction(name, field, bulk, searchType, domain_dimensions, codomain_dimensions)
//       {
//       }
      //MultipleFieldFunction() : FieldFunction() {}

      // addFields(std::vector<Field *> ffs) ...

#if 0
      computeAllFieldValues()
      {
        // loop over field functions
        // evaluate each to get the field values and store
      }

      operator() (in, out) = 0;

      // derived class can then use whatever field values it needs
      Derived::operator() (in, out)
      {
        // invoke helper:
        computeAllFieldValues();

        
        
      }

      // for l2 norm, can compose two functions
      // for StringFunction, can use these helpers also in evalFunctions
      
#endif      

    };

    

  }
}

#endif
