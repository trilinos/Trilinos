// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_MultipleFieldFunction_hpp
#define percept_MultipleFieldFunction_hpp

#include <cmath>
#include <math.h>
#include <map>

#include <typeinfo>

#include <percept/function/MDArray.hpp>

#include <percept/function/FunctionOperator.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/Function.hpp>
#include <percept/function/internal/HasValue.hpp>
#include <percept/function/StringFunction.hpp>
#include <percept/function/ElementOp.hpp>
#include <percept/function/BucketOp.hpp>

#include <percept/PerceptMesh.hpp>

#include <percept/norm/IntrepidManager.hpp>

  namespace percept
  {

    /** This class compbines several Fields into a single function, e.g., density_field*temperature_field/pressure_field
     *  It uses Intrepid2 to compute basis functions then evaluates all the required fields and stores them in MDArray's 
     *  ready for the operator() to compute the actual function.
     */
    class MultipleFieldFunction : public FieldFunction
    {
    public:
//       MultipleFieldFunction(const char *name, stk::mesh::FieldBase *field, stk::mesh::BulkData *bulk, SearchType searchType = SIMPLE_SEARCH,
//                              Dimensions domain_dimensions = Dimensions(),
//                              Dimensions codomain_dimensions = Dimensions(),
//                              unsigned integration_order = 0) : FieldFunction(name, field, bulk, searchType, domain_dimensions, codomain_dimensions)
//       {
//       }
      //MultipleFieldFunction() : FieldFunction() {}

      // add_fields(std::vector<Field *> ffs) ...

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

#endif
