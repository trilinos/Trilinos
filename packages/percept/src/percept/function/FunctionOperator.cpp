// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cmath>
#include <math.h>
#include <map>
#include <vector>

#include <typeinfo>

#include <percept/function/MDArray.hpp>

#include <percept/function/FunctionOperator.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/Function.hpp>
#include <percept/function/CompositeFunction.hpp>
#include <percept/function/StringFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/function/ElementOp.hpp>
#include <percept/function/BucketOp.hpp>
#include <percept/function/MultipleFieldFunction.hpp>

#include <percept/function/internal/HasValue.hpp>
#include <percept/function/internal/IntegratedOp.hpp>

#include <percept/PerceptMesh.hpp>

//#include <percept/norm/IntrepidManager.hpp>

  namespace percept
  {


#if 0
    void FunctionOperator::integrate( Function& integrand, Function& result)
    {
      std::cout << "type= " << typeid(integrand).name() << " " << typeid(FieldFunction).name() << std::endl;
      if (typeid(integrand) == typeid(FieldFunction) ||
          typeid(integrand) == typeid(StringFunction))
        {
          const FieldFunction *field_function_const = dynamic_cast<const FieldFunction *>(&integrand);
          FieldFunction *field_function = const_cast<FieldFunction *>(field_function_const);
          PerceptMesh meshUtil(const_cast<stk::mesh::MetaData *>(&MetaData::get(field_function->get_bulk_data())), field_function->get_bulk_data());
#if 0
          NoOpScalar no_op(integrand, field_function->get_field());
          meshUtil.elementOpLoop(no_op , field_function->get_field());
#endif
        }
    }
#endif

#if 0
    void FunctionOperator::integrate(FunctionWithIntrepidRequest& integrand, Function& result)
    {
      const FieldFunction *field_function_const = dynamic_cast<const FieldFunction *>(&integrand);
      FieldFunction *field_function = const_cast<FieldFunction *>(field_function_const);
      PerceptMesh meshUtil(const_cast<stk::mesh::MetaData *>(&MetaData::get(field_function->get_bulk_data())), field_function->get_bulk_data());
      IntegratedOp general_integrand(integrand, TURBO_NONE);
      meshUtil.elementOpLoop(general_integrand , field_function->get_field());
    }
#endif


  }//namespace percept
