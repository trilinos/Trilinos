#include <cmath>
#include <math.h>
#include <map>
#include <vector>

#include <typeinfo>

#include <stk_percept/function/MDArray.hpp>

#include <stk_percept/function/FunctionOperator.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/Function.hpp>
#include <stk_percept/function/CompositeFunction.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/function/ElementOp.hpp>
#include <stk_percept/function/BucketOp.hpp>
#include <stk_percept/function/MultipleFieldFunction.hpp>

#include <stk_percept/function/internal/HasValue.hpp>
#include <stk_percept/function/internal/IntegratedOp.hpp>

#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/norm/IntrepidManager.hpp>

namespace stk_classic
{
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
          PerceptMesh meshUtil(const_cast<mesh::fem::FEMMetaData *>(&MetaData::get(field_function->get_bulk_data())), field_function->get_bulk_data());
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
      PerceptMesh meshUtil(const_cast<mesh::fem::FEMMetaData *>(&MetaData::get(field_function->get_bulk_data())), field_function->get_bulk_data());
      IntegratedOp general_integrand(integrand, TURBO_NONE);
      meshUtil.elementOpLoop(general_integrand , field_function->get_field());
    }
#endif


  }//namespace percept
}//namespace stk_classic

