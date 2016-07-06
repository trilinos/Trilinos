/*
 * ReducedDependencyTransferBase.hpp
 *
 *  Created on: Apr 28, 2016
 *      Author: tcfishe
 */

#ifndef STK_STK_UNIT_TESTS_STK_TRANSFER_REDUCEDDEPENDENCYTRANSFERBASE_HPP_
#define STK_STK_UNIT_TESTS_STK_TRANSFER_REDUCEDDEPENDENCYTRANSFERBASE_HPP_

namespace stk {
namespace transfer {

class ReducedDependecyTransferBase {
public :
  ReducedDependecyTransferBase(){};
  virtual ~ReducedDependecyTransferBase(){};
  void initialize() {
    coarse_search();
    communication();
    local_search();
  }
  virtual void coarse_search() = 0;
  virtual void communication() = 0;
  virtual void communicate_destination_points() = 0;
  virtual void local_search()  = 0;
  virtual void apply()         = 0;
};
}
}




#endif /* STK_STK_UNIT_TESTS_STK_TRANSFER_REDUCEDDEPENDENCYTRANSFERBASE_HPP_ */
