#pragma once
#ifndef __COO_HPP__
#define __COO_HPP__

namespace Example { 
  
  using namespace std;
  
  template<typename ValueType,
           typename OrdinalType>
  class Coo {
  public:
    OrdinalType _i,_j;
    ValueType _val;
    
    Coo() {}
    Coo(const Coo& b)
      : _i(b._i),
        _j(b._j),
        _val(b._val) {}

    Coo<ValueType,OrdinalType>& operator=(const Coo<ValueType,OrdinalType> &y) {
      this->_i = y._i;
      this->_j = y._j;
      this->_val = y._val;

      return *this;
    }

    bool operator<(const Coo<ValueType,OrdinalType> &y) const {
      OrdinalType r_val = (this->_i - y._i);
      return (r_val == 0 ? this->_j < y._j : r_val < 0);
    }  
    
    bool operator==(const Coo<ValueType,OrdinalType> &y) const {
      return !(bool)((this->_i - y._i) + (this->_j - y._j));
    }  
    
    bool operator!=(const Coo<ValueType,OrdinalType> &y) const {
      return !(*this == y);
    }  
  };
  
}
#endif
