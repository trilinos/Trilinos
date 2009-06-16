#ifndef KOKKOS_SERIALNODE_HPP_
#define KOKKOS_SERIALNODE_HPP_

#include <Kokkos_StandardNodeMemoryModel.hpp>

class SerialNode : public StandardMemoryModel {
  public:
    template <class WDP>
    void execute1D(unsigned int length, WDP wd) {
      wd(0,length);
    }

    template <class WDP>
    void reduce1D(unsigned int length, WDP &wd) {
      for (unsigned int i=0; i<length; ++i) {
        wd.result = wd.reduce( wd.result, wd.generate(i) );
      }
    }

};

#endif
