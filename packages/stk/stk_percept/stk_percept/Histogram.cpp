#include <stk_percept/Histogram.hpp>

namespace stk {
  namespace percept {

    template<typename T>
    void example_histogram()
    {
       Histogram<T> h;
       T data[] = {1,2,3,5,8,13};
       h.assign(data,data+6);
       h.compute_ranges();
       h.compute_uniform_buckets(2);
       h.print(std::cout);
       h.compute_uniform_buckets(2, true); // use log scale
       h.print(std::cout);
       T ranges[] = {1,4,9,13};
       std::vector<T> vranges(ranges,ranges+4);
       h.set_ranges(vranges);
       h.print(std::cout);
    }

    template void example_histogram<double>();
    template<> void example_histogram<int>();
  }
}
