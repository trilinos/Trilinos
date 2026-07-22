#ifndef KRINO_KRINO_KRINO_LIB_AKRI_DISTRIBUTEDVECTOR_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_DISTRIBUTEDVECTOR_HPP_
#include <cstddef>
#include <initializer_list>
#include <vector>
#include <stk_util/parallel/Parallel.hpp>

namespace krino {

class DistributedVector
{
private:
    std::vector<double> myData;
    size_t myLocalSize;
    stk::ParallelMachine myComm;

public:
    DistributedVector() : myLocalSize(0), myComm(stk::parallel_machine_null()) {}

    // Constructor for purely local data with size=localSize, defaults to 0.
    DistributedVector(const size_t size, const double initialValue = 0.0)
        : myData(size, initialValue), myLocalSize(size), myComm(stk::parallel_machine_null()) {}

    // Constructor for purely local data from initializer_list
    DistributedVector(std::initializer_list<double> initList)
        : myData(initList), myLocalSize(initList.size()), myComm(stk::parallel_machine_null()) {}

    // Constructor for mix of local and nonlocal data, defaults to 0.
    DistributedVector(const stk::ParallelMachine comm, const size_t size, const size_t localSize, const double initialValue = 0.0 )
        : myData(size, initialValue), myLocalSize(localSize), myComm(comm) {}

    DistributedVector(const DistributedVector& other)
        : myData(other.myData), myLocalSize(other.myLocalSize), myComm(other.myComm) {}

    DistributedVector& operator=(const DistributedVector& other)
    {
        if (this != &other)
        {
            myData = other.myData;
            myLocalSize = other.myLocalSize;
            myComm = other.myComm;
        }
        return *this;
    }

    DistributedVector& operator=(std::initializer_list<double> initList)
    {
      myData = initList;
      myLocalSize = initList.size();
      myComm = stk::parallel_machine_null();
      return *this;
    }

    double& operator[](size_t index) { return myData[index]; }
    const double& operator[](size_t index) const { return myData[index]; }

    double * data() { return myData.data(); }
    const double * data() const { return myData.data(); }

    std::vector<double>::iterator begin() { return myData.begin(); }
    std::vector<double>::const_iterator begin() const { return myData.begin(); }
    std::vector<double>::iterator end() { return myData.end(); }
    std::vector<double>::const_iterator end() const { return myData.end(); }

    size_t size() const  { return myData.size(); }
    size_t local_size() const  { return myLocalSize; }
    std::pair<size_t,size_t> sizes() const { return std::make_pair(myData.size(), myLocalSize); }

    void assign(const stk::ParallelMachine comm, const size_t size, const size_t localSize, const double value);
    void resize(const size_t newSize);
    void resize(const std::pair<size_t,size_t> newSizes);

    std::vector<double> & get() { return myData; }
    const std::vector<double> & get() const { return myData; }

    stk::ParallelMachine comm() const { return myComm; }
};

DistributedVector xpby(const DistributedVector & x, const double b, const DistributedVector & y);

DistributedVector scalar_times_vector(const double a, const DistributedVector & x);

double Dot(const DistributedVector & x, const DistributedVector & y);

DistributedVector vectorSubtract(const DistributedVector& x, const DistributedVector& y);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_DISTRIBUTEDVECTOR_HPP_ */
