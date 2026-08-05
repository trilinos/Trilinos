#ifndef KRINO_KRINO_KRINO_LIB_AKRI_DISTRIBUTEDVECTOR_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_DISTRIBUTEDVECTOR_HPP_
#include <cstddef>
#include <initializer_list>
#include <vector>
#include <stk_util/parallel/Parallel.hpp>

namespace krino {

template<typename T>
class DistributedQuantities
{
private:
    std::vector<T> myData;
    size_t myLocalSize;
    stk::ParallelMachine myComm;

public:
    DistributedQuantities() : myLocalSize(0), myComm(stk::parallel_machine_null()) {}

    // Constructor for purely local data with size=localSize, default constructs.
    DistributedQuantities(const size_t size, const T initialValue = T())
        : myData(size, initialValue), myLocalSize(size), myComm(stk::parallel_machine_null()) {}

    // Constructor for purely local data from initializer_list
    DistributedQuantities(std::initializer_list<T> initList)
        : myData(initList), myLocalSize(initList.size()), myComm(stk::parallel_machine_null()) {}

    // Constructor for mix of local and nonlocal data, default constructs.
    DistributedQuantities(const stk::ParallelMachine comm, const size_t size, const size_t localSize, const T initialValue = T())
        : myData(size, initialValue), myLocalSize(localSize), myComm(comm) {}

    // Constructor for local and nonlocal vectors
    DistributedQuantities(const stk::ParallelMachine comm, const std::vector<T> & localVals, const std::vector<T> & remoteVals)
    : myLocalSize(localVals.size()), myComm(comm)
    {
      myData = localVals;
      myData.insert(myData.end(), remoteVals.begin(), remoteVals.end());
    }

    DistributedQuantities(const DistributedQuantities& other)
        : myData(other.myData), myLocalSize(other.myLocalSize), myComm(other.myComm) {}

    DistributedQuantities(DistributedQuantities&&) noexcept = default;
    DistributedQuantities& operator=(DistributedQuantities&&) noexcept = default;

    DistributedQuantities& operator=(const DistributedQuantities& other)
    {
        if (this != &other)
        {
            myData = other.myData;
            myLocalSize = other.myLocalSize;
            myComm = other.myComm;
        }
        return *this;
    }

    DistributedQuantities& operator=(std::initializer_list<T> initList)
    {
      myData = initList;
      myLocalSize = initList.size();
      myComm = stk::parallel_machine_null();
      return *this;
    }

    T& operator[](size_t index) { return myData[index]; }
    const T& operator[](size_t index) const { return myData[index]; }

    T * data() { return myData.data(); }
    const T * data() const { return myData.data(); }

    typename std::vector<T>::iterator begin() { return myData.begin(); }
    typename std::vector<T>::const_iterator begin() const { return myData.begin(); }
    typename std::vector<T>::iterator end() { return myData.end(); }
    typename std::vector<T>::const_iterator end() const { return myData.end(); }
    typename std::vector<T>::iterator begin_local() { return myData.begin(); }
    typename std::vector<T>::const_iterator begin_local() const { return myData.begin(); }
    typename std::vector<T>::iterator end_local() { return myData.begin() + myLocalSize; }
    typename std::vector<T>::const_iterator end_local() const { return myData.begin() + myLocalSize; }
    typename std::vector<T>::iterator begin_remote() { return myData.begin() + myLocalSize; }
    typename std::vector<T>::const_iterator begin_remote() const { return myData.begin() + myLocalSize; }
    typename std::vector<T>::iterator end_remote() { return myData.end(); }
    typename std::vector<T>::const_iterator end_remote() const { return myData.end(); }

    size_t size() const  { return myData.size(); }
    size_t local_size() const  { return myLocalSize; }
    std::pair<size_t,size_t> sizes() const { return std::make_pair(myData.size(), myLocalSize); }

    void assign(const stk::ParallelMachine comm, const size_t size, const size_t localSize, const T value);
    void resize(const size_t newSize);
    void resize(const std::pair<size_t,size_t> newSizes);

    std::vector<T> & get() { return myData; }
    const std::vector<T> & get() const { return myData; }

    stk::ParallelMachine comm() const { return myComm; }
};

template<typename T>
void DistributedQuantities<T>::assign(const stk::ParallelMachine comm, const size_t size, const size_t localSize, const T value)
{
  myComm = comm;
  myLocalSize = localSize;
  myData.assign(size, value);
}

template<typename T>
void DistributedQuantities<T>::resize(const size_t newSize)
{
  STK_ThrowAssert(myData.size() == myLocalSize);
  myData.resize(newSize);
  myLocalSize = newSize;
}

template<typename T>
void DistributedQuantities<T>::resize(const std::pair<size_t,size_t> newSizes)
{
  myData.resize(newSizes.first);
  myLocalSize = newSizes.second;
}

using DistributedVector = DistributedQuantities<double>;

DistributedVector xpby(const DistributedVector & x, const double b, const DistributedVector & y);

DistributedVector scalar_times_vector(const double a, const DistributedVector & x);

double Dot(const DistributedVector & x, const DistributedVector & y);

DistributedVector vectorSubtract(const DistributedVector& x, const DistributedVector& y);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_DISTRIBUTEDVECTOR_HPP_ */
