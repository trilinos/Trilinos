#include "stk_util/stk_config.h"
#include "ObjectLifetimeSpy.hpp"
#include <iostream>

namespace stk::unit_test_util {

static int s_objectLifetimeSpy_numConstructions = 0;
static int s_objectLifetimeSpy_numCopyConstructions = 0;
static int s_objectLifetimeSpy_numMoveConstructions = 0;
static int s_objectLifetimeSpy_numDestructions = 0;
static int s_objectLifetimeSpy_numCopyAssignments = 0;
static int s_objectLifetimeSpy_numMoveAssignments = 0;

void objectLifetimeSpy_clearCounts()
{
  s_objectLifetimeSpy_numConstructions = 0;
  s_objectLifetimeSpy_numCopyConstructions = 0;
  s_objectLifetimeSpy_numMoveConstructions = 0;
  s_objectLifetimeSpy_numDestructions = 0;
  s_objectLifetimeSpy_numCopyAssignments = 0;
  s_objectLifetimeSpy_numMoveAssignments = 0;
}

int objectLifetimeSpy_getNumConstructions()
{
  return s_objectLifetimeSpy_numConstructions;
}

int objectLifetimeSpy_getNumCopyConstructions()
{
  return s_objectLifetimeSpy_numCopyConstructions;
}

int objectLifetimeSpy_getNumMoveConstructions()
{
  return s_objectLifetimeSpy_numMoveConstructions;
}

int objectLifetimeSpy_getNumDestructions()
{
  return s_objectLifetimeSpy_numDestructions;
}

int objectLifetimeSpy_getNumCopyAssignments()
{
  return s_objectLifetimeSpy_numCopyAssignments;
}

int objectLifetimeSpy_getNumMoveAssignments()
{
  return s_objectLifetimeSpy_numMoveAssignments;
}

bool objectLifetimeSpy_checkBalancedConstructionsDestructions()
{
  return ( (s_objectLifetimeSpy_numConstructions +
            s_objectLifetimeSpy_numCopyConstructions +
            s_objectLifetimeSpy_numMoveConstructions) == s_objectLifetimeSpy_numDestructions);
}


ObjectLifetimeSpy::ObjectLifetimeSpy()
  : m_i(0)
{
#ifdef STK_USE_OPENMP
#pragma omp critical
#endif
  ++s_objectLifetimeSpy_numConstructions;
  if constexpr (m_printEvents) {
    std::cout << "ObjectLifetimeSpy(): Default Constructor " << m_i << std::endl;
  }
}

ObjectLifetimeSpy::ObjectLifetimeSpy(int i) : m_i(i)
{
#ifdef STK_USE_OPENMP
#pragma omp critical
#endif
  ++s_objectLifetimeSpy_numConstructions;
  if constexpr (m_printEvents) {
    std::cout << "ObjectLifetimeSpy(): Constructor " << m_i << std::endl;
  }
}

ObjectLifetimeSpy::~ObjectLifetimeSpy()
{
#ifdef STK_USE_OPENMP
#pragma omp critical
#endif
  ++s_objectLifetimeSpy_numDestructions;
  if constexpr (m_printEvents) {
    std::cout << "ObjectLifetimeSpy(): Destructor " << m_i << std::endl;
  }
}

ObjectLifetimeSpy::ObjectLifetimeSpy(const ObjectLifetimeSpy& rhs)
{
#ifdef STK_USE_OPENMP
#pragma omp critical
#endif
  ++s_objectLifetimeSpy_numCopyConstructions;
  m_i = rhs.m_i;
  if constexpr (m_printEvents) {
    std::cout << "ObjectLifetimeSpy(): Copy Constructor " << m_i << std::endl;
  }
}

ObjectLifetimeSpy::ObjectLifetimeSpy(ObjectLifetimeSpy&& rhs) noexcept
{
#ifdef STK_USE_OPENMP
#pragma omp critical
#endif
  ++s_objectLifetimeSpy_numMoveConstructions;
  m_i = rhs.m_i;
  if constexpr (m_printEvents) {
    std::cout << "ObjectLifetimeSpy(): Move Constructor " << m_i << std::endl;
  }
}

ObjectLifetimeSpy&
ObjectLifetimeSpy::operator=(const ObjectLifetimeSpy& rhs)
{
#ifdef STK_USE_OPENMP
#pragma omp critical
#endif
  ++s_objectLifetimeSpy_numCopyAssignments;
  m_i = rhs.m_i;
  if constexpr (m_printEvents) {
    std::cout << "ObjectLifetimeSpy(): Copy assignment operator " << m_i << std::endl;
  }
  return *this;
}

ObjectLifetimeSpy&
ObjectLifetimeSpy::operator=(ObjectLifetimeSpy&& rhs)
{
#ifdef STK_USE_OPENMP
#pragma omp critical
#endif
  ++s_objectLifetimeSpy_numMoveAssignments;
  m_i = rhs.m_i;
  if constexpr (m_printEvents) {
    std::cout << "ObjectLifetimeSpy(): Move assignment operator " << m_i << std::endl;
  }
  return *this;
}

}
