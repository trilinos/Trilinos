// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// An example that demonstrates usage of Sacado::TemplateContainer<> for storing
// objects templated on the scalar type

#include "Sacado_TemplateContainer.hpp"

#include <iostream>
#include "Sacado_mpl_vector.hpp"
#include "Sacado_mpl_placeholders.hpp"
#include "Sacado.hpp"

// A templated class that will be instantiated for several types T
template <class T>
struct MyClass {
  T x;
  MyClass() : x() {}  // to avoid uninitialized value warning below
};

// A functor to initialize a container of objects of type MyClass<T>
template <class Container>
struct SetFunctor {
  Container& container;
  double val;
  SetFunctor(Container& c, double v) : container(c), val(v) {}
  template <typename T> void operator()(T) const {
    container.template get<T>().x = val;
  }
};

int main() {
  using Sacado::mpl::placeholders::_;

  // Our scalar types
  typedef Sacado::Fad::DFad<double> FadType;
  typedef Sacado::mpl::vector<double,FadType> MyTypes;

  // Container to hold objects of type MyClass<T> for each T in MyTypes
  typedef Sacado::TemplateContainer< MyTypes,MyClass<_> > MyObjects;
  MyObjects myObjects;

  // Print out their initial values
  std::cout << myObjects.get<double>().x << std::endl;
  std::cout << myObjects.get<FadType>().x << std::endl << std::endl;

  // Set the values to 1.0 using container_for_each
  Sacado::container_for_each_no_kokkos( myObjects, SetFunctor<MyObjects>(myObjects,1) );
  std::cout << myObjects.get<double>().x << std::endl;
  std::cout << myObjects.get<FadType>().x << std::endl << std::endl;

  // Set the values to 2.0 using mpl::for_each
  Sacado::mpl::for_each_no_kokkos<MyTypes>( SetFunctor<MyObjects>(myObjects,2) );
  std::cout << myObjects.get<double>().x << std::endl;
  std::cout << myObjects.get<FadType>().x << std::endl << std::endl;

  // Set the values to 3.0 using mpl::for_each
  Sacado::mpl::for_each_no_kokkos<MyObjects>( SetFunctor<MyObjects>(myObjects,3) );
  std::cout << myObjects.get<double>().x << std::endl;
  std::cout << myObjects.get<FadType>().x << std::endl << std::endl;

  // Test
  bool passed = ( myObjects.get<double>().x == 3.0 &&
                  myObjects.get<FadType>().x.val() == 3.0 );
  if (passed)
    std::cout << "Test Passed!" << std::endl;
  else
    std::cout << "Test Failed!" << std::endl;

  /*
  // This requires C++14 to use generalized lambdas
  Sacado::mpl::for_each_no_kokkos<MyObjects>( [&](auto x) {
       typedef decltype(x) T;
       myObjects.get<T>().x = 4;
    });
  Sacado::mpl::for_each_no_kokkos<MyObjects>( [&](auto x) {
       typedef decltype(x) T;
       std::cout << myObjects.get<T>().x << std::endl;
    });
  */

  return passed ? 0 : 1;
}
