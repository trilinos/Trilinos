/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_Version.hpp"

class A { 
 public: 
   A() {}
   virtual ~A(){} 
   virtual void f(){} 
};   
class B1 : virtual public A {};
class B2 : virtual public A {};
class C : public B1, public B2 {};

using namespace Teuchos;

int main(int argc, char* argv[])
{

  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  // Create some reference-counted pointers.
  // Create a reference-counted NULL pointer of type A.
  RCP<A>	           a_null_ptr;
  // Create a reference-counted pointer of non-const type A.
  RCP<A>             a_ptr   = rcp(new A);
  // Create a reference-counted pointer of const type A.
  RCP<const A>       ca_ptr  = rcp(new A);
  // Create a const reference-counted pointer of non-const type A.
  const RCP<A>       a_cptr  = rcp(new A); 
  // Create a const reference-counted pointer of const type A.
  const RCP<const A> ca_cptr = rcp(new A); 

  // Perform implicit conversions between a derived class and its base class.
  RCP<B1> b1_ptr  = rcp(new B1);
  RCP<A> a_ptr1 = b1_ptr;

  /* Other non-implicit type conversions like static, dynamic, or const casts
     can be taken care of by non-member template functions.
  */
  RCP<const C>  c_ptr     = rcp(new C);
  // Implicit cast from C to B2.
  RCP<const B2> b2_ptr    = c_ptr;                              
  // Safe cast, type-checked, from C to A.
  RCP<const A>  ca_ptr1   = rcp_dynamic_cast<const A>(c_ptr); 
  // Unsafe cast, non-type-checked, from C to A.
  RCP<const A>  ca_ptr2   = rcp_static_cast<const A>(c_ptr);  
  // Cast away const from B2.
  RCP<B2>       nc_b2_ptr = rcp_const_cast<B2>(b2_ptr);           

  /* Using a reference-counted pointer is very similar to using a raw C++ pointer.  Some
     of the operations that are common to both are:
  */
  RCP<A>
    a_ptr2 = rcp(new A),  // Initialize reference-counted pointers.
    a_ptr3 = rcp(new A);  // ""
  A  *ra_ptr2 = new A,    // Initialize non-reference counted pointers.
    *ra_ptr3 = new A;     // ""
  a_ptr2 = rcp(ra_ptr3);  // Assign from a raw pointer (only do this once!)
  a_ptr3 = a_ptr1;        // Assign one smart pointer to another.
  a_ptr2 = rcp(ra_ptr2);  // Assign from a raw pointer (only do this once!)
  a_ptr2->f();            // Access a member of A using ->
  ra_ptr2->f();           // ""
  *a_ptr2 = *a_ptr3;      // Dereference the objects and assign.
  *ra_ptr2 = *ra_ptr3;    // ""

  // Get the raw C++ pointer.
  A* true_ptr = 0;
  true_ptr = a_ptr1.get();

  return 0;

}
