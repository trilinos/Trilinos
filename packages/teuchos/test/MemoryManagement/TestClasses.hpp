// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_TEST_CLASSES_HPP
#define TEUCHOS_TEST_CLASSES_HPP


#include "Teuchos_RCP.hpp"


// Return constants from class functions
const int A_g_return  = 1;
const int A_f_return  = 2;
const int B1_g_return = 3;
const int B1_f_return = 4;
const int B2_g_return = 5;
const int B2_f_return = 6;
const int C_g_return  = 7;
const int C_f_return  = 8;
const int D_g_return  = 9;
const int D_f_return  = 10;
const int E_g_return  = 11;
const int E_f_return  = 12;


/*

 Polymorphic multiple inheritance example

            -----
           |  A  |
            -----
             /|\
              | 
         ------------
        |            |
      -----        ------
     |  B1 |      |  B2  |
      -----        ------
       /|\          /|\
        |            |
         ------------
              |
            -----
           |  C  |
            -----

*/


class C;


class A {
	int A_g_, A_f_;
public:
	A() : A_g_(A_g_return), A_f_(A_f_return) {}
	virtual ~A(); // See below
	virtual int A_g() { return A_g_; }
	virtual int A_f() const { return A_f_; }
  int call_C_f();
private:
  Teuchos::RCP<C> c_;
public:
  void set_C(const Teuchos::RCP<C> &c ) { c_ = c; }
};


class B1 : virtual public A {
	int B1_g_, B1_f_;
public:
	B1() : B1_g_(B1_g_return), B1_f_(B1_f_return) {}
	~B1() { B1_g_ = -1; B1_f_ = -1; }
	virtual int B1_g() { return B1_g_; }
	virtual int B1_f() const { return B1_f_; }
};


class B2 : virtual public A {
	int B2_g_, B2_f_;
public:
	B2() : B2_g_(B2_g_return), B2_f_(B2_f_return) {}
	~B2() { B2_g_ = -1; B2_f_ = -1; }
	virtual int B2_g() { return B2_g_; }
	virtual int B2_f() const { return B2_f_; }
};


class C : virtual public B1, virtual public B2
{
	int C_g_, C_f_;
public:
	C() : C_g_(C_g_return), C_f_(C_f_return), call_A_on_delete_(false)
    {
      A_g_on_delete_ = -2;
    }
	~C()
    {
      C_g_ = -1; C_f_ = -1;
      if (call_A_on_delete_) {
        // VERY BAD THING TO DO!
        A_g_on_delete_ = call_A_g();
        // NOTE: If a_ is a weak pointer and the underlying 'A' object has
        // already been deleted, then this destructor will throw an exception.
        // This is *never* a good thing to do in production code.  However, I
        // am allowing this destructor to throw an exception so I can write a
        // unit test to detect this.
      }
    }
	virtual int C_g() { return C_g_; }
	virtual int C_f() const { return C_f_; }
  void call_A_on_delete(bool call_A_on_delete_in)
    { call_A_on_delete_ = call_A_on_delete_in; }
  int call_A_g() { return a_->A_g(); }
  static int get_A_g_on_delete() { return A_g_on_delete_; }
private:
  Teuchos::RCP<A> a_;
  bool call_A_on_delete_;
  static int A_g_on_delete_;
public:
  void set_A(const Teuchos::RCP<A> &a ) { a_ = a; }
  Teuchos::RCP<A> get_A() { return a_; }
};


// Need to put these here if we have circular references

inline
A::~A() { A_g_ = -1; A_f_ = -1; }


inline
int A::call_C_f() { return c_->C_f(); }


class Get_A_f_return {
  const A *a_;
  int *a_f_return_;
  Get_A_f_return();
public:
  Get_A_f_return( const A *a, int *a_f_return ) : a_(a), a_f_return_(a_f_return) {}
  ~Get_A_f_return() { *a_f_return_ = a_->A_f(); }
};


void deallocA(A* ptr);


void deallocHandleA(A** handle);


/*

 Non-polymophic classes hiearchy examlpe

            -----
           |  D  |
            -----
             /|\
              | 
            -----
           |  E  |
            -----

*/


class D 
{
	int D_g_, D_f_;
public:
	D() : D_g_(D_g_return), D_f_(D_f_return) {}
	int D_g() { return D_g_; }
	int D_f() const { return D_f_; }
};


class E : public D
{
	int E_g_, E_f_;
public:
	E() : E_g_(E_g_return), E_f_(E_f_return) {}
	int E_g() { return E_g_; }
	int E_f() const { return E_f_; }
};


/*

Typedef to pointer for undefined struct as an opaque object type without a
specialization of TypeNameTraits.

This simulates what happens with a lot of MPI implementations.

*/

struct UndefinedType; // Forward declared but never defined!
typedef UndefinedType* Opaque_handle;
const Opaque_handle OPAQUE_HANDLE_NULL = 0;
Opaque_handle createOpaque();
const int getOpaqueValue_return = 5;
int getOpaqueValue( Opaque_handle opaque );
void destroyOpaque( Opaque_handle * opaque );


/*

Typedef to pointer for an undefiend struct as an opaque object type out a
specialization of TypeNameTraits of the actually type.

This allows it to be stored in an RCP object itself.

*/

struct UndefinedType2; // Forward declared but never defined!
typedef UndefinedType2* Opaque2_handle;
const Opaque2_handle OPAQUE2_HANDLE_NULL = 0;
Opaque2_handle createOpaque2();
const int getOpaque2Value_return = 8;
int getOpaque2Value( Opaque2_handle opaque );
void destroyOpaque2( Opaque2_handle * opaque );


namespace Teuchos {


// Here we define the traits for the underlying type itself.
template<>
class TypeNameTraits<UndefinedType2> {
public:
  static std::string name() { return "UndefinedType2"; }
  static std::string concreteName(const UndefinedType2&)
    { return name(); }
};


} // namespace Teuchos


/*

Typedef to pointer for an undefiend struct as an opaque object type out a
specialization of TypeNameTraits of the actually type.

This allows handles to the type be used with Array, ArrayRCP, and ArrayView.
However, this type can *not* be used with RCP since it does not define a
TypeNameTraits specialization for the underlying undefined type.

This simulates what can happen with MPI implementations.

*/

struct UndefinedType3; // Forward declared but never defined!
typedef UndefinedType3* Opaque3_handle;
const Opaque3_handle OPAQUE3_HANDLE_NULL = 0;


namespace Teuchos {

// Here we only define the traits class for the handle type and we don't even
// need to worry about what the underlying type is (unless we already have a
// speicalization defined for it).
template<>
class TypeNameTraits<Opaque3_handle> {
public:
  static std::string name() { return "Opaque3_handle"; }
  static std::string concreteName(Opaque3_handle)
    { return name(); }
};


} // namespace Teuchos


#endif // TEUCHOS_TEST_CLASSES_HPP
