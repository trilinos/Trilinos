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

#ifndef TEUCHOS_CONST_NONCONST_OBJECT_CONTAINER_HPP
#define TEUCHOS_CONST_NONCONST_OBJECT_CONTAINER_HPP

#include "Teuchos_RCP.hpp"


namespace Teuchos {


/** \brief Simple class supporting the "runtime protection of const" idiom.
 *
 * This is a foundational class for supporting the "runtime protection of
 * const" idiom.  The problem this class is designed to help solve is the
 * general issue of const protection and const handling for "held" objects
 * inside of "container" objects.  The use case this class is designed to
 * support involves having the client create the "held" object, give it to the
 * "container" object, and where the "container" object has functions to give
 * the "held" object back again.  In this case, there are to specific roles
 * the "container" object is performing.  One role, the primary role, is the
 * primary function the "container" object was designed to perform where it
 * needs functionality of the the "held" object where only the const interface
 * of the "held" object is required.  If this primary role were the only
 * consideration, we could just write the "container" class as:

 \code

  // Basic "container" implementation that does not consider general
  // const/non-const issues.
  class Container {
  public:
    setHeld(const RCP<const Held> &held)
      { held_ = held; }
    RCP<const Held> getHeld() const
      { return held_; }
    void doSomething() // The primary role!
      { stuff = held_->computeSomething(...); } // const interface of Held!
  private:
    RCP<const Held> held_;
  };

 \endcode

 * The problem with this design of the "container" class is that it does not
 * well support the second natural role of any such "container" object, and
 * that is to act as a general object container that can be used to store and
 * extract the "held" object.  The difficulty that occurs is when the client
 * has a non-const reference to the "held" object, gives it to the "container"
 * object and then needs to get back a non-const reference to the "held"
 * object later.  With the current design, the client code must do a const
 * cast such as in:

 \code

  void setUpContainer( const Ptr<Container> &container )
  {

    // A non-const version of Held
    RCP<Held> myHeld = createNewHeld(...);

    // Give my non-const RCP to Held as a const RCP to held to container
    container->setHeld(myHeld);

  }


  void updateContainer( const Ptr<Container> &container )
  {

    // Get back a non-const version of Held (WARNING: const_cast!)
    RCP<Held> myHeld = rcp_const_cast<Held>(container->getHeld());

    // Change Held
    myHeld->changeSomething(...);

    // Put back Held
    container->setHeld(myHeld);

  }

 \endcode

 * Code like shown above if very common and exposes the core problem.  The
 * client should not have to const cast to get back a non-const verison of the
 * "held" object that it put in the "container" object in the first place.
 * The "container" object should know that it was given a non-const version of
 * "held" object and it should be able to give back a non-const version of the
 * "held" object.  As much as possible, const casting should be eliminated
 * from the code, especially user code.  Const casting is a source of defects
 * in C++ programs and violates the flow of C++ programming (See Item 94
 * "Avoid casting away const" in the book "C++ Coding Standards").
 *
 * The design of the "container" class using this class
 * ConstNonconstObjectContainer that resolves the problem is:

 \code

  // Implementation of container that uses the "runtime protection of const"
  // to hold and give up the "held" object.
  class Container {
  public:
    setNonconstHeld(const RCP<Held> &held)
      { held_ = held; }
    setHeld(const RCP<const Held> &held)
      { held_ = held; }
    RCP<const Held> getNonconstHeld()
      { return held_.getNonconstObj(); }
    RCP<const Held> getHeld() const
      { return held_.getConstObj(); }
    void doSomething() // The primary role!
      { stuff = held_->computeSomething(...); } // const interface of Held
  private:
    ConstNonconstObjectContainer<Held> held_;
  };

 \endcode

 * Now the client code can be written with no const casting as:

 \code

  void setUpContainer( const Ptr<Container> &container )
  {

    // A non-const version of Held
    RCP<Held> myHeld = createNewHeld(...);

    // Give my non-const RCP to Held now stored as a non-const object
    container->setNonconstHeld(myHeld);

  }


  void updateContainer( const Ptr<Container> &container )
  {

    // Get back a non-const version of Held (No const cating!)
    RCP<Held> myHeld = container->getNonconstHeld();

    // Change Held
    myHeld->changeSomething(...);

    // Put back Held
    container->setNonconstHeld(myHeld);

  }

 \endcode

 * The "runtime protection of const" idiom allows you to write a single
 * "container" class that can hold both non-const and const forms of a "held"
 * object, protects the const of objects being set as const, and can give back
 * non-const references to objects set as non-const.  The price one pays for
 * this is that the typical compile-time const protection provided by C++ is
 * instead replaced with a runtime check.  For example, the following code
 * with thrown a <tt>NonconstAccessError</tt> exception object:

 \code

  void fooThatThrows(const Ptr<Container> &container)
  {
 
    // A non-const version of Held
    RCP<Held> myHeld = createNewHeld(...);
  
    // Accidentally set a const version of Held
    container->setHeld(myHeld);
    
    // Try to get back a non-const version of Held
    RCP<Held> myHeldAgain = container->getNonconstHeld(); // Throws NonconstAccessError!

  }

 \endcode

 * These types of exceptions can be confuing to developers if they don't
 * understand the idiom.
 *
 * The alternative to the "runtime protection of const" idiom is to use
 * compile-time protection.  However, using compile-time const protection
 * would require two different versions of a the "container" class: a
 * "Container" class and a "ConstContainer" class.  I will not go into detail
 * about what these classes look like but this is ugly, more confusing, and
 * hard to maintain.
 *
 * Note that classes like RCP and boost:shared_ptr provide for compile-time
 * protection of const with just one (template) class definition.  RCP objects
 * of type RCP<Held> allow non-const access while RCP objects of type
 * RCP<const Held> only allow const access and protect const at compile time.
 * How can one class like RCP protect const at compile-time while a class like
 * Container shown above can't?  The reason of course is that RCP<Held> and
 * RCP<const Held> are realy *two* different C++ classes.  The template
 * mechanism in C++ made it easy to create these two different class types but
 * they are two seperate types none the less.
 *
 * Note that the "runtime protection of const" idiom using this
 * ConstNonconstObjectContainer is not necessary when the "container" object
 * needs a non-const "held" object to do its primary work.  In this case, a
 * client can't give a "container" object a non-const version of the "held"
 * object because it could not even do its primary role.  In cases where a
 * non-const version of "held" is needed for the primary role, the "container"
 * class can be written more simply without ConstNonconstObjectContainer as:

 \code

  // Simpler implementation of "container" where a non-const version of the
  // "held" object is needed to perform the primary role.
  class Container {
  public:
    setHeld(const RCP<Held> &held)
      { held_ = held; }
    RCP<const Held> getNonconstHeld()
      { return held_; }
    RCP<const Held> getHeld() const
      { return held_.; }
    void doSomething() // The primary role!
      { held_->changeSomething(...); } // non-const interface of Held
  private:
    RCP<Held> held_;
  };

 \endcode

 * NOTE: The default copy constructor and assignment operator functions are
 * allowed and result in shallow copy (i.e. just the RCP objects are copied).
 * However, the protection of const will be maintained in the copied/assigned
 * objects correctly.
 *
 * NOTE: Assignment for an RCP<const ObjType> is also supported due to the
 * implicit conversion from RCP<const ObjType> to
 * ConstNonconstObjectContainer<ObjType> that this class supports through its
 * constructor.
 */
template<class ObjType>
class ConstNonconstObjectContainer {
public:
  /** \brief. Constructs to uninitialized */
  ConstNonconstObjectContainer()
    :constObj_(null),isConst_(true) {}
  /** \brief. Calls <tt>initialize()</tt> with a non-const object. */
  ConstNonconstObjectContainer( const RCP<ObjType> &obj )
    { initialize(obj); }
  /** \brief. Calls <tt>initialize()</tt> with a const object. */
  ConstNonconstObjectContainer( const RCP<const ObjType> &obj )
    { initialize(obj); }
  /** \brief. Initialize using a non-const object.
   * Allows both const and non-const access to the contained object. */
  void initialize( const RCP<ObjType> &obj )
    {
      TEUCHOS_TEST_FOR_EXCEPTION(is_null(obj), NullReferenceError, "Error!");
      constObj_ = obj;
      isConst_ = false;
    }
  /** \brief. Initialize using a const object.
   * Allows only const access enforced with a runtime check. */
  void initialize( const RCP<const ObjType> &obj )
    {
      TEUCHOS_TEST_FOR_EXCEPTION(is_null(obj), NullReferenceError, "Error!");
      constObj_ = obj; 
      isConst_ = true;
    }
  /** \biref Uninitialize. */
  void uninitialize()
    { constObj_=null; isConst_=true; }
  /** \breif Assign to null. */
  ConstNonconstObjectContainer<ObjType>& operator=(ENull)
    { uninitialize(); return *this; }
  /** \brief Returns true if const-only access to the object is allowed. */
  bool isConst() const
    { return isConst_; }
  /** \brief Get an RCP to the non-const contained object.
   *
   * <b>Preconditions:</b>
   * <ul>
   * <li> [<tt>getConstObj().get()!=NULL</tt>] <tt>isConst()==false</tt>
   *      (throws <tt>NonconstAccessError</tt>)
   * </ul>
   *
   * <b>Postconditions:</b>
   * <ul>
   * <li>[<tt>getConstObj().get()==NULL</tt>] <tt>return.get()==NULL</tt>
   * <li>[<tt>getConstObj().get()!=NULL</tt>] <tt>return.get()!=NULL</tt>
   * </ul>
   */
  RCP<ObjType> getNonconstObj() const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        constObj_.get() && isConst_, NonconstAccessError,
        "Error, the object of reference type \""<<TypeNameTraits<ObjType>::name()
        <<"\" was given as a const-only object and non-const access is not allowed."
        );
      return rcp_const_cast<ObjType>(constObj_);
    }
  /** \brief Get an RCP to the const contained object.
   *
   * If <tt>return.get()==NULL</tt>, then this means that no object was given
   * to <tt>*this</tt> data container object.
   */
  RCP<const ObjType> getConstObj() const
    { return constObj_; }
  /** \brief Perform shorthand for <tt>getConstObj(). */
  RCP<const ObjType> operator()() const
    { return getConstObj(); }
  /** \brief Pointer (<tt>-></tt>) access to underlying const object.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
   * </ul>
   */
  const ObjType* operator->() const
    { return &*getConstObj(); } // Does assert also!
  /** \brief Dereference the underlying object.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
   * </ul>
   */
  const ObjType& operator*() const
    { return *getConstObj(); }
  /** \brief Perform an implicit conversion to an RCP<const ObjType>. */
  operator RCP<const ObjType>() const
    { return getConstObj(); }
  /** \brief Return the internal count. */
  int count() const
    { return constObj_.count(); }

private:
  RCP<const ObjType> constObj_;
  bool isConst_;
};


/** \brief Returns true if <tt>p.get()==NULL</tt>.
 *
 * \relates ConstNonconstObjectContainer
 */
template<class T>
bool is_null( const ConstNonconstObjectContainer<T> &p )
{ return is_null(p.getConstObj()); }


/** \brief Returns true if <tt>p.get()!=NULL</tt>.
 *
 * \relates ConstNonconstObjectContainer
 */
template<class T>
bool nonnull( const ConstNonconstObjectContainer<T> &p )
{ return nonnull(p.getConstObj()); }


} // namespace Teuchos


#endif // TEUCHOS_CONST_NONCONST_OBJECT_CONTAINER_HPP
