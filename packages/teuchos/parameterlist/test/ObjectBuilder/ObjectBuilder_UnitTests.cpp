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

#include "Teuchos_ObjectBuilder.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"

namespace Teuchos {

const std::string ObjectType_name = "Foo Type";

class Foo : virtual public ParameterListAcceptor {
  public:
    Foo() {}
    virtual ~Foo() {}
    virtual std::string getString() const =0;
    virtual void setDefaults() =0;
    void setParameterList(const RCP<ParameterList> & paramList) {
      if (!is_null(paramList)) {
        paramList->validateParameters(*this->getValidParameters());
        paramList_ = paramList;
      }
      setDefaults();
    }
    RCP<ParameterList> getNonconstParameterList() {
      return paramList_;
    }
    RCP<ParameterList> unsetParameterList() {
      RCP<ParameterList> pl = paramList_;
      paramList_ = null;
      return pl;
    }
    RCP<const ParameterList> getParameterList() const {
      return paramList_;
    }
  private:
    RCP<ParameterList> paramList_;
};
class FooA : virtual public Foo {
  public:
    FooA() {
      setDefaults();
    }
    virtual ~FooA() {}
    std::string getString() const {
      return foo_;
    }
    void setDefaults() {
      RCP<ParameterList> pl = this->getNonconstParameterList();
      if (is_null(pl)) {
        foo_ = "A";
      } else {
        foo_ = pl->get("String",foo_);
      }
    }
    RCP<const ParameterList> getValidParameters() const {
      static RCP<ParameterList> validPL;
      if (is_null(validPL)) {
        RCP<ParameterList> pl = parameterList();
        pl->set( "String", foo_ );
        validPL = pl;
      }
      return validPL;
    }
  private:
    std::string foo_;
};
class FooB : virtual public Foo {
  public:
    FooB() {
      setDefaults();
    }
    virtual ~FooB() {}
    std::string getString() const {
      return foo_;
    }
    void setDefaults() {
      RCP<ParameterList> pl = this->getNonconstParameterList();
      if (is_null(pl)) {
        foo_ = "B";
      } else {
        foo_ = pl->get("String",foo_);
      }
    }
    RCP<const ParameterList> getValidParameters() const {
      static RCP<ParameterList> validPL;
      if (is_null(validPL)) {
        RCP<ParameterList> pl = parameterList();
        pl->set( "String", foo_ );
        validPL = pl;
      }
      return validPL;
    }
  private:
    std::string foo_;
};
class FooC : virtual public Foo {
  public:
    FooC() {
      setDefaults();
    }
    virtual ~FooC() {}
    std::string getString() const {
      return foo_;
    }
    void setDefaults() {
      RCP<ParameterList> pl = this->getNonconstParameterList();
      if (is_null(pl)) {
        foo_ = "C";
      } else {
        foo_ = pl->get("String",foo_);
      }
    }
    RCP<const ParameterList> getValidParameters() const {
      static RCP<ParameterList> validPL;
      if (is_null(validPL)) {
        RCP<ParameterList> pl = parameterList();
        pl->set( "String", foo_ );
        validPL = pl;
      }
      return validPL;
    }
  private:
    std::string foo_;
};

// The following happens at construction:
// 1.  initializeDefaults_ is called
//     a)  object_name_ = "Object"
//     b)  objectType_name_ = "Object Type"
//     c)  defaultObject_ = "None"
//     d)  validObjectNames_ just has "None"
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, constructor) {
  RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>();
  TEST_EQUALITY_CONST( ob->getObjectName(), "None" );
  TEST_EQUALITY_CONST( ob->create(), null );
  RCP<const ParameterList> pl;
  TEST_NOTHROW( pl = ob->getValidParameters() );
  TEST_EQUALITY_CONST( pl->get<std::string>("Object Type"), "None" );
  TEST_NOTHROW( ob = null );
}

// Tests setObjectName and setObectTypeName
// Note:  it should throw an exception if the string is ""
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, setNames) {
  {
    const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>();
    TEST_THROW( ob->setObjectName(""), std::logic_error );
    TEST_THROW( ob->setObjectTypeName(""), std::logic_error );
  }
  {
    RCP<ObjectBuilder<Foo> > ob;
    TEST_THROW( ob = objectBuilder<Foo>("","Foo Type"), std::logic_error );
    TEST_THROW( ob = objectBuilder<Foo>("Foo",""), std::logic_error );
    TEST_THROW( ob = objectBuilder<Foo>("",""), std::logic_error );
  }
  {
    const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>();
    ob->setObjectName("Foo");
    ob->setObjectTypeName("Foo Type");
    const RCP<const ParameterList> validpl = ob->getValidParameters();
    // Now we check that the parameterlist is correct
    TEST_EQUALITY_CONST( validpl->get<std::string>("Foo Type"), "None" );
    const ParameterEntry pe = validpl->getEntry("Foo Type");
    TEST_EQUALITY_CONST( pe.docString(), 
        "Determines the type of Foo object that will be built.\nThe parameters for each Foo Type are specified in this sublist"
        );
  }
  {
    const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>("Foo","Foo Type");
    const RCP<const ParameterList> validpl = ob->getValidParameters();
    // Now we check that the parameterlist is correct
    TEST_EQUALITY_CONST( validpl->get<std::string>("Foo Type"), "None" );
    const ParameterEntry pe = validpl->getEntry("Foo Type");
    TEST_EQUALITY_CONST( pe.docString(), 
        "Determines the type of Foo object that will be built.\nThe parameters for each Foo Type are specified in this sublist"
        );
  }
}

// setObjectFactory does four things:
// 1.  adds a new object name
// 1a.  if object name is "" it throws an exception
// 2.  adds a new object factory
// 3.  sets defaultObject_
// 4.  deletes the validParamList_
//
// Notes about how to sense the changes:
// 1.  The new object name is appended to the list of valid names and shows up in the valid parameter list
// 2.  The new object factory is appended to the list of factories and is only accessible through create
// 3.  The default Object is accessible through both getObjectName and the valid parameter list.
// 4.  The validParameterList is deleted and this can only be sensed through calling getValidParameters
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, setObjectFactory) {
  const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>("Foo","Foo Type");
  TEST_EQUALITY_CONST( ob->getObjectName(), "None" );
  ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
  TEST_EQUALITY_CONST( ob->getObjectName(), "Foo A" );  // 3.
  RCP<const ParameterList> pl = ob->getValidParameters();
  TEST_EQUALITY_CONST( pl->get<std::string>("Foo Type"), "Foo A" ); // 1.
  TEST_EQUALITY_CONST( pl->sublist("Foo A").get<std::string>("String"), "A" ); // 1.
  const RCP<Foo> foo = ob->create();
  const RCP<FooA> fooA = rcp_dynamic_cast<FooA>(foo,false);
  TEST_EQUALITY_CONST( is_null(fooA), false ); // 2.
  ob->setObjectFactory(abstractFactoryStd<Foo,FooB>(),"Foo B");
  pl = ob->getValidParameters();
  TEST_EQUALITY_CONST( pl->get<std::string>("Foo Type"), "Foo B" ); // 4.
  TEST_THROW( ob->setObjectFactory(abstractFactoryStd<Foo,FooC>(),""), std::logic_error ); // 1a.
}

// We shouldn't be able to set two factories with the same name.
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, setObjectFactory_bad ) {
  {
    const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>("Foo","Foo Type");
    ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
    // ObjectBuilder will let you add the object, but will not throw until getValidParameters is called
#ifdef TEUCHOS_DEBUG
    TEST_THROW( ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A"), std::logic_error );
#else // TEUCHOS_DEBUG
    TEST_NOTHROW( ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A") );
    TEST_THROW( ob->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
  }
  {
    const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>("Foo","Foo Type");
    ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
    TEST_NOTHROW( ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"New Foo A") );
    TEST_NOTHROW( ob->getValidParameters() );
  }
}

// getObjectName returns the default in the parameter list (if given), or the
// default in the valid parameter list (if no parameter list is given)
// 1.  no parameter list is given, uses default in valid parameter list.
// 2.  parameter list is given, and uses its default
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, getObjectName) {
  const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>("Foo", "Foo Type");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooB>(),"Foo B");
  const RCP<ParameterList> pl = parameterList();
  pl->setParameters(*ob->getValidParameters()); // copy parameters 
  pl->set("Foo Type", "Foo A"); // change default
  // 1.
  TEST_EQUALITY_CONST( ob->getObjectName(), "Foo B" );
  // 2.
  ob->setParameterList(pl);
  TEST_EQUALITY_CONST( ob->getObjectName(), "Foo A" );
}

// create has many cases
// 1.  It should return a null RCP if no factories are set
// 2.  It should return a null RCP if "Object Type" is set to "None" in the provided parameterList
// 3.  It should return the correct object consistent with the "Object Type" setting in the parameterList if no string is passed
// 3a.  It should return the correct object consistent with the "Object Type"
// setting in the valid parameterList if no string is passed and no
// parameterList is provided.
// 4.  It should return the correct object consistent with the input string regardless of the parameterLists
// 4a.  It should throw an exception if an invalid input string is provided
// 5.  If no parameter list is provided, then it will use the valid parameter list to set parameters on the object
// 5a.  If a parameter list is provided, then it will use that parameter list to set parameters on the object
// 6.  It will throw an exception with a nice message if the factory creates a null RCP
//     Under what conditions could this happen?
// 7.  [03/05/09 tscoffe: found bug]  create() uses objectValidator_, so
// getValidParameters must be valid at the beginning to avoid a null
// dereference of the objectValidator_ pointer in the case that we ask for an
// object by name and the validParamList_ has not been set up yet.
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, create) {
  const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>("Foo", "Foo Type");
  TEST_EQUALITY_CONST( ob->create("None"), null ); // 7.
  TEST_EQUALITY_CONST( ob->create(), null ); // 1.
  ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooB>(),"Foo B");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooC>(),"Foo C");
  out << "op.getValidParamters():\n";
  printValidParameters(*ob, out);
  const RCP<ParameterList> pl = parameterList();
  pl->setParameters(*ob->getValidParameters());
  pl->set("Foo Type","None");
  ob->setParameterList(pl);
  TEST_EQUALITY_CONST( ob->create(), null ); // 2.
  pl->set("Foo Type", "Foo B");
  pl->sublist("Foo B").set("String","BB");
  pl->sublist("Foo C").set("String","CC");
  {
    const RCP<Foo> foo = ob->create();
    const RCP<FooB> fooB = rcp_dynamic_cast<FooB>(foo,false);
    TEST_EQUALITY_CONST( is_null(fooB), false ); // 3.
    TEST_EQUALITY_CONST( foo->getString(), "BB" ); // 5a.
  }
  ob->unsetParameterList();
  {
    const RCP<Foo> foo = ob->create();
    const RCP<FooC> fooC = rcp_dynamic_cast<FooC>(foo,false);
    TEST_EQUALITY_CONST( is_null(fooC), false ); // 3a.
    TEST_EQUALITY_CONST( foo->getString(), "C" ); // 5.
  }
  {
    const RCP<Foo> foo = ob->create("Foo A");
    const RCP<FooA> fooA = rcp_dynamic_cast<FooA>(foo,false);
    TEST_EQUALITY_CONST( is_null(fooA), false ); // 4.
  }
  ob->setParameterList(pl);
  {
    const RCP<Foo> foo = ob->create("Foo A");
    const RCP<FooA> fooA = rcp_dynamic_cast<FooA>(foo,false);
    TEST_EQUALITY_CONST( is_null(fooA), false ); // 4.
  }
  {
    RCP<Foo> foo;
    TEST_THROW( foo = ob->create("Foo D"), std::logic_error ); // 4a.
  }
  // 6. ???
}

// There are many places that the parameter list is validated to ensure that we
// catch invalid parameter lists before we use them.  This is particularly
// important because we're storing a pointer to the parameter list and the user
// can change it without ObjectBuilder knowing about it.
// The parameter list is validated in four places:
// 1. setParameterList 
// 2. unsetParameterList (only in debug mode)
// 3. create (only in debug mode)
// 4. destructor (only in debug mode)
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, setParameterList) {
    RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>();
    ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
    RCP<ParameterList> pl = null;
    TEST_NOTHROW( ob->setParameterList(pl) );
    pl = parameterList();
    TEST_NOTHROW( ob->setParameterList(pl) );
    pl->set("Hello","World");
    TEST_THROW( ob->setParameterList(pl), std::logic_error ); // 1.
#ifdef TEUCHOS_DEBUG
    TEST_THROW( ob->unsetParameterList(), std::logic_error ); // 2.
    TEST_THROW( ob->create(), std::logic_error ); // 3.
    TEST_THROW( ob = null, std::logic_error ); // 4.
#else // TEUCHOS_DEBUG
    TEST_NOTHROW( ob->unsetParameterList() ); 
    RCP<Foo> foo;
    TEST_NOTHROW( foo = ob->create() );
    const RCP<FooA> fooA = rcp_dynamic_cast<FooA>(foo,false);
    TEST_EQUALITY_CONST( is_null(fooA), false );
    TEST_NOTHROW( ob = null );
#endif // TEUCHOS_DEBUG
}

// Here we test 
// 1.  That it returns a null RCP before we give it a parameter list.
// 2.  That we can set up a valid parameter list, give it to the ObjectBuilder, and get it back out.
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, getParameterList) {
  const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>();
  ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
  const RCP<const ParameterList> pl = ob->getParameterList();
  TEST_EQUALITY_CONST( is_null(pl), true ); // 1.
  const RCP<ParameterList> nonconstPL = parameterList();
  nonconstPL->set("Object Type","None");
  TEST_NOTHROW( ob->setParameterList(nonconstPL) );
  {
    const RCP<const ParameterList> newPL = ob->getParameterList();
    TEST_EQUALITY_CONST( nonconstPL.get(), newPL.get() ); // 2.
  }
}

// Same as getParameterList
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, getNonconstParameterList) {
  const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>();
  ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
  RCP<ParameterList> pl = ob->getNonconstParameterList();
  TEST_EQUALITY_CONST( is_null(pl), true );
  pl = parameterList();
  pl->set("Object Type","None");
  TEST_NOTHROW( ob->setParameterList(pl) );
  {
    RCP<ParameterList> newPL = null;
    newPL = ob->getNonconstParameterList();
    TEST_EQUALITY_CONST( pl.get(), newPL.get() );
  }
}

// Here we're checking:
// 1.  That we can set a parameter list on it and it uses it and then we can
// unset it and it goes back to using the valid parameter list.
// 1a.  We get back the same parameter list we set
// 2.  In debug mode, the parameter list is validated when unsetParameterList
// is called.
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, unsetParameterList) {
  RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>();
  ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
  const RCP<ParameterList> pl = parameterList();
  pl->set("Object Type","None");
  ob->setParameterList(pl);
  RCP<Foo> foo = ob->create();
  TEST_EQUALITY_CONST( is_null(foo), true );
  RCP<ParameterList> newPL = ob->unsetParameterList();
  TEST_EQUALITY_CONST( pl.get(), newPL.get() ); // 1a.
  foo = ob->create();
  const RCP<FooA> fooA = rcp_dynamic_cast<FooA>(foo,false);
  TEST_EQUALITY_CONST( is_null(fooA), false ); // 1.
  ob->setParameterList(pl);
  pl->set("Hello","World");
  newPL = null;
#ifdef TEUCHOS_DEBUG
  TEST_THROW( newPL = ob->unsetParameterList(), std::logic_error ); // 2.
  TEST_EQUALITY_CONST( is_null(newPL), true );
  TEST_THROW( ob = null, std::logic_error );
#else // TEUCHOS_DEBUG
  TEST_NOTHROW( newPL = ob->unsetParameterList() );
  TEST_EQUALITY_CONST( pl.get(), newPL.get() ); // 1a.
  TEST_NOTHROW( ob = null );
#endif // TEUCHOS_DEBUG
}

// This function does several things.
// 1.  It creates the validParameterList whenever it is deleted [already tested in setObjectFactory]
// 2.  It creates the objectValidator 
// 3.  It adds a docstring to the "Object Type" parameter in the parameter list [already tested in setNames]
// 4.  It fills the parameter list out with the valid parameteres for each object it can create
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, getValidParameters) {
  {
    const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>();
    ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
    const RCP<ParameterList> pl = parameterList();
    pl->set("Object Type","Foo B");
    TEST_THROW( ob->setParameterList(pl), std::logic_error ); // 2.
  }
  {
    const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>();
    ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
    ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo B");
    ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo C");
    const RCP<ParameterList> validPL = parameterList();
    validPL->set("Object Type","Foo C");
    validPL->sublist("Foo A").set("String","A");
    validPL->sublist("Foo B").set("String","B");
    validPL->sublist("Foo C").set("String","C");
    Array<std::string> validObjectNames;
    validObjectNames.push_back("None");
    validObjectNames.push_back("Foo A");
    validObjectNames.push_back("Foo B");
    validObjectNames.push_back("Foo C");
    const RCP<const StringToIntegralParameterEntryValidator<int> > 
      objectValidator = rcp(
        new StringToIntegralParameterEntryValidator<int>(
          validObjectNames,"Object Type"
          )
        );
    validPL->set(
      "Object Type","Foo C"
      ,(std::string("Determines the type of Object object that will be built.\n")
        + "The parameters for each Object Type are specified in this sublist" 
        ).c_str()
      ,objectValidator
      );
    const RCP<const ParameterList> pl = ob->getValidParameters();
    TEST_NOTHROW( pl->validateParameters(*validPL) ); // 4.
    validPL->set("Object Type","Foo A");
    TEST_NOTHROW( pl->validateParameters(*validPL) ); // 4.
    validPL->set("Object Type","Foo B");
    TEST_NOTHROW( pl->validateParameters(*validPL) ); // 4.
    validPL->set("Object Type","None");
    TEST_NOTHROW( pl->validateParameters(*validPL) ); // 4.
  }
}

// Now we verify that the parameter lists are coming out with Used parameters in the correct state
// 1.  Pass in empty parameter list and create an object.  We should get a
// sublist and used parameters on the sublist for the object we created, but no
// other sublists.
// 2.  Pass in a full parameter list and create an object.  We should get 
// used parameters for only the sublist of the object we created.
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, usedParameters) {
  const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>("Foo","Foo Type");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooB>(),"Foo B");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooC>(),"Foo C");
  {
    const RCP<ParameterList> pl = parameterList();
    ob->setParameterList(pl);
    const RCP<Foo> foo = ob->create("Foo A");
    TEST_EQUALITY_CONST( foo->getString(), "A" );
    TEST_EQUALITY_CONST( pl->isSublist("Foo A"), true ); // 1.
    TEST_EQUALITY_CONST( pl->sublist("Foo A").isParameter("String"), true ); // 1.
    const ParameterEntry& pe = pl->sublist("Foo A").getEntry("String");
    TEST_EQUALITY_CONST( pe.isUsed(), true ); // 1.
    TEST_EQUALITY_CONST( pe.isDefault(), true ); // 1.
    // verify the other sublists are missing
    TEST_EQUALITY_CONST( pl->isSublist("Foo B"), false ); // 1.
    TEST_EQUALITY_CONST( pl->isSublist("Foo C"), false ); // 1.
    ob->unsetParameterList();
  }
  {
    RCP<ParameterList> pl = parameterList();
    pl->setParameters(*ob->getValidParameters());
    pl->sublist("Foo A").set("String","AA");
    ob->setParameterList(pl);
    pl = null;
    const RCP<Foo> foo = ob->create("Foo A");
    TEST_EQUALITY_CONST( foo->getString(), "AA" );
    const RCP<const ParameterList> outPL = ob->getParameterList();
    TEST_EQUALITY_CONST( outPL->isSublist("Foo A"), true );
    TEST_EQUALITY_CONST( outPL->sublist("Foo A").isParameter("String"), true );
    const ParameterEntry& pe = outPL->sublist("Foo A").getEntry("String");
    TEST_EQUALITY_CONST( pe.isUsed(), true ); // 2.
    TEST_EQUALITY_CONST( pe.isDefault(), false ); // 2.
    // verify the other sublists are unused
    TEST_EQUALITY_CONST( outPL->sublist("Foo B").getEntry("String").isUsed(), false ); // 2. 
    TEST_EQUALITY_CONST( outPL->sublist("Foo C").getEntry("String").isUsed(), false ); // 2. 
    ob->unsetParameterList();
  }
}  

TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, setDefaultObject_withOneUsePL ) {
  const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>("Foo","Foo Type");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooB>(),"Foo B");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooC>(),"Foo C");
  {
    const RCP<ParameterList> pl = parameterList();
    ob->setParameterList(pl);
    const RCP<Foo> foo = ob->create();
    RCP<FooC> fooC = Teuchos::rcp_dynamic_cast<FooC>(foo,false);
    TEST_ASSERT( !is_null(fooC) );
  }
  {
    const RCP<ParameterList> pl = parameterList();
    ob->setParameterList(pl);
    ob->setDefaultObject("Foo A");
    const RCP<Foo> foo = ob->create();
    RCP<FooA> fooA = Teuchos::rcp_dynamic_cast<FooA>(foo,false);
    TEST_ASSERT( !is_null(fooA) );
  }
  {
    const RCP<ParameterList> pl = parameterList();
    ob->setParameterList(pl);
    ob->setDefaultObject("None");
    const RCP<Foo> foo = ob->create();
    TEST_ASSERT( is_null(foo) );
  }
  {
#ifdef TEUCHOS_DEBUG
    TEST_THROW(ob->setDefaultObject("Foo D"), std::logic_error);
#else
    ob->setDefaultObject("Foo D");
    TEST_THROW(ob->getValidParameters(), std::logic_error);
#endif // TEUCHOS_DEBUG
  }
}
TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, setDefaultObject_withMultipleUsePL ) {
  const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>("Foo","Foo Type");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooB>(),"Foo B");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooC>(),"Foo C");
  const RCP<ParameterList> pl = parameterList();
  ob->setParameterList(pl);
  {
    const RCP<Foo> foo = ob->create();
    RCP<FooC> fooC = Teuchos::rcp_dynamic_cast<FooC>(foo,false);
    TEST_ASSERT( !is_null(fooC) );
    // Note:  At this point, pl contains "Foo Type = Foo C"
    // And this pl was set on the ObjectBuilder, so defaultObject does no good.
  }
  {
    ob->setDefaultObject("Foo A");
    const RCP<Foo> foo = ob->create();
    RCP<FooA> fooA = Teuchos::rcp_dynamic_cast<FooA>(foo,false);
    TEST_ASSERT( is_null(fooA) );
  }
  {
    ob->setDefaultObject("None");
    const RCP<Foo> foo = ob->create();
    TEST_ASSERT( !is_null(foo) );
  }
}

TEUCHOS_UNIT_TEST( Teuchos_ObjectBuilder, setDefaultObject_withoutPL ) {
  const RCP<ObjectBuilder<Foo> > ob = objectBuilder<Foo>("Foo","Foo Type");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooA>(),"Foo A");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooB>(),"Foo B");
  ob->setObjectFactory(abstractFactoryStd<Foo,FooC>(),"Foo C");
  {
    const RCP<Foo> foo = ob->create();
    RCP<FooC> fooC = Teuchos::rcp_dynamic_cast<FooC>(foo,false);
    TEST_ASSERT( !is_null(fooC) );
  }
  {
    ob->setDefaultObject("Foo A");
    const RCP<Foo> foo = ob->create();
    RCP<FooA> fooA = Teuchos::rcp_dynamic_cast<FooA>(foo,false);
    TEST_ASSERT( !is_null(fooA) );
  }
  {
    ob->setDefaultObject("None");
    const RCP<Foo> foo = ob->create();
    TEST_ASSERT( is_null(foo) );
  }
}

} // namespace Teuchos



