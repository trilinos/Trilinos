//@HEADER
// ***********************************************************************
//
//                           Teuchos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Teuchos_OBJECT_BUILDER_H
#define Teuchos_OBJECT_BUILDER_H

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

// 03/03/09 Todd Coffey <tscoffe@sandia.gov>:
// This is a generic builder class that provides a validated parameter list and
// can build anything that can be constructed with a default constructor and
// accepts a parameter list through setParameterList (e.g. it derives from
// ParameterListAcceptor).
// Note the following:
// * The default object name is "Object" (this can be changed through setObjectName)
// * The default object type name is "Object Type" (this can be changed through setObjectTypeName)
// * The valid parameter list has a parameter named "Object Type" with a default value of "None"
// * The builder will create a null RCP if no factories have been set on it with setObjectFactory
// * A parameter list need not be set on the builder to call create, it will
// simply create the default factory which is either "None" if no factories
// have been set or it will be the last factory that was set.  
// * Setting a parameter list on the builder allows you to specify which object
// will be created by default.


namespace Teuchos {

template<class ObjectType>
class ObjectBuilder : virtual public ParameterListAcceptor
{
public:

  /** \brief . */
  ObjectBuilder();

  /** \brief . */
  ~ObjectBuilder();

  /** \brief Set the name of the object this will be a builder for, e.g. "Object". */
  void setObjectName(
      const std::string &objectName
      );

  /** \brief Set the name of the parameterlist selector, e.g. "Object Type". */
  void setObjectTypeName(
      const std::string &objectTypeName
      );

  /** \brief Set a new Object factory object. */
  void setObjectFactory(
    const RCP<const AbstractFactory<ObjectType> > &objectFactory,
    const std::string &objectFactoryName
    );

  /** \brief Get the name of the Object that will be created
   * on the next call to <tt>this->create()</tt>.
   */
  std::string getObjectName() const;

  /** \brief . */
  RCP<ObjectType> create(
    const std::string &objectName = ""
    ) const;
  
  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(const RCP<ParameterList> & paramList);
  
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  
  /** \brief. */
  RCP<const ParameterList> getParameterList() const;

  /** \brief. */
  RCP<const ParameterList> getValidParameters() const;
 
  //@}
  
private:

  // //////////////////////////////////////
  // Private types

  typedef RCP<const AbstractFactory<ObjectType > > object_fcty_t;

  // //////////////////////////////////////
  // Private data members

  RCP<ParameterList> paramList_;
  mutable RCP<const ParameterList> validParamList_;
  mutable RCP<const StringToIntegralParameterEntryValidator<int> > objectValidator_;

  std::string object_name_;
  std::string objectType_name_;

  Array<std::string> validObjectNames_;
  Array<object_fcty_t> objectArray_;
  std::string defaultObject_;

  // //////////////////////////////////////
  // Private member functions

  void initializeDefaults_();

};


// Nonmember constructors


template<class ObjectType>
RCP<ObjectBuilder<ObjectType> > objectBuilder()
{
  RCP<ObjectBuilder<ObjectType> > ob = rcp(new ObjectBuilder<ObjectType>() );
  return ob;
}


template<class ObjectType>
RCP<ObjectBuilder<ObjectType> >
objectBuilder(const std::string& objectName, const std::string& objectTypeName)
{
  RCP<ObjectBuilder<ObjectType> > ob = rcp(new ObjectBuilder<ObjectType>() );
  ob->setObjectName(objectName);
  ob->setObjectTypeName(objectTypeName);
  return ob;
}


template<class ObjectType>
ObjectBuilder<ObjectType>::ObjectBuilder()
{
  this->initializeDefaults_();
}


template<class ObjectType>
ObjectBuilder<ObjectType>::~ObjectBuilder()
{
#ifdef TEUCHOS_DEBUG
  // Validate that we read the parameters correctly!
  if(!is_null(paramList_)) {
    paramList_->validateParameters(*this->getValidParameters());
  }
#endif    
}


template<class ObjectType>
void ObjectBuilder<ObjectType>::setObjectFactory(
  const RCP<const AbstractFactory<ObjectType > > &objectFactory,
  const std::string &objectName
  )
{
  TEST_FOR_EXCEPT( objectName.length() == 0 );
  validObjectNames_.push_back(objectName);
  objectArray_.push_back(objectFactory);
  defaultObject_ = objectName;
  validParamList_ = null;
}


template<class ObjectType>
std::string
ObjectBuilder<ObjectType>::getObjectName() const
{
  if(is_null(validParamList_)) {
    this->getValidParameters();
  }
  // If the user has not specified a ParameterList, then use the ValidParameterList.
  RCP<ParameterList> pl = null;
  if (!is_null(paramList_)) {
    pl = paramList_;
  } else {
    pl = parameterList();
    pl->setParameters(*this->getValidParameters());
  }
  return objectValidator_->getStringValue(*pl, objectType_name_, defaultObject_);
}


template<class ObjectType>
void ObjectBuilder<ObjectType>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  if (!is_null(paramList)) {
    paramList->validateParameters(*this->getValidParameters());
    paramList_ = paramList;
  }
}


template<class ObjectType>
RCP<ParameterList>
ObjectBuilder<ObjectType>::getNonconstParameterList()
{
  return paramList_;
}


template<class ObjectType>
RCP<ParameterList>
ObjectBuilder<ObjectType>::unsetParameterList()
{
#ifdef TEUCHOS_DEBUG
  // Validate that we read the parameters correctly!
  if(!is_null(paramList_))
    paramList_->validateParameters(*this->getValidParameters());
#endif    
  RCP<ParameterList> _paramList = paramList_;
  paramList_ = null;
  return _paramList;
}


template<class ObjectType>
RCP<const ParameterList>
ObjectBuilder<ObjectType>::getParameterList() const
{
  return paramList_;
}


template<class ObjectType>
RCP<const ParameterList>
ObjectBuilder<ObjectType>::getValidParameters() const
{
  if(!validParamList_.get()) {
    RCP<ParameterList> validParamList = parameterList();
    // Object Types
    objectValidator_ = rcp(
      new StringToIntegralParameterEntryValidator<int>(
        validObjectNames_, objectType_name_
        )
      );
    validParamList->set(
      objectType_name_, defaultObject_,
      (std::string("Determines the type of " + object_name_ + " object that will be built.\n")
        + "The parameters for each " + objectType_name_ + " are specified in this sublist" 
        ).c_str(),
      objectValidator_
      );
    for( int i = 0; i < static_cast<int>(objectArray_.size()); ++i ) {
      const std::string
        &sname = validObjectNames_[i+1];
      const RCP<ObjectType >
        object = objectArray_[i]->create();
      validParamList->sublist(sname).setParameters(
        *object->getValidParameters()).disableRecursiveValidation();
    }
    validParamList_ = validParamList;
  }
  return validParamList_;
}


template<class ObjectType>
RCP<ObjectType >
ObjectBuilder<ObjectType>::create(
  const std::string &objectName
  ) const
{
  const std::string
    sname = ( objectName.length()
             ? objectName
             : this->getObjectName() );
  RCP<ObjectType> object = null; 
  // Get the index of this object factory (this will validate!)
  const int
    s_idx = objectValidator_->getIntegralValue(sname, objectType_name_);
  if (s_idx != 0) {
    // Create the uninitialized object
    object = objectArray_[s_idx-1]->create(); 
    TEST_FOR_EXCEPTION( is_null(object), std::logic_error,
        (std::string("Error!  ObjectBuilder attempted to create an object of type ") 
         + validObjectNames_[s_idx] + " and it came back as a null RCP!").c_str()
        );
    // Allow the user to not set a parameterlist (this requires copying the
    // parameters in the valid parameter list into a new parameter list:
    RCP<ParameterList> pl = null;
    if (is_null(paramList_)) {
      pl = parameterList();
      pl->setParameters(this->getValidParameters()->sublist(sname));
    } else {
#ifdef TEUCHOS_DEBUG
      // We're validating the parameter list here again because we're storing a
      // pointer to it and the user could have changed it.
      paramList_->validateParameters(*this->getValidParameters());
#endif // TEUCHOS_DEBUG
      pl = sublist(paramList_,sname);
    }
    // Now set the parameters for the object 
    object->setParameterList(pl);
  }
  return object;
}


template<class ObjectType>
void ObjectBuilder<ObjectType>::setObjectName(
    const std::string &objectName
    )
{
  TEST_FOR_EXCEPT(objectName.length() == 0);
  object_name_ = objectName;
  validParamList_ = null;
}


template<class ObjectType>
void ObjectBuilder<ObjectType>::setObjectTypeName(
    const std::string &objectTypeName
    )
{
  TEST_FOR_EXCEPT(objectTypeName.length() == 0);
  objectType_name_ = objectTypeName;
  validParamList_ = null;
}


template<class ObjectType>
void ObjectBuilder<ObjectType>::initializeDefaults_()
{

  object_name_ = "Object";
  objectType_name_ = "Object Type";

  defaultObject_ = "None";
  validObjectNames_.resize(0);
  validObjectNames_.push_back(defaultObject_);

}


} // namespace Teuchos


#endif //Teuchos_OBJECT_BUILDER_H
