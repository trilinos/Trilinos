#ifndef __Clone_Factory_hpp__
#define __Clone_Factory_hpp__

#include "Teuchos_RCP.hpp"

namespace PB {

/** Base class for cloneable objects */
class Cloneable {
public:
   /** Function that clones this object.  This
     * is not neccessarily a copy, but it is an
     * object of the same dynamic type.
     *
     * \returns A pointer to an object of the same
     *          type is returned
     */
   virtual Teuchos::RCP<Cloneable> clone() const = 0;
};

/** Class that provides an easy way to create a cloneable
  * class. All that is required is that the user implements
  * a default constructor. This also serves as a model for
  * creating other AutoClone-like interfaces. 
  *
  * \note Use of this class assumes the BaseType class has
  *       a default constructor and does not override
  *       the clone method.
  */
template <class BaseType>
class AutoClone : public Cloneable, public BaseType {
public:
   /** Simple default constructor, calls the default
     * constructor of BaseType 
     */
   AutoClone() 
      : BaseType()
   { }
   
   /** Function that clones this object.  This
     * is not neccessarily a copy, but it is an
     * object of the same dynamic type.
     *
     * \returns A pointer to an object of the same
     *          type is returned
     */
   virtual Teuchos::RCP<Cloneable> clone() const
   { return Teuchos::rcp(new AutoClone<BaseType>()); }

private:
   // hidden
   AutoClone(const AutoClone<BaseType> & );
};

/** This class eases the construction of clone factories.  
  * It takes any Cloneable object and associates it with
  * a string. It will then perform the dynamic cast to
  * whatever type the user specifies using the CloneBaseType
  * template parameter. 
  *
  * \note A word of warning. This class does not provide compile
  *       time type safety. 
  */
template <class CloneBaseType>
class CloneFactory {
public:
   //! Default constructor
   CloneFactory() {}

   //! Copy constructor
   CloneFactory(const CloneFactory<CloneBaseType> & cf) : parentClones_(cf.parentClones_) {}

   /** Build a clone of the object associated with the string. This
     * object is automatically cast to the desired base type. This will
     * permit the easy use of the AutoClone class. 
     *
     * \note This method will throw an exception if the clone is not the right
     *       type (i.e. the dynamic cast to CloneBaseType fails). Also, this method
     *       returns null if the clone associated with the string can't be
     *       found.
     *
     * \param[in] str String associated with object to build.
     *
     * \returns A cloned object, correctly casted to CloneBaseType. If the
     *          string cannot be found then null is returned.
     */
   virtual Teuchos::RCP<CloneBaseType> build(const std::string & str) const
   { 
      std::map<std::string,Teuchos::RCP<const Cloneable> >::const_iterator itr 
            = parentClones_.find(str);
      if(itr==parentClones_.end()) return Teuchos::null;
      return Teuchos::rcp_dynamic_cast<CloneBaseType>(itr->second->clone(),true); 
   }

   /** Add a string associated clone to the factory. This object can be used
     * later to build a clone of itself. If this method is called twice with
     * the same string, the later clone will be maintained.
     *
     * \note The object to be cloned is stored until the CloneFactory is destroyed.
     *
     * \param[in] str String associated with this object.
     * \param[in] clone Object to be cloned.
     */
   virtual void addClone(const std::string & str,const Teuchos::RCP<Cloneable> & clone)
   { parentClones_[str] = clone; }

   //! Return the number of clones stored in this factory
   virtual int cloneCount() const
   { return parentClones_.size(); }

protected:
   //! stores the clonable objects
   std::map<std::string,Teuchos::RCP<const Cloneable> > parentClones_;
};

};

#endif
