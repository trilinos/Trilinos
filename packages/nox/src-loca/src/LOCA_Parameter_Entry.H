// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_PARAMETER_ENTRY_H
#define LOCA_PARAMETER_ENTRY_H

namespace LOCA {

  namespace Parameter {

    //! Abstract interface for all entries in LOCA::Parameter::Library
    /*!
     * This class doesn't have much of an interface and really only serves
     * the purpose of having a common parent class for parameter entries
     * of all value types.
     */
    class AbstractEntry {
    public:

      //! Default contructor
      AbstractEntry() {}

      //! Destructor
      virtual ~AbstractEntry() {}
    };

    //! Parameter entry interface class templated on \em ValueType
    /*!
     * This class provides the interface that all parameter entry classes
     * should implement.  It is templated on the \em ValueType, which is the
     * type that the underlying parameter is stored as.
     */
    template <class ValueType>
    class Entry : public AbstractEntry {
    public:

      //! Default constructor
      Entry() {}

      //! Destructor
      virtual ~Entry() {}

      //! Set parameter this object represents to \em value
      virtual void setValue(const ValueType& value) = 0;

      //! Get parameter value this object represents
      virtual ValueType getValue() const = 0;

      //! Informs entry that it is now stored in the library
      /*!
       * This is used primarily for informing the entry on how to delete
       * itself when deleting the library.
       */
      virtual void setIsInLibrary() = 0;
    };

    //! Standard parameter entry class using a function object
    /*!
     * This is the standard parameter entry class that uses a function object
     * to actually set/retrieve parameter values.  The nice thing about using
     * a function object is it allows one to set parameters that don't actually
     * exist in the code, for example, setting a dimensionless group value
     * by modifiying a number of physical parameters.  By supplying an
     * appropriate
     * function object, this class should suffice for setting/retrieving
     * parameter values in nearly all cases.
     *
     * The constructor takes a pointer to the supplied function object.  It is
     * assumed that this class then owns that pointer, and in particular, calls
     * delete in the destructor if the entry is successfully added to the
     * library.  It does not delete the function object otherwise.
     */
    template <class FunctorType, class ValueType>
    class StandardEntry :
      public Entry<ValueType> {

    public:

      //! Constructor
      StandardEntry(FunctorType* fctr) : deleteFctr(false), func(fctr) {}

      //! Destructor
      virtual ~StandardEntry() { if (deleteFctr) delete func; }

      //! Set parameter this object represents to \em value
      virtual void setValue(const ValueType& value) { func->set(value); }

      //! Get parameter value this object represents
      virtual ValueType getValue() const { return func->get(); }

      //! Informs entry that it is now stored in the library
      /*!
       * This is used primarily for informing the entry on how to delete
       * itself when deleting the library.
       */
      virtual void setIsInLibrary() { deleteFctr = true; }

    private:

      //! Private to prohibit copying
      StandardEntry(const StandardEntry&);

      //! Private to prohibit copying
      StandardEntry& operator=(const StandardEntry&);

    protected:

      //! Flag indicating whether to delete function object in destructor
      bool deleteFctr;

      //! Stores pointer to function object
      FunctorType* func;
    };

    /*!
     * \brief Default function object for setting a single parameter in a
     * single object using a data member pointer
     */
    /*!
     * The constructor takes a reference to an object \em object of type\
     * \em ObjectType and a pointer \em object_val_ptr to a data member of
     * class of \em ObjectType of type \em ValueType.  The parameter is set to
     * \em value via
     * \code
     *    object.*object_val_ptr = value
     * \endcode
     */
    template <class ObjectType, class ValueType>
    class DefaultFunctor {
    public:

      //! Constructor
      /*!
       * \em object is a reference to the object to set the parameter in, and
       * \em object_val_ptr is a pointer to a data member of type
       * \em ValueType of that class.
       */
      DefaultFunctor(ObjectType& object,
             ValueType ObjectType::* object_val_ptr) :
    obj(object), obj_val_ptr(object_val_ptr) {}

      //! Destructor
      virtual ~DefaultFunctor() {}

      //! Set parameter using object and data member pointer
      virtual void set(const ValueType& value) { obj.*obj_val_ptr = value; }

       //! Get parameter value this object represents
      virtual ValueType get() const { return obj.*obj_val_ptr; }

    private:

      //! Private to prohibit copying
      DefaultFunctor(const DefaultFunctor&);

      //! Private to prohibit copying
      DefaultFunctor& operator=(const DefaultFunctor&);

    protected:

      //! Reference to object to set parameter it
      ObjectType& obj;

      //! Pointer to data member in \em ObjectType of type \em ValueType
      ValueType ObjectType::* obj_val_ptr;
    };
  }
}

#endif
