// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

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

      //! Informs entry that it is now stored in the library
      /*!
       * This is used primarily for informing the entry on how to delete
       * itself when deleting the library.
       */
      virtual void setIsInLibrary() { deleteFctr = true; }

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
      ~DefaultFunctor() {}

      //! Set parameter using object and data member pointer
      void set(const ValueType& value) { obj.*obj_val_ptr = value; }

    protected:

      //! Reference to object to set parameter it
      ObjectType& obj;

      //! Pointer to data member in \em ObjectType of type \em ValueType
      ValueType ObjectType::* obj_val_ptr;
    };
  }
}

#endif
