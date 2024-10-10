// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_OBJECT_H
#define GALERI_OBJECT_H

/*! 
 * \file Galeri_core_Object.h
 *
 * \class core::Object
 *
 * \brief Basic object for Galeri/pfem.
 *
 * \author Marzio Sala, ETH.
 *
 * \date Last modified on Aug-06.
 */
namespace Galeri {

namespace core {

/*! \brief Basic class for all Galeri/pfem objects.
 */
class Object
{
  public:
    //! @{ ctor, dtor, and operator =
    //! Constructor with specified label and ID.
    Object(const std::string& Label = "Galeri::core::Object", const int ID = 0) 
    {
      setLabel(Label);
      setID(ID);
    }

    //! Copy constructor.
    Object(const Object& rhs) 
    {
      setLabel(rhs.getLabel());
      setID(rhs.getID());
    }

    //! Copies the object from \c rhs.
    Object& operator=(const Object& rhs)
    {
      setLabel(rhs.getLabel());
      setID(rhs.getID());
      return(*this);
    }

    //! Virtual dtor.
    virtual ~Object() {}

    //! @}
    //! @{ Set and Get Methods
    
    //! Sets the label associated with \c this object.
    virtual void setLabel(const std::string& label)
    {
      label_ = label;
    }

    //! Gets the label associated with \c this object.
    virtual std::string getLabel() const
    {
      return(label_);
    }

    //! Sets the ID associated with \c this object.
    virtual void setID(const int& ID)
    {
      ID_ = ID;
    }

    //! Gets the ID associated with \c this object.
    virtual int getID() const
    {
      return(ID_);
    }

    //! @}
    //! @{ Other Methods
    
    /*! \brief Prints Object to the specified output stream.
     *
     * To be customized by derived classes. Operator << uses
     * this method. Default implementation prints the object label
     * to \c os.
     */
    virtual void print(ostream & os) const
    {
      os << endl;
      os << "** Object label: " << getLabel() << endl;
      os << "** Object ID: " << getID() << endl;
      os << endl;
      return;
    }

    //! @}
    //! @{ Private Variables and Methods.
    
  private:
    //! Label associated to \c this object.
    std::string label_;
    //! Integer-typed ID associated to \c this object.
    int ID_;

    //! @}
}; // class Object

inline ostream& operator<<(ostream& os, const Object& obj)
{
  obj.print(os);

  return(os);
}

} // namespace core

}; // namespace Galeri

#endif
