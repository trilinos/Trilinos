// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
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
    Object(const string& Label = "Galeri::core::Object", const int ID = 0) 
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
    virtual void setLabel(const string& label)
    {
      label_ = label;
    }

    //! Gets the label associated with \c this object.
    virtual string getLabel() const
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
    string label_;
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

