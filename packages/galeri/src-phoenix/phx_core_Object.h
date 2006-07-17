// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef PHX_OBJECT_H
#define PHX_OBJECT_H

namespace phx {

namespace core {

class Object
{
  public:
    Object() 
    {
      label_ = "";
      ID_ = 0;
    }

    Object(const string& Label, const int ID = 0) 
    {
      setLabel(Label);
      setID(ID);
    }

    Object(const Object& rhs) 
    {
      setLabel(rhs.getLabel());
      setID(rhs.getID());
    }

    Object& operator=(const Object& rhs)
    {
      setLabel(rhs.getLabel());
      setID(rhs.getID());
      return(*this);
    }

    virtual ~Object() {}

    virtual string getLabel() const
    {
      return(label_);
    }

    virtual void setLabel(const string& label)
    {
      label_ = label;
    }

    virtual int getID() const
    {
      return(ID_);
    }

    virtual void setID(const int& ID)
    {
      ID_ = ID;
    }

    //! Print Object to an output stream
    virtual void print(ostream & os) const
    {
      return;
    }

  private:
    string label_;
    int ID_;

}; // class Object

inline ostream& operator<<(ostream& os, const Object& obj)
{
  os << obj.getLabel() << endl;
  obj.print(os);

  return(os);
}

} // namespace core

}; // namespace phx

#endif

