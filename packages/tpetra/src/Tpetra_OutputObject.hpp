// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef TPETRA_OUTPUTOBJECT_HPP
#define TPETRA_OUTPUTOBJECT_HPP

#include "Tpetra_ConfigDefs.hpp" // for iostream and string
#include "Tpetra_OutputManager.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace Tpetra
{

//! Tpetra::OutputObject:  The base Tpetra class.
/*! The OutputObject class provides capabilities common to all Tpetra objects,
    such as a label that identifies an object instance, constant definitions,
    enum types.
*/

class OutputObject
{
  public:
  //@{ \name Constructors/destructor.
  //! OutputObject Constructor.
  /*! OutputObject is the base class any class in Tpetra that wants to overload ostream<<.  This class is seldom
      used explictly.
  */
  OutputObject() 
    : om_(Teuchos::null) {
    om_ = Teuchos::rcp(new OutputManager());
  }

  //! OutputObject Constructor passing in an existing output manager.
  OutputObject(const Teuchos::RefCountPtr<OutputManager> &om)
    : om_(om){}

  //! OutputObject Copy Constructor.
  /*! Makes an exact copy of an existing OutputObject instance.
  */
  OutputObject(OutputObject const& obj)
    : om_(obj.om_){}

  //! OutputObject Destructor.
  /*! Completely deletes an OutputObject object.  
  */
  virtual ~OutputObject(){}
	
  //@}
  
  //@{ \name Attribute set/get methods.

  //! Set the output manager we should use.
  /*! Defines the label used to describe the \e this object.  
  */
  virtual void setOutputManager(const Teuchos::RefCountPtr<OutputManager> &om){om_ = om;}

  //! Set the output manager to that of another OutputObject.
  void setOutputManager(const OutputObject &outputObject) {om_= outputObject.getOutputManager(); return;}
  
  //! Set the output manager to default: no manager.
  void unsetOutputManager() {om_=Teuchos::null; return;}
  
  //! OutputObject Label access funtion.
  /*! Returns the string used to define this object.  
  */
  Teuchos::RefCountPtr<OutputManager> getOutputManager() const {return(om_);}

	
  //@}

  //@{ \name Virtual print methods.  These should be overloaded by the derived class.

  //! Print signature of an object to an ostream.
  /*! Prints a presumably small amount of information to an ostream.  This information
      should be sufficient to uniquely identify an object or set of identical objects.
      \param os (InOut) The ostream to which information will be printed.
  */
  virtual void printSignature(ostream& os) const {os << "Replace this version of printSignature with your own." << endl;} // Default is to print nothing
	
  //! Print summary contents of an object to an ostream.
  /*! Prints the contents of the attributes of an object.  The output from this method
      should be a fairly complete printout of object contents in a human-readable form.
      \param os (InOut) The ostream to which information will be printed.
  */
  virtual void printSummary(ostream& os) const {os << "Replace this version of printSummary with your own." << endl;} // Default is to print nothing

  //! Print watch results of an object to an ostream.
  /*! Prints information about an object when the state of an object changes.
      This method should print out only the information about what is changing
      about the state of the object.
      \param os (InOut) The ostream to which information will be printed.
      \warning Calls to this method should normally be surrounded by #ifdef HAVE_DEBUG...#endif since performance will be potential be impacted.
  */
  virtual void printWatch(ostream& os) const {} // Default is to print nothing

  //! New print method that is aware of Output manager settings.
  /*! Prints information about an object based on the settings of the Output Manager associated with this object.
      \param os (InOut) The ostream to which information will be printed.
  */
  virtual void newPrint(ostream& os) const {

    Teuchos::RefCountPtr<OutputManager> om = getOutputManager();
    
    if (om->isVerbosityAndPrint(Tpetra::Signature)) printSignature(os);
    if (om->isVerbosityAndPrint(Tpetra::Summary)) printSummary(os);
    return;
  }

private:
	
  
  //! Assignment operator (declared but not defined, do not use)
  OutputObject& operator = (OutputObject const& source);
  
  Teuchos::RefCountPtr<OutputManager> om_;

}; // class OutputObject

#ifdef HAVE_NEW_OUTPUT
inline ostream& operator<<(ostream& os, Tpetra::OutputObject const& outputObject) {
  
  Teuchos::RefCountPtr<OutputManager> om = outputObject.getOutputManager();
  
  if (om->isVerbosityAndPrint(Tpetra::Signature)) outputObject.printSignature(os);
  if (om->isVerbosityAndPrint(Tpetra::Summary)) outputObject.printSummary(os);
  return(os);
}
#endif
} // namespace Tpetra


#endif // TPETRA_OUTPUTOBJECT_HPP
