#ifndef MLAPI_BASEOBJECT_H
#define MLAPI_BASEOBJECT_H

//! MLAPI: Default namespace for all ML API classes.
namespace MLAPI {

static int count = 0;

/*!
 * \class BaseObject
 *
 * \brief Basic class for MLAPI objects
 *
 * BaseObject is the basic class for all MLAPI objects. Currently, it 
 * contains the name of the object.
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Last updated on 07-Jan-05
 */
class BaseObject {

public:
  //! Constructor with empty name.
  BaseObject() 
  {
    Name_ = "obj_" + count;
    ++count;
  }

  //! Constructor with given name.
  BaseObject(const string& Name)
  {
    Name_ = Name;
  }

  //! Sets the name of this object to \c Name.
  void SetName(const string& Name)
  {
    Name_ = Name;
  }

  //! Returns the name of this object.
  const string& Name() const
  {
    return(Name_);
  }

private:
  //! Name of this object.
  string Name_;

};

} // namespace MLAPI

#endif
