#ifndef MLAPI_BASEOBJECT_H
#define MLAPI_BASEOBJECT_H

namespace MLAPI {

/*!
 * \class BaseObject
 *
 * \brief Basic class for MLAPI objects
 *
 * BaseObject is the basic class for all MLAPI objects. Currently, it 
 * contains the label of the object.
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Last updated on Feb-05.
 */
class BaseObject {

public:
  //! Constructor with empty label.
  BaseObject() 
  {
    Label_ = "obj_" + GetString(count);
    ++count;
  }

  //! Constructor with given Label.
  BaseObject(const string& Label)
  {
    Label_ = Label;
  }

  //! Destructor.
  virtual ~BaseObject() {};

  //! Sets the Label of this object to \c Label.
  void SetLabel(const string& Label)
  {
    Label_ = Label;
  }

  //! Returns the Label of this object.
  const string& GetLabel() const
  {
    return(Label_);
  }

  //! Prints information on stream.
  virtual std::ostream& Print(std::ostream& os, 
                              const bool Verbose = true) const = 0;

private:
  //! Label of this object.
  string Label_;

  static int count;
};

int BaseObject::count = 0;

std::ostream& operator<< (std::ostream& os, const BaseObject& obj)
{
  return(obj.Print(os));
}

} // namespace MLAPI

#endif
