#ifndef MLAPI_BASEOBJECT_H
#define MLAPI_BASEOBJECT_H

#define ML_THROW(str,val) { \
  std::cerr << "ERROR: In function/method " << __func__ << "()" << endl; \
  std::cerr << "ERROR: File " << __FILE__ << ", line " << __LINE__ << endl; \
  std::cerr << "ERROR: " << str << endl; \
  throw(val); \
  }

//! MLAPI: Default namespace for all ML API classes.
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
 * \date Last updated on 07-Jan-05
 */
class BaseObject {

public:
  //! Constructor with empty label.
  BaseObject() 
  {
    Label_ = "not-set";
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

  string toString(const int& x) const {
    char s[100];
    sprintf(s, "%d", x);
    return string(s);
  }

  string toString(const double& x) const {
    char s[100];
    sprintf(s, "%g", x);
    return string(s);
  }

private:
  //! Label of this object.
  string Label_;

};

std::ostream& operator << (std::ostream& os, const BaseObject& obj)
{
  return(obj.Print(os));
}

} // namespace MLAPI

#endif
