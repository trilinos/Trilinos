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

