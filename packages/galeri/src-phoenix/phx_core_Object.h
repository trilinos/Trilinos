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
      tag_ = 0;
    }

    Object(const string& Label, const int Tag = 0) 
    {
      setLabel(Label);
      setTag(Tag);
    }

    Object(const Object& rhs) 
    {
      setLabel(rhs.getLabel());
      setTag(rhs.getTag());
    }

    Object& operator=(const Object& rhs)
    {
      setLabel(rhs.getLabel());
      setTag(rhs.getTag());
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

    virtual int getTag() const
    {
      return(tag_);
    }

    virtual void setTag(const int& tag)
    {
      tag_ = tag;
    }

    //! Print Object to an output stream
    virtual void print(ostream & os) const
    {
      return;
    }

  private:
    string label_;
    int tag_;

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

