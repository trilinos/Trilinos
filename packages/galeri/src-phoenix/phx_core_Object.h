#ifndef PHX_OBJECT_H
#define PHX_OBJECT_H

namespace phx {

namespace core {

class Object
{
  public:
    Object() 
    {
      Label_ = "";
      Tag_ = 0;
      Modifiable_ = true;
    }

    Object(const string& Label, const int Tag = 0) 
    {
      setLabel(Label);
      setTag(Tag);
      setModifiable(true);
    }

    Object(const Object& rhs) 
    {
      setLabel(rhs.getLabel());
      setTag(rhs.getTag());
      setModifiable(rhs.getModifiable());
    }

    Object& operator=(const Object& rhs)
    {
      setLabel(rhs.getLabel());
      setTag(rhs.getTag());
      setModifiable(rhs.getModifiable());
      return(*this);
    }

    virtual ~Object() {}

    virtual string getLabel() const
    {
      return(Label_);
    }

    virtual void setLabel(const string& Label)
    {
      Label_ = Label;
    }

    virtual int getTag() const
    {
      return(Tag_);
    }

    virtual void setTag(const int& Tag)
    {
      Tag_ = Tag;
    }

    virtual bool getModifiable() const
    {
      return(Modifiable_);
    }

    virtual void setModifiable(const bool& Modifiable)
    {
      Modifiable_ = Modifiable;
    }

    //! Print Object to an output stream
    virtual void Print(ostream & os) const
    {
      return;
    }

  private:
    string Label_;
    int Tag_;
    bool Modifiable_;

}; // class Object

inline ostream& operator<<(ostream& os, const Object& obj)
{
  os << obj.getLabel() << endl;
  obj.Print(os);

  return(os);
}

} // namespace core

}; // namespace phx

#endif

