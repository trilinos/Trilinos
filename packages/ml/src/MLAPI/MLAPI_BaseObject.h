#ifndef MLAPI_BASEOBJECT_H
#define MLAPI_BASEOBJECT_H

namespace MLAPI {

static int count = 0;

class BaseObject {

public:
  BaseObject() 
  {
    Name_ = "obj_" + count;
    ++count;
  }

  BaseObject(const string& Name)
  {
    Name_ = Name;
  }

  void SetName(const string& Name)
  {
    Name_ = Name;
  }

  const string& Name() const
  {
    return(Name_);
  }

private:
  string Name_;

};

} // namespace MLAPI

#endif
