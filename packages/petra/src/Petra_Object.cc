
#include "Petra_Object.h"


//=============================================================================
Petra_Object::Petra_Object(void) 
  : Label_(0)
{
  SetLabel("Petra::Object");
}
//=============================================================================
Petra_Object::Petra_Object(const char * const Label) 
  : Label_(0)
{
  SetLabel(Label);
}
//=============================================================================
Petra_Object::Petra_Object(const Petra_Object& Object)
  : Label_(0)
{
  SetLabel(Object.Label_);
}
//=============================================================================
void Petra_Object::Print(ostream & os) const {
  // os << Label_; // No need to print label, since ostream does it already
  return;
}
//=============================================================================
Petra_Object::~Petra_Object(void)  
{
  if (Label_!=0) delete Label_;
}
//=============================================================================
char * Petra_Object::Label(void) const {
  return(Label_);
}
//=============================================================================
void Petra_Object::SetLabel(const char * const Label) { 
  if (Label_!=0) delete Label_;
  Label_ = new char[strlen(Label)+1];
  strcpy(Label_,Label);
  return;
}
