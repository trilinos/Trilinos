#include "Tpetra_Object.hpp"

namespace Tpetra
{
//=============================================================================
Object::Object(int tracebackModeIn) : label_(0)
{
  setLabel("Tpetra::Object");
  tracebackMode = (tracebackModeIn != -1) ? tracebackModeIn : tracebackMode;
}
//=============================================================================
Object::Object(const char* const label, int tracebackModeIn) : label_(0)
{
  setLabel(label);
  tracebackMode = (tracebackModeIn != -1) ? tracebackModeIn : tracebackMode;
}
//=============================================================================
Object::Object(const Object& Obj) : label_(0)
{
  setLabel(Obj.label());
}
// Set TracebackMode value to default
int Object::tracebackMode(-1);

void Object::setTracebackMode(int tracebackModeValue)
{
  if (tracebackModeValue < 0)
    tracebackModeValue = 0;
  Object tempObject(tracebackModeValue);
}

int Object::getTracebackMode()
{
  int temp = Object::tracebackMode;
  if (temp == -1)
    temp = Tpetra_DefaultTracebackMode;
  return(temp);
}
//=============================================================================
void Object::print(ostream& os) const
{
  // os << label_; // No need to print label, since ostream does it already
}
//=============================================================================
Object::~Object()  
{
  if (label_!=0) {
    delete [] label_;
		label_ = 0;
	}
}
//=============================================================================
char* Object::label() const
{
  return(label_);
}
//=============================================================================
void Object::setLabel(const char* const label)
{ 
  if (label_ != 0)
    delete [] label_;
  label_ = new char[strlen(label) + 1];
  strcpy(label_, label);
}
} // namespace Tpetra
