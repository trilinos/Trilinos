/*Paul
27-May-2002 General cleanup. Changed method names to fit namingConvention. Moved some code in from Tpetra::Object.h. Commented out reportError method.
09-June-2002 Switched comment format.
06-August-2002 Changed to images (nothing changed). Also touched up a few things.
*/

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
int Tpetra::Object::reportError(const string message, int errorCode) const
{
#ifndef TPETRA_NO_ERROR_REPORTS
  // NOTE:  We are extracting a C-style string from Message because 
  //        the SGI compiler does not have a real string class with 
  //        the << operator.  Some day we should get rid of ".c_str()"
  cerr << endl << "Error in Tpetra Object with label:  " << label_ << endl
       << "Tpetra Error:  " << message.c_str() << "  Error Code:  " << errorCode << endl;
#endif
  return(errorCode);
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
