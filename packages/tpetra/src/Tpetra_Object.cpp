/*Paul
27-May-2002 General cleanup. Changed method names to fit namingConvention. Moved some code in from Tpetra::Object.h. Commented out reportError method.
09-June-2002 Switched comment format.
*/

namespace Tpetra
{
//=============================================================================
Object::Object(int TracebackModeIn) : Label_(0)
{
  setLabel("Tpetra::Object");
  TracebackMode = (TracebackModeIn != -1) ? TracebackModeIn : TracebackMode;
}
//=============================================================================
Object::Object(const char* const Label, int TracebackModeIn) : Label_(0)
{
  setLabel(Label);
  TracebackMode = (TracebackModeIn != -1) ? TracebackModeIn : TracebackMode;
}
//=============================================================================
Object::Object(const Object& Obj) : Label_(0)
{
  setLabel(Obj.Label_);
}
// Set TracebackMode value to default
int Object::TracebackMode(-1);

void Object::setTracebackMode(int TracebackModeValue)
{
  if (TracebackModeValue < 0)
    TracebackModeValue = 0;
  Object TempObject(TracebackModeValue);
}

int Object::getTracebackMode()
{
  int temp = Object::TracebackMode;
  if (temp==-1)
    temp = Tpetra_DefaultTracebackMode;
  return(temp);
}
//=============================================================================
void Object::print(ostream& os) const
{
  // os << Label_; // No need to print label, since ostream does it already
  return;
}
//=============================================================================
int Tpetra::Object::reportError(const string Message, int ErrorCode) const
{
#ifndef TPETRA_NO_ERROR_REPORTS
  // NOTE:  We are extracting a C-style string from Message because 
  //        the SGI compiler does not have a real string class with 
  //        the << operator.  Some day we should get rid of ".c_str()"
  cerr << endl << "Error in Tpetra Object with label:  " << Label_ << endl
       << "Tpetra Error:  " << Message.c_str() << "  Error Code:  " << ErrorCode << endl;
#endif
  return(ErrorCode);
}
//=============================================================================
Object::~Object()  
{
  if (Label_!=0)
    delete [] Label_;
}
//=============================================================================
char * Object::label() const
{
  return(Label_);
}
//=============================================================================
void Object::setLabel(const char* const Label)
{ 
  if (Label_!=0)
    delete [] Label_;
  Label_ = new char[strlen(Label)+1];
  strcpy(Label_,Label);
  return;
}
} // namespace Tpetra
