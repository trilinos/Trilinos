/*Paul
27-May-2002 General cleanup. Checked for newNamingConvention (already done). Moved some code here from Tpetra_CompObject.h
06-August-2002 Changed to images (nothing changed).
*/

// Constructor


namespace Tpetra
{
CompObject::CompObject()
{
  flopCounter_ = 0;
}

// Copy Constructor

CompObject::CompObject(const CompObject& source) : flopCounter_(source.flopCounter_)
{
}

// Destructor

CompObject::~CompObject()
{
  flopCounter_ = 0;
}

void CompObject::setFlopCounter(const Flops &flopCounter)
{
  flopCounter_= (Flops *)& flopCounter; 
  return;
}
  
void CompObject::setFlopCounter(const CompObject &compObject) 
{
  flopCounter_= (Flops *) (compObject.getFlopCounter()); 
  return;
}

void CompObject::unsetFlopCounter() 
{
  flopCounter_= 0; 
  return;
}
  
Flops * CompObject::getFlopCounter() const 
{
  return(flopCounter_);
}

void CompObject::resetFlops() const 
{
  if (flopCounter_!=0) 
    flopCounter_->resetFlops(); 
  return;
}

double CompObject::flops() const 
{
  if (flopCounter_!=0) 
    return(flopCounter_->flops()); 
  else return(0.0);
}

void CompObject::updateFlops(int flops) const 
{
  if (flopCounter_!=0) 
    flopCounter_->updateFlops(flops); 
  return;
}

void CompObject::updateFlops(long int flops) const 
{
  if (flopCounter_!=0) 
    flopCounter_->updateFlops(flops); 
  return;
}

void CompObject::updateFlops(double flops) const 
{
  if (flopCounter_!=0) 
    flopCounter_->updateFlops(flops); 
  return;
}

void CompObject::updateFlops(float flops) const 
{
  if (flopCounter_!=0) 
    flopCounter_->updateFlops(flops); 
  return;
}
} // namespace Tpetra
