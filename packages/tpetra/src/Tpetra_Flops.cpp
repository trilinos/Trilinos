/*
27-May-2002 General cleanup. Changed method names to fit namingConvention(already done). Moved some code in from Tpetra_Flops.h
06-August-2002 Changed to images.
*/

namespace Tpetra
{
Flops::Flops(void) : flopCounter_(0.0)
{
}

Flops::Flops(const Flops& flops) : flopCounter_(0.0)
{
}

Flops::~Flops(void)  
{
  flopCounter_ = 0.0;
}

double Flops::flops() const 
{
  double tmp = flopCounter_; 
  flopCounter_ = 0.0; 
  return(tmp);
}

void Flops::resetFlops() 
{
  flopCounter_=0.0;
}

void Flops::updateFlops(int flops) const 
{
  flopCounter_ += (double) flops;
}

void Flops::updateFlops(long int flops) const 
{
  flopCounter_ += (double) flops;
}

void Flops::updateFlops(double flops) const
{
  flopCounter_ += flops;
}

void Flops::updateFlops(float flops) const 
{
  flopCounter_ +=(double) flops;
}

}  // namespace Tpetra
