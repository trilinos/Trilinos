#include "NLS_PetraVector.H"

NLS_PetraVector::NLS_PetraVector(const Petra_RDP_Vector& copyFrom, bool doCopyEntries)
{
  if (doCopyEntries) {
    // deep copy
    petraVec = new Petra_RDP_Vector(copyFrom); 
  }
  else {
    // copy map and fill with zeros
    petraVec = new Petra_RDP_Vector(copyFrom.Map()); 
  }

  doDeletePetraVec = true;
}

NLS_PetraVector::NLS_PetraVector(Petra_RDP_Vector& pointTo)
{
  petraVec = &pointTo;		// copy pointer only
  doDeletePetraVec = false;	// do not delete when this is deleted
}

NLS_PetraVector::~NLS_PetraVector()
{
  if (doDeletePetraVec)
    delete petraVec;
  petraVec = NULL;
}

NLS_PetraVector& NLS_PetraVector::operator=(const NLS_PetraVector& copyFrom)
{
  if (petraVec == NULL) {
    // If petraVec is empty, fill it...
    petraVec = new Petra_RDP_Vector(*(copyFrom.petraVec)); // deep copy
    doDeletePetraVec = true;
  }

  else {
    // Otherwise, copy into existing petraVec
    int errcode = petraVec->Update(1.0, *(copyFrom.petraVec), 0.0);
    if (errcode != 0) 
      cerr << "Error in NLS_Petra_Vec::operator=!" << endl;
  }
  
  return *this;
}

void NLS_PetraVector::random()
{
  petraVec->Random();
}

void NLS_PetraVector::putScalar(double alpha)
{
  petraVec->PutScalar(alpha);
}

void NLS_PetraVector::update(double scalarA, 
			     const NLS_Vector& y, 
			     double scalar)
{
  throw;
}

void NLS_PetraVector::update(double scalarA, 
			     const NLS_PetraVector& y, 
			     double scalar)
{
  petraVec->Update(scalarA, *(y.petraVec), scalar);
}

void NLS_PetraVector::scale(double alpha)
{
  petraVec->Scale(alpha);
}

double NLS_PetraVector::dot(const NLS_Vector& y) const
{
  throw;
}

double NLS_PetraVector::dot(const NLS_PetraVector& y) const
{
  double dot;
  petraVec->Dot(*(y.petraVec), &dot);
  return dot;
}

double NLS_PetraVector::normInf() const
{
  double norm;
  petraVec->NormInf(&norm);
  return norm;
}

double NLS_PetraVector::norm1() const
{
  double norm;
  petraVec->Norm1(&norm);
  return norm;
}

double NLS_PetraVector::norm2() const
{
  double norm;
  petraVec->Norm2(&norm);
  return norm;
}

double NLS_PetraVector::norm2(const NLS_Vector& weights) const
{
  throw;
}

double NLS_PetraVector::norm2(const NLS_PetraVector& weights) const
{
  double norm;
  petraVec->NormWeighted(*(weights.petraVec), &norm);
  return norm;
}

void NLS_PetraVector::Abs(const NLS_Vector& base)
{
  throw; 
}

void NLS_PetraVector::Abs(const NLS_PetraVector& base)
{
  petraVec->Abs(*(base.petraVec)); 
}

double NLS_PetraVector::minValue() const
{
  double minVal;
  petraVec->MinValue(&minVal);
  return minVal;
}

double NLS_PetraVector::maxValue() const
{
  double maxVal;
  petraVec->MaxValue(&maxVal);
  return maxVal;
}

double NLS_PetraVector::meanValue() const
{
  double avgVal;
  petraVec->MeanValue(&avgVal);
  return avgVal;
}



