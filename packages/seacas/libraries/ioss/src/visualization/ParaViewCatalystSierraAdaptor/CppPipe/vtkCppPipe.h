#ifndef VTKCPPPIPE_H
#define VTKCPPPIPE_H

#include <vtkCPPipeline.h>
#include <string>

class vtkCPDataDescription;

class vtkCppPipe : public vtkCPPipeline
{
public:
  static vtkCppPipe* New();
  vtkTypeMacro(vtkCppPipe,vtkCPPipeline);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  virtual int Initialize();

  virtual int RequestDataDescription(vtkCPDataDescription* dataDescription);

  virtual int CoProcess(vtkCPDataDescription* dataDescription);

protected:
  vtkCppPipe();
  virtual ~vtkCppPipe();

private:
  vtkCppPipe(const vtkCppPipe&); // Not implemented
  void operator=(const vtkCppPipe&); // Not implemented
};
#endif
