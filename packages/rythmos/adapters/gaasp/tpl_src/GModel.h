#ifndef INC_Lorenz_GModel_h
#define INC_Lorenz_GModel_h

#include "GModelBase.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <limits>
#include <cmath>

namespace GAASP {

class Lorenz_GModel:virtual public GModelBase
{
public:
	// Constructor and destructor
	Lorenz_GModel();								//	Default constructor
	Lorenz_GModel(std::vector<double>);			//	Constructor with parameters passed in
	~Lorenz_GModel();								//	Destructor

	// Functions
	void Clear();
	int calcDerivs(double *yin, double *yout, double t);
	void updateParameters(std::vector<double>);
	int getDim();

	// Data
	int sysDim;

private:
	// Data

	// Parameters used by the model
	std::vector<double> params;

	// Functions
	void Initialize();

};

} // namespace GAASP

#endif // INC_Lorenz_GModel_h
