#include "GModel.h"

namespace GAASP {

using namespace std;

// Default constructor
Lorenz_GModel::Lorenz_GModel()
{
	Initialize();
}
// Constructor for parameters
Lorenz_GModel::Lorenz_GModel(std::vector<double> p)
{
	for(int i=0; i<p.size(); i++)
	{
		params.push_back(p[i]);
	}
	sysDim = 3;
}

Lorenz_GModel::~Lorenz_GModel()
{
	Clear();
}

void Lorenz_GModel::Clear()
{
}

int Lorenz_GModel::calcDerivs(double *yin, double *yout, double t)
{
    // Parameters
    // params[0]	(Jeff's value = 10.0)
    // params[1]	(Jeff's value =	28.0)
    // params[2]	(Jeff's value =	8.0/3.0)

    yout[0] = -params[0] * yin[0] + params[0] * yin[1];
    yout[1] = params[1] * yin[0] - yin[1] - yin[0] * yin[2];
    yout[2] = -params[2] * yin[2] + yin[0] * yin[1];

	return(0);
}

void Lorenz_GModel::Initialize()
{
	params.push_back(10.0);
	params.push_back(28.0);
	params.push_back(8.0/3.0);
	sysDim = 3;
	return;
}

void Lorenz_GModel::updateParameters(std::vector<double> p)
{
	params.clear();
	for(int i=0; i<p.size(); i++)
	{
		params.push_back(p[i]);
	}
}

int Lorenz_GModel::getDim()
{
	return sysDim;
}

} // namespace GAASP

