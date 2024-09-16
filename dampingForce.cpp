#include "dampingForce.h"
#include <iostream>

dampingForce::dampingForce(elasticRod &m_rod, timeStepper &m_stepper, 
	double m_viscosity, double m_eta_per, double m_eta_par)
{
	rod = &m_rod;
	stepper = &m_stepper;
	viscosity = m_viscosity;
	dt = rod->dt;
	eta_per = m_eta_per;
	eta_par = m_eta_par;
		
	Id3<<1,0,0,
		0,1,0,
        0,0,1;

    forceVec = VectorXd::Zero(rod->ndof);
}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd(double m_currentTime)
{
	forceVec = VectorXd::Zero(rod->ndof);

	for (int i = 0; i < rod->ndof; i++)
	{
		double localForce = - viscosity * rod->massArray[i] * (rod->x(i) - rod->x0(i)) / rod->dt;

		//double localForce = - viscosity * (rod->x(i) - rod->x0(i)) / rod->dt;

		if (m_currentTime < 45.0)
		{
			localForce = - 1.0 * rod->massArray[i] * (rod->x(i) - rod->x0(i)) / rod->dt;
		}

		forceVec(i) = forceVec(i) + localForce;

		stepper->addForce(i, - localForce); // subtracting external force
	}

}

void dampingForce::computeJd(double m_currentTime)
{
	for (int i = 0; i < rod->ndof; i++)
	{
		double localJacobian = - viscosity * rod->massArray[i] / rod->dt;

		//double localJacobian = - viscosity / rod->dt;

		if (m_currentTime < 45.0)
		{
			localJacobian = - 1.0 * rod->massArray[i] / rod->dt;
		}

		stepper->addJacobian(i, i, - localJacobian);
	}

}
