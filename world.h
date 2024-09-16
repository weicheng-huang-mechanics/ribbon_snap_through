#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include elastic rod class
#include "elasticRod.h"

// include force classes
#include "elasticStretchingForce.h"
#include "elasticRibbonForce.h"
#include "externalGravityForce.h"
#include "inertialForce.h"

// include external force
#include "dampingForce.h"

// include time stepper
#include "timeStepper.h"

// include input file and option
#include "setInput.h"

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	void setRodStepper();
	void updateTimeStep();
	int simulationRunning();
	int numPoints();
	double getScaledCoordinate(int i);
	double getBoundaryCoordination_left(int i);
	double getBoundaryCoordination_right(int i);
	double getCurrentTime();
	double getTotalTime();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);
		
private:

	// Physical parameters
	double RodLength;
	double helixradius, helixpitch;
	double rodRadius;
	int numVertices;
	double youngM;
	double Poisson;
	double shearM;
	double deltaTime;
	double totalTime;
	double density;
	Vector3d gVector;
	double viscosity;
	double speed;
	double width;
	double inputShear;
	double inputCompress;
	double inputAngle;

	double inputDeltaTime;

	double totalCompress;
	double totalShear;
	double totalAngle;

	double plateLength;
	double plateWidth;
	double plateDeltaEdge;

	// Viscous drag coefficients
    double eta_per, eta_par;
    
	double tol, stol;
	int maxIter; // maximum number of iterations
	double characteristicForce;
	double forceTol;
	
	// Geometry
	MatrixXd vertices;
	VectorXd theta;
	
	// Rod
	elasticRod *rod;
	
	// set up the time stepper
	timeStepper *stepper;
	double *totalForce;
	double currentTime;
	
	// declare the forces
	elasticStretchingForce *m_stretchForce;
	elasticRibbonForce *m_elasticRibbonForce;
	inertialForce *m_inertialForce;	
	externalGravityForce *m_gravityForce;
	dampingForce *m_dampingForce;
	
	int Nstep;
	int timeStep;
	int iter;

	void rodGeometry();
	void rodBoundaryCondition();
    
	bool render; // should the OpenGL rendering be included?
	bool saveData; // should data be written to a file?

	VectorXd reactionForce;
	void computeReactionForce();

	int ifResetMap;

	int dataStep;
};

#endif
