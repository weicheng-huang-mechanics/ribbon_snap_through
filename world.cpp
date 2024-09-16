#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean

	// Physical parameters
	RodLength = m_inputData.GetScalarOpt("RodLength");      // meter
    helixradius = m_inputData.GetScalarOpt("helixradius");  // meter
    gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
    maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
    helixpitch = m_inputData.GetScalarOpt("helixpitch");    // meter
	rodRadius = m_inputData.GetScalarOpt("rodRadius");      // meter
	numVertices = m_inputData.GetIntOpt("numVertices");     // int_num
	youngM = m_inputData.GetScalarOpt("youngM");            // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");          // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");          // kg/m^3
	viscosity = m_inputData.GetScalarOpt("viscosity");      // viscosity in Pa-s
	speed = m_inputData.GetScalarOpt("speed");
	width = m_inputData.GetScalarOpt("width");
	inputShear = m_inputData.GetScalarOpt("inputShear");
	inputCompress = m_inputData.GetScalarOpt("inputCompress");
	inputAngle= m_inputData.GetScalarOpt("inputAngle");
	dataStep = m_inputData.GetIntOpt("dataStep");

	inputDeltaTime = deltaTime;
	
	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus
	
	// Viscous drag coefficients using Resistive Force Theory
	eta_per = 4.0*M_PI*viscosity/( log(2*helixpitch/rodRadius) + 0.5);
    eta_par = 2.0*M_PI*viscosity/( log(2*helixpitch/rodRadius) - 0.5);

    totalShear = 0.0;

    totalCompress = 0.0;

    inputShear = inputShear * RodLength;

    inputCompress = inputCompress * RodLength;

    totalAngle = 0.0;

    ifResetMap = 0;
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	time_t current_time = time(0);

	// Open an input file named after the current time
	ostringstream name;
	name.precision(6);
	name << fixed;

    name << "datafiles/simDER_1";
    name << "_compress_" << inputCompress / RodLength;
    name << "_shear_" << inputShear / RodLength;
    name << "_angle_" << inputAngle;
    name << "_speed_" << speed;
    name << "_viscosity_" << viscosity;
    name << ".txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
		return;

	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
		return;

	if ( timeStep % dataStep != 0)
	{
		return;
	}

	if (currentTime >= 50.0)
	//if (timeStep == Nstep)
	{
		Vector3d xMid = rod->getVertex( (rod->nv - 1) / 2);

		outfile << totalAngle << " " << xMid(2) << endl;

		//cout << totalAngle << " " << xMid(2) << endl;

		//outfile << xMid(2) << endl;
	}
}

void world::setRodStepper()
{
	// Set up geometry
	rodGeometry();	

	// Create the rod 
	rod = new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
		youngM, shearM, RodLength, theta, width);

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = pow(rodRadius, 3) * width / 12.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;
	
	// Set up boundary condition
	rodBoundaryCondition();
	
	// setup the rod so that all the relevant variables are populated
	rod->setup();
	// End of rod setup
	
	// set up the time stepper
	stepper = new timeStepper(*rod);
	totalForce = stepper->getForce();

	// declare the forces
	m_stretchForce = new elasticStretchingForce(*rod, *stepper);
	m_elasticRibbonForce = new elasticRibbonForce(*rod, *stepper);
	m_inertialForce = new inertialForce(*rod, *stepper);	
	m_gravityForce = new externalGravityForce(*rod, *stepper, gVector);
	m_dampingForce = new dampingForce(*rod, *stepper, viscosity, eta_per, eta_par);

	Nstep = totalTime/deltaTime;

	// Allocate every thing to prepare for the first iteration
	rod->updateTimeStep();
	
	timeStep = 0;
	currentTime = 0.0;

	reactionForce = VectorXd::Zero(rod->ndof);
}

// Setup geometry
void world::rodGeometry()
{
    double deltaLength = RodLength / (numVertices - 1);

    vertices = MatrixXd(numVertices, 3);

    for (int i = 0; i < numVertices; i++)
    {
    	vertices(i, 0) = 0.0;
    	vertices(i, 1) = deltaLength * i - RodLength / 2;
    	vertices(i, 2) = 0.0;
    } 

    // initial theta should be zeros
    theta = VectorXd::Zero(numVertices - 1);
}

void world::rodBoundaryCondition()
{
	// Apply boundary condition
	rod->setVertexBoundaryCondition(rod->getVertex(0), 0);
	rod->setThetaBoundaryCondition(rod->getTheta(0), 0);
	rod->setVertexBoundaryCondition(rod->getVertex(1), 1);

	rod->setVertexBoundaryCondition(rod->getVertex(rod->nv - 2), rod->nv - 2);
	rod->setThetaBoundaryCondition(rod->getTheta(rod->nv - 2), rod->nv - 2);
	rod->setVertexBoundaryCondition(rod->getVertex(rod->nv - 1), rod->nv - 1);
}
	

void world::updateTimeStep()
{
	if (currentTime < 20.0)
	{
		deltaTime = 5e-3;
		rod->dt = 5e-3;
	}
	else if (currentTime < 50.0)
	{
		deltaTime = 5e-3;
		rod->dt = 5e-3;
	}
	else
	{
		deltaTime = inputDeltaTime;
		rod->dt = inputDeltaTime;
	}

	rod->setThetaBoundaryCondition(rod->getTheta(0), 0);
	rod->setThetaBoundaryCondition(rod->getTheta(rod->nv - 2), rod->nv - 2);

	if (currentTime >= 1.0 && currentTime < 20.0 && totalCompress <= inputCompress)
	{
		Vector3d x_0 = rod->getVertex(0);
		Vector3d x_1 = rod->getVertex(1);

		x_0(1) = x_0(1) + deltaTime * 0.02;
		x_1(1) = x_1(1) + deltaTime * 0.02;

		rod->setVertexBoundaryCondition(x_0, 0);
		rod->setVertexBoundaryCondition(x_1, 1);

		Vector3d x_2 = rod->getVertex(rod->nv - 2);
		Vector3d x_3 = rod->getVertex(rod->nv - 1);

		x_2(1) = x_2(1) - deltaTime * 0.02;
		x_3(1) = x_3(1) - deltaTime * 0.02;	

		rod->setVertexBoundaryCondition(x_2, rod->nv - 2);
		rod->setVertexBoundaryCondition(x_3, rod->nv - 1);	

		totalCompress = totalCompress + 2 * deltaTime * 0.02;
	}

	if (currentTime >= 20.0 && totalShear <= inputShear)
	{
		Vector3d x_0 = rod->getVertex(0);
		Vector3d x_1 = rod->getVertex(1);

		x_0(0) = x_0(0) + deltaTime * 0.02;
		x_1(0) = x_1(0) + deltaTime * 0.02;

		rod->setVertexBoundaryCondition(x_0, 0);
		rod->setVertexBoundaryCondition(x_1, 1);

		Vector3d x_2 = rod->getVertex(rod->nv - 2);
		Vector3d x_3 = rod->getVertex(rod->nv - 1);

		x_2(0) = x_2(0) - deltaTime * 0.02;
		x_3(0) = x_3(0) - deltaTime * 0.02;	

		rod->setVertexBoundaryCondition(x_2, rod->nv - 2);
		rod->setVertexBoundaryCondition(x_3, rod->nv - 1);	

		totalShear = totalShear + 2 * deltaTime * 0.02;
	}

	if (currentTime >= 50.0 && totalAngle <= inputAngle)
	{
		Vector3d x_0 = rod->getVertex(0);
		Vector3d x_1 = rod->getVertex(1);

		double edgeLength = (x_1 - x_0).norm();

		Vector3d tangentVec;

		tangentVec(0) = 0.0;
		tangentVec(1) = 1.0;
		tangentVec(2) = 0.0;

		tangentVec = edgeLength * tangentVec;

		Vector3d xNew;

		xNew = x_1;

		xNew(0) = x_0(0);
		xNew(1) = x_0(1) + tangentVec(1) * cos(totalAngle + deltaTime * speed) + tangentVec(2) * sin(totalAngle + deltaTime * speed);
		xNew(2) = x_0(2) - tangentVec(1) * sin(totalAngle + deltaTime * speed) + tangentVec(2) * cos(totalAngle + deltaTime * speed);

		//cout << x_1.transpose() << " " << xNew.transpose() << endl;

		totalAngle = totalAngle + deltaTime * speed;

		rod->setVertexBoundaryCondition(x_0, 0);
		rod->setVertexBoundaryCondition(xNew, 1);

		rod->setThetaBoundaryCondition(rod->getTheta(0), 0);

		rod->setVertexBoundaryCondition(rod->getVertex(rod->nv - 2), rod->nv - 2);
		rod->setThetaBoundaryCondition(rod->getTheta(rod->nv - 2), rod->nv - 2);
		rod->setVertexBoundaryCondition(rod->getVertex(rod->nv - 1), rod->nv - 1);
	}

	/*

	if (currentTime >= 30.0)
	{
		Vector3d x_2 = rod->getVertex(rod->nv - 2);
		Vector3d x_3 = rod->getVertex(rod->nv - 1);

		x_2(1) = x_2(1) + deltaTime * speed;
		x_3(1) = x_3(1) + deltaTime * speed;	

		rod->setVertexBoundaryCondition(x_2, rod->nv - 2);
		rod->setVertexBoundaryCondition(x_3, rod->nv - 1);
	}

	*/

	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	iter = 0;

	// Start with a trial solution for our solution x
	rod->updateGuess(); // x = x0 + u * dt
		
	while (solved == false)
	{
		rod->prepareForIteration();
		
		stepper->setZero();

		// Compute the forces and the jacobians
			
		m_stretchForce->computeFs();
		m_stretchForce->computeJs();

		m_elasticRibbonForce->computeFribbon();
		m_elasticRibbonForce->computeJribbon();

		if (currentTime < 15.0)
		{
			m_gravityForce->gVector(0) = 0.0000;
			m_gravityForce->gVector(1) = 0.0000;
			m_gravityForce->gVector(2) = 100.0000;

			m_gravityForce->setGravity();

			m_gravityForce->computeFg();
			m_gravityForce->computeJg();
		}

		if (currentTime > 15.0 && currentTime < 40.0)
		{
			m_gravityForce->gVector(0) = 10.0000;
			m_gravityForce->gVector(1) = 0.0000;
			m_gravityForce->gVector(2) = 0.0000;

			m_gravityForce->setGravity();

			m_gravityForce->computeFg();
			m_gravityForce->computeJg();
		}

		if (currentTime > 45.0)
		{
			m_gravityForce->gVector(0) = 0.0000;
			m_gravityForce->gVector(1) = 0.0000;
			m_gravityForce->gVector(2) = 0.0000;

			m_gravityForce->setGravity();

			m_gravityForce->computeFg();
			m_gravityForce->computeJg();
		}

		//if (currentTime < 49.0)
		{
			m_inertialForce->computeFi();
			m_inertialForce->computeJi();
		
			m_dampingForce->computeFd(currentTime);
			m_dampingForce->computeJd(currentTime);
		}

		// Compute norm of the force equations.
		normf = 0;
		for (int i=0; i < rod->uncons; i++)
		{
			normf += totalForce[i] * totalForce[i];
		}
		normf = sqrt(normf);
		if (iter == 0) 
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}
		
		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			rod->updateNewtonX(totalForce); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}
	
	rod->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << " iter=" << iter << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;
	
	if (solved == false)
	{
		timeStep = Nstep; // we are exiting
	}
}

int world::simulationRunning()
{
	//if (totalAngle <= inputAngle * 0.95 || timeStep<Nstep)
	if (timeStep<Nstep) 
		return 1;
	else 
	{
		return -1;
	}
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int i)
{
	return rod->x[i] / RodLength;
}

double world::getBoundaryCoordination_left(int i)
{
	int j = i / 4;
	int k = i - 4 * j;

	return	(rod->x[i] + 0.5 * rod->m1(j, k) * width) / RodLength;
}

double world::getBoundaryCoordination_right(int i)
{
	int j = i / 4;
	int k = i - 4 * j;

	return	(rod->x[i] - 0.5 * rod->m1(j,k) * width) / RodLength;
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}

void world::computeReactionForce()
{
	reactionForce = VectorXd::Zero(rod->ndof);

	//reactionForce = m_inertialForce->forceVec - m_dampingForce->forceVec - m_elasticRibbonForce->forceVec - m_stretchForce->forceVec;
	reactionForce = - m_elasticRibbonForce->forceVec - m_stretchForce->forceVec;
}