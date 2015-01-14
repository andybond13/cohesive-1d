/*******************************************************************************

  Class <CartRing> 

  Author: Martin Hautefeuille
          Duke Computational Mechanics Lab (DCML)
          Duke University - Durham, NC (USA)
  E-mail: mh186@duke.edu

  Copyright (c) 2010 Martin Hautefeuille. All rights reserved. No warranty. No
  liability.

*******************************************************************************/
#include "cartRing.h"
#include <iostream>
#include <math.h>
#include "matPropGen.h"
#include "mpi.h"
#include "parallelCombiner.h"

using namespace std;

int main (int argc, char* argv[]) {

    // ring specification
    double L   = .05;   		// circumference, in m
    double A   = 1.0e-6;   	// cross-sectional area, in m^2
    double rho = 2.75e+3;  	// density, in kg/m3
    double E   = 2.75e+11; 	// young's modulus, in Pa
    int nx     = 100;		// number of elements
	
	//path to results directory
//    std::string path = "~/andrew/Duke/results";
	std::string path = "/Volumes/Stershic-HybridHD/Users/andrewstershic/Code/flaming-shame/results";
    std::string multiResultsPath = "/Volumes/Stershic-HybridHS/Users/andrewstershic/Code/flaming-shame/results";

	//create ring object
    CartRing ring( L, A, rho, E, nx, path );

	//apply boundary conditions
	//ring.applyForc( "CONS", "THETA", 6.0e-1 );
    //ring.initVel( "RADIA", 2.5*1.592e+1 );
    //ring.applyForc( "LINE", "RADIA", 1.592e+1 );
//    ring.applyForc( 6.0e-1 );
    ring.applyVel( "CONST_ENDSR", 0.004*SR);	//apply strain rate = 1000 1/s;

    // cohesive law
    /**/
    std::vector<std::vector<double> > cohPar( 2 );
    /*cohPar[0].assign( nx, 4.0e+12 ); //4.0e+12 (assign all strengths manually)
    cohPar[0][nx-100] = 4.4e+9; //(with exceptions)
    cohPar[0][nx-125] = 4.37e+9;
    cohPar[0][nx-150] = 4.37e+9;
    cohPar[0][nx-525] = 4.4e+9;*/

	//assign strengths via distribution
    std::vector<unsigned> cohNums(4);
    cohNums[0] = 41; cohNums[1] = 92; cohNums[2] = 99; cohNums[3] = 107;
    std::string distrib = "Normal";
    double param1 = 3e+8; 	//4.37785e+8 Pa;		//Theoretically stress_max = v0*sqrt(E*rho)
    double param2 = 0;
    double param3 = 0;

    MatPropGen properties(distrib, param1, param2, param3, nx);
    properties.assign(cohPar[0], cohNums); //assign randomly generated values to cohesive strength: cohPar[0]
    cohPar[1].assign( 1, 100 );	//Gc (in N/m = Pa*m) value for cohesive zone: cohPar[1]
    ring.setCohLaw( "SQRTSG", cohPar );	//set cohesive zone parameters
    ring.plotCohLaw( cohNums );
    

    ring.display( 1, 1 );//print vtk files every x timesteps

    unsigned nodes[] = { 0, 1 };
    std::vector<unsigned> node2plot ( nodes, nodes+2 );
    ring.plotAtNodes( node2plot );

    unsigned elm[] = { 0, 1 };
    std::vector<unsigned> elm2plot ( elm, elm+2 );
    ring.plotAtElms( elm2plot );

    ring.plotEnergies();
    ring.plotFrags();
    ring.plotSTheta();
    ring.defectLimit(0.000);		//minimum distance between opening locations



    double totalTime = 5.0e-06;		//end time of simulation
    unsigned printFreq = 1;			//print table/graph data every x timesteps
	bool allowPlateauEnd = true; 	//allow the program to exit early when cohesive/fracture energy plateaus
	bool checkEnergy = true;			//report on energy balance, end when too much energy generated
    double refine = 0.01;				//refine the timestep by a factor of x per level of refinement (x or x^2)
	ring.solve( totalTime, printFreq, refine, allowPlateauEnd, checkEnergy );
    ring.printHisto();

	//grab results for parallel combining
	double runTime = 0.0; unsigned numFrag = 32768; unsigned nIter = 0; double Wcoh0; double Wsum; double Wmax; double WsprD;
	std::vector<double> fragLength; double meanFragLength; std::vector<unsigned> fHisto;
	std::vector<std::vector<double> >fragInvCDF;
	ring.grabInfo( runTime, numFrag, nIter, Wcoh0, Wsum, Wmax, fragLength, meanFragLength, WsprD, fHisto, fragInvCDF);

	if (nIter>0) {	//should eliminate _myid >0
		//std::cout << "runtime: " << runTime << "   numFrag: " << numFrag << " nIter: " << nIter << std::endl;

		if (argc > 1) {
			//collect results off compute nodes before combining
			string command = "python distributedCollector.py " + path + " " + multiResultsPath + " ";
			for (int i = 1; i < argc; ++i) {
				command += (string)argv[i] + " ";
			}
			system( command.c_str() );
		}

		ParallelCombiner pc; pc.run(1,multiResultsPath,nx,true);
	}
    return 0;
}

