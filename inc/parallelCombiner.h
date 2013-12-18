/*******************************************************************************

  Class <parallelCombiner> 

  Author: Andrew Stershic
          Duke Computational Mechanics Lab (DCML)
          Duke University - Durham, NC (USA)
  E-mail: ajs84@duke.edu

  Copyright (c) 2013 Andrew Stershic. All rights reserved. No warranty. No
  liability.

*******************************************************************************/
#include <vector>
#include <string>
#include <algorithm>

#ifndef __PARALLELCOMBINER_H__
#define __PARALLELCOMBINER_H__

class ParallelCombiner {

public:
    //! Constructor
    /*! 
        \brief Loads values for distribution parameters and the number of elements
        \param string : distribution name
        \param parameter1 : first parameter
        \param parameter2 : second parameter
        \param nx : number of nodal points
    */
    //! CTL Empty constructor
    ParallelCombiner ();

    //! Standard destructor
    ~ParallelCombiner (); 

    //! Assign
    /*! 
        \brief Switches to the probability distribution methods and writes them to cohPar
        \param cohPar : cohesive law parameter list, either sigma or delta
    */

	void run(const int procs, const std::string resultsPath, int Nx, bool deleteFlag);

private:

	std::vector<std::vector<std::string> > _fileSet;

	void findSets(std::string inPath);
	void combineSet(std::vector<std::string> inSet, int Nx);
	bool stringsMatch(std::string in1, std::string in2);

	void combineCohLaw(std::string inPath, bool deleteFlag);
	void combineSTheta(std::string inPath, bool deleteFlag);

};

#endif//__PARALLELCOMBINER_H__
