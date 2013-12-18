/*******************************************************************************

  Class <matPropGen> 

  Author: Martin Hautefeuille
          Duke Computational Mechanics Lab (DCML)
          Duke University - Durham, NC (USA)
  E-mail: mh186@duke.edu

  Copyright (c) 2010 Martin Hautefeuille. All rights reserved. No warranty. No
  liability.

*******************************************************************************/
#include "matPropGen.h"
#include "mpi.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
//#include "str_ops.hpp"
#include <stdlib.h>
#include <time.h>
#include "boost/random.hpp"


/*******************************************************************************

                M E T H O D S    F O R    R E S O L U T I O N 

*******************************************************************************/

/*------------------------------- P U B L I C --------------------------------*/

MatPropGen::MatPropGen ( const std::string distrib, const double parameter1,
                     const double parameter2, const double parameter3, const int nx) {

    // assign input values to attributes

    _distrib = distrib;
    _nx = nx;

    //Gamma KT - gamma distribution defined by k & theta (unused, k - shape, theta - scale)
    if (distrib == "GammaKT") {
        _k = parameter2;
        _theta = parameter3;
    }

    //GammaAB - gamma distribution defined by alpha & beta (unused, alpha - shape, beta - rate)
    if (distrib == "GammaAB") {
        _alpha = parameter2;
        _beta = parameter3;
    }

    //Gamma - gamma distribution defined by mean and variance (mean, variance, unused)
    if (distrib == "Gamma") {
        _mean = parameter1;
        _variance = parameter2;
    }

    //LogNorm - log-normal distribution defined by mean and standard deviation (location, scale, unused)
    if (distrib == "LogNorm") {
        _mean = parameter1;
        _stdev = parameter2;
    }

    //Normal - normal distribution defined by mean and standard deviation (mean, stdev, unused)
    if (distrib == "Normal") {
        _mean = parameter1;
        _stdev = parameter2;
    }

    //Uniform - uniform distribution defined by minimum and maximum (maximum, minimum, unused)
    if (distrib == "Uniform") {
        _min = parameter2;
        _max = parameter1;
    }

    //Weibull - Weibull distribution defined by lambda and k (lambda - shape, k - scale, unused)
    if (distrib == "Weibull") {
        _lambda = parameter1;
        _k = parameter2;
    }

    //Poisson Process - Weak Points (one-node length) (normal, weak, lambda)
    if (distrib == "Poisson") {
        _lambda = parameter3 / _nx;	//Expected number of failures in ring
		_normal = parameter1;   //Strength of non-affected regions
		_weak = parameter2;	//Strength of affected regions
    }    

}

void MatPropGen::assign ( std::vector<double>& cohPar, std::vector<unsigned>& cohNums) {

    cohPar.resize(_nx);

    if (_distrib == "Normal") {
	std::vector<double> output;
	Normal(output);
	for (int i = 0; i < _nx; i++) {
	    cohPar[i] = output[i];
	}
    }
    if (_distrib == "GammaKT") {
	std::vector<double> output;
	GammaKT(output);
	for (int i = 0; i < _nx; i++) {
	    cohPar[i] = output[i];
	}
    }
    if (_distrib == "GammaAB") {
	std::vector<double> output;
	GammaAB(output);
	for (int i = 0; i < _nx; i++) {
	    cohPar[i] = output[i];
	}
    }
    if (_distrib == "Gamma") {
	std::vector<double> output;
	Gamma(output);
	for (int i = 0; i < _nx; i++) {
	    cohPar[i] = output[i];
	}
    }
    if (_distrib == "LogNorm") {
	std::vector<double> output;
	LogNorm(output);
	for (int i = 0; i < _nx; i++) {
	    cohPar[i] = output[i];
	}
    }
    if (_distrib == "Uniform") {
	std::vector<double> output;
	Uniform(output);
	for (int i = 0; i < _nx; i++) {
	    cohPar[i] = output[i];
	}
    }
    if (_distrib == "Weibull") {
	std::vector<double> output;
	Weibull(output);
	for (int i = 0; i < _nx; i++) {
	    cohPar[i] = output[i];
	}
    }

    if (_distrib == "Poisson") {
	std::vector<double> output;
	Poisson(output);
	for (int i = 0; i < _nx; i++) {
	    cohPar[i] = output[i];
	}

	//Seeds poisson weak spots to cohNums up to weak spot count and original size
	int last = -1;
	unsigned count = 0;
	for (int i = 0; (i < _nx) && (count < cohNums.size()); i++) {
	    if ((output[i] == _weak) && (i > last)) {
		last = i;
		cohNums[count] = i;
		count++;
	    }
	}
    }
}

MatPropGen::MatPropGen () {}
MatPropGen::~MatPropGen () { }



/*------------------------------ P R I V A T E -------------------------------*/


void MatPropGen::Normal (std::vector<double>& output) {
    //Create a normal distrib'd number
    boost::mt19937 rng(static_cast<unsigned> (time(0)));

    boost::normal_distribution<> nd(_mean, _stdev);

    boost::variate_generator<boost::mt19937, boost::normal_distribution<> > var_nor(rng, nd);

    output.resize(_nx);
    for (int i = 0; i < _nx; i++){
	output[i] = var_nor();
    }
}

void MatPropGen::LogNorm (std::vector<double>& output) {
    //Create a lognormal distrib'd number
    boost::mt19937 rng(static_cast<unsigned> (time(0))); 

    boost::lognormal_distribution<> lnd(_mean, _stdev);

    boost::variate_generator<boost::mt19937, boost::lognormal_distribution<> > var_lognor(rng, lnd);

    output.resize(_nx);
    for (int i = 0; i < _nx; i++){
	output[i] = var_lognor();
    }
}

void MatPropGen::Uniform (std::vector<double>& output) {
    //Create a uniform distrib'd number
    boost::mt19937 rng(static_cast<unsigned> (time(0))); 

    boost::uniform_real<> ud(_min, _max);

    boost::variate_generator<boost::mt19937, boost::uniform_real<> > var_unif(rng, ud);

    output.resize(_nx);
    for (int i = 0; i < _nx; i++){
	output[i] = var_unif();
    }
}

void MatPropGen::GammaAB (std::vector<double>& output) {
    //Create a gamma distrib'd number
    boost::mt19937 rng(static_cast<unsigned> (time(0)));

    boost::gamma_distribution<> gd(_alpha);

    boost::variate_generator<boost::mt19937, boost::gamma_distribution<> > var_gamma(rng, gd);

    output.resize(_nx);
    for (int i = 0; i < _nx; i++){
	output[i] = 1 / _beta * var_gamma();
    }
}

void MatPropGen::GammaKT (std::vector<double>& output) {
    //Create a gamma distrib'd number
    boost::mt19937 rng(static_cast<unsigned> (time(0)));

    boost::gamma_distribution<> gd(_k);

    boost::variate_generator<boost::mt19937, boost::gamma_distribution<> > var_gamma(rng, gd);

    output.resize(_nx);
    for (int i = 0; i < _nx; i++){
	output[i] = _theta * var_gamma();
    }
}

void MatPropGen::Gamma (std::vector<double>& output) {
    //Calculate input parameters
    _k = _mean * _mean / _variance;
    _theta = _variance / _mean;

    //Create a gamma distrib'd number
    boost::mt19937 rng(static_cast<unsigned> (time(0)));

    boost::gamma_distribution<> gd(_k);

    boost::variate_generator<boost::mt19937, boost::gamma_distribution<> > var_gamma(rng, gd);

    output.resize(_nx);
    for (int i = 0; i < _nx; i++){
	output[i] = _theta * var_gamma();
    }
}

void MatPropGen::Weibull (std::vector<double>& output) {
    //Create a uniform distrib'd number
    boost::mt19937 rng(static_cast<unsigned> (time(0)));

    boost::uniform_real<> wbd(0, 32767);

    boost::variate_generator<boost::mt19937, boost::uniform_real<> > var_wbd(rng, wbd);

    output.resize(_nx);
    for (int i = 0; i < _nx; i++){
        //Create a Weibull distrib'd number from the uniform distrib'd number
	output[i] = pow( -log(1 - (var_wbd() / 32768)), 1/_k) * _lambda;
    }
}


void MatPropGen::Poisson (std::vector<double>& output) {
    //Create a uniform distrib'd number
    boost::mt19937 rng(static_cast<unsigned> (time(0))); 

    boost::uniform_real<> ud(0, 32768);

    boost::variate_generator<boost::mt19937, boost::uniform_real<> > var_unif(rng, ud);

    output.resize(_nx);
    double prob = 0;
    unsigned count = 0;
    unsigned count2 = 0;
    double x = 0;
    double var;

    //Generates an exponentially distributed random variable based
    //distance from the last weak spot and compares it to a uniform 0-1 R.V.
    //to decide whether or not this is a weak spot.
    for (int i = 0; i < _nx; i++){
        //Create an exponentially distrib'd number from the uniform distrib'd number
	x = (double)(count + 1)/(double)(_nx);
	prob = _lambda * exp(-_lambda * x);
	var = 1 - (var_unif() / 32768);
	if (var < prob) {
	    output[i] = _weak;
	    count = 0;
	    count2++;
	} else {
	    output[i] = _normal;
	    count++;
	}
    }
    std::cout << "-------" << std::endl;
    std::cout << "Poisson Process generated weak spots:  " << count2 << std::endl;
    std::cout << "-------" << std::endl;
}


