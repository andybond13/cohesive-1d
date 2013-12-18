/*******************************************************************************

  Class <matPropGen> 

  Author: Martin Hautefeuille
          Duke Computational Mechanics Lab (DCML)
          Duke University - Durham, NC (USA)
  E-mail: mh186@duke.edu

  Copyright (c) 2010 Martin Hautefeuille. All rights reserved. No warranty. No
  liability.

*******************************************************************************/
#include <vector>
#include <string>
#include <algorithm>

#ifndef __MATPROPGEN_H__
#define __MATPROPGEN_H__

class MatPropGen {

public:
    //! Constructor
    /*! 
        \brief Loads values for distribution parameters and the number of elements
        \param string : distribution name
        \param parameter1 : first parameter
        \param parameter2 : second parameter
        \param nx : number of nodal points
    */
    MatPropGen ( const std::string distrib, const double parameter1,
                     const double parameter2, const double parameter3, const int nx);

    //! CTL Empty constructor
    MatPropGen ();

    //! Standard destructor
    ~MatPropGen (); 

    //! Assign
    /*! 
        \brief Switches to the probability distribution methods and writes them to cohPar
        \param cohPar : cohesive law parameter list, either sigma or delta
    */
    void assign (std::vector<double>& cohPar, std::vector<unsigned>& cohNums);

private:

    //!
    //!  M E T H O D S    F O R    R E S O L U T I O N
    //!

    //! Distribution Parameters
    /*! 
        \brief Input parameters for the individual probability distributions
    */
    double _k;
    double _theta;
    double _alpha;
    double _beta;
    double _lambda;
    double _mean;
    double _stdev;
    double _variance;
    double _min;
    double _max;
    double _normal;
    double _weak;

    //! Other inputs
    /*! 
        \brief Which distribution is specified and the size of cohPar.
    */
    std::string _distrib;
    int _nx;

    //! Distribution Methods
    /*! 
        \brief Methods to generate random variables according to the distributions named
			for use in defining material properties - e.g. _SigCoh.
    */
    void GammaKT (std::vector<double>& output);
    void GammaAB (std::vector<double>& output);
    void Gamma (std::vector<double>& output);
    void Normal (std::vector<double>& output);
    void LogNorm (std::vector<double>& output);
    void Uniform (std::vector<double>& output);
    void Weibull (std::vector<double>& output);
    void Poisson (std::vector<double>& output);

};

#endif//__MATPROPGEN_H__
