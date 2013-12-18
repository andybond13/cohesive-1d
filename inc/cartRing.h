/*******************************************************************************

  Class <cartRing> 

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
#include <ctime> 

#ifndef __CARTRING_H__
#define __CARTRING_H__

class CartRing {

public:
    //! Constructor
    /*! 
        \brief Constructor used for deterministic cases where no random fields
               are involved for material parameters.
        \param length length of the ring
        \param crossSec cross section of the ring
        \param density density of the constitutive material
        \param YoungMod Young modulus of the constitutive material
        \param numNod number of nodal point
    */
    CartRing ( const double length, const double crossSec, const double density,
               const double YoungMod, const int nodNum, const std::string path);

    //! CTL Empty constructor
    CartRing ();

    //! Standard destructor
    ~CartRing (); 

    //! set the cohesive law
    /*!
        \brief Method to set the type of cohesive law in the cohesive links as
               well as the corresponding parameters.
        \param lawTyp : LINSD -> LINear defined by Sigma and Delta
                        LINSG -> LINear defined by Sigma and Gc (frac energy)
        \param param : vector of vector valued parameters that defined the
                       cohesive law. The size of the vector for each value is
                       either 1 if the value is constant in each cohesive link
                       or the number of cohesive links.
    */
    void setCohLaw ( const std::string& lawTyp,
                     const std::vector<std::vector<double> >& param );

    //! Method to apply forces on the mesh
    /*!
        \brief Method to apply body force or torque depending on the load
               direction
        \param loadTyp : CONS, LINE, ...
        \param loadDir : RADIAL, THETA
        \param loadVal : value of the loading in s.i.u
    */
    void applyForc ( const std::string& loadTyp, const std::string& loadDir,
                     const double loadVal );

    //! Method to apply constant velocity / strain rate on the mesh
    /*!
        \brief Method to apply velocity depending on the load direction
        \param loadDir : RADIAL, THETA
        \param loadVal : value of the loading in s.i.u
    */
    void applyVel (  const std::string& velDir, const double velVal  );

    //! Method to set initial conditions
    /*!
        \brief Method to set initial velocity
        \param velDir : RADIAL, THETA
        \param velVal : initial value of the prescribed velocity
    */
    void initVel ( const std::string& velDir, const double velVal );

    //! Solve the problem
    /*!
        \brief Solve the mechanical problem until a given time T
        \param endTime: time T in s.i.u.
	  \param printFrequency : prints all data every x iterations
	  \param refine : time step refinement factor
    */
    void solve ( const double endTime, const unsigned printFrequency, const double refine, const bool allowPlateauEnd, const bool checkEnergy );

    //! Plot nodal velocities, elem stress and global energies
    /*!
        \brief Prescribe a list of nodal Ids where velocities should
               be plotted
               Prescribe a list of elem Ids where stress should be plotted
		   Method to print fragmentation histogram to console
    */
    void plotAtNodes ( const std::vector<unsigned>& NodalIds );
    void plotAtElms ( const std::vector<unsigned>& ElmIds );
    void plotEnergies ();
    void plotCohLaw ( const std::vector<unsigned>& cohNum);
    void plotFrags ();
    void printHisto ();
    void plotSTheta ();


    //! Set defection limit
    /*!
        \brief +/- from opening on either side
        \param defectLimit : distance away from current crack opening up to which addition
					openings are prohibited. If = 0, then there is no limit.
    */
    void defectLimit ( const double& defectRange);


    //! Print Vtk files every prescribed time step number
    /*!
        \brief Method to prescribe the time step at which a vtk file must be
               printed during both the elastic and the softening regimes.
        \param timStepNumElas number of elastic time step between which vtk
               files must be printed
        \param timStepNumFrac number of time step after the first failure
               occurs between which vtk files must be printed
    */
    void display ( const unsigned timStepNumElas,
                   const unsigned timStepNumFrac );


    //! Return information after run to post-processor; full and scalar only
    /*!
        \brief Method to return completed-run information to post-processor
        \param runTime _T
        \param numFrag _numFrag
	  \param nIter _Nt
        \param Wcoh0 _Wcoh[0]
        \param Wsum _Wsum
	  \param Wmax _Wmax
	  \param fragLength _fragLength
	  \param fHisto fragmentation histogram
	  \param fragInvCDF fragmentation CDF
    */
    void grabInfo ( double& runTime, unsigned& numFrag, unsigned& nIter, double& Wcoh0, double& Wsum,
			double& Wmax, std::vector<double>& fragLength, double& meanFragLength, double& WsprD,
			std::vector<unsigned>& fHisto, std::vector<std::vector<double> >& fragInvCDF);

    void grabInfo ( double& runTime, unsigned& numFrag, unsigned& nIter, double& Wcoh0, double& Wsum,
			double& Wmax, double& meanFragLength, double& WsprD);



private:

    //!
    //!  M E T H O D S    F O R    R E S O L U T I O N
    //!

	//! Perform domain decomposition
	/*!
		\brief divide mesh into different processors; assign _local & _owned to reveal if nodes are locally owned and if not, who the owner is
	*/	
	void domainDecomposition();


	//convert integer to string
	std::string convertInt ( int in) const;

    //! Method to build the mesh
    void buildDiscretization ();

    //! Compute sprVec & cohVec
    /*!
        \brief sprVec : Method to compute the vector which direct a given
               spring element from its first node to its second as
               prescribed in _SprCon.
               sprVecPred : the same but with the prediction of displacements
               added.
               cohVecPred : Method to compute the vector which direct a given
               cohesive link from its first node to its second as set up in
               _CohCon. It includes the predictions of displacements
        \param sprNum given spring element
        \param cohNum given cohesive link
        \return vector which direct the spring element or the cohesive link
    */
    std::vector<double> sprVec ( const unsigned sprNum ) const;
    std::vector<double> sprVecPred ( const unsigned sprNum ) const;
    std::vector<double> cohVecPred ( const unsigned cohNum ) const;

    //! Compute cosTheta at a node
    /*!
        \brief Method to compute the cosine between the current nodal position
               and the X-axis
        \param nodNum node
        \return cosine
    */
    double cosTheta ( const unsigned nodNum ) const;
    double cosThetaPred ( const unsigned nodNum ) const;

    //! Compute sinTheta at a node
    /*!
        \brief Method to compute the sine between the current nodal position
               and the Y-axis
        \param nodNum node
        \return sine
    */
    double sinTheta ( const unsigned nodNum ) const;
    double sinThetaPred ( const unsigned nodNum ) const;

    //! Methods for second-order explicit Newmark integration scheme
    /*!
        \brief Method to predict velocities and displacements for 2nd order
               explicit Newmark scheme. ( beta = 0 and gamma = 0.5 )
               Method to solve the accelerations for 2nd order explicit Newmark
               scheme.
               Method to correct both the velocities and the displacements for
               2nd order explicit Newmark scheme.
    */
    void NewmarkPred ();
    void NewmarkReso ();
    void NewmarkCorr ();

    //! Compute spring forces for a given spring
    /*!
        \brief sprForc Method to compute the elastic spring force
        \param sprNum spring element
        \return elastic energy for spring element
    */
    double sprForc ( const unsigned sprNum );

    //! Compute cohesive forces for a given cohesive element
    /*!
        \brief cohForc Method to compute the cohesive force
        \param cohNum cohesive element
        \return both dissipated and elastic energy for the cohesive element
    */
    std::vector<double> cohForc ( const unsigned nodNum );

    //! Compute external force for a node 
    /*!
        \brief extForc Method to compute the external force
        \param nodNum node
        \return external energy at each node
    */
    double extForc ( const unsigned nodNum );

    //! Compute required force for a node to maintain constant velocity
    /*!
        \brief calcVelForc Method to compute the external force
        \param nodNum node
    */
    void calcVelForc ( const unsigned i );

    //! Compute stress at a given spring element
    /*!
        \brief Method to compute the stress at a given spring element
        \param sprNum given spring element
        \return stress value
    */
    double stress ( const unsigned sprNum );

    //! Compute stress at a given cohesive link
    /*!
        \brief Method to compute the stress at a given cohesive element. This contains the cohesive zone law formula.
        \param cohNum given cohesive link
        \return stress value
    */
    void cohStr ( const unsigned cohNum);

    //! Check the balance of energy
    /*!
        \brief Method to check the balance of energy at each time step
               as described in ( Belytschko pp. 315 )
    */
    void energBalance ();

    //! Count the number of fragments
    /*!
        \brief Method to count the number of fragments at each time step
    */
    void fragCount ();

    //! Check the stability
    /*!
        \brief Method to check stability and decide whether or not to enable time step decrease
    */
    void checkStable ();

    //! Refine the time step, if necessary
    /*!
        \brief Method to affect time step refinement, if necessary
    */
    void timeStepRefine (const double refine);

    //! Update all the time dependent value
    /*!
        \brief Method to update all the time dependant values. In short, it
               amounts to replace the values of the variables computed at time
               t_n by the one computed at time t_n+1
    */
    void update ();

    //! Use MPI to exchange the boundary nodes for calculating spring force and spring force
    /*!
        \brief Method to define the boundary nodes
			   Method to exchange the boundary nodes for calculating spring force and spring force, using MPI
    */
    void defineBoundaryNodes ();
    void exchangeBoundaryNodes ();
	void exchangeSprForc ();

    //! Use MPI to exchange the fragment number and location info
    /*!
        \brief Method to define the boundary nodes
			   Method to xchange the fragment number and location info, using MPI
    */
	void exchangeFragInfo ();

    //!
    //!  M E T H O D S    T O    P R I N T    S T U F F
    //!

    //! Print Spring and Cohesive zone connectivities
    /*!
        \brief Method to print on screen the connectivities of the Springs and 
               the cohesive zones
    */
    void printConnec() const;

    //! Print available pieces of information at time t_n
    /*!
        \brief Method to print available pieces of information at time t_n
        \param t_n current time-step
    */
    void printVtk ( const unsigned timStepNum ) const;

    //! Print stuff in the vtk file format
    /*!
        \brief Method to print in a new file the usual vtk header
               Method to print in a file the mesh of the ring
               according to the vtk file format. All the nodes are written
               twice to handle failure and fragmentation.
               Method to print in a file the nodal values of the displacements
               and of the velocities
               Method to print in a file the elementary values of the stress
        \alert This method use the values of variables at time t_n
        \param vtkFile vtk file
    */
    void printHeader ( const std::string& vtkFile ) const;
    void printMesh ( const std::string& vtkFile ) const;
    void printPointData ( const std::string& vtkFile ) const;
    void printCellData ( const std::string& vtkFile ) const;

    //! Print elementary, nodal, global, and cohesive link pieces of information
    /*!
        \brief Method to print nodal pieces of information for a prescribed
               node in a file
               Method to print elm pieces of information for a prescribed
               elm in a file
               Method to print global pieces of information
		   Method to print cohesive link information
		   Method to print fragmentation information
        \alert Both these methods should placed after the update
    */
    void printNodalInfo () const;
    void printElmInfo () const;
    void printGlobalInfo () const;
    void printCohLaw () const;
    void printFrags () const;
    void printSTheta () const;
    void printClean () const;
    void plotHisto ();

    //!
    //!  A T T R I B U T E S
    //!

    //! Scalar attributes of the class 
    /*!
        \brief _L   : Total length of the ring
               _A   : Cross section of the ring
               _rho : Density of the constitutive material 
               _E   : Young's modulus of the consitutive material 
               _R0  : Radius of the ring
               _Nx  : Number of elements / Number of nodes
               _c   : Wave celerity in the ring  
               _Dx  : Space interval
               _m   : nodal mass
		   _path: output path (results folder name)
		   _logPath: log file name and path
		   _defectRange: crack opening limit range, +/- from opening on either side
		
    */
	unsigned _myid;
	unsigned _numprocs;
    double _L;
    double _A;
    double _rho;
    double _E;
    double _R0;
    unsigned _Nx;
	std::vector<unsigned> _local;
	std::vector<unsigned> _owner;
	unsigned _begin;
	unsigned _end;
    double _c;
    double _Dx;
    double _m;
    std::string _path;
    std::string _logPath;
    double _defectRange;
	std::clock_t _start;
	bool _allowPlateauEnd;
	bool _checkEnergy;
	std::string _lawTyp;

	//boundary node and owner/neighbor list
	std::vector<unsigned> _nodeList;
	std::vector<unsigned> _originList;
	std::vector<unsigned> _destList;
//	std::vector<int> _dirList;	//1 for owner 1->owner 2. 2 for owner2->owner1.

    //! Attributes describing the mesh
    /*!
        \brief _NodPos : nodal position -- (x,y) = (_NodPos[i][0],_NodPos[i][0])
	         _NodPosOrig : orig. nodal position -- (x,y) =(_NodPos[i][0],_NodPos[i][0])
               _SprCon : spring connectivity -- _SprCon[i].first = 1st node
                         & _SprCon[i].second = 2nd node
               _CohCon : cohesive connectivity -- _CohCon[i].first = 1st node
                         & _CohCon[i].second = 2nd node
               _NodCon : _NodCon[i].first = spring connected to node i
                         & _NodCon[i].second = cohesive joint connected to i
    */
    std::vector<std::vector<double> > _NodPos;
    std::vector<std::vector<double> > _NodPosOrig;
    std::vector<std::pair<unsigned,unsigned> > _SprCon;
    std::vector<std::pair<unsigned,unsigned> > _CohCon;
    std::vector<std::pair<unsigned,unsigned> > _NodCon;

    //! Attributes describing the kinematics
    /*!
        \brief _Dis : _Dis[i][0] = displacement of node i at time t_n
                    & _Dis[i][1] = prediction of disp of node i at time t_p
                    & _Dis[i][2] = prediction of disp of node i at time t_(n+1)
               _Vel : _Vel[i][0] = velocity of node i at time t_n
                    & _Vel[i][1] = prediction of velocity at time t_p
                    & _Vel[i][1] = velocity of node i at time t_(n+1)
               _Acc : _Acc[i][0] = acceleration of node i at previous time-step
                    & _Vel[i][1] = acceleration of node i at next time-step
    */
    std::vector<std::vector<std::vector<double> > > _Dis;
    std::vector<std::vector<std::vector<double> > > _Vel;
    std::vector<std::vector<std::vector<double> > > _Acc;

    //! Forces
    /*!
        \brief _Fspr : _Fspr[i] = [Fx, Fy] spring force at node i
               _Fcoh : _Fcoh[i] = [Fx, Fy] cohesive force at node i
               _Fext : _Fext[i][0] = [Fx, Fy] external force at node i at t_n
                       _Fext[i][1] = [Fx, Fy] external force at node i at t_n+1
    */
    std::vector<std::vector<double> > _Fspr;
    std::vector<std::vector<double> > _Fcoh;
    std::vector<std::vector<std::vector<double> > > _Fext;

    //! _Stress
    /*!
        \brief _Stress[i]: Stress value at element i
    */
    std::vector<double> _Stress;

    //! Attributes describing the Boundary and Initial Conditions
    /*!
        \brief _NodForcBC[i] : R (i=0), T (i=1)
                             -> 1 (CONS), 2 (LINE: 0 -> Val)
               _ValForcBC[i] : R (i=0), T (i=1)
                             -> Value of force or torque (i=1)
               _ValVelBC[i] : R = 1, T = 2 (i=0) -> Value of velocity or rotation (i=1)
		   _ConstSRFlag : 0 if not, 1 if so
		   _VelForcReq[i][j] : the force required to maintain constant velocity for
						at node i. X-direction at j=0, y-direction at j = 1;
    */
    std::vector<unsigned> _NodForcBC;
    std::vector<double> _ValForcBC;
    std::vector<double> _ValVelBC;
    std::vector<std::vector<double> > _VelForcReq;
    unsigned _ConstSRFlag;

    //! Attributes dealing with time evolution
    /*!
        _T    : current time
        _Nt   : current time-step number
        _Dt   : current time-step value 
        _Dt_c : current critical time-step value
	  _tFlag: flag to signal tighter time-step, _tFlag[0]=past, _tFlag[1]=current
	  _deactive : indicator of stage of progress of time-step relaxing 0:0.0001:1.0
	  _deactive2 : indicator of stage of progress of time-step relaxing 0:0.0001:1.0
    */
    double _T;
    unsigned _Nt;
    double _Dt;
    double _Dt_c;
    std::vector<unsigned> _tFlag;
    double _deactive;
    double _deactive2;
    bool _stopFlag;

    //! Attributes to plot velocities, stresses and energies into files
    /*!
        _NodesToPlot : list of nodal Ids where velo should be plotted
        _NodeFiles   : list of file names, one per node
        _ElmsToPlot  : list of elem Ids where stress should be plotted
        _ElmFiles    : list of file names, one per elm
        _EnrgFile    : file name where to print energies
	  _CohLawFile  : file name where to print Cohesion information
        _FragFile    : file name where to print fragmentation information
	  _HistoFile   : file name where to print fragmentation histogram information
	  _SThetaFile  : file name where to print stress vs. theta should be plotted
    */
    std::vector<unsigned> _NodesToPlot;
    std::vector<std::string> _NodeFiles;
    std::vector<unsigned> _ElmsToPlot;
    std::vector<std::string> _ElmFiles;
    std::string _EnrgFile;
    std::string _CohLawFile;
    std::string _FragFile;
    std::string _HistoFile;
    std::string _SThetaFile;

    //! Attributes to display meshes of the ring at different time-step
    /*!
        \brief true if vtk files are desired
               time intervale between two vtk printings in elastic regime
               time intervale between two vtk printings after the first failure
    */
    bool _DisplayFlag;
    unsigned _DtPrintElas;
    unsigned _DtPrintFrac;

    //! Attributes describing the cohesive law
    /*!
        \brief _SigC     : Ultimate stress
               _DelC     : Critical opening of cracks
               _CohParam : Other params...
    */
    std::vector<double> _SigC;
    std::vector<double> _DelC;

    //! Attributes for internal variables
    /*!
        \brief _ActivCoh : if true the cohesive zone is activated
               _D        : _delta / _DelC -- _D[0] = former & _D[1] = future
		   _sigCoh : cohesive link stress
		   _delta : cohesive link separation/extension
		   _cLaw : list of cohesive links to plot _sigCoh vs. _delta
		   _sprDamage : damage of springs allowed to have perfect plastic damaging

    */
    std::vector<int> _ActivCoh;
    std::vector<std::vector<double> > _D;
    std::vector<double> _sigCoh;
    std::vector<double> _delta;
    std::vector<unsigned> _cLaw;
    std::vector<double> _sprDamage;

    //! Attributes for fragmentation information
    /*!
	 \brief  _numFrag : number of fragments (sum of _D values >=1)
		  _fragLength : length of each fragment
		  _fragLoc : location of each fragment
		  _DSum : approximate fractional fragmentation (sum of _D values)
		  _fMean : fragment mean length
\		  _fMed : fragment median length
		  _fMax : fragment maxmum length
		  _fMin : fragment minimum length
		  _fStDev : fragment length std. deviation (Rayleigh)
		  _fRange : fragment length range
		  _fSkew : fragment length skew
		  _fExKurtosis : fragment length excess kurtosis
		  _fHisto : array of histogram bins, by 10%'s.
		  _fragInvCDF : inverse cdf, size vs. count (fragsize>size)
    */
	// 
    unsigned _numFrag;
    std::vector<int> _fragLoc;
    std::vector<double> _fragLength;
    double _DSum;
    double _fMean;
    double _fMed;
    double _fMax;
    double _fMin;
    double _fStDev;
    double _fRange;
    double _fSkew;
    double _fExKurtosis;
    std::vector<unsigned> _fHisto;
    std::vector<std::vector<double> > _fragInvCDF;

    //! Attributes for energies
    /*!
        \brief _Wkin : Total kinematic energy 
               _Wext : Total external energy 
               _WextT: Total external energy, global
               _Wspr : Total spring energy
               _Wcoh : _Wcoh[0] = Total dissipated energy
                       _Wcoh[1] = Total elastic energy
		   _Wsum : Total system energy
		   _Wmax : Maximum energy component
		   _dWcoh : [difference in Wcoh0, previous value of Wcoh0, previous _dWcoh, d2Wcoh, time of negative]
		   _Wcoh100 : [Wcoh value 100 iterations previous, time of 100 iterations previous, tangent slope, time of negative d2) 
		   _WsprD : energy of active cohesive links in present and past
    */
    double _Wkin;
    double _Wext;
    double _WextT;
    double _Wspr;
    std::vector<double> _Wcoh;
    double _Wsum;
    double _Wmax;
    std::vector<double> _dWcoh;
    std::vector<double> _Wcoh100;
    double _WsprD;
};

#endif//__CARTRING_H__
