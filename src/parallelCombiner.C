/*******************************************************************************

  Class <parallelCombiner> 

  Author: Andrew Stershic
          Duke Computational Mechanics Lab (DCML)
          Duke University - Durham, NC (USA)
  E-mail: ajs84@duke.edu

  Copyright (c) 2013 Andrew Stershic. All rights reserved. No warranty. No
  liability.

*******************************************************************************/
#include "parallelCombiner.h"
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "alphanum.h"

using namespace std;
using namespace boost::filesystem; 

/*******************************************************************************

                M E T H O D S    F O R    R E S O L U T I O N 

*******************************************************************************/

/*------------------------------- P U B L I C --------------------------------*/

ParallelCombiner::ParallelCombiner ()  {}

ParallelCombiner::~ParallelCombiner ()  {}

void ParallelCombiner::run(const int procs, const std::string resultsPath, int Nx, bool deleteFlag) {

	/*
	 CURRENTLY COMBINES VTK AND cohLaw.dat AND stressTheta.dat FILES ONLY!!!	
	 -takes results files from each sub-domain and writes them into one unified file
	*/

	string vtkPath = resultsPath + "/vtkFiles/";
	string datPath = resultsPath + "/datFiles/";


	//VTK-----------------------------------------------
	//check vtk files for sets
	findSets(vtkPath);

	//combine sets
	for (unsigned i = 0; i < _fileSet.size(); ++i) {
		combineSet(_fileSet[i], Nx);

		//remove unnecessary files
		if (deleteFlag) {												
			for (unsigned j = 0; j < _fileSet[i].size(); ++j) {
				string in2 = _fileSet[i][j];
	//			cout << (exists(in2) ? "Found: " : "Not found: ") << in2 << endl;
				remove(in2);
	//			cout << (exists(in2) ? "Found: " : "Not found: ") << in2 << endl;
	//			cout << endl;
			}
		}
	}

	cout << "*** vtk files combined" << endl;

	//cohLaw.dat----------------------------------------
	//combine cohLaw.dat's
	combineCohLaw(datPath,deleteFlag);

	//stresstheta.dat-----------------------------------
	//combine stresstheta.dat's
	combineSTheta(datPath,deleteFlag);

	return;
}

/*------------------------------ P R I V A T E -------------------------------*/

void ParallelCombiner::findSets(string inPath) {

	string p = inPath;
	vector<string> allFiles(0);

	//get file listing
	if (is_directory(p)) {
		for (directory_iterator itr(p); itr!=directory_iterator(); ++itr) {
//			cout << itr->path().filename() << " "; // display filename only
			if (is_regular_file(itr->status())) {
				allFiles.push_back(itr->path().string());
			}
		}
	} else {
		//cout << (exists(p) ? "Found: " : "Not found: ") << p << endl;
	}

	std::sort(allFiles.begin(), allFiles.end(), doj::alphanum_less<std::string>());

	//sort into sets
	_fileSet.resize(0);
	for (unsigned i = 0; i < allFiles.size()-1; ++i) {
		if (i >= allFiles.size()) return;
		vector<string> files(0);
		files.push_back(allFiles[i]);
		for (unsigned j = i+1; j < allFiles.size(); ++j) {
			if (stringsMatch(allFiles[i],allFiles[j])) {
				files.push_back(allFiles[j]);
				allFiles.erase(allFiles.begin()+j);
				j--;
			}	
		}

		if (files.size() >= 2) {
			_fileSet.push_back(files);
			allFiles.erase(allFiles.begin()+i);
			i--;
		}
	}
	
	return;
}

void ParallelCombiner::combineSet(std::vector<std::string> inSet, int nx) {
	
	//0: find end of points dataset
	//other: find beginning and end of points data set
	//copy in sequentially
	//repeat for cells,displacement,velocity,stress
	//add together numbers


	//open target file
	assert(exists(inSet[0]));
	ofstream outFile;
	string outFileName = inSet[0].substr(0,inSet[0].length()-4) + "_multi.vtk";

	std::sort(inSet.begin(), inSet.end(),doj::alphanum_less<std::string>());

	outFile.open(outFileName.c_str());	

	//open all files in set
	string line;
	string oldLine;
	int linecount = 0;
	std::ifstream files[inSet.size()];
	for(unsigned i = 0; i < inSet.size(); i++) {
		assert(exists(inSet[i]));
		string filename = inSet[i];
		files[i].open(filename.c_str());
	}
	
	//***find POINTS, write them to outFile
	for (unsigned i = 0; i < inSet.size(); ++i) {
		while ( getline( files[i] , line ) ) {
			unsigned found = line.find("POINTS");
			if (found == 0) {
				if (i == 0) outFile << "POINTS " << 2*nx << " float" << endl;
				break;
			}
			if (i == 0) outFile << line << endl;
			if (i == 0) linecount++ ;
		}
	}	
	while ( getline( files[inSet.size()-1] , line ) ) {
			if (line.length() == 0) {
				outFile << oldLine << endl;
				break;
			}
			oldLine = line;
	}
	files[inSet.size()-1].close();
	files[inSet.size()-1].open(inSet[inSet.size()-1].c_str());
	while ( getline( files[inSet.size()-1] , line ) ) {
		unsigned found = line.find("POINTS");
		if (found == 0) {
			break;
		}
	}
	oldLine = "";
	for (unsigned i = 0; i < inSet.size(); ++i) {
		while ( getline( files[i] , line ) ) {
			if (line.length() == 0) break;
			if (i == inSet.size()-1) outFile << oldLine;
			else outFile << line << '\n';
			if (i == 0) linecount++ ;
			if (i == inSet.size()-1) oldLine = line+'\n';
		}
	}
	outFile << endl;

	//find CELLS, write them to outFile
	for (unsigned i = 0; i < inSet.size(); ++i) {
		while ( getline( files[i] , line ) ) {
			unsigned found = line.find("CELLS");
			if (found == 0) {
				if (i == 0) outFile << "CELLS " << nx << " " << 3*nx << endl;
				break;
			}
			if (i == 0) outFile << line << endl;
			if (i == 0) linecount++ ;
		}
		while ( getline( files[i] , line ) ) {
			if (line.length() == 0) break;
			outFile << line << endl;
			if (i == 0) linecount++ ;
		}
	}
	outFile << endl;

	//find CELL_TYPES, write them to outFile
	for (unsigned i = 0; i < inSet.size(); ++i) {
		while ( getline( files[i] , line ) ) {
			unsigned found = line.find("CELL_TYPES");
			if (found == 0) {
				if (i == 0) outFile << "CELL_TYPES " << nx << endl;
				break;
			}
			if (i == 0) outFile << line << endl;
			if (i == 0) linecount++ ;
		}
		while ( getline( files[i] , line ) ) {
			if (line.length() == 0) break;
			outFile << line << endl;
			if (i == 0) linecount++ ;
		}
	}
	outFile << endl;

	//***find POINT_DATA, write them to outFile
	for (unsigned i = 0; i < inSet.size(); ++i) {
		while ( getline( files[i] , line ) ) {
			unsigned found = line.find("POINT_DATA");
			if (found == 0) {
				if (i == 0) outFile << "POINT_DATA " << 2*nx << endl;
				if (i == 0) outFile << "VECTORS displacements float" << endl;
				getline( files[i] , line );
				break;
			}
			if (i == 0) outFile << line << endl;
			if (i == 0) linecount++ ;
		}
	}
	while ( getline( files[inSet.size()-1] , line ) ) {
			if (line.length() == 0) {
				outFile << oldLine << endl;
				break;
			}
			oldLine = line;
	}
	files[inSet.size()-1].close();
	files[inSet.size()-1].open(inSet[inSet.size()-1].c_str());
	while ( getline( files[inSet.size()-1] , line ) ) {
		unsigned found = line.find("VECTORS displacements float");
		if (found == 0) {
			break;
		}
	}
	oldLine = "";
	for (unsigned i = 0; i < inSet.size(); ++i) {
		while ( getline( files[i] , line ) ) {
			if (line.length() == 0) break;
			if (i == inSet.size()-1) outFile << oldLine;
			else outFile << line << '\n';
			if (i == 0) linecount++ ;
			if (i == inSet.size()-1) oldLine = line+'\n';
		}
	}
	outFile << endl;

	//***find VECTORS velocities float, write them to outFile
	for (unsigned i = 0; i < inSet.size(); ++i) {
		while ( getline( files[i] , line ) ) {
			unsigned found = line.find("VECTORS velocities float");
			if (found == 0) {
				if (i == 0) outFile << "VECTORS velocities float" <<  endl;
				break;
			}
			if (i == 0) outFile << line << endl;
			if (i == 0) linecount++ ;
		}
	}
	while ( getline( files[inSet.size()-1] , line ) ) {
			if (line.length() == 0) {
				outFile << oldLine << endl;
				break;
			}
			oldLine = line;
	}
	files[inSet.size()-1].close();
	files[inSet.size()-1].open(inSet[inSet.size()-1].c_str());
	while ( getline( files[inSet.size()-1] , line ) ) {
		unsigned found = line.find("VECTORS velocities float");
		if (found == 0) {
			break;
		}
	}
	oldLine = "";
	for (unsigned i = 0; i < inSet.size(); ++i) {
		while ( getline( files[i] , line ) ) {
			if (line.length() == 0) break;
			if (i == inSet.size()-1) outFile << oldLine;
			else outFile << line << '\n';
			if (i == 0) linecount++ ;
			if (i == inSet.size()-1) oldLine = line+'\n';
		}
	}
	outFile << endl;

	//find CELL_DATA velocities float, write them to outFile
	for (unsigned i = 0; i < inSet.size(); ++i) {
		while ( getline( files[i] , line ) ) {
			unsigned found = line.find("CELL_DATA");
			if (found == 0) {
				if (i == 0) outFile << "CELL_DATA " << nx << endl;
				if (i == 0) outFile << "SCALARS stress float" << endl;
				if (i == 0) outFile << "LOOKUP_TABLE default" << endl;
				getline( files[i] , line );
				getline( files[i] , line );
				break;
			}
			if (i == 0) outFile << line << endl;
			if (i == 0) linecount++ ;
		}
		while ( getline( files[i] , line ) ) {
			if (line.length() == 0) break;
			outFile << line << endl;
			if (i == 0) linecount++ ;
		}
	}
	outFile << endl;


	for(unsigned i = 0; i < inSet.size(); i++) {
		files[i].close();
	}

	outFile.close();

	return;
}

bool ParallelCombiner::stringsMatch(std::string in1, std::string in2) {

	//make sure both start with ring and then delete it
	unsigned find1 = in1.find("ring_");
	if (find1 < in1.length()) in1.erase(in1.begin(),in1.begin()+find1+5);
	else return false;
	unsigned find2 = in2.find("ring_");
	if (find2 < in2.length()) in2.erase(in2.begin(),in2.begin()+find2+5);
	else return false;

	//make sure both end with /vtk and then delete it
	find1 = in1.find(".vtk");
	if (find1 < in1.length()) in1.erase(in1.end()-4,in1.end());
	else return false;
	find2 = in2.find(".vtk");
	if (find2 < in2.length()) in2.erase(in2.end()-4,in2.end());
	else return false;
	
	//delete processor number
	find1 = in1.find("_");
	if (find1 < in1.length()) in1.erase(in1.begin(),in1.begin()+find1+1);
	else return false;
	find2 = in2.find("_");
	if (find2 < in2.length()) in2.erase(in2.begin(),in2.begin()+find2+1);
	else return false;

	if (in1==in2) return true;

	return false;
}

void ParallelCombiner::combineCohLaw(std::string inPath, bool deleteFlag) {

	string p = inPath;
	vector<string> inSet(0);

	//get file listing
	if (is_directory(p)) {
		for (directory_iterator itr(p); itr!=directory_iterator(); ++itr) {
//			cout << itr->path().filename() << " "; // display filename only
			if (is_regular_file(itr->status()) && itr->path().string().find("cohLaw_")!=std::string::npos ) {
				inSet.push_back(itr->path().string());
			}
		}
	} else {
		//cout << (exists(p) ? "Found: " : "Not found: ") << p << endl;
	}

	if (inSet.size() == 0) return;

	//open target file
	assert(exists(inSet[0]));
	ofstream outFile;

	string outFileName = inSet[0].substr(0,inSet[0].length()-6) + ".dat";
	outFile.open(outFileName.c_str());	

	//open all files in set
	string line;
	int linecount = 0;
	std::ifstream files[inSet.size()];
	for(unsigned i = 0; i < inSet.size(); i++) {
		assert(exists(inSet[i]));
		string filename = inSet[i];
		files[i].open(filename.c_str());
	}

	//***find last header line 
	while ( getline( files[0] , line ) ) {
		unsigned found = line.find("time");
		outFile << line << endl;
		linecount++ ;
		if (found == 8) 	break;
	}

	while (1==1) {
		for(unsigned i = 0; i < inSet.size(); i++) {
			getline( files[i] , line );
			size_t found = line.find("N/A");
			if (line.length() == 0) goto OutOfWhile;

			if (i > 0) outFile << " ";
			if (i == 0) 	outFile << line.substr(0,12);

			if (i == 0 && found!=std::string::npos) {
				//do nothing
			} else if ( i != 0 && found!=std::string::npos) {
				//do nothing
			} else {
				outFile << line.substr(12,line.length());
			}
			linecount++ ;
		}
		outFile << endl;
	}

	OutOfWhile:
	//close files
	for(unsigned i = 0; i < inSet.size(); i++) {
		files[i].close();
	}
	outFile.close();
	cout << "*** cohLaw files combined" << endl;

		//remove unnecessary files
	if (deleteFlag) {												
		for (unsigned j = 0; j < inSet.size(); ++j) {
			string in2 = inSet[j];
//			cout << (exists(in2) ? "Found: " : "Not found: ") << in2 << endl;
			remove(in2);
//			cout << (exists(in2) ? "Found: " : "Not found: ") << in2 << endl;
//			cout << endl;
		}
	}

	return;
}

void ParallelCombiner::combineSTheta(std::string inPath, bool deleteFlag) {

	string p = inPath;
	vector<string> inSet(0);

	//get file listing
	if (is_directory(p)) {
		for (directory_iterator itr(p); itr!=directory_iterator(); ++itr) {
//			cout << itr->path().filename() << " " << endl; // display filename only
			if (is_regular_file(itr->status()) && itr->path().string().find("stresstheta_")!=std::string::npos ) {
				inSet.push_back(itr->path().string());
			}
		}
	} else {
		//cout << (exists(p) ? "Found: " : "Not found: ") << p << endl;
	}

	if (inSet.size() == 0) return;
	if (inSet.size() == 1) {
		remove(inSet[0]);
		cout << "*** no stresstheta data found" << endl;
		return;
	}

	//open target file
	assert(exists(inSet[0]));
	ofstream outFile;

	string outFileName = inSet[0].substr(0,inSet[0].length()-6) + ".dat";
	outFile.open(outFileName.c_str());	

	//open all files in set
	string line;
	string oldLine[inSet.size()];
	int linecount = 0;
	std::ifstream files[inSet.size()];
	for(unsigned i = 0; i < inSet.size(); i++) {
		assert(exists(inSet[i]));
		string filename = inSet[i];
		files[i].open(filename.c_str());
	}

	//***find last header line 
	while ( getline( files[0] , line ) ) {
		unsigned found = line.find("#Time");
		outFile << line << endl;
		linecount++ ;
		if (found == 0) 	break;
	}

	while (1==1) {
		for(unsigned i = 0; i < inSet.size(); i++) {
			outFile << oldLine[i];
			while (getline( files[i] , line )) {

				if (line.length() == 0) {
					getline( files[i] , line );
					if (line.length() == 0 && i == inSet.size()-1) { goto OutOfWhile2;}
					oldLine[i] = line + '\n';
					break;
				}
				oldLine[i] = "";
				outFile << line << endl;
				linecount++ ;
			}
			//NextWhile:
		}
		outFile << endl;		outFile << endl;
	}

	OutOfWhile2:
	//close files
	for(unsigned i = 0; i < inSet.size(); i++) {
		files[i].close();
	}
	outFile.close();
	cout << "*** stresstheta files combined" << endl;

		//remove unnecessary files
	if (deleteFlag) {												
		for (unsigned j = 0; j < inSet.size(); ++j) {
			string in2 = inSet[j];
//			cout << (exists(in2) ? "Found: " : "Not found: ") << in2 << endl;
			remove(in2);
//			cout << (exists(in2) ? "Found: " : "Not found: ") << in2 << endl;
//			cout << endl;
		}
	}

	return;
}

