
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/cstdint.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include "omp.h"
#include "empfuneval_TGA.hpp"
#include "psopoint_ublas_2.hpp"
#include <sys/stat.h> // Used for checking the results folder

int main() {
  std::cout << "Sorbentfit is running." << std::endl;
  double dummy;
  time_t tim; // Used for creating the time string used to differentiate b/w simulations
  tim=time(NULL);
  std::stringstream Conv;
  Conv << tim;
  std::string sTime = Conv.str();
  time_t Start;
  time_t End;
  Start=time(NULL); // Used to create a timer to tell how long the code has been running
  //Checking files
  std::string exdata;
  std::ifstream filestream ("filelist/filelist.txt");
  if (filestream.is_open()) {
    for (int i=0; i<3; ++i)
      filestream >> dummy;
    while (filestream >> exdata) {
      std::ifstream datstream;
      datstream.open(exdata.c_str(), std::ios::in);
      if(!datstream.is_open()) {
	std::cerr << "The following data file: (" << exdata << ") could not be found. Please review the product manual as to how all files need to be formatted and organized." << std::endl;
	std::cout << "To exit, press enter or return." << std::endl;
	std::cin.ignore();
	return 1;
      }
      datstream.close();
    }
    filestream.close();
  }
  else {
    std::cerr << "There is an issue with the filelist file. Please review the product manual as to how all files need to be formatted and organized." << std::endl;
    std::cout << "To exit, press enter or return." << std::endl;
    std::cin.ignore();
    return 1;
  }

  // Checking if config file is present
  std::ifstream cgfile ("config/config_TGA_empfit.txt");
  if (!cgfile.is_open()) {
    std::cerr << "There is an issue with the config file. Please review the product manual as to how all files need to be formatted and organized."<<std::endl;
    std::cout << "To exit, press enter or return." << std::endl;
    std::cin.ignore();
    return 1;
  }

  //Checking if results folder is present
  struct stat SBS;
  if (stat("./results/", &SBS) != 0) {
    if (0 != system("mkdir results")) {
      std::cerr << "There is an issue with the results folder. Please review the product manual as to how all files need to be formatted and organized." << std::endl;
      std::cout << "To exit, press enter or return." << std::endl;
      std::cin.ignore();
      return 1;
    } }

  // Code chunk dealing with copying the config file used, timestamping it, and booting it to the results folder
  //std::stringstream cstream; // Creates a string to be used in new config file
  std::string cstream ("results/config_TGA_");
  std::ofstream ocfile; // Creates an output file stream to be used in new config file
  std::ifstream icfile; // Used for the old config file
  std::string midman; // String that operates as a middleman to copy the buffer over
  //cstream.str("");  // Opens string
  //cstream << "results/config"<<tim<< std::endl;
  cstream += sTime;
  cstream += ".txt";

  icfile.open("config/config_TGA_empfit.txt");
  ocfile.open(cstream.c_str(),std::ios::out);

  while (!icfile.eof()) { // Reads until the end of the file
    getline(icfile, midman, '\n');
    ocfile << midman << std::endl;
  }

  icfile.close();
  ocfile.close();

  // Code chunk dealing with copying the filelist file used, stamping it, and booting it to the results folder
  std::string flstream ("results/filelist_TGA_"); // Creates a string to be used for the new filelist file
  std::ofstream oflfile; // Creates an output file stream to be used for the new filelist file
  std::ifstream iflfile; // Creates an input file stream for the old filelist file
  //flstream.str("");  // Opens string
  //flstream << "results/filelist"<<tim<< std::endl;
  flstream += sTime;
  flstream += ".txt";


  iflfile.open("filelist/filelist.txt");
  oflfile.open(flstream.c_str(),std::ios::out);

  while (!iflfile.eof()) { // Reads until the end of the file
    getline(iflfile, midman, '\n');
    oflfile << midman << std::endl;
  }

  iflfile.close();
  oflfile.close();

  using namespace boost::numeric::ublas;
  typedef vector<double>::iterator uvi_t;
  typedef matrix<double>::iterator1 umi1_t;
  typedef matrix<double>::iterator2 umi2_t;
  typedef vector<bool>::iterator uvbi_t;

  uint32_t seed;
  boost::lagged_fibonacci44497 randfib(time(NULL));
  boost::uniform_int<uint32_t> uintdist(0, std::numeric_limits<uint32_t>::max());
  boost::variate_generator< boost::lagged_fibonacci44497&, boost::uniform_int<uint32_t> > getrand(randfib, uintdist);
  seed = getrand();
  boost::mt19937 randmt(seed);
  boost::uniform_01<> u01dist;
  boost::variate_generator< boost::mt19937&, boost::uniform_01<> > randstrm(randmt, u01dist);

  bool breakme = false;
  size_t bestnum, bestnum_it;
  long numit=0;
  double costmax, costmin;
  int numagent;
  int scale_mode = 0;

  empfuneval fun2("config/config_TGA_empfit.txt", "filelist/filelist.txt");
  scale_mode = fun2.scale_mode;
  int chec = 0;
  bool PSO=true;
  boost::numeric::ublas::vector<bool> logger(13);
  for (int i = 0; i != 13; i++) {
    if (fun2.override[i]) {
      chec++;}
    if (i == 3 || i == 8 || i == 12)
      logger(i) = 1;
    else
      logger(i) = 0;
  }

  if(chec==13)
    PSO=false;

  if(!PSO){
    std::cout << "The mode selected is a parameter evaluation. PSO will not be engaged! \n";
    fun2.printeval(fun2.params, logger, tim);
    std::cout << "To exit, press enter or return." << std::endl;
    std::cin.ignore();
    return 0;}
  //std::stringstream ostream; // Creates a string to be used in the iteration file
  std::string ostream ("results/optresults");
  std::ofstream ofile; // Creates an output file stream to be used in the iteration file
  //ostream.str("");  // This set of code changes the output to a text file located in the results
  //ostream << "results/optresults"<<tim<< std::endl; // Names the file with a time stamp that will be unique to each simulation
  ostream += sTime;
  ostream += ".txt";
  ofile.open(ostream.c_str(), std::ios::out); // Creates a file with the name of the of the above defined stream

#pragma omp parallel
  {
    empfuneval fun("config/config_TGA_empfit.txt", "filelist/filelist.txt");
    int size = omp_get_num_threads();
    int rank = omp_get_thread_num();
    int numagent = fun.nodesize*size;
    vector<double> pos(fun.numparams);
    vector<psopoint> agent(fun.nodesize);
	scale_mode = fun.scale_mode;

    if (rank == 0) {
      if (fun.error != 0) {
	std::cout << "Error with the config file. Check value set " << fun.error << std::endl;
	breakme = true;
      }
      std::cout << "Running on " << size << " Nodes" << std::endl; // Terminal Output
      std::cout << fun.paramorder << std::endl;

      ofile << "Running on " << size << " Nodes" << std::endl; // File Output
      ofile << fun.paramorder << std::endl;

      psopoint::costvalist.resize(numagent);
      psopoint::poslist.resize(numagent, pos.size());
      psopoint::posbestlist.resize(numagent, pos.size());
      psopoint::logten.resize(pos.size());
	  //for (int i = 0; i < agent.size(); i++) Caused parallel issues, so i removed it
		  //agent[i].scale_mode = scale_mode;
      for (uvbi_t logten_it=psopoint::logten.begin();
	   logten_it != psopoint::logten.end();
	   ++logten_it)
	*logten_it = 0;

      int k = 0;
      for (int i = 0; i != 13; i++) {
	if (fun.override[i]) {
	  k++;
	}
	else if (i == 3 || i == 8 || i == 12)
	  psopoint::logten(i-k) = 1;
      }
    }
    //std::cout << psopoint::logten << std::endl;
#pragma omp barrier
    for (int i = 0; i < agent.size(); i++) {
      int k = 0;
      for (int j = 0; j != 13; j++) {
	if (!fun.override(j)) {
	  pos[k] = fun.lowbounds[j] + randstrm()*(fun.highbounds[j] - fun.lowbounds[j]);
	  k++;
	}
      }
      agent[i].init(pos, rank*fun.nodesize+i, numagent+size, fun);
    }

    while (1) {
#pragma omp barrier

      if (rank == 0) {

	bestnum_it = 0;
	bestnum = 0;
	costmax = 0;
	costmin = *psopoint::costvalist.begin();

	for (uvi_t costvalist_it=psopoint::costvalist.begin();
	     costvalist_it != psopoint::costvalist.end();
	     ++costvalist_it, ++bestnum_it) {

	  if (*costvalist_it > costmax)
	    costmax = *costvalist_it;

	  if (*costvalist_it < costmin) {
	    costmin = *costvalist_it;
	    bestnum = bestnum_it;
	  }
	}

	std::cout << numit << " " << costmin << " " << costmax << "\n" << project(psopoint::posbestlist, range(bestnum, bestnum+1), range(0,fun.numparams)) << "\n"; // Terminal Output
	// ++numit;

	ofile << numit << " " << costmin << " " << costmax << "\n" << project(psopoint::posbestlist, range(bestnum, bestnum+1), range(0,pos.size())) << std::endl; // File Output
	++numit;
	// Note if you have both terminal & file out active, keep only the last ++numit uncommented

	if ((costmax - costmin)/costmin < 0.01)
	  breakme = true;
      }

      if (breakme)
	break;

      for (size_t i=0; i<fun.nodesize; i++)
	agent[i].move(fun, randstrm);
    }
    if (rank == 0) {
      End=time(NULL); // Used to create a timer to tell how long the code has been running
      double Dur=(End-Start);

      std::cout << "The Final Parameters are as Follows : \n";
      std::cout << project(psopoint::posbestlist, range(bestnum, bestnum+1), range(0,pos.size())) << "\n";

      ofile << "The Final Parameters are as Follows:" << std::endl;
      ofile << project(psopoint::posbestlist, range(bestnum, bestnum+1), range(0,pos.size())) << "\n"; // Writes final parameters into text file

      if (Dur>60 && Dur<3600) { // Checks if the run time is in minutes
	Dur = Dur/60;
	std::cout << "Sorbentfit has been running for : " << Dur << " minutes." << std::endl;
	ofile << "Sorbentfit has been running for : " << Dur << " minutes." << std::endl;
      }
      else if (Dur>3600) { // Run time is in hours
	Dur = Dur/3600;
	std::cout << "Sorbentfit has been running for : " << Dur << " hours." << std::endl;
	ofile << "Sorbentfit has been running for : " << Dur << " hours." << std::endl;
      }
      ofile.close();

      fun.printeval(row(psopoint::posbestlist, bestnum), psopoint::logten, tim);
    }
  }
  std::cout << "To exit, press enter or return." << std::endl;
  std::cin.ignore();
  return 0;
}
