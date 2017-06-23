
#define USE_MPI

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
#include <mpi.h>
#include "empfuneval.hpp"
#include "psopoint_ublas_2.hpp"
#include <sys/stat.h> // Used for checking the results folder

int main() {
  int rank, size;
  int exiter=0;
  double dummy;
  MPI::Init();
  size = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();
  time_t tim; // Used for creating the time string used to differentiate b/w simulations
  tim=time(NULL);
  std::stringstream Conv;
  Conv << tim;
  std::string sTime = Conv.str();
  time_t Start;
  time_t End;
  Start=time(NULL); // Used to create a timer to tell how long the code has been running
  
  if (rank == 0) {
    std::cout << "Sorbentfit is running." << std::endl;
    //Checking rest of files
    std::string exdata;
    std::ifstream filestream ("filelist/filelist");
    if (filestream.is_open()) 
      {
	for (int i=0; i<2; ++i)
	  filestream >> dummy;
	while (filestream >> exdata) 
	  {
	      
	    std::ifstream datstream;
	    datstream.open(exdata.c_str(), std::ios::in);
	    if(!datstream.is_open())
	      {
		std::cerr << "The following data file: (" << exdata << ") could not be found. Please review the product manual as to how all files need to be formatted and organized." << std::endl;
		exiter=1;
	      }
	    datstream.close();
	  }
	filestream.close();
      }
    else
      {
	std::cerr << "There is an issue with the filelist file. Please review the product manual as to how all files need to be formatted and organized." << std::endl;
	exiter=1;
      }

    // Checking if config file is present
    std::ifstream cgfile ("config/config");
    if (!cgfile.is_open()) {
      std::cerr << "There is an issue with the config file. Please review the product manual as to how all files need to be formatted and organized."<<std::endl;
      exiter=1;
    }
      
    //Checking if results folder is present
    struct stat SBS;
    if (stat("./results/", &SBS) != 0) {
      if (0 != system("mkdir results")) {
	std::cerr << "There is an issue with the results folder. Please review the product manual as to how all files need to be formatted and organized." << std::endl;
	exiter=1;
      } }
      
    for (int i=1; i<size; ++i) 
      MPI::COMM_WORLD.Send(&exiter, 1, MPI::UNSIGNED_LONG, i, 0);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
 
  for (int i=1; i<size; ++i) 
    {
      if (rank == i)
	MPI::COMM_WORLD.Recv(&exiter, 1, MPI::UNSIGNED_LONG, 0, 0);
    }  
  if (exiter==1)
    {
      MPI::Finalize();
      return 0;
    }

  if (rank==0){
    // Code chunk dealing with copying the config file used, timestamping it, and booting it to the results folder
    //std::stringstream cstream; // Creates a string to be used in new config file
    std::string cstream ("results/config");
    std::ofstream ocfile; // Creates an output file stream to be used in new config file
    std::ifstream icfile; // Used for the old config file
    std::string midman; // String that operates as a middleman to copy the buffer over
    //cstream.str("");  // Opens string
    //cstream << "results/config"<<tim<< std::endl;
    cstream += sTime;
    cstream += ".txt";

    icfile.open("config/config");
    ocfile.open(cstream.c_str(),std::ios::out); 

    while (!icfile.eof()) { // Reads until the end of the file
      getline(icfile, midman, '\n');
      ocfile << midman << std::endl;
    }

    icfile.close();
    ocfile.close();

    // Code chunk dealing with copying the filelist file used, stamping it, and booting it to the results folder
    std::string flstream ("results/filelist"); // Creates a string to be used for the new filelist file
    std::ofstream oflfile; // Creates an output file stream to be used for the new filelist file
    std::ifstream iflfile; // Creates an input file stream for the old filelist file
    //flstream.str("");  // Opens string
    //flstream << "results/filelist"<<tim<< std::endl;
    flstream += sTime;
    flstream += ".txt";

    iflfile.open("filelist/filelist");
    oflfile.open(flstream.c_str(),std::ios::out);

    while (!iflfile.eof()) { // Reads until the end of the file
      getline(iflfile, midman, '\n');
      oflfile << midman << std::endl;
    }

    iflfile.close();
    oflfile.close();
  }

  using namespace boost::numeric::ublas;
  typedef vector<double>::iterator uvi_t;
  typedef matrix<double>::iterator1 umi1_t;
  typedef matrix<double>::iterator2 umi2_t;
  typedef vector<bool>::iterator uvbi_t;

  uint32_t seed;
  if (rank == 0) {
    boost::lagged_fibonacci44497 randfib(time(NULL));
    boost::uniform_int<uint32_t> uintdist(0, std::numeric_limits<uint32_t>::max());
    boost::variate_generator< boost::lagged_fibonacci44497&, boost::uniform_int<uint32_t> > getrand(randfib, uintdist);
    for (int i=1; i<size; ++i) {
      seed = getrand();
      MPI::COMM_WORLD.Send(&seed, 1, MPI::UNSIGNED_LONG, i, 0);
    }
    seed = getrand();
  }
  for (int i=1; i<size; ++i) {
    if (rank == i)
      MPI::COMM_WORLD.Recv(&seed, 1, MPI::UNSIGNED_LONG, 0, 0);
  }
    
  bool breakme = false;

  empfuneval fun("config/config", "filelist/filelist");
  int numagent = fun.nodesize*size;
  vector<double> pos(fun.numparams);
  vector<psopoint> agent(fun.nodesize);

  std::string ostream ("results/optresults");
  std::ofstream ofile; // Creates an output file stream to be used in the iteration file
  //ostream.str("");  // This set of code changes the output to a text file located in the results
  //ostream << "results/optresults"<<tim<< std::endl; // Names the file with a time stamp that will be unique to each simulation
  ostream += sTime;
  ostream += ".txt";
  
  if (rank == 0) {
    if (fun.error != 0) {
      std::cout << "Error with the config file. Check value set" << fun.error << std::endl;
      breakme = true;
    }
    if (!breakme) {
      std::cout << "Running on " << size << " Nodes" << std::endl;
      std::cout << fun.paramorder << std::endl;

      ofile.open(ostream.c_str(),std::ios::out); // Creates a file with the name of the of the above defined stream
      ofile << "Running on " << size << " Nodes" << std::endl;
      ofile << fun.paramorder << std::endl;
      //ofile.close();
    }
  }
    
  psopoint::costvalist.resize(numagent, pos.size());
  psopoint::poslist.resize(numagent, pos.size());
  psopoint::posbestlist.resize(numagent, pos.size());
  psopoint::logten.resize(pos.size());
    
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

  boost::mt19937 randmt(seed);
  boost::uniform_01<> u01dist;
  boost::variate_generator< boost::mt19937&, boost::uniform_01<> > randstrm(randmt, u01dist);
	
  for (int i = 0; i < agent.size(); i++) {
    int k = 0;
    for (int j = 0; j != 13; j++) {
      if (!fun.override(j)) {
	pos[k] = fun.lowbounds[j] + randstrm()*(fun.highbounds[j] - fun.lowbounds[j]);
	k++;
      }
    }
    agent[i].init(pos, rank*fun.nodesize+i, numagent+size, fun); // Trouble b/w lines 169-179. High probab 177?
  }

  size_t bestnum, bestnum_it;
  long numit=0;
  double costmax, costmin;

  //ofile.open((ostream.str()+".txt").c_str(), std::ios::app);

  while (1)  {

    for (int i=0; i<size; ++i) {

      MPI::COMM_WORLD.Bcast((psopoint::costvalist.data().begin()+i*fun.nodesize), fun.nodesize, MPI::DOUBLE, i);
      MPI::COMM_WORLD.Bcast((psopoint::poslist.data().begin()+i*fun.nodesize*pos.size()), fun.nodesize*pos.size(), MPI::DOUBLE, i);
      MPI::COMM_WORLD.Bcast((psopoint::posbestlist.data().begin()+i*fun.nodesize*pos.size()), fun.nodesize*pos.size(), MPI::DOUBLE, i);
    }

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

      std::cout << numit << " " << costmin << " " << costmax << "\n" << project(psopoint::posbestlist, range(bestnum, bestnum+1), range(0,pos.size())) << "\n";
      // ++numit;  // These lines change the iteration output to the terminal

      // ofile.open((ostream.str()).c_str(), std::ios::app);
      ofile << numit << " " << costmin << " " << costmax << "\n" << project(psopoint::posbestlist, range(bestnum, bestnum+1), range(0,pos.size())) << std::endl;
      //ofile.close(); // Check here, may be able to speed code up?
      ++numit;
      // Note if you have both terminal & file out active, keep only the last ++numit uncommented

      if ((costmax - costmin)/costmin < 0.01)
	breakme = true;
    }

    MPI::COMM_WORLD.Bcast(&breakme, 1, MPI::BOOL, 0);
    if (breakme)
      break;

    for (size_t i=0; i<fun.nodesize; i++)
      agent[i].move(fun, randstrm);
  }

  if (rank == 0) {

    std::cout << "The Final Parameters are as follows: \n";
    std::cout << project(psopoint::posbestlist, range(bestnum, bestnum+1), range(0,pos.size())) << "\n"; // Writes final parameters into the terminal

    //ofile.open((ostream.str()+".txt").c_str(), std::ios::app);
    ofile << "The Final Parameters are as follows:" << std::endl;
    ofile << project(psopoint::posbestlist, range(bestnum, bestnum+1), range(0,pos.size())) << "\n"; // Writes final parameters into text file

    End=time(NULL);
    double Dur= (End-Start);
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
  MPI::Finalize();
  return 0;
}
