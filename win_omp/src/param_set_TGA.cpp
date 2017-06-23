
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
#include "empfuneval_TGA.hpp"
#include "boost/algorithm/string.hpp"
#include <boost/lexical_cast.hpp>
#include <sys/stat.h> // Used for checking the results folder

int main() {
  std::cout << "Sorbentfit ParamSet is running." << std::endl;
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
  std::string cstream ("results/config_TGA_paramset");
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
  std::string flstream ("results/filelist_TGA_paramset"); // Creates a string to be used for the new filelist file
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

  //Random Stuff
  uint32_t seed;
  boost::lagged_fibonacci44497 randfib(time(NULL));
  boost::uniform_int<uint32_t> uintdist(0, std::numeric_limits<uint32_t>::max());
  boost::variate_generator< boost::lagged_fibonacci44497&, boost::uniform_int<uint32_t> > getrand(randfib, uintdist);
  seed = getrand();
  boost::mt19937 randmt(seed);
  boost::uniform_01<> u01dist;
  boost::variate_generator< boost::mt19937&, boost::uniform_01<> > randstrm(randmt, u01dist);
  //Break variable to know when to exit
  bool breakme = false;

  empfuneval fun("config/config_TGA_empfit.txt", "filelist/filelist.txt");
  fun.bays=true;
  boost::numeric::ublas::vector<bool> logger(13);
  boost::numeric::ublas::vector<double> params(13);
  for (int i = 0; i != 13; i++) {
    if (i == 3 || i == 8 || i == 12)
      logger(i) = 1;
    else
      logger(i) = 0;}
  //Get Parameters and run printeval
  int error = 0;
  std::ifstream constream;
  constream.open("data/Paramset_TGA.txt", std::ios::in);
  int i = 0;
  std::string line;
  while (!constream.eof() ) {
    getline(constream, line, '\n');
    std::vector<std::string> str;
    if (line[0] != '#' && line[0] != ' ' && line.length() > 0) {
      boost::split(str, line, boost::is_any_of(" "));
      if (str.size() == 13) {
	for(int m = 0; m<13; m++)
	  params[m] = boost::lexical_cast<double>(str[m]);
	fun.printeval(params, logger, i);
	std::cout << "Evaluating Paramset: " << i << std::endl;
      }
      else
	std::cout << "The Set " << i << " had an issue..." << std::endl;
      i++;
    }
    if(i==0)
      std::cout << "There was an issue reading Paramset file " << std::endl;
  }
  std::cout << "ParamSet is finished. Thank You!" << std::endl;
  std::cout << "To exit, press enter or return." << std::endl;
  std::cin.ignore();
  
  return 0;
}
