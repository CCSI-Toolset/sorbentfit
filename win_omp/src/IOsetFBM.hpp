// This is the file that handles  the file I/O for the the MCMC routine of sorbentfit

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/stat.h> // Used for checking the results folder
#include "boost/algorithm/string.hpp"
#include <boost/lexical_cast.hpp>

class IOset{

public:

  IOset(std::string const &confile, std::string const &listfile, time_t &time) {filecheck(confile, listfile, time); filelistset(listfile); configset(confile);};
  void configset(std::string const &confile);
  void filelistset(std::string const &listfile);
  void filecheck(std::string const &confile, std::string const &listfile, time_t &time);
  boost::numeric::ublas::vector<double> lowbounds;
  boost::numeric::ublas::vector<double> highbounds;
  boost::numeric::ublas::vector<double> means;
  boost::numeric::ublas::vector<double> sdev;
  boost::numeric::ublas::vector<int> numsamp;
  boost::numeric::ublas::vector<double> vara;
  boost::numeric::ublas::vector<double> varb;
  std::vector<std::string> sfilelistnames;
  int mcmcsteps, numfiles, rcount;
  int error1, error2, xpoints, ypoints;
  //double vara, varb;

private:
  double delt, taur, taua;
  double Patm, Tatm, rho, por;
  std::vector<std::string> filelist;
};

void IOset::configset(std::string const &confile) {
	means.resize(15);
	lowbounds.resize(15);
	highbounds.resize(15);
	sdev.resize(15);
	vara.resize(3);
	varb.resize(3);
	error1 = 0;
	std::ifstream constream;
	constream.open(confile.c_str(), std::ios::in);
	int i = 0;
	std::string line;
	while (!constream.eof() ) {
		getline(constream, line, '\n');
		//std::cout << line << " " << i << std::endl;
		std::vector<std::string> str;
		if (line[0] != '#' && line[0] != ' ' && line.length() > 0) {
			boost::split(str, line, boost::is_any_of(" "));
			if (i == 0) {
				if (str.size() == 3) {
					delt = boost::lexical_cast<double>(str[0]);
					taur = boost::lexical_cast<double>(str[1]);
					taua = boost::lexical_cast<double>(str[2]);
				}
				else error1 = i+1;
			}
			else if (i == 1) {
				if (str.size() == 4) {
					Patm = boost::lexical_cast<double>(str[0]);
					Tatm = boost::lexical_cast<double>(str[1]);
					rho = boost::lexical_cast<double>(str[2]);
					por = boost::lexical_cast<double>(str[3]);
				}
				else {error1 = i+1;}
			}
			else if (i == 2) {
				if (str.size() == 3) {
					double rad = boost::lexical_cast<double>(str[0]);
					double dl = boost::lexical_cast<double>(str[1]);
					double L = boost::lexical_cast<double>(str[2]);
				}
				else { error1 = i + 1; }
			}
			else if (i == 3) {
				if (str.size() == 3) {
					bool Pdrop = boost::lexical_cast<int>(str[0]);
					bool CO2 = boost::lexical_cast<int>(str[1]);
					bool H2O = boost::lexical_cast<int>(str[2]);
				}
				else { error1 = i + 1; }
			}
			else if (i > 3 && i < 19) {
				if (str.size() == 4) {
					lowbounds[i-4] = boost::lexical_cast<double>(str[0]);
					highbounds[i-4] = boost::lexical_cast<double>(str[1]);
					means[i-4] = boost::lexical_cast<double>(str[2]);
					sdev[i-4] = boost::lexical_cast<double>(str[3]);
				}
				else {error1 =i+1;}
			}
			else if (i == 19) {
				if (str.size() == 3) {
					vara[0] = boost::lexical_cast<double>(str[0]);
					vara[1] = boost::lexical_cast<double>(str[1]);
					vara[2] = boost::lexical_cast<double>(str[2]);
				}
				else { error1 = i + 1; }
			}
			else if (i == 20) {
				if (str.size() == 3) {
					varb[0] = boost::lexical_cast<double>(str[0]);
					varb[1] = boost::lexical_cast<double>(str[1]);
					varb[2] = boost::lexical_cast<double>(str[2]);
				}
				else { error1 = i + 1; }
			}
			else if (i == 21) {
				if (str.size() == 1)
					mcmcsteps = boost::lexical_cast<int>(str[0]);
				else {error1 = i+1;}
			}
			else if (i == 22) {
				if (str.size() == 1)
					rcount = boost::lexical_cast<int>(str[0]);
				else {error1 = i+1;}
			}
			
			i++;
		}
	}

	if (error1 != 0)
		std::cout << "Error with the config_bays file. Check value set " << error1 << std::endl;
}

void IOset::filelistset(std::string const &listfile){
  //Getting data file names
  std::string exdata;
  std::ifstream filestream (listfile.c_str());
  std::vector<std::string> filenames;
  if (filestream.is_open()) {
    filestream >> numfiles;
    filestream >> xpoints >> ypoints;
    filenames.resize(numfiles);
    for (int i = 0; i<numfiles; i++) {
      filestream >> exdata;
      filenames[i]=exdata;
    }
    filestream.close();
  }
  //Now writing subfilelists "sfilelist_X" X being an integer
  //Creating input file name vector
  for (int i = 0; i<filenames.size(); i++){
	  std::string num1;
      std::stringstream convert1;
	  convert1 << i;
	  num1 = convert1.str();
       //Building File Name
	  std::string inputfile;
	  inputfile= "filelist/.sfilelist_";
	  inputfile.append(num1);
	  inputfile.append(".txt");
	  sfilelistnames.push_back(inputfile);
	  
	  //Write subfilelist files
 	  std::ofstream oflfile; // Creates an output file stream to be used for the new filelist file
 	  std::string temp = sfilelistnames[i];
      oflfile.open(temp.c_str());
      oflfile << "1" << std::endl << xpoints << " " << ypoints << std::endl <<  filenames[i] << std::endl;
      oflfile.close();
  }
}
  

void IOset::filecheck(std::string const &confile, std::string const &listfile, time_t &time){
  int numfiles, xpoints;
  int numsamptemp;
  numsamptemp=0;
  error2=0;
  //Checking files------------------------------------------------------------------------------------------
  std::string exdata;
  std::ifstream filestream (listfile.c_str());
  if (filestream.is_open()) {
    filestream >> numfiles;
    filestream >> xpoints >> ypoints;
	numsamp.resize(numfiles);
    for (int i = 0; i<numfiles; i++) {
      filestream >> exdata;
      //Check for file to open
      std::ifstream datstream;
      datstream.open(exdata.c_str(), std::ios::in);
      if(!datstream.is_open()) {
	std::cerr << "The following data file: (" << exdata << ") could not be found. Please review the product manual as to how all files need to be formated and organized." << std::endl;
	error2=1;
      }
      if (error2!=1){
    	  datstream >> numsamptemp;
    	  numsamp[i] = numsamptemp;}
      datstream.close();
    }
    filestream.close();
  }
  else {
    std::cerr << "There is an issue with the filelist_bays file. Please review the product manual as to how all files need to be formated and organized." << std::endl;
    error2=1;
  }
  //Checking if results folder is present
  struct stat SBS;
  if (stat("./results/", &SBS) != 0) 
    system("mkdir results");

  // Code chunk dealing with copying the config file used, timestamping it, and booting it to the results folder
  std::stringstream cstream; // Creates a string to be used in new config file
  std::ofstream ocfile; // Creates an output file stream to be used in new config file
  std::ifstream icfile; // Used for the old config file
  std::string midman; // String that operates as a middleman to copy the buffer over
  cstream.str("");  // Opens string
  cstream << "results/config_FBM_bayes"<<time << ".txt";

  icfile.open(confile.c_str());
  ocfile.open((cstream.str()).c_str()); 

  while (!icfile.eof()) { // Reads until the end of the file
    getline(icfile, midman, '\n');
    ocfile << midman << std::endl;
  }

  icfile.close();
  ocfile.close();

  // Code chunk dealing with copying the filelist file used, stamping it, and booting it to the results folder
  std::stringstream flstream; // Creates a string to be used for the new filelist file
  std::ofstream oflfile; // Creates an output file stream to be used for the new filelist file
  std::ifstream iflfile; // Creates an input file stream for the old filelist file
  flstream.str("");  // Opens string
  flstream << "results/filelist_FBM_bayes"<<time << ".txt";

  iflfile.open(listfile.c_str());
  oflfile.open((flstream.str()).c_str());

  while (!iflfile.eof()) { // Reads until the end of the file
    getline(iflfile, midman, '\n');
    oflfile << midman << std::endl;
  }

  iflfile.close();
  oflfile.close();
}
