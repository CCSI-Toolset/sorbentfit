

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
#include "empfuneval_FBM.hpp"
#include "boost/algorithm/string.hpp"
#include <boost/lexical_cast.hpp>
//#include "omp.h"
#include <sys/stat.h> // Used for checking the results folder

int main() {
	std::cout << "Sorbentfit ParamSet is running." << std::endl;
	double dummy;
	time_t tim; // Used for creating the time string used to differentiate b/w simulations
	tim = time(NULL);
	std::stringstream Conv;
	Conv << tim;
	std::string sTime = Conv.str();
	time_t Start;
	time_t End;
	Start = time(NULL); // Used to create a timer to tell how long the code has been running
	//Checking files
	std::string exdata;
	std::ifstream filestream("filelist/filelist.txt");
	if (filestream.is_open()) {
		for (int i = 0; i < 3; ++i)
			filestream >> dummy;
		while (filestream >> exdata) {
			std::ifstream datstream;
			datstream.open(exdata.c_str(), std::ios::in);
			if (!datstream.is_open()) {
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
	std::ifstream cgfile("config/config_FBM_empfit.txt");
	if (!cgfile.is_open()) {
		std::cerr << "There is an issue with the config file. Please review the product manual as to how all files need to be formatted and organized." << std::endl;
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
		}
	}

	// Code chunk dealing with copying the config file used, timestamping it, and booting it to the results folder
	//std::stringstream cstream; // Creates a string to be used in new config file
	std::string cstream("results/config_FBM_paramset");
	std::ofstream ocfile; // Creates an output file stream to be used in new config file
	std::ifstream icfile; // Used for the old config file
	std::string midman; // String that operates as a middleman to copy the buffer over
	//cstream.str("");  // Opens string
	//cstream << "results/config"<<tim<< std::endl;
	cstream += sTime;
	cstream += ".txt";

	icfile.open("config/config_FBM_empfit.txt");
	ocfile.open(cstream.c_str(), std::ios::out);

	while (!icfile.eof()) { // Reads until the end of the file
		getline(icfile, midman, '\n');
		ocfile << midman << std::endl;
	}

	icfile.close();
	ocfile.close();

	// Code chunk dealing with copying the filelist file used, stamping it, and booting it to the results folder
	std::string flstream("results/filelist_FBM_paramset"); // Creates a string to be used for the new filelist file
	std::ofstream oflfile; // Creates an output file stream to be used for the new filelist file
	std::ifstream iflfile; // Creates an input file stream for the old filelist file
	//flstream.str("");  // Opens string
	//flstream << "results/filelist"<<tim<< std::endl;
	flstream += sTime;
	flstream += ".txt";


	iflfile.open("filelist/filelist.txt");
	oflfile.open(flstream.c_str(), std::ios::out);

	while (!iflfile.eof()) { // Reads until the end of the file
		getline(iflfile, midman, '\n');
		oflfile << midman << std::endl;
	}

	iflfile.close();
	oflfile.close();

	//Break variable to know when to exit
	bool breakme = false;
	int numset;
	empfuneval fun("config/config_FBM_empfit.txt", "filelist/filelist.txt");
	fun.bays = true;
	boost::numeric::ublas::vector<bool> logger(15);
	boost::numeric::ublas::vector<double> params(15);
	boost::numeric::ublas::matrix<double>param_mat;
	for (int i = 0; i != 15; i++) {
		if (i == 3 || i == 8 || i == 13)
			logger(i) = 1;
		else
			logger(i) = 0;
	}
	//Get Parameters and run printeval
	int error = 0;
	std::ifstream constream;
	constream.open("data/Paramset_FBM.txt", std::ios::in);
	int i = 0;
	std::string line;
	getline(constream, line, '\n');
	std::vector<std::string> str;
	boost::split(str, line, boost::is_any_of(" "));
	if (str.size() == 1) {
		numset = boost::lexical_cast<int>(str[0]);
	}
	else {
		std::cerr << "Paramset_FBM not formated properly" << std::endl;
		std::cin.get();
		return 1;
	}
	param_mat.resize(numset, 15);
	while (!constream.eof()) {
		getline(constream, line, '\n');
		if (line[0] != '#' && line[0] != ' ' && line.length() > 0) {
			boost::split(str, line, boost::is_any_of(" "));
            //std::cout << str.size() << std::endl;
			if (str.size() == 15) {
				for (int m = 0; m < 15; m++)
					params[m] = boost::lexical_cast<double>(str[m]);
				row(param_mat, i) = params;
			}
			else
				std::cout << "The Set " << i << " had an issue..." << std::endl;
			i++;
		}
		if (i == 0)
			std::cout << "There was an issue reading Paramset file " << std::endl;
	}
	for (int j = 0; j < numset; j++) {
		std::cout << "Evaluating Paramset: " << j << std::endl;
		params = row(param_mat, j);
		fun.printeval(params, logger, j);	
	}
	std::cout << "ParamSet is finished. Thank You!" << std::endl;
	std::cout << "To exit, press enter or return." << std::endl;
	std::cin.ignore();

	return 0;
}
