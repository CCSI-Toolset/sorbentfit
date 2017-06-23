
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
#include <omp.h>
#include <stdlib.h>
#include <random>
#include <vector>
#include "empfuneval_FBM.hpp"
#include "IOsetFBM.hpp"



int main() {
	//Initialize Program (Send files to IOset for parsing and organization. Create funeval objects.
	bool first = true;
	time_t tim; // Used for creating the time string used to differentiate b/w simulations
	tim = time(NULL);
	std::string confile = "config/config_FBM_bayes.txt"; 
	std::string filelistfile = "filelist/filelist.txt";
	IOset IO(confile, filelistfile, tim);
	//std::cout << "Program Beginning" << std::endl;
	if (IO.error1 > 0 || IO.error2 > 0) {
		std::cout << "There was a problem with file IO" << std::endl;
		std::cin.get();
		return 1;
	}
	int filenum = IO.numfiles;
	int ypoints = IO.ypoints;
	int numparams = IO.means.size(); //Number of parameters
	boost::numeric::ublas::vector<double> nsamp = IO.numsamp;
	int rcount = IO.rcount;
	//std::cout << nsamp << std::endl;
	//Make funeval objects***********************************************************************************************************************************
	std::string sfilelist;
	std::vector<empfuneval> fun(filenum);
	for (int i = 0; i < fun.size(); i++) {
		sfilelist = IO.sfilelistnames[i]; //list of subfilelists for different fun objects
		empfuneval temp(confile, sfilelist);
		temp.bays = true;
		fun[i] = temp;
	}
	//Create stringsd for file output *****************************************************************************************************************************
	//Create name for Bays stat results
	std::stringstream baystream1;
	baystream1.str("");  // Opens string
	baystream1 << "results/BayesParameterResults_" << tim << ".txt";
	//Create name for Bays stat results
	std::stringstream baystream2;
	baystream2.str("");  // Opens string
	baystream2 << "results/BayesAcceptanceStats_" << tim << ".txt";
	//Create name for Bays Psi results
	std::stringstream baystream3;
	baystream3.str("");  // Opens string
	baystream3 << "results/BayesPsiResults_" << tim << ".txt";
	//Retrieve values from IO ****************************************************************************************************************************************
	boost::numeric::ublas::vector<double> lowbounds = IO.lowbounds; //Get low bounds from IO
	boost::numeric::ublas::vector<double> highbounds = IO.highbounds; //Get high bounds from IO
	boost::numeric::ublas::vector<double> means = IO.means; //Get means from IO
	boost::numeric::ublas::vector<double> sdev = IO.sdev; //Get standard deviation from IO
	boost::numeric::ublas::vector<double> nuvect = IO.vara; //Get nu vector from IO
	boost::numeric::ublas::vector<double> tauvect = IO.varb; //Get tau vector from IO
	//Create Bins************************************************************************************************************************************************************
	boost::numeric::ublas::vector<double> draws(numparams); //Vector for draws
	boost::numeric::ublas::vector<double> Pset(numparams); //Vector of curent values of parameters
	boost::numeric::ublas::vector<double> Ptemp(numparams); //Vector of proposed parameters
	boost::numeric::ublas::matrix<double> Rpsimat(rcount, filenum*ypoints);
	boost::numeric::ublas::vector<double> likevect(rcount);
	boost::numeric::ublas::matrix<double> psi(filenum,ypoints);
	boost::numeric::ublas::matrix<double> Rnumat(rcount, filenum);
	boost::numeric::ublas::matrix<double> Rtaumat(rcount, filenum);
	boost::numeric::ublas::matrix<double> numat(filenum, nuvect.size());
	boost::numeric::ublas::matrix<double> taumat(filenum, tauvect.size());
	boost::numeric::ublas::matrix<double> mod_result(filenum, ypoints);
	std::vector<double> Paccepts(numparams);
	std::vector<double> Prates(numparams);
	std::vector<double> PE(filenum);
	std::vector<std::normal_distribution<double>> Pdist(numparams); //Vector of distributions
	std::default_random_engine generator(time(NULL));
	boost::numeric::ublas::matrix<double> results(rcount, numparams); //Matrix to store results
	boost::numeric::ublas::matrix<double> results2(rcount, numparams); //Matrix to store results
	std::vector<double> likely(filenum);
	std::vector<double> likebar(filenum);
	double completelikebar, ratetemp;
	double psitemp;
	double tau, nu;
	bool inbounds, failure, breaker;
	breaker = false;
	int accepter = 0;
	double totlike = 0;
	double completelike = 0;
	double tempdraw, draw;
	boost::numeric::ublas::vector<bool>logten(numparams);
	//CREATING RANDOM NUMBER***********************************************************************************************************************************************************
	std::uniform_real_distribution<double> unif(0, 1);
	std::default_random_engine re(time(NULL));
	double rand = unif(re);
	//Initialze Bins*******************************************************************************************************************************************************************
	//Initialize rate bins to 0
	for (int i = 0; i < numparams; i++) {
		Paccepts[i] = 0;
		Prates[i] = 0;
	}
	//Logten fill in
	for (int i = 0; i < numparams; i++) {
		if (i == 3 || i == 8 || i == 12)
			logten[i] = true;
		else
			logten[i] = false;
	}
	//Create Distributions
	for (int i = 0; i < numparams; i++) {
		std::normal_distribution<double> n(means[i], sdev[i]);
		Pdist[i] = n;
	}

	//Make initial draws************************************************************************************************************************************************************
	for (int i = 0; i < numparams; i++) {
		inbounds = false;
		while (!inbounds) {
			double temp = Pdist[i](generator);
			draws[i] = temp;
			if (lowbounds[i] <= draws[i] && draws[i] <= highbounds[i]) 
				inbounds = true;
		}
	}
	Pset = means; //Set initial set to means
	//Now setup up initial Likely hood calc before paralyzing************************************************************************************************************************
#pragma omp parallel for
	for (int k = 0; k < filenum; k++) {
		row(mod_result,k) = fun[k].objeval(means, logten, 0);
		if (!fun[k].converge) {
			std::cout << "Initialize of Parameters failed. Please choose acceptable starting point. Thank you! " << std::endl;
			breaker = true;
		}
	}
	if (breaker)
		return 1;
	//Calculate initial psi's
	for (int i = 0; i < filenum; i++) {
		for (int j = 0; j < ypoints; j++) {
			nu = (2 * nuvect[j] + nsamp[i]) / 2;
			tau = ((2 * tauvect[j]) + mod_result(i,j)) / 2;
			std::gamma_distribution<double> PSIdist(nu, 1 / tau);
			psitemp = PSIdist(generator);
			psi(i,j)= 1 / psitemp;}
	}
	//Calculate total likelyhood
	for (int k = 0; k < filenum; k++) 
		for (int l = 0; l < ypoints; l++) 
			totlike = totlike + (-nsamp[k] / 2)*log(psi(k, l))+(-mod_result(k,l) / (2 * psi(k,l)));

	completelikebar = totlike;
	//Write initial to file*********************************************************************************************************
	std::cout << "Bayesian Calibration running \n"; // Terminal Output
	std::ofstream ofile1;
	ofile1.open((baystream1.str()).c_str(), std::ios::out);
	if (!ofile1.is_open()) {
		std::cout << "Error writing param result file." << std::endl;
		return 1;
	}
	ofile1 << "Bayesian Calibration running \n"; // Psi File Output
	ofile1.close();

	std::ofstream ofile2;
	ofile2.open((baystream2.str()).c_str(), std::ios::out);
	if (!ofile2.is_open()) {
		std::cout << "Error writing stat result file." << std::endl;
		return 1;
	}
	ofile2 << "Bayesian Calibration running \n"; // Psi File Output
	ofile2.close();

	std::ofstream ofile3;
	ofile3.open((baystream3.str()).c_str(), std::ios::out);
	if (!ofile3.is_open()) {
		std::cout << "Error writing stat result file." << std::endl;
		return 1;
	}
	ofile3 << "Bayesian Calibration running \n"; // Psi File Output
	ofile3.close();

	int count = 0;
	int failcount, rejecter;
	Ptemp = means;
	//Loop for MCMC******************************************************************************************************************
	for (int i = 0; i < IO.mcmcsteps; i++) {
		failcount = 0;
		rejecter = 0;
		Ptemp = Pset;
		for (int j = 0; j < numparams; j++) {
			totlike = 0;
			//Update Jth param
			Ptemp[j] = draws[j];
			//Likelyhood calculation
			failure = false;
#pragma omp parallel for
			for (int k = 0; k < filenum; k++) {
				row(mod_result,k)= fun[k].objeval(Ptemp, logten,0);
				if (!fun[k].converge) {
					failure = true;
				}
			}
			//Calculate total likelyhood
			totlike = 0;
			for (int k = 0; k < filenum; k++)
				for (int l = 0; l < ypoints; l++)
					totlike = totlike + (-nsamp[k] / 2)*log(psi(k, l)) + (-mod_result(k, l) / (2 * psi(k, l)));
			completelike = totlike;
			//std::cout<< totlike <<std::endl;
			//Acceptance Criterion
			if (!failure) {
				if (completelike >= completelikebar) {
					completelikebar = completelike;
					likebar = likely;
					Pset[j] = Ptemp[j];
					Paccepts[j] = Paccepts[j] + 1;
					Prates[j] = Paccepts[j] / (i + 1);
					//std::cout << ratetemp << std::endl;
					accepter++;
				}
				else {
					double test = unif(re);
					double crit = log(test);
					//std::cout << (completelike-completelikebar) << " " << crit << std::endl;
					if ((completelike - completelikebar) > crit) {
						likebar = likely;
						completelikebar = completelike;
						Pset[j] = Ptemp[j];
						accepter++;
						rejecter++;
						Paccepts[j] = Paccepts[j] + 1;
						Prates[j] = Paccepts[j] / (i + 1);
						//std::cout << ratetemp << std::endl;
					}
					//std::cout << "REJECT ACCEPTED!!!!" << std::endl;}
				}
			}
			else { failcount++; }
			Ptemp = Pset;
		}
		for (int k = 0; k < filenum; k++) {
			for (int o = 0; o < ypoints; o++)
				Rpsimat(count, (k*ypoints)+o) = psi(k,o);
		}
		likevect[count] = completelikebar;
		//std::cout << likevect[count] << std::endl;
		//Calculate new psi
		for (int k = 0; k < filenum; k++) {
			for (int o = 0; o < ypoints; o++) {
				nu = (2 * nuvect[o] + nsamp[k]) / 2;
				tau = ((2 * tauvect[o]) + mod_result(k, o)) / 2;
				std::gamma_distribution<double> PSIdist(nu, 1 / tau);
				psitemp = PSIdist(generator);
				psi(k, o) = 1 / psitemp;
			}
		}

		//Write results to matrix
		for (int h = 0; h < Pset.size(); h++) {
			results(count, h) = Pset[h];
			results2(count, h) = Prates[h];
		}
		//Obtain new draws
		for (int k = 0; k < numparams; k++) {
			inbounds = false;
			while (!inbounds) {
				tempdraw = Pdist[k](generator);
				if (tempdraw >= means[k])
					draw = tempdraw + (Pset[k] - means[k]);
				if (tempdraw <= means[k])
					draw = tempdraw - (means[k] - Pset[k]);
				if (lowbounds[k] <= draw && draw <= highbounds[k]) {
					draws[k] = draw;
					inbounds = true;
				}
			}
		}
		//Check if file IO for results is needed
		count++;
		if ((i + 1) % rcount == 0) {
			double Arate = (double)accepter / (numparams * (i + 1));
			std::cout << "Current Step: " << (i + 1) << std::endl << "Accept Rate: " << Arate << std::endl;
			for (int k = 0; k < filenum; k++) {
				PE[k] = fun[k].printevalwithPE(Pset, logten, tim, i + 1, k + 1);
				std::cout << "File Number " << (k + 1) << " Percent Error: " << PE[k] << "%" << std::endl;
			}
			//write param results
			ofile1.open((baystream1.str()).c_str(), std::ios::app);
			for (int m = 0; m < results.size1(); m++) {
				//Write Parameters
				for (int n = 0; n < results.size2(); n++) {
					ofile1 << results(m, n) << " ";
				}//Write Likelyhood
				ofile1 << " " << likevect[m];
				for (int n = 0; n < PE.size(); n++) {
					ofile1 << " " << PE[n];
				}
				ofile1 << std::endl;
			}
			ofile1.close();
			//Write Stat Results
			ofile2.open((baystream2.str()).c_str(), std::ios::app);
			for (int m = 0; m < results2.size1(); m++) {
				for (int n = 0; n < results2.size2(); n++) {
					ofile2 << results2(m, n) << " ";
				}
				ofile2 << std::endl;
			}
			ofile2.close();
			//Writ PSI results
			ofile3.open((baystream3.str()).c_str(), std::ios::app);
			for (int m = 0; m < Rpsimat.size1(); m++) {
				for (int n = 0; n < Rpsimat.size2(); n++) {
					ofile3 << Rpsimat(m, n) << " ";
				}
				ofile3 << std::endl;
			}
			ofile3.close();
			count = 0;
		}
	}
	std::cout << "Calibration Finished! Thank you for using our product!" << std::endl;
	return 0;
}
