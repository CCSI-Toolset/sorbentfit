

#include "funeval_base.hpp" 
#include "boost/algorithm/string.hpp"
#include <boost/lexical_cast.hpp>
#include "psopoint_ublas_2.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>


class empfuneval : public funeval_base {

public:

	empfuneval(std::string const &confile, std::string const &datfiles) {configset(confile); datfile(datfiles); counter = 0; };
	empfuneval(std::string const &confile) {configset(confile); datload = false; counter = 0; };
	empfuneval() { confload = false; datload = false; counter = 0; };
	void configset(std::string const &confile);

	boost::numeric::ublas::vector<bool> override;
	boost::numeric::ublas::vector<double> lowbounds;
	boost::numeric::ublas::vector<double> highbounds;
	boost::numeric::ublas::vector<double> params;
	std::string paramorder;
	int numparams;
	int nodesize, scale_mode;
	int error;
	bool bays;

protected:

	bool function(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::matrix<double> &ypass);

private:
	double delt, taur, taua;
	double nv, nv_l, nv_h;
	double nh, nh_l, nh_h;
	double dl, L, por, mob, Vv, Vs,Vt, area; //Variables for FBM
	int numcv;
	double rho;
	double dh_kapc, ds_kapc, dh_kc, zeta_kc; //Carbamate Rate Constants
	double dh_kaph, ds_kaph, dh_kh, zeta_kh; //Water Rate Constants
	double dh_kapb, ds_kapb, dh_kb, zeta_kb; //Bicabonate Rate Constants

	// Equilibrium constants and rate constants
	double kapc, kapb, kc, kb, kaph, kh;
	
	// Switches for outputs
	bool Pdrop, H2O, CO2;

	// The temperature and pressures and inputs
	double T, Patm, Tatm, Fc1, Fh1, Fn1, Fin1, Fc2, Fh2, Fn2, Fin2;
	double P = 101325;
	// Gas constant
	double R;
	
	// Computational variables
	double nv_inv,nh_inv,st1,st2;
	double Rx1, Rx2, Ra1, Ra2, Rb1, Rb2;
	boost::numeric::ublas::vector<double> X1, A1, B1, C1, H1, P1, Conc;

	// Jacobian
	int jacsize;
	boost::numeric::ublas::vector<double> A;
	boost::numeric::ublas::vector<int> IA, JA;

	void f(boost::numeric::ublas::vector<double> const &X2, boost::numeric::ublas::vector<double> const &A2, boost::numeric::ublas::vector<double> const &B2, boost::numeric::ublas::vector<double> const &C2, boost::numeric::ublas::vector<double> const &H2, boost::numeric::ublas::vector<double> const &P2, boost::numeric::ublas::vector<double> &res);

	void jacsolve(boost::numeric::ublas::vector<double> const &X2, boost::numeric::ublas::vector<double> const &A2, boost::numeric::ublas::vector<double> const &B2, boost::numeric::ublas::vector<double> const &C2, boost::numeric::ublas::vector<double> const &H2, boost::numeric::ublas::vector<double> const &P2, boost::numeric::ublas::vector<double> const &bb, boost::numeric::ublas::vector<double> &s);


};

void empfuneval::configset(std::string const &confile) {
	override.resize(15);
	lowbounds.resize(15);
	highbounds.resize(15);
	params.resize(15);
	error = 0;
	std::ifstream constream;
	constream.open(confile.c_str(), std::ios::in);
	int i = 0;
	std::string line;
	while (!constream.eof()) {
		getline(constream, line, '\n');
		std::vector<std::string> str;
		if (line[0] != '#' && line[0] != ' ' && line.length() > 0) {
			boost::split(str, line, boost::is_any_of(" "));
			if (i == 0) {
				if (str.size() == 3) {
					delt = boost::lexical_cast<double>(str[0]);
					taur = boost::lexical_cast<double>(str[1]);
					taua = boost::lexical_cast<double>(str[2]);
				}
				else error = i + 1;
			}
			else if (i == 1) {
				if (str.size() == 4) {
					Patm = boost::lexical_cast<double>(str[0]);
					Tatm = boost::lexical_cast<double>(str[1]);
					rho = boost::lexical_cast<double>(str[2]);
					por = boost::lexical_cast<double>(str[3]);
				}
				else { error = i + 1; }
			}
			else if (i == 2) {
				if (str.size() == 3) {
					area = boost::lexical_cast<double>(str[0]);
					dl = boost::lexical_cast<double>(str[1]);
					L = boost::lexical_cast<double>(str[2]);
				}
				else { error = i + 1; }
			}
			else if (i == 3) {
				if (str.size() == 3) {
					Pdrop = boost::lexical_cast<int>(str[0]);
					CO2 = boost::lexical_cast<int>(str[1]);
					H2O = boost::lexical_cast<int>(str[2]);
				}
				else { error = i + 1; }
			}
			else if (i > 3 && i < 19) {
				if (str.size() == 1) {
					params[i - 4] = boost::lexical_cast<double>(str[0]);
					override[i - 4] = true;
				}
				else if (str.size() == 2) {
					lowbounds[i - 4] = boost::lexical_cast<double>(str[0]);
					highbounds[i - 4] = boost::lexical_cast<double>(str[1]);
				}
				else { error = i + 1; }
			}
			else if (i == 19) {
				if (str.size() == 1)
					nodesize = boost::lexical_cast<int>(str[0]);
				else { error = i + 1; }
			}
			else if (i == 20) {
				if (str.size() == 1)
					scale_mode = boost::lexical_cast<int>(str[0]);
				else { error = i + 1; }
			}
			i++;
		}
	}
	constream.close();
	numparams = 0;
	bool first = true;
	for (int i = 0; i != 15; i++) {
		if (!override[i]) {
			numparams++;
			if (!first) {
				paramorder += " | ";
			}
			else first = false;
			switch (i) {
			case 0:
				paramorder += "dh_kapc";
				break;
			case 1:
				paramorder += "ds_kapc";
				break;
			case 2:
				paramorder += "dh_kc";
				break;
			case 3:
				paramorder += "zeta_kc";
				break;
			case 4:
				paramorder += "nv";
				break;
			case 5:
				paramorder += "dh_kaph";
				break;
			case 6:
				paramorder += "ds_kaph";
				break;
			case 7:
				paramorder += "dh_kh";
				break;
			case 8:
				paramorder += "zeta_kh";
				break;
			case 9:
				paramorder += "nh";
				break;
			case 10:
				paramorder += "dh_kapb";
				break;
			case 11:
				paramorder += "ds_kapb";
				break;
			case 12:
				paramorder += "dh_kb";
				break;
			case 13:
				paramorder += "zeta_kb";
				break;
			case 14:
				paramorder += "mobility";
				break;
			}
		}
		else {
			switch (i) {
			case 0:
				dh_kapc = params[i];
				break;
			case 1:
				ds_kapc = params[i];
				break;
			case 2:
				dh_kc = params[i];
				break;
			case 3:
				zeta_kc = pow(10.0, params[i]);
				break;
			case 4:
				nv = params[i];
				break;
			case 5:
				dh_kaph = params[i];
				break;
			case 6:
				ds_kaph = params[i];
				break;
			case 7:
				dh_kh = params[i];
				break;
			case 8:
				zeta_kh = pow(10.0, params[i]);
				break;
			case 9:
				nh = params[i];
				break;
			case 10:
				dh_kapb = params[i];
				break;
			case 11:
				ds_kapb = params[i];
				break;
			case 12:
				dh_kb = params[i];
				break;
			case 13:
				zeta_kb = pow(10.0, params[i]);
				break;
			case 14:
				mob = params[i];
				break;
			}
		}
	}
	
	// Calculate Number of CV's
	double temp = (L / dl);
	numcv = floor(temp);
	//std::cout << "NumCV = " << numcv << std::endl;
	// Resize Global vectors to correct size
	//Creates Jacobian matrix, then passes to gmres solver to obtain search direction
	jacsize = 6 * numcv * 6 * numcv;
	JA.resize(jacsize);
	IA.resize(jacsize);
	A.resize(jacsize);
	X1.resize(numcv);
	A1.resize(numcv);
	B1.resize(numcv);
	C1.resize(numcv);
	H1.resize(numcv);
	P1.resize(numcv);
	Conc.resize(numcv);
	confload = true;
	bays = false;
}

bool empfuneval::function(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::matrix<double> &ypass) {
	int j = 0;
	bool inbounds = true;
	for (int i = 0; i != 15; i++) {
		if (!bays) {
			if (!override[i]) {
				if ((par[j] < lowbounds[i] || par[j] > highbounds[i]) && !psopoint::logten[j]) {
					//std::cout << lowbounds[i] << " " << par[j] << " " << highbounds[i] << std::endl;
					inbounds = false;
					break;
				}
				else if ((par[j] < pow(10.0, lowbounds[i]) || par[j] > pow(10.0, highbounds[i])) && psopoint::logten[j]) {
					//std::cout << lowbounds[i] << " " << i << " " << par[j] << " " << j << " " << highbounds[i] << std::endl;
					//std::cout << psopoint::logten << std::endl;
					inbounds = false;
					break;
				}

				switch (i) {
				case 0:
					dh_kapc = par[j];
					break;
				case 1:
					ds_kapc = par[j];
					break;
				case 2:
					dh_kc = par[j];
					break;
				case 3:
					zeta_kc = par[j];
					break;
				case 4:
					nv = par[j];
					break;
				case 5:
					dh_kaph = par[j];
					break;
				case 6:
					ds_kaph = par[j];
					break;
				case 7:
					dh_kh = par[j];
					break;
				case 8:
					zeta_kh = par[j];
					break;
				case 9:
					nh = par[j];
					break;
				case 10:
					dh_kapb = par[j];
					break;
				case 11:
					ds_kapb = par[j];
					break;
				case 12:
					dh_kb = par[j];
					break;
				case 13:
					zeta_kb = par[j];
					break;
				case 14:
					mob = par[j];
					break;
				}
				j++;
			}
		}
		if (bays) {
			switch (i) {
			case 0:
				dh_kapc = par[i];
				break;
			case 1:
				ds_kapc = par[i];
				break;
			case 2:
				dh_kc = par[i];
				break;
			case 3:
				zeta_kc = par[i];
				break;
			case 4:
				nv = par[i];
				break;
			case 5:
				dh_kaph = par[i];
				break;
			case 6:
				ds_kaph = par[i];
				break;
			case 7:
				dh_kh = par[i];
				break;
			case 8:
				zeta_kh = par[i];
				break;
			case 9:
				nh = par[i];
				break;
			case 10:
				dh_kapb = par[i];
				break;
			case 11:
				ds_kapb = par[i];
				break;
			case 12:
				dh_kb = par[i];
				break;
			case 13:
				zeta_kb = par[i];
				break;
			case 14:
				mob = par[i];
				break;
			}
		}
	}
	if (!inbounds)
		return 0;

	using namespace boost::numeric::ublas;
	typedef vector<double> uvec_t;
	typedef vector<double>::iterator uvi_t;
	typedef matrix<double> gvvec_t;
	typedef matrix<double>::iterator1 gvi1_t;
	typedef matrix<double>::iterator2 gvi2_t;
	typedef matrix_column<matrix<double> > umc_t;

	//some physical constants
	double kB = 1.3806503e-23;   //Boltzmann's constant
	double Na = 6.0221415e23;    //Avogadro's number
	R = kB*Na;            //universal gas constant
	double pi = 3.14159265;
	double M = 43.989829244e-3;  //molar mass of CO2
	double H = 18.015e-3;        //molar mass of H2O
	double F = 96485.3335;  // Faraday's constant
	//Calculation of Volumes
	Vt = (area*dl);
	Vv = por*Vt;
	Vs = (1 - por)*Vt;
	double Conc_atm = (Patm*Vv) / (R*Tatm);
	//CHANGE DESCRIPTION
	// The parameters in _xpass are, for each row, the CO2
	// concentration in Pa, the temperature in deg.-K, and the time.
	// The state of the sorbent at time zero is assumed to be free of
	// CO2 under a purge gas, at a temperature equal to the first
	// temperature in the data set.

	double TT;
	long NT;
	//for(int j=0; j<xpass.size1(); j++){
	//std::cout << xpass(j,0) << " " << xpass(j,1) << " " << xpass(j,2) << " " << xpass(j,3) << " " << xpass(j,4) << std::endl; }

	TT = *(xpass.rbegin1()).rbegin() + 5.0;
	NT = (long)TT / delt;
	nh_inv = 1.0 / nh;
	nv_inv = 1.0 / nv;
	// Vectors for timestep plus 1 and test in newton solver
	uvec_t C2, H2, X2, A2, B2, P2, Ct(numcv), Ht(numcv), Xt(numcv), At(numcv), Bt(numcv), Pt(numcv);
	// Initialize field variables as 0.0 and Patm for pressure
	for (int d = 0; d < numcv; d++) {
		C1[d] = 0.0;
		H1[d] = 0.0;
		X1[d] = 0.0;
		A1[d] = 0.0;
		B1[d] = 0.0;
		P1[d] = Patm;
	}
	C2 = C1; //Set Variables "2" to state "1"
	H2 = H1;
	X2 = X1;
	A2 = A1;
	B2 = B1;
	P2 = P1;
	Fc1 = 0.0; //Set inputs for time 0 equal to 0.0
	Fh1 = 0.0;
	Fn1 = 0.0;
	Fin1 = 0.0;
	// Vectors and constants to use for the Newton iterations
	uvec_t res(6 * numcv), s(6 * numcv), restest(6 * numcv), diff(6 * numcv),Tres(6*numcv),H2Oper(numcv),CO2per(numcv);
	double nres, nrestest =40, lambda, r0,crit1, crit2;
	long newcount, newcount2; //Counters for max itteration in newton solver
	bool ypasset = false;
	bool cut;
	int cutter = 0;
	long timestep = 0;
	gvi1_t xpass_it1 = xpass.begin1();
	gvi1_t ypass_it = ypass.begin1();
	do {
		//Advance the iterator when a point in _xpass is passed.
		if (xpass_it1 != xpass.end1()) {
			if ((timestep + 1)*delt >= *xpass_it1.rbegin()) {
				++xpass_it1;
				ypasset = true;
			}
		}

		//If we still haven't passed the first point in _xpass, set the
		//temperature to that first point and interpolate between the
		//first pressure in _xpass and the purge state (zero pressure).
		if (xpass_it1 == xpass.begin1()) {
			T = *(xpass_it1.begin() + 3);// + 273.15;
			Fc2 = *xpass_it1.begin()*(timestep + 1)*delt / (*xpass_it1.rbegin());
			Fh2 = (*(xpass_it1.begin() + 1)*(timestep + 1)*delt / (*xpass_it1.rbegin()));
			Fn2 = *(xpass_it1.begin() + 2)*(timestep + 1)*delt / (*xpass_it1.rbegin());
		}
		//If we've passed the end point in _xpass, maintain the
		//temperature and pressure at their values in the last point of
		//_xpass (in other words, do nothing).  But if not, interpolate.
		else if (xpass_it1 != xpass.end1()) {
			T =    (*((xpass_it1 - 1).begin() + 3)*(*xpass_it1.rbegin() - (timestep + 1)*delt) + *(xpass_it1.begin() + 3)*((timestep + 1)*delt - *(xpass_it1 - 1).rbegin())) / (*xpass_it1.rbegin() - *(xpass_it1 - 1).rbegin());
			Fc2 = ((*(xpass_it1 - 1).begin()*(*xpass_it1.rbegin() - (timestep + 1)*delt) + *xpass_it1.begin()*((timestep + 1)*delt - *(xpass_it1 - 1).rbegin())) / (*xpass_it1.rbegin() - *(xpass_it1 - 1).rbegin()));
			Fh2 = ((*((xpass_it1 - 1).begin() + 1)*(*xpass_it1.rbegin() - (timestep + 1)*delt) + *(xpass_it1.begin() + 1)*((timestep + 1)*delt - *(xpass_it1 - 1).rbegin())) / (*xpass_it1.rbegin() - *(xpass_it1 - 1).rbegin()));
			Fn2 = ((*((xpass_it1 - 1).begin() + 2)*(*xpass_it1.rbegin() - (timestep + 1)*delt) + *(xpass_it1.begin() + 2)*((timestep + 1)*delt - *(xpass_it1 - 1).rbegin())) / (*xpass_it1.rbegin() - *(xpass_it1 - 1).rbegin()));
		}
		//if (timestep == 10612)
			//std::cout << "I am here" << std::endl;
		//Caclulate New Concentrations
		Conc = (P1*Vv) / (R*T);
		//Calculate Constants for timestep
		Fin2 = Fc2 + Fh2 + Fn2;
		kapc = exp(ds_kapc / R)*exp(-dh_kapc / (R*T)) / P;
		kapb = exp(ds_kapb / R)*exp(-dh_kapb / (R*T)) / P;
		kc = zeta_kc*T*exp(-dh_kc / (R*T));
		kb = zeta_kb*T*exp(-dh_kb / (R*T));
		kaph = exp(ds_kaph / R)*exp(-dh_kaph / (R*T)) / P;
		kh = zeta_kh*T*exp(-dh_kh / (R*T));
		//Newton Itterations-------------------------------------------------------------------------------------------------------------------------------
		f(X2, A2, B2, C2, H2, P2, res); //Initail residual calculations
		nres = norm_2(res); //Norm of residual
		r0 = nres;
		newcount = 0;
		crit1 = (taur*r0 + taua);
		//While the residual is larger than the tolerance
		while (nres > crit1) {
			jacsolve(X2, A2, B2, C2, H2, P2, -res, s); //Solve jacobian for search direction
			if (norm_2(s) != norm_2(s) || nres != nres)
 				return 0;
			lambda = 1.0; //Set step size factor to 1 for a whole step
			newcount2 = 0;
			do {
				//Add step for all field variables
				int i = 0;
				while (i < numcv) {
					Xt[i] = X2[i] + (lambda*s[i]);
					//if (Xt[i] < 0.0)
					//	Xt[i] = 0.0;
					At[i] = A2[i] + (lambda*s[numcv + i]);
					//if (At[i] < 0.0)
					//	At[i] = 0.0;
					Bt[i] = B2[i] + (lambda*s[(2 * numcv) + i]);
					//if (Bt[i] < 0.0)
					//	Bt[i] = 0.0;
					Ct[i] = C2[i] + (lambda*s[(3 * numcv) + i]);
					//if (Ct[i] < 0.0)
					//	Ct[i] = 0.0;
					Ht[i] = H2[i] + (lambda*s[(4 * numcv) + i]);
					//if (Ht[i] < 0.0)
					//	Ht[i] = 0.0;
					Pt[i] = P2[i] + (lambda*s[(5 * numcv) + i]);
					i++;
				}
				f(Xt, At, Bt, Ct, Ht, Pt, restest); //Evaluate step results
				nrestest = norm_2(restest);
				//for (int z = 0; z < res.size(); z++)
					//Tres[z] = abs(res[z]) - abs(restest[z]);
				if (nrestest != nrestest)
					return 0;
				lambda /= 2; //Cut Step in half
				newcount2++;
				if (newcount2 > 1000) {
					//std::cout << "No Convergence at: " << timestep << "\n" << "Residual: " << nres << std::endl;
					//std::cin.get();
					return 0;
				}
				crit2 = nres*(1.0 - 1e-4/lambda);
			} while (nrestest > crit2);
			//Set guess for second state equal to the test step
			X2 = Xt;
			A2 = At;
			B2 = Bt;
			C2 = Ct;
			H2 = Ht;
			P2 = Pt;
			res = restest;
			nres = nrestest;
			++newcount;
			cutter = 0;
			if (newcount > 1000) {
				//std::cout << "No Convergence at: " << timestep << "\n" << "Residual: " << nres << std::endl;
				//std::cin.get();
				return 0;
			}
		}
		//Newton Itterations-------------------------------------------------------------------------------------------------------------------------------
		//If a point in xpass has been stepped through in the last time
		//iteration, write to ypass.
		if (ypasset) {
			/*for (int k = 0; k < numcv; k++) {
				H2Oper[k] = (H2[k] / Conc[k]) * 100;
				CO2per[k] = (C2[k] / Conc[k]) * 100;
			}*/
			H2Oper[numcv-1] = (H2[numcv-1] / Conc_atm) * 100;
			CO2per[numcv-1] = (C2[numcv-1] / Conc_atm) * 100;
			if (Pdrop)
				*ypass_it.begin() = (P2[0]-Patm)/P;
			else
				*ypass_it.begin() = 0.0;
			if (CO2)
				*(ypass_it.begin() + 1) = CO2per[numcv-1];
			else
				*(ypass_it.begin()+1) = 0.0;
			if (H2O)
				*(ypass_it.begin() + 2) = H2Oper[numcv-1];
			else
				*(ypass_it.begin()+2) = 0.0;
			++ypass_it;
			ypasset = false;
		}
		++timestep; // Increase timestep counter
		
		//std::cout << timestep*delt << std::endl;
		X1 = X2;
		A1 = A2;
		B1 = B2;
		C1 = C2;
		H1 = H2;
		P1 = P2;
		Fc1 = Fc2;
		Fh1 = Fh2;
		Fn1 = Fn2;
		Fin1 = Fin2;
		
	} while (timestep <= NT);

	return 1;
}

void empfuneval::f(boost::numeric::ublas::vector<double> const &X2, boost::numeric::ublas::vector<double> const &A2, boost::numeric::ublas::vector<double> const &B2, boost::numeric::ublas::vector<double> const &C2, boost::numeric::ublas::vector<double> const &H2, boost::numeric::ublas::vector<double> const &P2, boost::numeric::ublas::vector<double> &res) {
	//Caclulate New Concentrations
	Conc = (P1*Vv) / (R*T);
	//Calculate all the residuals for all CV's 
	for (int k = 0; k < numcv; k++) {
		//First obtain Available site constants
		st1 = 1.0 - 2 * X1[k] - B1[k];
		st2 = 1.0 - 2 * X2[k] - B2[k];

		//Calculate formation rates for x, a, b for states 1 and 2
		Rx1 = (kc*((st1*st1*P1[k] * (C1[k] / Conc[k])) - (X1[k] * X1[k] / kapc)));
		Rx2 = (kc*((st2*st2*P1[k] * (C2[k] / Conc[k])) - (X2[k] * X2[k]/ kapc)));
		Ra1 = (kh*((P1[k] * (H1[k] / Conc[k])*(1-A1[k])) - (A1[k] / kaph)));
		Ra2 = (kh*((P1[k] * (H2[k] / Conc[k])*(1-A2[k])) - (A2[k] / kaph)));
		Rb1 = (kb*((st1*A1[k] * P1[k] * (C1[k] / Conc[k])) - (B1[k] * (1-A1[k]) / kapb)));
		Rb2 = (kb*((st2*A2[k] * P1[k] * (C2[k] / Conc[k])) - (B2[k] * (1-A2[k]) / kapb)));

		//Carbamate formation
		res[k] = X1[k] + (delt / 2)*(Rx1 + Rx2) - X2[k];

		//Water absorbtion
		res[numcv+k] = A1[k] + (delt / 2)*(Ra1+Ra2-Rb1-Rb2) - A2[k];

		//bicarbonate formation
		res[(2*numcv)+k] = B1[k] + (delt / 2)*(Rb1+Rb2) - B2[k];
		
		//Concentration of CO2 & H20 plus pressure
		if (k == 0) //Boundary Condition for begining of bed
		{
			//CO2
			res[(3 * numcv) + k] = C1[k] + (delt / 2)*((Fc1-(mob*Vv*((C1[k] * (P1[k] - P1[k + 1])) / (Conc[k] * dl*dl)))) + (Fc2-(mob*Vv*((C2[k] * (P1[k] - P1[k + 1])) / (Conc[k] * dl*dl)))) - (Vs)*(nv*(Rx1 + Rx2 + Rb1 + Rb2))) - C2[k];
			//H20
			res[(4 * numcv) + k] = H1[k] + (delt / 2)*((Fh1-(mob*Vv*((H1[k] * (P1[k] - P1[k + 1])) / (Conc[k] * dl*dl)))) + (Fh2-(mob*Vv*((H2[k] * (P1[k] - P1[k + 1])) / (Conc[k] * dl*dl)))) - (Vs)*(nh*(Ra1 + Ra2))) - H2[k];
			//Pressure
			res[(5 * numcv) + k] = (.5)*((Fin1-(mob*Vv*((P1[k] - P1[k + 1]) / (dl*dl)))) + (Fin2 - mob*Vv*(((P2[k] - P2[k + 1]) / (dl*dl)))) - (Vs)*(nv*(Rx1 + Rx2 + Rb1 + Rb2) + nh*(Ra1 + Ra2)));
		}
		else if (k == numcv - 1) //Boundary Condition for end of bed
		{
			//CO2
			res[(3 * numcv) + k] = C1[k] + (delt / 2)*(((mob*Vv*(((C1[k - 1] * (P1[k - 1] - P1[k])) / (Conc[k - 1] * dl*dl)) - ((C1[k] * (P1[k] - Patm)) / (Conc[k] * dl*dl)))) + (mob*Vv*(((C2[k - 1] * (P2[k - 1] - P2[k])) / (Conc[k - 1] * dl*dl)) - ((C2[k] * (P2[k] - Patm)) / (Conc[k] * dl*dl))))) - (Vs)*(nv*(Rx1 + Rx2 + Rb1 + Rb2))) - C2[k];
			//H20
			res[(4 * numcv) + k] = H1[k] + (delt / 2)*(((mob*Vv*(((H1[k - 1] * (P1[k - 1] - P1[k])) / (Conc[k - 1] * dl*dl)) - ((H1[k] * (P1[k] - Patm)) / (Conc[k] * dl*dl)))) + (mob*Vv*(((H2[k - 1] * (P2[k - 1] - P2[k])) / (Conc[k - 1] * dl*dl)) - ((H2[k] * (P2[k] - Patm)) / (Conc[k] * dl*dl))))) - (Vs)*(nh*(Ra1 + Ra2))) - H2[k];
			//Pressure
			res[(5 * numcv) + k] = (.5)*(((mob*Vv/(dl*dl))*(P1[k - 1] - (2 * P1[k]) + P2[k - 1] - (2 * P2[k]) + (2 * Patm))) - (Vs)*(nv*(Rx1 + Rx2 + Rb1 + Rb2) + nh*(Ra1 + Ra2)));
		}
		else {
			//CO2
			res[(3 * numcv) + k] = C1[k] + (delt / 2)*(((mob*Vv*(((C1[k - 1] * (P1[k - 1] - P1[k])) / (Conc[k - 1] * dl*dl)) - ((C1[k] * (P1[k] - P1[k + 1])) / (Conc[k] * dl*dl)))) + (mob*Vv*(((C2[k - 1] * (P2[k - 1] - P2[k])) / (Conc[k - 1] * dl*dl)) - ((C2[k] * (P2[k] - P2[k + 1])) / (Conc[k] * dl*dl))))) - (Vs)*(nv*(Rx1 + Rx2 + Rb1 + Rb2))) - C2[k];
			//H20
			res[(4 * numcv) + k] = H1[k] + (delt / 2)*(((mob*Vv*(((H1[k - 1] * (P1[k - 1] - P1[k])) / (Conc[k - 1] * dl*dl)) - ((H1[k] * (P1[k] - P1[k + 1])) / (Conc[k] * dl*dl)))) + (mob*Vv*(((H2[k - 1] * (P2[k - 1] - P2[k])) / (Conc[k - 1] * dl*dl)) - ((H2[k] * (P2[k] - P2[k + 1])) / (Conc[k] * dl*dl))))) - (Vs)*(nh*(Ra1 + Ra2))) - H2[k];
			//Pressure
			res[(5 * numcv) + k] = (.5)*(((mob*Vv / (dl*dl))*(P1[k - 1] - (2 * P1[k]) + P1[k + 1] + P2[k - 1] - (2 * P2[k]) + P2[k + 1])) - (Vs)*(nv*(Rx1 + Rx2 + Rb1 + Rb2) + nh*(Ra1 + Ra2)));
		}
	}
}

void empfuneval::jacsolve(boost::numeric::ublas::vector<double> const &X2, boost::numeric::ublas::vector<double> const &A2, boost::numeric::ublas::vector<double> const &B2, boost::numeric::ublas::vector<double> const &C2, boost::numeric::ublas::vector<double> const &H2, boost::numeric::ublas::vector<double> const &P2, boost::numeric::ublas::vector<double> const &bb, boost::numeric::ublas::vector<double> &s) {
	//Chemical Reaction Derivatives:
	int i = -1;
	for (int k = 0; k < numcv; k++) {
		st2 = 1.0 - 2 * X2[k] - B2[k];
		//Carbamate Formation ________________________________________________________________________________________________________________________________________________
		i++; IA[i] = k; JA[i] = k; //Equation 1 Respect to x (ind)
		A[i] = (delt / 2)*(kc*((-4+8*X2[k]+4*B2[k])*P1[k]*(C2[k]/Conc[k]) - (2*X2[k])/kapc)) - 1; //Equation 1 Respect to x

		i++; IA[i] = k; JA[i] = (2 * numcv) + k; //Equation 1 Respect to b (ind)
		A[i] = (delt / 2)*(kc*((-2 + 4 * X2[k] + 2 * B2[k])*P1[k] * (C2[k] / Conc[k]))); //Equation 1 Respect to b
		
		i++; IA[i] = k; JA[i] = (3 * numcv) + k; //Equation 1 Respect to c (ind)
		A[i] = (delt / 2)*(kc*(st2*st2*P1[k]/Conc[k])); //Equation 1 Respect to c
		
		//Water Absorbtion_______________________________________________________________________________________________________________________________________________________
		i++; IA[i] = (1 * numcv) + k; JA[i] = k; //Equation 2 Respect to x (ind)
		A[i] = (delt / 2)*(kb*((2*A2[k]*P1[k]*C2[k]/Conc[k]))); //Equation 2 Respect to x
		
		i++; IA[i] = (1 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 2 Respect to a (ind)
		A[i] = (delt / 2)*(-kh*(P1[k]*H2[k]/Conc[k] + 1/kaph)-kb*((st2*P1[k] * C2[k] / Conc[k]) + (B2[k]/kapb)))-1; //Equation 2 Respect to a
		
		i++; IA[i] = (1 * numcv) + k; JA[i] = (2 * numcv) + k; //Equation 2 Respect to b (ind)
		A[i] = (delt / 2)*(kb*((A2[k]*P1[k]*C2[k]/Conc[k])+(1-A2[k])/kapb)); //Equation 2 Respect to b
		
		i++; IA[i] = (1 * numcv) + k; JA[i] = (3 * numcv) + k; //Equation 2 Respect to c (ind)
		A[i] = (delt / 2)*(-kb*(st2*P1[k]*A2[k]/ Conc[k])); //Equation 2 Respect to c
		
		i++; IA[i] = (1 * numcv) + k; JA[i] = (4 * numcv) + k; //Equation 2 Respect to h (ind)
		A[i] = (delt / 2)*(kh*((P1[k]*(1-A2[k]))/Conc[k])); //Equation 2 Respect to h
		
		
		//Bicabonate Formation_________________________________________________________________________________________________________________________________________________
		i++; IA[i] = (2 * numcv) + k; JA[i] = k; //Equation 3 Respect to x (ind)
		A[i] = (delt / 2)*(kb*((-2 * A2[k] * P1[k] * C2[k] / Conc[k]))); //Equation 3 Respect to x

		i++; IA[i] = (2 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 3 Respect to a (ind)
		A[i] = (delt / 2)*(kb*(st2*P1[k] * C2[k] / Conc[k])+(B2[k]/kapb)); //Equation 3 Respect to a

		i++; IA[i] = (2 * numcv) + k; JA[i] = (2 * numcv) + k; //Equation 3 Respect to b (ind)
		A[i] = (delt / 2)*(-kb*((A2[k] * P1[k] * C2[k] / Conc[k])  + (1-A2[k]) / kapb))-1; //Equation 3 Respect to b

		i++; IA[i] = (2 * numcv) + k; JA[i] = (3 * numcv) + k; //Equation 3 Respect to c (ind)
		A[i] = (delt / 2)*(kb*(st2*P1[k] * A2[k] / Conc[k])); //Equation 3 Respect to c
	}
	//Concentration & Pressure Equation Derivatives
	//Begining of Bed+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	int k = 0;
	st2 = 1.0 - 2 * X2[k] - B2[k];
	//CO2 Concentration ________________________________________________________________________________________________________________________________________________
	i++; IA[i] = (3 * numcv) + k; JA[i] = k; //Equation 4 Respect to x (ind)
	A[i] = (delt / 2)*(-Vs)*(nv*(kc*((-4 + 8 * X2[k] + 4 * B2[k])*P1[k] * (C2[k] / Conc[k]) - (2 * X2[k]) / kapc) + (kb*((-2 * A2[k] * P1[k] * C2[k] / Conc[k]))))); //Equation 4 Respect to x

	i++; IA[i] = (3 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 4 Respect to a (ind)
	A[i] = (delt / 2)*(-Vs)*(nv*kb*(st2*P1[k] * C2[k] / Conc[k]) + (B2[k] / kapb)); //Equation 4 Respect to a

	i++; IA[i] = (3 * numcv) + k; JA[i] = (2 * numcv) + k; //Equation 4 Respect to b (ind)
	A[i] = (delt / 2)*(-Vs)*(nv*(kc*((-2 + 4 * X2[k] + 2 * B2[k])*P1[k] * (C2[k] / Conc[k])) - (kb*((A2[k] * P1[k] * C2[k] / Conc[k]) + (1 - A2[k]) / kapb)))); //Equation 4 Respect to b
    
	i++; IA[i] = (3 * numcv) + k; JA[i] = (3 * numcv) + k; //Equation 4 Respect to c (ind)
	A[i] = (delt / 2)*((mob*Vv*(P2[k + 1] - P2[k]) / (Conc[k] * dl*dl))-(Vs)*(nv*(kc*(st2*st2*P1[k] / Conc[k])) + (kb*(st2*P1[k] * A2[k] / Conc[k])))) - 1; //Equation 4 Respect to c

	i++; IA[i] = (3 * numcv) + k; JA[i] = (5 * numcv) + k; //Equation 4 Respect to P (ind)
	A[i] = (delt / 2)*(mob*Vv*(-(C2[k] / (Conc[k] * dl*dl)))); //Equation 4 Respect to P

	i++; IA[i] = (3 * numcv) + k; JA[i] = (5 * numcv) + k + 1; //Equation 4 Respect to P {i+1} (ind)
	A[i] = (delt / 2)*(mob*Vv*(C2[k]) / (Conc[k] * dl*dl)); //Equation 4 Respect to P {i+1}

	//H20 Concentration_______________________________________________________________________________________________________________________________________________
	i++; IA[i] = (4 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 5 Respect to a (ind)
	A[i] = (delt / 2)*(-Vs)*nh*(-kh*(P1[k] * H2[k] / Conc[k] + 1 / kaph)); //Equation 5 Respect to a
	
	i++; IA[i] = (4 * numcv) + k; JA[i] = (4 * numcv) + k; //Equation 5 Respect to h (ind)
	A[i] = (delt / 2)*((mob*Vv*(P2[k + 1] - P2[k]) / (Conc[k] * dl*dl))-(Vs)*(nh*(kh*((P1[k] * (1 - A2[k])) / Conc[k])))) - 1; //Equation 5 Respect to h
	
	i++; IA[i] = (4 * numcv) + k; JA[i] = (5 * numcv) + k; //Equation 5 Respect to P (ind)
	A[i] = (delt / 2)*(mob*Vv*(-(H2[k] / (Conc[k] * dl*dl)))); //Equation 5 Respect to P

	i++; IA[i] = (4 * numcv) + k; JA[i] = (5 * numcv) + k + 1; //Equation 5 Respect to P {i+1} (ind)
	A[i] = (delt / 2)*(mob*Vv*H2[k] / (Conc[k] * dl*dl)); //Equation 5 Respect to P {i+1}
	
	//Pressure Drop_________________________________________________________________________________________________________________________________________________
	i++; IA[i] = (5 * numcv) + k; JA[i] = k; //Equation 6 Respect to x (ind)
	A[i] = (.5)*(-Vs)*(nv*(kc*((-4 + 8 * X2[k] + 4 * B2[k])*P1[k] * (C2[k] / Conc[k]) - (2 * X2[k]) / kapc) + (kb*(-2 * A2[k] * P1[k] * C2[k] / Conc[k]))));//Equation 6 Respect to x

	i++; IA[i] = (5 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 6 Respect to a (ind)
	A[i] = (.5)*(-Vs)*(nh*(-kh*(P1[k] * H2[k] / Conc[k] + 1 / kaph)) + nv*(kb*(st2*P1[k] * C2[k] / Conc[k]) + (B2[k] / kapb))); //Equation 6 Respect to a

	i++; IA[i] = (5 * numcv) + k; JA[i] = (2 * numcv) + k; //Equation 6 Respect to b (ind)
	A[i] = (.5)*(-Vs)*(nv*(kc*((-2 + 4 * X2[k] + 2 * B2[k])*P1[k] * (C2[k] / Conc[k])) - (kb*((A2[k] * P1[k] * C2[k] / Conc[k]) + (1 - A2[k]) / kapb)))); //Equation 6 Respect to b
	
	i++; IA[i] = (5 * numcv) + k; JA[i] = (3 * numcv) + k; //Equation 6 Respect to c (ind)
	A[i] = (.5)*(-Vs)*(nv*(kc*(st2*st2*P1[k] / Conc[k])) + (kb*(st2*P1[k] * A2[k] / Conc[k]))); //Equation 6 Respect to c

	i++; IA[i] = (5 * numcv) + k; JA[i] = (4 * numcv) + k; //Equation 6 Respect to h (ind)
	A[i] = (.5)*(-Vs)*(nh*(kh*((P1[k] * (1 - A2[k])) / Conc[k]))); //Equation 6 Respect to h
	
	i++; IA[i] = (5 * numcv) + k; JA[i] = (5 * numcv) + k; //Equation 6 Respect to P (ind)
	A[i] = (.5)*((-mob*Vv / (dl*dl))); //Equation 6 Respect to P

	i++; IA[i] = (5 * numcv) + k; JA[i] = (5 * numcv) + k + 1; //Equation 6 Respect to P {i+1} (ind)
	A[i] = (.5)*(mob*Vv / (dl*dl)); //Equation 6 Respect to P {i+1}
	//Middle of Bed+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	for (int k = 1; k < numcv-1; k++) {
		st2 = 1.0 - 2 * X2[k] - B2[k];
		//CO2 Concentration ________________________________________________________________________________________________________________________________________________
		i++; IA[i] = (3 * numcv) + k; JA[i] = k; //Equation 4 Respect to x (ind)
		A[i] = (delt / 2)*(-Vs)*(nv*(kc*((-4 + 8 * X2[k] + 4 * B2[k])*P1[k] * (C2[k] / Conc[k]) - (2 * X2[k]) / kapc) + (kb*((-2 * A2[k] * P1[k] * C2[k] / Conc[k]))))); //Equation 4 Respect to x

		i++; IA[i] = (3 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 4 Respect to a (ind)
		A[i] = (delt / 2)*(-Vs)*(nv*kb*(st2*P1[k] * C2[k] / Conc[k]) + (B2[k] / kapb)); //Equation 4 Respect to a
		
		i++; IA[i] = (3 * numcv) + k; JA[i] = (2 * numcv) + k; //Equation 4 Respect to b (ind)
		A[i] = (delt / 2)*(-Vs)*(nv*(kc*((-2 + 4 * X2[k] + 2 * B2[k])*P1[k] * (C2[k] / Conc[k])) - (kb*((A2[k] * P1[k] * C2[k] / Conc[k]) + (1 - A2[k]) / kapb)))); //Equation 4 Respect to b
		
		i++; IA[i] = (3 * numcv) + k; JA[i] = (3 * numcv) + k; //Equation 4 Respect to c (ind)
		A[i] = (delt / 2)*((mob*Vv*(P2[k + 1] - P2[k]) / (Conc[k] * dl*dl))-(Vs)*(nv*(kc*(st2*st2*P1[k] / Conc[k])) + (kb*(st2*P1[k] * A2[k] / Conc[k])))) - 1; //Equation 4 Respect to c

		i++; IA[i] = (3 * numcv) + k; JA[i] = (3 * numcv) + k - 1; //Equation 4 Respect to c {i-1} (ind)
		A[i] = (delt / 2)*(mob*Vv*(P2[k - 1] - P2[k]) / (Conc[k - 1] * dl*dl)); //Equation 4 Respect to c {i-1}

		i++; IA[i] = (3 * numcv) + k; JA[i] = (5 * numcv) + k; //Equation 4 Respect to P (ind)
		A[i] = (delt / 2)*(mob*Vv*((-C2[k - 1] / (Conc[k - 1] * dl*dl)) - (C2[k] / (Conc[k] * dl*dl)))); //Equation 4 Respect to P

		i++; IA[i] = (3 * numcv) + k; JA[i] = (5 * numcv) + k - 1; //Equation 4 Respect to P {i-1} (ind)
		A[i] = (delt / 2)*(mob*Vv*(C2[k - 1])/(Conc[k - 1]*dl*dl)); //Equation 4 Respect to P {i-1}

		i++; IA[i] = (3 * numcv) + k; JA[i] = (5 * numcv) + k + 1; //Equation 4 Respect to P {i+1} (ind)
		A[i] = (delt / 2)*(mob*Vv*(C2[k]) / (Conc[k] * dl*dl)); //Equation 4 Respect to P {i+1}

		//H20 Concentration_______________________________________________________________________________________________________________________________________________
		i++; IA[i] = (4 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 5 Respect to a (ind)
		A[i] = (delt / 2)*(-Vs)*nh*(-kh*(P1[k] * H2[k] / Conc[k] + 1 / kaph)); //Equation 5 Respect to a
		
		i++; IA[i] = (4 * numcv) + k; JA[i] = (4 * numcv) + k; //Equation 5 Respect to h (ind)
		A[i] = (delt / 2)*((mob*Vv*(P2[k + 1] - P2[k]) / (Conc[k] * dl*dl))- (Vs)*(nh*(kh*((P1[k] * (1 - A2[k])) / Conc[k])))) - 1; //Equation 5 Respect to h

		i++; IA[i] = (4 * numcv) + k; JA[i] = (4 * numcv) + k - 1; //Equation 5 Respect to h {i-1} (ind)
		A[i] = (delt / 2)*(mob*Vv*(P2[k - 1] - P2[k]) / (Conc[k - 1] * dl*dl)); //Equation 5 Respect to h {i-1}

		i++; IA[i] = (4 * numcv) + k; JA[i] = (5 * numcv) + k; //Equation 5 Respect to P (ind)
		A[i] = (delt / 2)*(mob*Vv*((-H2[k - 1] / (Conc[k - 1] * dl*dl)) - (H2[k] / (Conc[k] * dl*dl)))); //Equation 5 Respect to P

		i++; IA[i] = (4 * numcv) + k; JA[i] = (5 * numcv) + k - 1; //Equation 5 Respect to P {i-1} (ind)
		A[i] = (delt / 2)*(mob*Vv*(H2[k - 1]) / (Conc[k - 1] * dl*dl)); //Equation 5 Respect to P {i-1}

		i++; IA[i] = (4 * numcv) + k; JA[i] = (5 * numcv) + k + 1; //Equation 5 Respect to P {i+1} (ind)
		A[i] = (delt / 2)*(mob*Vv*H2[k]/(Conc[k] * dl*dl)); //Equation 5 Respect to P {i+1}
		
		//Pressure Drop_________________________________________________________________________________________________________________________________________________
		i++; IA[i] = (5 * numcv) + k; JA[i] = k; //Equation 6 Respect to x (ind)
		A[i] = (.5)*(-Vs)*(nv*(kc*((-4 + 8 * X2[k] + 4 * B2[k])*P1[k] * (C2[k] / Conc[k]) - (2 * X2[k]) / kapc) + (kb*(-2 * A2[k] * P1[k] * C2[k] / Conc[k])))); //Equation 6 Respect to x

		i++; IA[i] = (5 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 6 Respect to a (ind)
		A[i] = (.5)*(-Vs)*(nh*(-kh*(P1[k] * H2[k] / Conc[k] + 1 / kaph)) + nv*(kb*(st2*P1[k] * C2[k] / Conc[k]) + (B2[k] / kapb))); //Equation 6 Respect to a

		i++; IA[i] = (5 * numcv) + k; JA[i] = (2 * numcv) + k; //Equation 6 Respect to b (ind)
		A[i] = (.5)*(-Vs)*(nv*(kc*((-2 + 4 * X2[k] + 2 * B2[k])*P1[k] * (C2[k] / Conc[k])) - (kb*((A2[k] * P1[k] * C2[k] / Conc[k]) + (1 - A2[k]) / kapb)))); //Equation 6 Respect to b
		
		i++; IA[i] = (5 * numcv) + k; JA[i] = (3 * numcv) + k; //Equation 6 Respect to c (ind)
		A[i] = (.5)*(-Vs)*(nv*(kc*(st2*st2*P1[k] / Conc[k])) + (kb*(st2*P1[k] * A2[k] / Conc[k]))); //Equation 6 Respect to c

		i++; IA[i] = (5 * numcv) + k; JA[i] = (4 * numcv) + k; //Equation 6 Respect to h (ind)
		A[i] = (.5)*(-Vs)*(nh*(kh*((P1[k] * (1 - A2[k])) / Conc[k]))); //Equation 6 Respect to h
		
		i++; IA[i] = (5 * numcv) + k; JA[i] = (5 * numcv) + k; //Equation 6 Respect to P (ind)
		A[i] = (.5)*((-2 * mob*Vv / (dl*dl))); //Equation 6 Respect to P

		i++; IA[i] = (5 * numcv) + k; JA[i] = (5 * numcv) + k - 1; //Equation 6 Respect to P {i-1} (ind)
		A[i] = (.5)*(mob*Vv/(dl*dl)); //Equation 6 Respect to P {i-1}

		i++; IA[i] = (5 * numcv) + k; JA[i] = (5 * numcv) + k + 1; //Equation 6 Respect to P {i+1} (ind)
		A[i] = (.5)*(mob*Vv/(dl*dl)); //Equation 6 Respect to P {i+1}
	}
	//End of Bed+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	k = numcv - 1;
	st2 = 1.0 - 2 * X2[k] - B2[k];
	//CO2 Concentration ________________________________________________________________________________________________________________________________________________
	i++; IA[i] = (3 * numcv) + k; JA[i] = k; //Equation 4 Respect to x (ind)
	A[i] = (delt / 2)*(-Vs)*(nv*(kc*((-4 + 8 * X2[k] + 4 * B2[k])*P1[k] * (C2[k] / Conc[k]) - (2 * X2[k]) / kapc) + (kb*((-2 * A2[k] * P1[k] * C2[k] / Conc[k]))))); //Equation 4 Respect to x

	i++; IA[i] = (3 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 4 Respect to a (ind)
	A[i] = (delt / 2)*(-Vs)*(nv*kb*(st2*P1[k] * C2[k] / Conc[k]) + (B2[k] / kapb)); //Equation 4 Respect to a

	i++; IA[i] = (3 * numcv) + k; JA[i] = (2 * numcv) + k; //Equation 4 Respect to b (ind)
	A[i] = (delt / 2)*(-Vs)*(nv*(kc*((-2 + 4 * X2[k] + 2 * B2[k])*P1[k] * (C2[k] / Conc[k])) - (kb*((A2[k] * P1[k] * C2[k] / Conc[k]) + (1 - A2[k]) / kapb)))); //Equation 4 Respect to b
	
	i++; IA[i] = (3 * numcv) + k; JA[i] = (3 * numcv) + k; //Equation 4 Respect to c (ind)
	A[i] = (delt / 2)*((mob*Vv*(Patm - P2[k]) / (Conc[k] * dl*dl))-(Vs)*(nv*(kc*(st2*st2*P1[k] / Conc[k])) + (kb*(st2*P1[k] * A2[k] / Conc[k])))) - 1; //Equation 4 Respect to c

	i++; IA[i] = (3 * numcv) + k; JA[i] = (3 * numcv) + k - 1; //Equation 4 Respect to c {i-1} (ind)
	A[i] = (delt / 2)*(mob*Vv*(P2[k - 1] - P2[k]) / (Conc[k - 1] * dl*dl)); //Equation 4 Respect to c {i-1}

	i++; IA[i] = (3 * numcv) + k; JA[i] = (5 * numcv) + k; //Equation 4 Respect to P (ind)
	A[i] = (delt / 2)*(mob*Vv*((-C2[k - 1] / (Conc[k - 1] * dl*dl)) - (C2[k] / (Conc[k] * dl*dl)))); //Equation 4 Respect to P

	i++; IA[i] = (3 * numcv) + k; JA[i] = (5 * numcv) + k - 1; //Equation 4 Respect to P {i-1} (ind)
	A[i] = (delt / 2)*(mob*Vv*(C2[k - 1]) / (Conc[k - 1] * dl*dl)); //Equation 4 Respect to P {i-1}

	//H20 Concentration_______________________________________________________________________________________________________________________________________________
	i++; IA[i] = (4 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 5 Respect to a (ind)
	A[i] = (delt / 2)*(-Vs)*nh*(-kh*(P1[k] * H2[k] / Conc[k] + 1 / kaph)); //Equation 5 Respect to a
	
	i++; IA[i] = (4 * numcv) + k; JA[i] = (4 * numcv) + k; //Equation 5 Respect to h (ind)
	A[i] = (delt / 2)*((mob*Vv*(Patm - P2[k]) / (Conc[k] * dl*dl))- (Vs)*(nh*(kh*((P1[k] * (1 - A2[k])) / Conc[k])))) - 1; //Equation 5 Respect to h

	i++; IA[i] = (4 * numcv) + k; JA[i] = (4 * numcv) + k - 1; //Equation 5 Respect to h {i-1} (ind)
	A[i] = (delt / 2)*(mob*Vv*(P2[k - 1] - P2[k]) / (Conc[k - 1] * dl*dl)); //Equation 5 Respect to h {i-1}

	i++; IA[i] = (4 * numcv) + k; JA[i] = (5 * numcv) + k; //Equation 5 Respect to P (ind)
	A[i] = (delt / 2)*(mob*Vv*((-H2[k - 1] / (Conc[k - 1] * dl*dl)) - (H2[k] / (Conc[k] * dl*dl)))); //Equation 5 Respect to P

	i++; IA[i] = (4 * numcv) + k; JA[i] = (5 * numcv) + k - 1; //Equation 5 Respect to P {i-1} (ind)
	A[i] = (delt / 2)*(mob*Vv*(H2[k - 1]) / (Conc[k - 1] * dl*dl)); //Equation 5 Respect to P {i-1}
	
	//Pressure Drop_________________________________________________________________________________________________________________________________________________
	i++; IA[i] = (5 * numcv) + k; JA[i] = k; //Equation 6 Respect to x (ind)
	A[i] = (.5)*(-Vs)*(nv*(kc*((-4 + 8 * X2[k] + 4 * B2[k])*P1[k] * (C2[k] / Conc[k]) - (2 * X2[k]) / kapc) + (kb*(-2 * A2[k] * P1[k] * C2[k] / Conc[k])))); //Equation 6 Respect to x

	i++; IA[i] = (5 * numcv) + k; JA[i] = (1 * numcv) + k; //Equation 6 Respect to a (ind)
	A[i] = (.5)*(-Vs)*(nh*(-kh*(P1[k] * H2[k] / Conc[k] + 1 / kaph)) + nv*(kb*(st2*P1[k] * C2[k] / Conc[k]) + (B2[k] / kapb))); //Equation 6 Respect to a

	i++; IA[i] = (5 * numcv) + k; JA[i] = (2 * numcv) + k; //Equation 6 Respect to b (ind)
	A[i] = (.5)*(-Vs)*(nv*(kc*((-2 + 4 * X2[k] + 2 * B2[k])*P1[k] * (C2[k] / Conc[k])) - (kb*((A2[k] * P1[k] * C2[k] / Conc[k]) + (1 - A2[k]) / kapb)))); //Equation 4 Respect to b

	i++; IA[i] = (5 * numcv) + k; JA[i] = (3 * numcv) + k; //Equation 6 Respect to c (ind)
	A[i] = (.5)*(-Vs)*(nv*(kc*(st2*st2*P1[k] / Conc[k])) + (kb*(st2*P1[k] * A2[k] / Conc[k]))); //Equation 6 Respect to c

	i++; IA[i] = (5 * numcv) + k; JA[i] = (4 * numcv) + k; //Equation 6 Respect to h (ind)
	A[i] = (.5)*(-Vs)*(nh*(kh*((P1[k] * (1 - A2[k])) / Conc[k]))); //Equation 6 Respect to h
	
	i++; IA[i] = (5 * numcv) + k; JA[i] = (5 * numcv) + k; //Equation 6 Respect to P (ind)
	A[i] = (.5)*((-2 * mob*Vv / (dl*dl))); //Equation 6 Respect to P

	i++; IA[i] = (5 * numcv) + k; JA[i] = (5 * numcv) + k - 1; //Equation 6 Respect to P {i-1} (ind)
	A[i] = (.5)*(mob*Vv / (dl*dl)); //Equation 6 Respect to P {i-1}
	i++;
	//Pass to GMRES------------------------------------------------------------------------------------------------------------------------------------------------------------
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(i);
	for (int m = 0; m < i; m++)
	{
		tripletList.push_back(T(IA[m], JA[m], A[m])); 
	}
	Eigen::SparseMatrix<double> Jac(6 * numcv, 6 * numcv);
	Eigen::VectorXd rhs(6 * numcv);
	Eigen::VectorXd x_estimate(6 * numcv);
	Jac.setFromTriplets(tripletList.begin(), tripletList.end());
	for (int i = 0; i < (6*numcv); i++) {
		rhs[i] = bb[i];}
	Eigen::GMRES<Eigen::SparseMatrix<double> >solver(Jac);
	solver.setTolerance(1e-8);
	x_estimate = solver.solve(rhs);
	for (int i = 0; i < (6 * numcv); i++) {
		s[i] = x_estimate[i];
	}
}
