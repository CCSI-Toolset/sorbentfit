/*
Variables from configset are named incorectly. Output supplies correct names, but for coding variables use following:
Dry: use variables with "1"
Wat: use variables with "h"
Humid: use variables with "2"

xb[0]=x
xb[1]=b
xb[2]=a

 */

#include "funeval_base.hpp" 
#include "boost/algorithm/string.hpp"
#include <boost/lexical_cast.hpp>
#include "psopoint_ublas_2.hpp"

class empfuneval : public funeval_base {

public:

	empfuneval(std::string const &confile, std::string const &datfiles) { configset(confile); datfile(datfiles); counter = 0; };
	empfuneval(std::string const &confile) { configset(confile); datload = false; counter = 0; };
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
	double delttrue, delt, taur, taua;
	double P, nv, nv_l, nv_h;
	double rho, wtpersite;
	double dh_kap1, dh_kap1_l, dh_kap1_h, ds_kap1, ds_kap1_l, ds_kap1_h, dh_k1, dh_k1_l, dh_k1_h, zeta_k1, lzeta_k1_l, lzeta_k1_h;
	double dh_kaph, dh_kaph_l, dh_kaph_h, ds_kaph, ds_kaph_l, ds_kaph_h, dh_kh, dh_kh_l, dh_kh_h, zeta_kh, lzeta_kh_l, lzeta_kh_h;
	double dh_kap2, dh_kap2_l, dh_kap2_h, ds_kap2, ds_kap2_l, ds_kap2_h, dh_k2, dh_k2_l, dh_k2_h, zeta_k2, lzeta_k2_l, lzeta_k2_h;

	// Equilibrium constants and rate constants
	double kap1, kap2, k1, k2, kaph, kh;

	// Interaction energy
	double gam;

	// The temperature and pressures
	double T, p1, p2;

	// Gas constant
	double R;

	// Computational variables
	double nv_inv, kap1_prime, st, w, stn, wn;
	boost::numeric::ublas::vector<double> up, diag, down, b;

	// Carbamate and bicarbonate fractions for the previous time step.
	boost::numeric::ublas::vector<double> xbn;

	// Jacobian
	boost::numeric::ublas::matrix<double> jac;

	void f(boost::numeric::ublas::vector<double> const &xb, boost::numeric::ublas::vector<double> &res);

	void jacsolve(boost::numeric::ublas::vector<double> const &xb, boost::numeric::ublas::vector<double> const &bb, boost::numeric::ublas::vector<double> &s);


};

void empfuneval::configset(std::string const &confile) {
	override.resize(13);
	lowbounds.resize(13);
	highbounds.resize(13);
	params.resize(13);
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
					delttrue = boost::lexical_cast<double>(str[0]);
					taur = boost::lexical_cast<double>(str[1]);
					taua = boost::lexical_cast<double>(str[2]);
				}
				else error = i + 1;
			}
			else if (i == 1) {
				if (str.size() == 2) {
					P = boost::lexical_cast<double>(str[0]);
					rho = boost::lexical_cast<double>(str[1]);
				}
				else { error = i + 1; }
			}
			else if (i > 1 && i < 15) {
				if (str.size() == 1) {
					params[i - 2] = boost::lexical_cast<double>(str[0]);
					override[i - 2] = true;
				}
				else if (str.size() == 2) {
					lowbounds[i - 2] = boost::lexical_cast<double>(str[0]);
					highbounds[i - 2] = boost::lexical_cast<double>(str[1]);
				}
				else { error = i + 1; }
			}
			else if (i == 15) {
				if (str.size() == 1)
					nodesize = boost::lexical_cast<int>(str[0]);
				else { error = i + 1; }
			}
			else if (i == 16) {
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
	for (int i = 0; i != 13; i++) {
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
				paramorder += "dh_kapb";
				break;
			case 10:
				paramorder += "ds_kapb";
				break;
			case 11:
				paramorder += "dh_kb";
				break;
			case 12:
				paramorder += "zeta_kb";
				break;
			}
		}
		else {
			switch (i) {
			case 0:
				dh_kap1 = params[i];
				break;
			case 1:
				ds_kap1 = params[i];
				break;
			case 2:
				dh_k1 = params[i];
				break;
			case 3:
				zeta_k1 = params[i];
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
				zeta_kh = params[i];
				break;
			case 9:
				dh_kap2 = params[i];
				break;
			case 10:
				ds_kap2 = params[i];
				break;
			case 11:
				dh_k2 = params[i];
				break;
			case 12:
				zeta_k2 = params[i];
				break;
			}
		}
	}

	jac.resize(3, 3);
	xbn.resize(3);
	up.resize(3);
	diag.resize(3);
	down.resize(3);
	b.resize(3);
	for (boost::numeric::ublas::vector<double>::iterator up_it = up.begin(), diag_it = diag.begin(), down_it = down.begin(), b_it = b.begin();
	up_it != up.end();
		++up_it, ++diag_it, ++down_it, ++b_it) {

		*up_it = 0.0;
		*diag_it = 0.0;
		*down_it = 0.0;
		*b_it = 0.0;
	}

	confload = true;
	bays = false;
}

bool empfuneval::function(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::matrix<double> &ypass) {
	int j = 0;
	bool inbounds = true;
	for (int i = 0; i != 13; i++) {
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
			}

			switch (i) {
			case 0:
				dh_kap1 = par[j];
				break;
			case 1:
				ds_kap1 = par[j];
				break;
			case 2:
				dh_k1 = par[j];
				break;
			case 3:
				zeta_k1 = par[j];
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
				dh_kap2 = par[j];
				break;
			case 10:
				ds_kap2 = par[j];
				break;
			case 11:
				dh_k2 = par[j];
				break;
			case 12:
				zeta_k2 = par[j];
				break;
			}
			j++;
		}
		if (bays) {
			switch (i) {
			case 0:
				dh_kap1 = par[i];
				break;
			case 1:
				ds_kap1 = par[i];
				break;
			case 2:
				dh_k1 = par[i];
				break;
			case 3:
				zeta_k1 = par[i];
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
				dh_kap2 = par[i];
				break;
			case 10:
				ds_kap2 = par[i];
				break;
			case 11:
				dh_k2 = par[i];
				break;
			case 12:
				zeta_k2 = par[i];
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

	// The parameters in _xpass are, for each row, the CO2
	// concentration in Pa, the temperature in deg.-K, and the time.
	// The state of the sorbent at time zero is assumed to be free of
	// CO2 under a purge gas, at a temperature equal to the first
	// temperature in the data set.

	double TT;
	long NT;

	TT = *(xpass.rbegin1()).rbegin() + 5.0;
	delt = delttrue;
	NT = (long)TT / delt;
	wtpersite = rho / nv;
	nv_inv = 1.0 / nv;

	// Carbamate and bicarbonate fractions -- they go in vectors, with
	// the first element of the vector the carbamate fraction.
	uvec_t xb(3);

	// Initialize
	xb[0] = 0.0;
	xb[1] = 0.0;
	xb[2] = 0.0;

	xbn = xb;

	// Vectors and constants to use for the Newton iterations
	uvec_t res(3), s(3), restest(3), xbtest(3);
	double nres, nrestest, lambda, r0;
	long newcount, newcount2;
	bool Tlogic = true;
	bool ypasset = false;
	long timestep = 0;
	gvi1_t xpass_it1 = xpass.begin1();
	gvi1_t ypass_it = ypass.begin1();

	do {
		if (timestep < 400 && Tlogic == true) {
			delt = .01;
		}
		else if (timestep == 400 && Tlogic == true) {
			timestep = (4 / delttrue);
			delt = delttrue;
			Tlogic = false;
		}

		/*  Set the parameters for this time step.  */


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

			T = *(xpass_it1.begin() + 2);// + 273.15;
			p1 = *xpass_it1.begin()*(timestep + 1)*delt / (*xpass_it1.rbegin());
			p2 = *(xpass_it1.begin() + 1)*(timestep + 1)*delt / (*xpass_it1.rbegin());

			//std ::cout << T << " " << p1 << " " << p2 << std::endl;
		}
		//If we've passed the end point in _xpass, maintain the
		//temperature and pressure at their values in the last point of
		//_xpass (in other words, do nothing).  But if not, interpolate.
		else if (xpass_it1 != xpass.end1()) {

			T = (*((xpass_it1 - 1).begin() + 2)*(*xpass_it1.rbegin() - (timestep + 1)*delt) + *(xpass_it1.begin() + 2)*((timestep + 1)*delt - *(xpass_it1 - 1).rbegin())) / (*xpass_it1.rbegin() - *(xpass_it1 - 1).rbegin());

			p1 = (*(xpass_it1 - 1).begin()*(*xpass_it1.rbegin() - (timestep + 1)*delt) + *xpass_it1.begin()*((timestep + 1)*delt - *(xpass_it1 - 1).rbegin())) / (*xpass_it1.rbegin() - *(xpass_it1 - 1).rbegin());

			p2 = (*((xpass_it1 - 1).begin() + 1)*(*xpass_it1.rbegin() - (timestep + 1)*delt) + *(xpass_it1.begin() + 1)*((timestep + 1)*delt - *(xpass_it1 - 1).rbegin())) / (*xpass_it1.rbegin() - *(xpass_it1 - 1).rbegin());
			//std ::cout << T << " " << p1 << " " << p2 << std::endl;
		}

		kap1 = exp(ds_kap1 / R)*exp(-dh_kap1 / (R*T)) / P;
		kap2 = exp(ds_kap2 / R)*exp(-dh_kap2 / (R*T)) / P;
		k1 = zeta_k1*T*exp(-dh_k1 / (R*T));
		k2 = zeta_k2*T*exp(-dh_k2 / (R*T));
		kaph = exp(ds_kaph / R)*exp(-dh_kaph / (R*T)) / P;
		kh = zeta_kh*T*exp(-dh_kh / (R*T));

		f(xb, res);
		nres = norm_2(res);
		r0 = nres;
		newcount = 0;

		while (nres > (taur*r0 + taua)) {

			jacsolve(xb, -res, s);
			if (norm_2(s) != norm_2(s) || nres != nres)
				return 0;

			xbtest = xb;
			lambda = 1.0;
			newcount2 = 0;

			do {
				xbtest = xb + lambda*s;
				f(xbtest, restest);
				nrestest = norm_2(restest);
				if (nrestest != nrestest)
					return 0;
				lambda /= 2;
				newcount2++;
				if (newcount2 > 10000) {
					//std::cout<< "No Convergence at: " << timestep << std::endl;
					return 0;
				}

			} while (nrestest > nres*(1.0 - 1e-4*lambda));
			xb = xbtest;
			res = restest;
			nres = nrestest;
			++newcount;
			if (newcount > 1000) {
				//std::cout<< "No Convergence at: " << timestep << std::endl;
				return 0;
			}
		}

		//If a point in xpass has been stepped through in the last time
		//iteration, write to ypass.
		if (ypasset) {

			*ypass_it.begin() = (((timestep + 1)*delt - *(xpass_it1 - 1).rbegin())*(xbn[0] * M*nv + xbn[1] * (M + H) + xbn[2] * H) + (*(xpass_it1 - 1).rbegin() - timestep*delt)*(xb[0] * M*nv + xb[1] * (M + H) + xb[2] * H)) / (wtpersite*nv*delt);

			++ypass_it;
			ypasset = false;
		}

		xbn = xb;
		++timestep;

	} while (timestep <= NT);

	return 1;
}

void empfuneval::f(boost::numeric::ublas::vector<double> const &xb, boost::numeric::ublas::vector<double> &res) {

	st = 1.0 - 2 * xb[0] - nv_inv*xb[1];
	w = xb[0] + nv_inv*xb[1];

	stn = 1.0 - 2 * xbn[0] - nv_inv*xbn[1];
	wn = xbn[0] + nv_inv*xbn[1];

	//res[0] = xbn[0] + k1*delt*(st*st*p1 - xb[0]*w/kap1) - xb[0];
	res[0] = xbn[0] + 0.5*k1*delt*(st*st*p1 - xb[0] * w / kap1 + stn*stn*p1 - xbn[0] * wn / kap1) - xb[0];

	//res[1] = xbn[1] + k2*delt*(st*xb[2]*p1 - w*xb[1]/kap2) - xb[1];
	res[1] = xbn[1] + 0.5*k2*delt*(st*xb[2] * p1 - w*xb[1] / kap2 + stn*xbn[2] * p1 - xbn[1] * wn / kap2) - xb[1];

	//res[2] = xbn[2] + kh*delt*(p2 - xb[2]/kaph) - xb[2];
	res[2] = xbn[2] + 0.5*delt*(kh*(2 * p2 - (xb[2] + xbn[2]) / kaph) - k2*(st*xb[2] * p1 - w*xb[1] / kap2 + stn*xbn[2] * p1 - xbn[1] * wn / kap2)) - xb[2];

}

void empfuneval::jacsolve(boost::numeric::ublas::vector<double> const &xb, boost::numeric::ublas::vector<double> const &bb, boost::numeric::ublas::vector<double> &s) {

	b = bb;

	*diag.begin() = 0.5*k1*delt*(-4 * st*p1 - (w + xb[0]) / kap1) - 1.0;

	*(diag.begin() + 1) = 0.5*k2*delt*(-nv_inv*xb[2] * p1 - (w + nv_inv*xb[1]) / kap2) - 1.0;

	*(diag.begin() + 2) = -0.5*delt*(kh / kaph + k2*st*p1) - 1.0;

	*up.begin() = -0.5*k1*delt*(2 * nv_inv*st*p1 + xb[0] * nv_inv / kap1);

	*(up.begin() + 1) = 0.5*k2*delt*st*p1;

	*down.begin() = 0.5*delt*k2*(2.0*xb[2] * p1 + xb[1] / kap2);

	*(down.begin() + 1) = -0.5*k2*delt*(2 * xb[2] * p1 + xb[1] / kap2);

	*(down.begin() + 2) = 0.5*delt*k2*(xb[2] * p1*nv_inv + (w + xb[1] * nv_inv) / kap2);

	typedef boost::numeric::ublas::vector<double>::iterator vi_t;
	typedef boost::numeric::ublas::vector<double>::reverse_iterator vri_t;

	// eliminate the bottom triangle
	*up.begin() /= *diag.begin();
	*b.begin() /= *diag.begin();
	for (vi_t iup = up.begin() + 1, idiag = diag.begin() + 1, idown = down.begin() + 1, ib = b.begin() + 1;
	iup != up.end() - 1;
		++iup, ++idiag, ++idown, ++ib) {
		*iup /= *idiag - *(iup - 1)*(*idown);
		*ib = (*ib - *(ib - 1)*(*idown)) / (*idiag - *(iup - 1)*(*idown));
	}
	*b.rbegin() = (*b.rbegin() - (*b.begin())*(*down.begin()) - *(b.rbegin() + 1)*(*down.rbegin() - (*up.begin())*(*down.begin()))) / (*diag.rbegin() - *(up.rbegin() + 1)*(*down.rbegin()));

	// back substitute
	vri_t rix = s.rbegin();
	*rix = *b.rbegin();
	for (vri_t rib = b.rbegin() + 1, riup = up.rbegin() + 1;
	rib != b.rend();
		++riup, ++rib) {
		*(rix + 1) = *rib - *riup*(*rix);
		++rix;
	}
}
