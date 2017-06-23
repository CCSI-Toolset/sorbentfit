
// Copyright 2011 David S. Mebane.  Use, modification and distribution
// are subject to the GNU General Public License, Version 3.  (A copy
// is available at <http://www.gnu.org/licenses/gpl.html>.)

// Function evaluation parent class for use with particle swarm
// optimizer.  Public member _objeval returns an evaluation of the
// objective function for the parameter space location _par.  The
// constructor binds the object to the configuration file for the
// model solver (_confile) and to a file containing the names of the
// files which contain the experimental data (_files).  The base class
// is written for the CCSI PEI-based sorbent.

// The parameters are given in the vector _par, with the vector
// _logten evaluating to 'true' for any parameter that is base-10
// logarithmically defined.  The file named in the string _datfiles
// contains a list of filenames, each of which contains a data set.
// Each set pertains to a single call of the associated model function
// _function.  Within each set, the data should be organized with
// independent variables (including any state variables) followed by a
// single dependent variable.  The file _files should begin with two
// integers, followed by the file names.  The numbers are the total
// number of data files and the total number of independent variables
// (assumed the same for each file).



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


class funeval_base {

public:

  funeval_base(std::string const &confile, std::string const &datfiles) {configset(confile); datfile(datfiles); counter = 0;};
  funeval_base(std::string const &confile) {configset(confile); datload = false; counter = 0;};
  funeval_base() {confload = false; datload = false; counter = 0;};

  bool eval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::matrix<double> &ypass);
  boost::numeric::ublas::vector<double> objeval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, int const &mode);
  void printeval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, time_t const &ti);
  double printevalwithPE(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, time_t const &ti, int const &itnum, int const &fnum);
  virtual void configset(std::string const &confile);
  void datfile(std::string const &datfiles);
  long getcounter() {return counter;};
  void setcounter(long count) {counter = count;};
  double baytest;
  bool converge;
protected:

  virtual bool function(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::matrix<double> &ypass);
  boost::numeric::ublas::matrix<double> xval;
  boost::numeric::ublas::matrix<double> yval;
  boost::numeric::ublas::matrix<double> ycalc;
  boost::numeric::ublas::vector<size_t> setnums;
  bool datload;
  bool printme;
  bool confload;
  int xpoints, ypoints;
  std::string ofilename;
  long counter;
  time_t ti;
  long datnum;


};

bool funeval_base::eval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::matrix<double> &ypass) {

  if (!confload) {
    std::cout << "\n\nno config file\n\n";
    return 0;
  }

  using namespace boost::numeric::ublas;
  typedef vector<double> uvec;
  uvec newpar(par);

  for (size_t i=0; i < newpar.size(); ++i) {

    if (logten[i])
      newpar[i] = pow(10.0,newpar[i]);
  }
  return function(newpar, xpass, ypass);

}

boost::numeric::ublas::vector<double> funeval_base::objeval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, int const &mode) {
	boost::numeric::ublas::vector<double> retvar(ypoints);
	for (int i = 0; i < retvar.size(); i++)
		retvar[i] = 0.0;
  if (!datload) {
    std::cout << "\n\nno data file\n\n";
    return retvar;
  }
  if (!confload) {
    std::cout << "\n\nno config file\n\n";
    return retvar;
  }
  converge = true;
  using namespace boost::numeric::ublas;
  typedef vector<double> uvec;
  typedef matrix<double> umat;
  typedef vector<size_t> uvsiz;
  uvec newpar(par);

  for (size_t i=0; i < newpar.size(); ++i) {

    if (logten[i])
      newpar[i] = pow(10.0,newpar[i]);
  }

  umat xproj;
  umat yproj,yvalproj;
  boost::numeric::ublas::matrix<double> diffm;
  boost::numeric::ublas::vector<double> diffv;
  boost::numeric::ublas::vector<double> temp;
  boost::numeric::ublas::vector<double> scale_factor;
  bool success;

  size_t total=0;
  for (uvsiz::iterator setnums_it=setnums.begin();
       setnums_it != setnums.end();
       ++setnums_it) {

    xproj.resize(*setnums_it, xval.size2());
    noalias(xproj) = project(xval, range(total,total + *setnums_it), range(0,xval.size2()));
    yproj.resize(*setnums_it, yval.size2());
	yvalproj.resize(*setnums_it, yval.size2());
	noalias(yvalproj) = project(yval, range(total, total + *setnums_it), range(0, yval.size2()));

    success = function(newpar, xproj, yproj);
    if (!success) {
		//for (umat::iterator1 yproj_it1 = yproj.begin1(); yproj_it1 != yproj.end1(); ++yproj_it1)
			//for (umat::iterator2 yproj_it2 = yproj_it1.begin(); yproj_it2 != yproj_it1.end(); ++yproj_it2)
				//*yproj_it2 = 0.0;
		converge = false;
    }

	if (converge) {
		diffm = yvalproj - yproj;
		scale_factor.resize(diffm.size2());
		//Make scale values
		//scale_factor[0] = 1; 
		//scale_factor[1] = 1;
		//scale_factor[2] = 1; 
		//Obtain Error
		diffv.resize(yval.size1());
		for (int i = 0; i < yval.size2(); i++) {
			diffv = column(diffm, i);
			if (mode == 0)
				retvar[i] = retvar[i] + norm_2(diffv);
			if (mode == 1)
				retvar[i] = retvar[i] + (norm_2(diffv));// / (scale_factor[i] + 1e-16));
		}
	}
    noalias(project(ycalc, range(total,total + *setnums_it),range(0,yval.size2()))) = yproj;
    total += *setnums_it;
  }
  if (!converge) {
	  for (int i = 0; i < yval.size2(); i++)
		  retvar[i] = 1e12;
  }
  return retvar;

}

void funeval_base::printeval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, time_t const &tii) {

  if (!datload) {
    std::cout << "\n\nno data file\n\n";
    return;
  }
  if (!confload) {
    std::cout << "\n\nno config file\n\n";
    return;
  }

  ti = tii;

  using namespace boost::numeric::ublas;
  typedef vector<double> uvec;
  typedef matrix<double> umat;
  typedef vector<size_t> uvsiz;
  uvec newpar(par);

  for (size_t i=0; i < newpar.size(); ++i) {

    if (logten[i])
      newpar[i] = pow(10.0,newpar[i]);
  }

  umat xproj;
  umat::iterator1 xproj_it1;
  umat::iterator2 xproj_it2;
  umat yproj;
  umat yvalproj;
  bool success;

  std::ofstream outfile;
  std::stringstream ofilestream;

  size_t total=0;
  datnum=0;
  for (uvsiz::iterator setnums_it=setnums.begin();
       setnums_it != setnums.end();
       ++setnums_it) {

    xproj.resize(*setnums_it, xval.size2());
    noalias(xproj) = project(xval, range(total,total + *setnums_it), range(0,xval.size2()));
    yvalproj.resize(*setnums_it,yval.size2());
    noalias(yvalproj) = project(yval, range(total, total + *setnums_it),range(0,yval.size2()));
    yproj.resize(*setnums_it,yval.size2());

    success = function(newpar, xproj, yproj);
    if (!success) {
      for (umat::iterator1 yproj_it1=yproj.begin1(); yproj_it1 != yproj.end1();++yproj_it1)
		  for (umat::iterator2 yproj_it2 = yproj_it1.begin(); yproj_it2 != yproj_it1.end(); ++yproj_it2)
			  *yproj_it2 = 0.0;
        //std::cout << "Paramset did not converge" << std::endl;
    }

    total += *setnums_it;
    //std::cout << "Writing Results" << std::endl;
    ofilestream.str("");
    ofilestream << "results/data" << ti << "_" << datnum;
    outfile.open((ofilestream.str()+".txt").c_str(),std::ios::out);
    xproj_it1=xproj.begin1();
    for (umat::iterator1 yproj_it1=yproj.begin1(), yvalproj_it1=yvalproj.begin1(); yproj_it1 != yproj.end1(); ++yproj_it1, ++yvalproj_it1, ++xproj_it1) {
		//Write Xpass
		for (xproj_it2=xproj_it1.begin(); xproj_it2 != xproj_it1.end(); ++xproj_it2) {
			outfile << *xproj_it2 << " ";
		}
		//Write Yval
		for (umat::iterator2 yvalproj_it2 = yvalproj_it1.begin(); yvalproj_it2 != yvalproj_it1.end(); ++yvalproj_it2) {
			outfile << *yvalproj_it2 << " ";
		}
		//Write Yproject
		for (umat::iterator2 yproj_it2 = yproj_it1.begin(); yproj_it2 != yproj_it1.end(); ++yproj_it2) {
			outfile << *yproj_it2 << " ";
		}
      outfile << "\n";
    }
    outfile.close();
    ++datnum;
  }
  ++counter;
}

void funeval_base::datfile(std::string const &files) {

  using namespace boost::numeric::ublas;

  std::ifstream filestream;
  filestream.open(files.c_str(), std::ios::in);

  size_t numfiles;
  filestream >> numfiles;
  filestream >> xpoints >> ypoints;
  setnums.resize(numfiles);

  vector<double> xval_temp(xpoints);
  vector<double> yval_temp(ypoints);
  std::string exdata;
  size_t n1, n2, n3, filen;

  n1=0;
  filen=0;
  while (filestream >> exdata) {

    std::ifstream datstream;
    datstream.open(exdata.c_str(), std::ios::in);
	if (!datstream.is_open()) {
		std::cout << "error" << std::endl;
	}
    datstream >> setnums[filen];

    xval.resize(xval.size1()+setnums[filen],xpoints,true);
    yval.resize(yval.size1()+setnums[filen],ypoints,true);

	//Read Xval values to xval temp
    n2=0;
    while (datstream >> xval_temp[0]) {
      //std::cout << xval_temp[0] << " " << n2;
      ++n2;
      while (n2 < xpoints) {
		  datstream >> xval_temp[n2];
		  //std::cout << " "  << xval_temp[n2] << " " << n2;
		  ++n2;
      }
	  //Write xvaltemp to xval matrix
      for(int m = 0; m < xval_temp.size(); m++){
		  xval(n1,m) = xval_temp[m];
		  //std::cout << xval(n1,m) << " ";
      }
	  //Read yval values to yval temp
	  n3 = 0;
	  while (n3 < ypoints) {
		  datstream >> yval_temp[n3];
		  //std::cout << " "  << xval_temp[n2] << " " << n2;
		  ++n3;
	  }
	  //Write yvaltemp to yval matrix
	  for (int m = 0; m < yval_temp.size(); m++) {
		  yval(n1, m) = yval_temp[m];
		  //std::cout << xval(n1,m) << " ";
	  }
	  n3 = 0;
      n2 = 0;
      ++n1;
    }
    datstream.close();
    ++filen;
  }
  filestream.close();

  ycalc.resize(yval.size1(),yval.size2());

  datload = true;

}

double funeval_base::printevalwithPE(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, time_t const &tii, int const &itnum, int const &fnum) {


	if (!datload) {
		std::cout << "\n\nno data file\n\n";
		return 0;
	}
	if (!confload) {
		std::cout << "\n\nno config file\n\n";
		return 0;
	}

	ti = tii;

	using namespace boost::numeric::ublas;
	typedef vector<double> uvec;
	typedef matrix<double> umat;
	typedef vector<size_t> uvsiz;
	uvec newpar(par);

	for (size_t i = 0; i < newpar.size(); ++i) {

		if (logten[i])
			newpar[i] = pow(10.0, newpar[i]);
	}

	umat xproj;
	umat::iterator1 xproj_it1;
	umat::iterator2 xproj_it2;
	umat yproj;
	umat yvalproj;
	bool success;

	std::ofstream outfile;
	std::stringstream ofilestream;

	size_t total = 0;
	datnum = 0;
	for (uvsiz::iterator setnums_it = setnums.begin();
	setnums_it != setnums.end();
		++setnums_it) {

		xproj.resize(*setnums_it, xval.size2());
		noalias(xproj) = project(xval, range(total, total + *setnums_it), range(0, xval.size2()));
		yvalproj.resize(*setnums_it, yval.size2());
		noalias(yvalproj) = project(yval, range(total, total + *setnums_it), range(0, yval.size2()));
		yproj.resize(*setnums_it, yval.size2());

		success = function(newpar, xproj, yproj);
		if (!success) {
			for (umat::iterator1 yproj_it1 = yproj.begin1(); yproj_it1 != yproj.end1(); ++yproj_it1)
				for (umat::iterator2 yproj_it2 = yproj_it1.begin(); yproj_it2 != yproj_it1.end(); ++yproj_it2)
					*yproj_it2 = 0.0;
		}

		total += *setnums_it;

		ofilestream.str("");
		ofilestream << "results/data" << ti << "_" << itnum << "_"<< fnum;
		outfile.open((ofilestream.str() + ".txt").c_str(), std::ios::out);
		xproj_it1 = xproj.begin1();
		for (umat::iterator1 yproj_it1 = yproj.begin1(), yvalproj_it1 = yvalproj.begin1(); yproj_it1 != yproj.end1(); ++yproj_it1, ++yvalproj_it1, ++xproj_it1) {
			//Write Xpass
			for (xproj_it2 = xproj_it1.begin(); xproj_it2 != xproj_it1.end(); ++xproj_it2) {
				outfile << *xproj_it2 << " ";
			}
			//Write Yval
			for (umat::iterator2 yvalproj_it2 = yvalproj_it1.begin(); yvalproj_it2 != yvalproj_it1.end(); ++yvalproj_it2) {
				outfile << *yvalproj_it2 << " ";
			}
			//Write Yproject
			for (umat::iterator2 yproj_it2 = yproj_it1.begin(); yproj_it2 != yproj_it1.end(); ++yproj_it2) {
				outfile << *yproj_it2 << " ";
			}
			outfile << "\n";
		}

		outfile.close();

		++datnum;

	}


	++counter;
    //Cacluate Percent error
	boost::numeric::ublas::matrix<double> diffm;
	boost::numeric::ublas::vector<double> diffv;
	double temp, retvar;
	double PE = 0;
	int n = 0;
	temp = 0;
	retvar = 0;
	diffm = yval - ycalc;
	diffv.resize(yval.size1());
	for (int i = 0; i < yval.size2(); i++) {
		diffv = column(diffm, i);
		for (int i = 0; i < diffv.size(); i++) {
			diffv[i] = sqrt(diffv[i] * diffv[i]);
		}
		for (int j = 0; j<yval.size1(); j++) {
			if (j>300) {
				if (yval(i,j)>0) {
					temp = diffv[j] / (yval(i,j));
					PE = PE + temp;
					n++;
				}
			}
		}
	}  
    //std::cout << PE << std::endl;
    PE=(PE/n)*100;
    return PE;
}

void funeval_base::configset(std::string const &confile) {

}

bool funeval_base::function(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::matrix<double> &ypass) {

  return 1;

}
