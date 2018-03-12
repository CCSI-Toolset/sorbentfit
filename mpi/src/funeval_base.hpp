// See LICENSE.md for license anc copyright details.

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

  bool eval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::vector<double> &ypass);
  double objeval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten);
  void printeval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, time_t const &ti);

  virtual void configset(std::string const &confile);
  void datfile(std::string const &datfiles);
  long getcounter() {return counter;};
  void setcounter(long count) {counter = count;};

protected:

  virtual bool function(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::vector<double> &ypass);
  boost::numeric::ublas::matrix<double> xval;
  boost::numeric::ublas::vector<double> yval;
  boost::numeric::ublas::vector<double> ycalc;
  boost::numeric::ublas::vector<size_t> setnums;
  bool datload;
  bool printme;
  bool confload;
  std::string ofilename;
  long counter;
  time_t ti;
  long datnum;


};

bool funeval_base::eval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::vector<double> &ypass) {

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

double funeval_base::objeval(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<bool> const &logten) {

  if (!datload) {
    std::cout << "\n\nno data file\n\n";
    return 0.0;
  }
  if (!confload) {
    std::cout << "\n\nno config file\n\n";
    return 0.0;
  }

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
  uvec yproj;
  bool success;

  size_t total=0;
  for (uvsiz::iterator setnums_it=setnums.begin();
       setnums_it != setnums.end();
       ++setnums_it) {

    xproj.resize(*setnums_it, xval.size2());
    noalias(xproj) = project(xval, range(total,total + *setnums_it), range(0,xval.size2()));
    yproj.resize(*setnums_it);

    success = function(newpar, xproj, yproj);
    if (!success) {
      for (uvec::iterator yproj_it=yproj.begin();
	   yproj_it != yproj.end();
	   ++yproj_it)
	*yproj_it = 0.0;
    }

    noalias(project(ycalc, range(total,total + *setnums_it))) = yproj;
    total += *setnums_it;
  }

  return norm_2(yval - ycalc);

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
  uvec yproj;
  uvec yvalproj;
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
    yvalproj.resize(*setnums_it);
    noalias(yvalproj) = project(yval, range(total, total + *setnums_it));
    yproj.resize(*setnums_it);

    success = function(newpar, xproj, yproj);
    if (!success) {
      for (uvec::iterator yproj_it=yproj.begin();
	   yproj_it != yproj.end();
	   ++yproj_it)
	*yproj_it = 0.0;
    }

    total += *setnums_it;

    ofilestream.str("");
    ofilestream << "results/data" << ti << "_" << datnum;
    outfile.open((ofilestream.str()+".txt").c_str(),std::ios::out);
    xproj_it1=xproj.begin1();
    for (uvec::iterator yproj_it=yproj.begin(), yvalproj_it=yvalproj.begin();
	 yproj_it != yproj.end();
	 ++yproj_it, ++yvalproj_it, ++xproj_it1) {

      for (xproj_it2=xproj_it1.begin();
	   xproj_it2 != xproj_it1.end();
	   ++xproj_it2) {

	outfile << *xproj_it2 << " ";
      }

      outfile << *yproj_it << " ";
      outfile << *yvalproj_it << "\n";
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

  size_t numfiles, xpoints;
  filestream >> numfiles;
  filestream >> xpoints;
  setnums.resize(numfiles);

  vector<double> xval_temp(xpoints);

  std::string exdata;
  size_t n1, n2, filen;

  n1=0;
  filen=0;
  while (filestream >> exdata) {

    std::ifstream datstream;
    datstream.open(exdata.c_str(), std::ios::in);
    datstream >> setnums[filen];

    xval.resize(xval.size1()+setnums[filen],xpoints,true);
    yval.resize(yval.size()+setnums[filen],true);

    n2=0;
    while (datstream >> xval_temp[0]) {

      ++n2;

      while (n2 < xpoints) {
	datstream >> xval_temp[n2];
	++n2;
      }

      row(xval, n1) = xval_temp;
      datstream >> yval[n1];

      n2=0;
      ++n1;
    }
    datstream.close();
    ++filen;
  }
  filestream.close();

  ycalc.resize(yval.size());

  datload = true;

}

void funeval_base::configset(std::string const &confile) {

}

bool funeval_base::function(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::matrix<double> &xpass, boost::numeric::ublas::vector<double> &ypass) {

  return 1;

}
