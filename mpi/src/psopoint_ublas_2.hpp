
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/cstdint.hpp>

#ifndef _PSOPOINT
#define _PSOPOINT

class psopoint {

public:

  psopoint() {isinit = false;};
  psopoint(boost::numeric::ublas::vector<double> const &pos, int const &num, int const &totsize, funeval_base &fun) {init(pos, num, totsize, fun);};

  bool init(boost::numeric::ublas::vector<double> const &pos, int const &numi, int const &totsize, funeval_base &fun);

  void move(funeval_base &fun, boost::variate_generator< boost::mt19937&, boost::uniform_01<> > &randstrm);
  double getcost() {return costval;};
  // boost::numeric::ublas::vector<double> getpos() {return pos;};
  // boost::numeric::ublas::vector<double> getposbest() {return posbest;};

  static boost::numeric::ublas::vector<bool> logten;
  static boost::numeric::ublas::vector<double> costvalist;
  static boost::numeric::ublas::matrix<double> poslist;
  static boost::numeric::ublas::matrix<double> posbestlist;

protected:

  double distance(size_t i);

  double costval;
  // boost::numeric::ublas::vector<double> posbest;
  // boost::numeric::ublas::vector<double> pos;


  boost::numeric::ublas::vector<double> velocity;
  boost::numeric::ublas::vector<double> Pm;
  boost::numeric::ublas::vector<double> dist;
  int num;
  int neigh1;
  int neigh2;
  double dist1;
  double dist2;
  double phi1;
  double phi2;

  //  double dist_calc;

  bool isinit;

};

boost::numeric::ublas::vector<bool> psopoint::logten;
boost::numeric::ublas::vector<double> psopoint::costvalist;
boost::numeric::ublas::matrix<double> psopoint::poslist;
boost::numeric::ublas::matrix<double> psopoint::posbestlist;

bool psopoint::init(boost::numeric::ublas::vector<double> const &posi, int const &numi, int const &totsize, funeval_base &fun) {
  using namespace boost::numeric::ublas;
  typedef vector<double> uvec;

  if (logten.size() != posi.size() || poslist.size2() != posi.size() || costvalist.size() == 0 || posbestlist.size2() != posi.size()) {
    std::cout << "\n\npsopoint init fail\n\n" << poslist.size2();
    return 0;
  }

  num = numi;

  // pos.resize(posi.size());
  // pos = posi;
  // posbest.resize(posi.size());
  // posbest = posi;
  velocity.resize(posi.size());
  Pm.resize(posi.size());
  dist.resize(costvalist.size()-1);
  for (uvec::iterator velocity_it=velocity.begin();
       velocity_it != velocity.end();
       ++velocity_it)
    *velocity_it = 0.0;
  costval = fun.objeval(posi, logten);
  costvalist[num] = costval;

  row(poslist, num) = posi;
  row(posbestlist, num) = posi;
	
  isinit = true;

  return 1;

}

void psopoint::move(funeval_base &fun, boost::variate_generator< boost::mt19937&, boost::uniform_01<> > &randstrm)	{

  using namespace boost::numeric::ublas;

  for (int i=0; i<costvalist.size(); i++)	{

    if (i == num)
      continue;
    else if (i < num)
      dist[i] = distance(i);
    else
      dist[i-1] = distance(i);
  }

  dist1 = dist[0];
  dist2 = dist[1];
  neigh1 = 0;
  neigh2 = 1;
  for (int i=2; i<dist.size(); i++)	{

    if (dist[i] < dist1 || dist[i] < dist2)	{

      if (dist1 < dist2)	{

	neigh2 = i;
	dist2 = dist[i];
      }
      else	{

	neigh1 = i;
	dist1 = dist[i];
      }
    }
  }

  if (neigh1 >= num)
    ++neigh1;

  if (neigh2 >= num)
    ++neigh2;
	
  phi1 = (randstrm())*4.1/2;
  phi2 = (randstrm())*4.1/2;
  // for (size_t i=0; i<velocity.size(); i++)	{
  //   Pm[i] = (phi1*posbestlist(neigh1,i) + phi2*posbestlist(neigh2,i))/(phi1 + phi2);
  //   velocity[i] = 0.7298*(velocity[i] + (phi1 + phi2)*(Pm[i] - poslist(num,i)));
  //   poslist(num,i) += velocity[i];
  // }

  Pm = (phi1*row(posbestlist, neigh1) + phi2*row(posbestlist, neigh2))/(phi1 + phi2);
  velocity = 0.7298*(velocity + (phi1 + phi2)*(Pm - row(poslist, num)));
  row(poslist, num) += velocity;
	
  //vector<double> vec(poslist.size2());
  //vec = row(poslist, num);
  costval = fun.objeval(row(poslist,num), logten);
	
  if (costval < costvalist[num])	{
    costvalist[num] = costval;
    // for (size_t i=0; i<velocity.size(); i++)
    //   posbestlist(num,i) = poslist(num,i);

    row(posbestlist, num) = row(poslist, num);
  }

}

double psopoint::distance(size_t i)	{

  //dist_calc = 0.0;
  // for (size_t j=0; j<velocity.size(); j++)
  //   dist_calc += (poslist(i,j) - poslist(num,j))*(poslist(i,j) - poslist(num,j));

  // return sqrt(dist_calc);

  return norm_2(row(poslist, i) - row(poslist, num));

}

#endif
