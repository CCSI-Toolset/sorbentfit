# sorbentfit
Sorbent kinetic/equilibrium model code for parameter estimation for the 1st generation sorbent model. This contains the Sorbentfit suite of routines that fit combined thermodynamic and kinetic models of amine-based CO2 solid sorbents to laboratory and bench-scale data.  There are two levels of models included with the code: a lumped kinetic model and a reaction-diffusion model with zwitterion chemistry. The lumped kinetic model is faster, while the reaction-diffusion model is more accurate. Both models can be fit to thermogravimetry and fixed-bed data. Fits can be either point estimates of parameters, achieved using a particle swarm optimizer, or Bayesian calibration, which results in a distribution of parameter values. Sorbentfit does not require any third-party software; it runs in Windows, Mac and UNIX/Linux operating systems.



## Development Practices

* Code development will be peformed in a forked copy of the repo. Commits will not be 
  made directly to the ngt-archive repo. Developers will submit a pull 
  request that is then merged by another team member, if another team member is available.
* Each pull request should contain only related modifications to a feature or bug fix.  
* Sensitive information (secret keys, usernames etc) and configuration data 
  (e.g database host port) should not be checked in to the repo.
* A practice of rebasing with the main repo should be used rather that merge commmits.

## Getting Started

TBD

## Authors

* Keenan Kocan
* David Mebane
* Brian Logsdon
* Conor Pyles

See also the list of [contributors](https://github.com/CCSI-Toolset/sorbentfit/contributors) who participated in this project.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, 
see the [tags on this repository](https://github.com/sorbentfit/tags). 

## License

See [LICENSE.md](LICENSE.md) file for details

## Copyright Notice

TBD
