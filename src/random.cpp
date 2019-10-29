#include "random.h"

RandomNumbers::RandomNumbers (unsigned long int s)
{
	if (s==0){ 
		std::random_device rd;
		seed = rd();
	}else{
		seed=s;
	}
	rng = std::mt19937(seed);
}

void RandomNumbers::uniform_double(std::vector<double>& tab, double lower, double upper)
{
	for (int n=0; n<tab.size(); ++n){
		tab[n]= uniform_double(lower, upper);
	}
}

double RandomNumbers::uniform_double(double lower, double upper)
{
	std::uniform_real_distribution<> dunif (lower, upper);
	return dunif(rng);
}

void RandomNumbers::normal(std::vector<double>& tab, double mean, double sd)
{	
	for (int n=0; n<tab.size(); ++n){
		tab[n]= normal(mean, sd);
	}
}

double RandomNumbers::normal(double mean, double sd)
{
	std::normal_distribution<double> dnorm (mean, sd);
	return dnorm(rng);
}

void RandomNumbers::poisson(std::vector<int>& tab, double mean)
{
	for (int n=0; n<tab.size(); ++n){
		tab[n]= poisson(mean);
	}
}

int RandomNumbers::poisson(double mean)
{
	std::poisson_distribution<> dpoiss(mean);
	return dpoiss(rng);
}


