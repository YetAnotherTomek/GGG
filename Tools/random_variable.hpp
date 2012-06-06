// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// r a n d o m _ v a r i b l e . h p p : 
//
// This contains an array which accumulates moments of a random
// variable during simulation.
//
// These then can be used at the end of the sim. to calculate
// the mean, variance and skewness.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef _RANDOM_VARIABLE_HPP
#define _RANDOM_VARIABLE_HPP

#include<cstdlib>
#include<cmath>
#include<algorithm>

class random_var
{
        public:
                random_var(void);
                random_var(int size);
                ~random_var(void) {delete [] moments;};
                double moment(const int);
                double mean(void);
                double variance(void);
                double std_dev(void);
                double skewness(void);
                void rnd_add(const double,const double);
                void operator=(random_var &);
                bool add_everything(random_var&);
                bool operator+=(random_var&);
        public:
                unsigned int max_moments;
                unsigned int events;
                double weight;
                double* moments;
};

#endif // _RANDOM_VARIABLE_HPP
