// ~~~~~~~~~~~~~~~~~~~~~~~
// r e a c t i o n . h p p
// ~~~~~~~~~~~~~~~~~~~~~~~

#ifndef _REACTION_HPP
#define _REACTION_HPP

#include <cfloat>
#include <cstdlib>
#include <iostream>


using namespace std;

enum reac_type { PRODUCTION, DECAY, BINDING, UNBINDING, DIMERIZATION };

class reaction_class
{
        
    public:

        reaction_class(void);
        double prop(double *conf);
        void fire(double *conf);
        void init(int, reac_type, double, int, int, int, int);
        void info(void);
        
        reaction_class& operator=(const reaction_class&);
        
        int number;
        reac_type type;
        int partner[2];
        int product[2];
        double rate;
};


#endif // _REACTION_HPP

