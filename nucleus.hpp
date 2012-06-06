//  ~~~~~~~~~~~~~~~~~~~~~
//  n u c l e u s . h p p
//  ~~~~~~~~~~~~~~~~~~~~~

#ifndef _NUCLEUS_HPP
#define _NUCLEUS_HPP

#include <cfloat>
#include <cstdlib>

#include "./Tools/array2D.hpp"
#include "./Tools/mersenne_twister.hpp"
#include "reaction.hpp"
#include "observable.hpp"

using namespace std;

class nucleus
{
    // Methods
    public:
        nucleus();
        ~nucleus(){ delete[] reaction; delete[] diff_para; delete[] configuration; delete[] obs;
                    delete[] diff_prop; delete[] reac_prop; delete[] my_neig; };
        void init(const int nr, const int ns, const int np);
        void init_reac_para(const double* reac_parameter);
        void init_reactions(const reaction_class* reaction_ext);
        void init_diff_para(const double* diff_parameter);
        void init_observables(const observable* obs_ext);
        void init_config(const double *cnfg);
        void init_neighbours(const unsigned int no_neighbours);
        void set_neighbour(const unsigned int index,const unsigned int position);
        void init_propensities(MTRand & mtrand);
        
        int  gillespie_step(MTRand & mtrand,unsigned int & neighbour,unsigned int & neig_spec);
        int  add_molecule(const unsigned int this_spec,const double time,MTRand &mtrand);
        void update_observables(void);        
        
        inline unsigned int get_no_reac(void){return no_reac;};
        inline unsigned int get_no_spec(void){return no_spec;};
        inline unsigned int get_no_obs(void){return no_obs;};
        inline double get_config(const unsigned int spec){  return configuration[spec]; };
        inline double get_observable(const unsigned int s){ return obs[s].value; };
        inline double get_time(void){return next_reac_time;};
        inline double set_time(double t){next_reac_time = t;};
        inline double get_prop(const unsigned int idx){ return next_reac_time;};
        inline double get_neig(const unsigned int idx){ return no_neig;};

    // Variables
    public:
        int last_reac;
        
    private:
      
        static const double ARTIFICIAL_INFINITY = 1.0e+300; // NOTE was DBL_MAX before, but this very rarely caused problems at temporary outputs
        
        unsigned int no_reac;
        unsigned int no_spec;
        unsigned int no_obs;
        
        // Parameters
        reaction_class* reaction;
        double* diff_para;
        double* configuration;
        observable* obs;
        
        // Neighbour table
        unsigned int no_neig;
        unsigned int* my_neig;
        
        // Propensities
        double* diff_prop;
        double* reac_prop;
        double sum_diff_prop;
        double sum_reac_prop;
        
        double next_reac_time;        
        
};

#endif // _NUCLEUS_HPP
