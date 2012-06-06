// ~~~~~~~~~~~~~~~~~~~
// c o r t e x . h p p
// ~~~~~~~~~~~~~~~~~~~

#ifndef _CORTEX_HPP
#define _CORTEX_HPP

#include "./Tools/mersenne_twister.hpp"
#include "./Tools/binary_tree.hpp"

#include "parameters.hpp"
#include "nucleus.hpp"

using namespace std;

class cortex
{
    // Methods
    public:
        cortex(void){};
        ~cortex(void){delete[] nuclear_array;};
        void init_cortex(MTRand & mtrand,const run_parameter input);
        void gillespie_step(MTRand & mtrand);
        inline double get_config(const unsigned int r_nucl,const unsigned int z_nucl,const unsigned int spec);
        inline double get_observable(const unsigned int r_nucl,const unsigned int z_nucl,const unsigned int obs);
        inline double get_time(const unsigned int r_nucl,const unsigned int z_nucl);
        inline double set_time(const unsigned int r_nucl,const unsigned int z_nucl, double t);
        
    private:
        int assign_neighbours(const run_parameter input);
        int assign_reactions(MTRand & mtrand,const run_parameter input);
        int assign_diff_para(MTRand & mtrand,const run_parameter input);
        int assign_observables(MTRand & mtrand,const run_parameter input);
        int assign_configuration(MTRand & mtrand,const run_parameter input);
        int init_nuclear_queue(void);
    
    // Variables
    public:
        nucleus* nuclear_array;
        binary_tree nuclear_queue;
        
    private:
        unsigned int no_nucl;
        unsigned int no_r_nucl;
        unsigned int no_z_nucl;
        unsigned int no_neig;
};

inline double cortex::get_config(const unsigned int r_nucl,const unsigned int z_nucl,const unsigned int spec)
{
    return nuclear_array[z_nucl+r_nucl*no_z_nucl].get_config(spec);
}

inline double cortex::get_observable(const unsigned int r_nucl,const unsigned int z_nucl,const unsigned int obs)
{
    return nuclear_array[z_nucl+r_nucl*no_z_nucl].get_observable(obs);
}

inline double cortex::get_time(const unsigned int r_nucl,const unsigned int z_nucl)
{
    return nuclear_array[z_nucl+r_nucl*no_z_nucl].get_time();
}

inline double cortex::set_time(const unsigned int r_nucl,const unsigned int z_nucl, double t)
{
    return nuclear_array[z_nucl+r_nucl*no_z_nucl].set_time(t);
}

#endif // _CORTEX_HPP

