//  ~~~~~~~~~~~~~~~~~~~~~~~~~~
//  p a r a m e t e rs . h p p
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#include<ctime>
#include<fstream>
#include<iostream>
#include<cstdlib>
#include<stdio.h>
#include<string.h>
#include<assert.h>

#include "./Tools/random_variable.hpp"
#include "./Tools/histo1D.hpp"
#include "./Tools/mersenne_twister.hpp"
#include "reaction.hpp"
#include "observable.hpp"

static double sqrt3 = 1.732050807568877;

using namespace std;

class filenames
{

  public:
    
        char* input_file;
        char* new_input_file;
        char* reactions_file;
        char* observables_file;
        char* temp_file;
        char* current_temp_file;
        char* traject_file;
        char* average_file;
        char* bnd_dst_file;
        
        filenames(){
         
            input_file        = new char[500];
            new_input_file    = new char[500];
            reactions_file    = new char[500];
            observables_file  = new char[500];
            temp_file         = new char[500];
            current_temp_file = new char[500];
            traject_file      = new char[500];
            average_file      = new char[500];
            bnd_dst_file      = new char[500];
        };        
};

class parameter
{
  
  public:
    
    parameter() {};
    
    // Simulation properties and RNG
    unsigned int what_to_do;
    unsigned int no_trajec;
    unsigned int max_interrupts;
    unsigned int interrupts;
    unsigned int instant_outputs;
    bool output_all;
    
    unsigned int init_random_seed;
    unsigned int random_seed;
    MTRand::uint32* generator_state;
    unsigned int generator_state_len;
    // NOTE: init_random_seed is the standard seed passed via
    // the input_file. However it can be overriden by random_seed
    // which is read from temp_file. This can be used to reset the
    // seed at interrupts, for example.

    // Time constants
    double rlx_time;
    double meas_time;
    double main_time;
    unsigned long rlx_steps;
    unsigned long meas_steps;
    unsigned long main_steps;
    
    // Time variables
    double last_time;
    double next_time;
    unsigned long ac_run_step;
    double ac_run_time;
    unsigned long ac_meas_step;
    double bgn_meas;
    double end_meas;
    
    // Next-event tree
    int tree_size;
    double* tree_value;
    int* tree_index;
    double* event_time;

    // Geometry    
    double subv_size;

    unsigned int no_nucl;
    unsigned int no_z_nucl;
    unsigned int no_r_nucl;
    unsigned int no_reac;
    unsigned int no_para;
    unsigned int no_spec;
    unsigned int no_obs;
    unsigned int no_neig;    

    // Parameters and configuration
    reaction_class* reaction;
    double* reac_parameter;
    double* diff_rate_const;
    double* configuration;

    // Observables
    observable* obs;
    double* threshold;
    
    random_var** obs_config;
    random_var** avg_obs_config;
    histo1D* obs_boundary_pos;
    histo1D* obs_boundary_no;    
    
    char dump[200];

};

class run_parameter : public parameter, public filenames
{
    public:
        run_parameter() : filenames(), parameter() {};
        ~run_parameter() {};
                
        int read_filenames(int argc, char *argv[]);
        int read_parameter(void);
        int read_temp_file(void);
        int open_output(ofstream& fout,const char * filename, bool header);
        int open_output(ofstream& fout,const char * filename, bool header, int prec);
        int open_output(ofstream& fout,const char * filename, ios_base::openmode, int prec);
        int write_header(ofstream& fout,const char * filename);
        int close_output(ofstream& fout);
};
#endif // _PARAMETERS_HPP
