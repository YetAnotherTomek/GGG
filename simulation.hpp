// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
// s i m u l a t i o n . h p p
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef _SIMULATION_HPP
#define _SIMULATION_HPP

#include "./Tools/mersenne_twister.hpp"
#include "./Tools/random_variable.hpp"
#include "./Tools/histo1D.hpp"

#include "parameters.hpp"
#include "cortex.hpp"

using namespace std;


enum sim_mode { TRAJECTORY, AVERAGES, ALL };

class simulation
{
  public:
      simulation(run_parameter*);
      ~simulation();
      
      void init(run_parameter*);
      void run(run_parameter*);
      
      void output_temp_file(run_parameter*);
      void output_input_file(run_parameter*);
      void output_averages(run_parameter*);
      
      void calculate_thresholds(run_parameter*);
      
      
  private:
      void measure(void);      
      void output_instant(run_parameter*);
      
  public:
    
      // Siumulation mode
      sim_mode mode;
      
      // The simulated system
      cortex fly;
      int no_r_nucl;
      int no_z_nucl;
      double subv_size;
      
      // Mersenne Twister RNG
      MTRand mtrand;
 
      // Output files and related
      ofstream traf;
      ofstream trav;
      ofstream tempf;
      char tmp_string[500];
      char trav_name[1000];
      
      // Observables / Measurement
      unsigned int no_obs;              // no. of observables
      double* threshold;
      random_var** obs_config;          // random variable array for configurations
      random_var** avg_obs_config;      // random variable array for averaged configurations
      histo1D* obs_boundary_pos;        // histogram for bnd. positions
      histo1D* obs_boundary_no;         // histogram for no. of boundaries detected
      unsigned int* no_obs_boundary;    // array for measured number of boundary crossings, needed in the production run
      
      double* obs_int;                  // system-wide integral of observables
      double* CoM;                      // arrays for concentration CoM measurement
      double* avg_obs_temp;             // temporary averaging array
      
      unsigned long N_measurements;     // monitors the number of performed measurements
      int interrupts;                   // contains no. of previous temporary outputs
      int frame;			// monitors the no. of instantaneous outputs
      
       // Variables needed for time control
      double last_time, next_time, meas_time;
      double bgn_meas, end_meas;
      double weight_time;
      double ac_run_time;
      unsigned long ac_run_step;
      unsigned long ac_meas_step;
      unsigned long this_sim_run_step;     
      unsigned long break_steps;        // number of simulation steps between two temporary outputs
      unsigned long instant_steps;      // number of simulation steps between two instant outputs
      
      // Flags indicating different states of the simulation
      // The flags are always set to false initially, also after an interrupted output
      // At the begin of every simulation the state is set in the first loop automatically
      // before any measurement etc. is performed
      bool relaxed_flag;    // whether the system is relaxed
      bool measure_flag;    // whether a measurement should be performed in the next loop
      bool instant_flag;    // whether to output instantaneous data in the next loop
                            
      bool break_now;       // The main loop is broken if this is set to true      
              
};


#endif // _SIMULATION_HPP
