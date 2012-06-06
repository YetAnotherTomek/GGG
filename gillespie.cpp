//  ~~~~~~~~~~~~~~~~~~~~~~~~~
//  g i l l e s p i e . c p p
//  ~~~~~~~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//        T.R. Sokolowski    
//        T.   Erdmann
//        2009 - 2012
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <iostream>

#include "parameters.hpp"
#include "simulation.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 2) 
    {
        cerr << "You have to pass arguments to function main!" << endl;
        exit(EXIT_FAILURE);
    }
    
    // Create input parameter class and read the input.
    // In case of previous simulation interrupts, also read the tempfile.
    run_parameter input;    
    input.read_filenames(argc, argv);
    input.read_parameter();
    input.read_temp_file();
    
    simulation sim(&input);

    switch(input.what_to_do)
    {
        case 0:
          
            /* In this mode, only a time-trajectory of the observables
               and an input and temp file for restarts is generated.
            */
            sim.mode = TRAJECTORY;
            sim.init(&input);
            
            sim.run(&input);
            
            sim.output_input_file(&input);
            sim.output_temp_file(&input);
            
            break;
            
        case 1:
          
            /* In this mode, also the accumulation of the observable moments
               is performed, processed and written out at the end of the sim.
            */
            sim.mode = AVERAGES;
            sim.init(&input);
            
            sim.run(&input);
            
            sim.output_input_file(&input);
            sim.output_temp_file(&input);
            
            // If this is the last interrupt, run the averaging.
            // Note that input.max_interrupts initially contains the number
            // of previous outputs.
            // This number is increased by sim.output_temp_file()
            if(input.interrupts>=input.max_interrupts){
              
                sim.calculate_thresholds(&input);  
                sim.output_averages(&input);
            }                
            
            break;
            
        case 2:

            /* In this mode, everything is calculated and written out,
               i.e. a time-trajectory of the observables is generated,
               the observable moments are accumulated and processed at
               the end and also the value of the FFS reaction coordinate
               is calculated at the end.
            */
            sim.mode = ALL;
            sim.init(&input);
            
            sim.run(&input);
            
            sim.output_input_file(&input);
            sim.output_temp_file(&input);
            
            // If this is the last interrupt, run the averaging.
            // Note that input.max_interrupts initially contains the number
            // of previous outputs.
            // This number is increased by sim.output_temp_file()
            if(input.interrupts>=input.max_interrupts){
              
                sim.calculate_thresholds(&input);  
                sim.output_averages(&input);
            }
            
            cout.flush();
            
            break;

        default :
            cerr << "Don't know what_to_do. Select [0-2]!\n";
            
    }
    
    return(EXIT_SUCCESS);
}
