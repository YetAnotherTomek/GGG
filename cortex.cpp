// ~~~~~~~~~~~~~~~~~~~
// c o r t e x . c p p
// ~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~
//   T.   Erdmann
//   T.R. Sokolowski    
//   2009 - 2012
// ~~~~~~~~~~~~~~~~~~~

#include "cortex.hpp"

using namespace std;

void cortex::init_cortex(MTRand & mtrand,const run_parameter input)
{
    // Declare array of nuclei
    no_r_nucl = input.no_r_nucl;
    no_z_nucl = input.no_z_nucl;

    no_nucl = no_r_nucl*no_z_nucl;

    // Create spatial array of nuclei
    nuclear_array = new nucleus[no_r_nucl*no_z_nucl];
    for(int n=0; n<no_r_nucl*no_z_nucl; n++)
        nuclear_array[n].init(input.no_reac,input.no_spec,input.no_obs);

    // Find the neighbours and assign the parameters
    assign_neighbours(input);
    assign_reactions(mtrand,input);
    assign_diff_para(mtrand,input);
    assign_observables(mtrand,input);
    assign_configuration(mtrand,input);

    // Calculate propensities and next reaction time
    for(unsigned int nucl = 0; nucl < no_nucl; nucl++)
    {
        nuclear_array[nucl].init_propensities(mtrand);
    }
    
    // Init the next-event-time scheduler
    init_nuclear_queue();
}

int cortex::assign_neighbours(const run_parameter input)
{
    no_neig = input.no_neig;
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            nuclear_array[z_nucl+r_nucl*no_z_nucl].init_neighbours(no_neig);
        }
    }
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            unsigned int this_nucl = z_nucl+r_nucl*no_z_nucl;
            unsigned int neig_nucl;

            unsigned int r_neig_nucl;
            unsigned int z_neig_nucl;

            r_neig_nucl = r_nucl;
            if(z_nucl == 0)
                z_neig_nucl = 0;
            else
                z_neig_nucl = z_nucl-1;
            neig_nucl = z_neig_nucl+r_neig_nucl*no_z_nucl;
            nuclear_array[this_nucl].set_neighbour(0,neig_nucl);

            if(r_nucl == 0)
                r_neig_nucl = no_r_nucl-1;
            else
                r_neig_nucl = r_nucl-1;
            z_neig_nucl = z_nucl;
            neig_nucl = z_neig_nucl+r_neig_nucl*no_z_nucl;
            nuclear_array[this_nucl].set_neighbour(1,neig_nucl);

            r_neig_nucl = r_nucl;
            if(z_nucl == no_z_nucl-1)
                z_neig_nucl = no_z_nucl-1;
            else
                z_neig_nucl = z_nucl+1;
            neig_nucl = z_neig_nucl+r_neig_nucl*no_z_nucl;
            nuclear_array[this_nucl].set_neighbour(2,neig_nucl);

            if(r_nucl == no_r_nucl-1)
                r_neig_nucl = 0;
            else
                r_neig_nucl = r_nucl+1;
            z_neig_nucl = z_nucl;
            neig_nucl = z_neig_nucl+r_neig_nucl*no_z_nucl;
            nuclear_array[this_nucl].set_neighbour(3,neig_nucl);
        }
    }
    return(EXIT_SUCCESS);
}

int cortex::assign_reactions(MTRand & mtrand,const run_parameter input)
{
    unsigned int no_reac = input.no_reac;
    if(no_reac != nuclear_array[mtrand.randInt(no_nucl-1)].get_no_reac())
    {
        cerr << "Parameter input.no_reac does not match no_reac in constructed nucleus!" << endl;
        exit(EXIT_FAILURE);
    }
    
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            nuclear_array[z_nucl+r_nucl*no_z_nucl].init_reactions(input.reaction);            
        }
    }    

    return(EXIT_SUCCESS);
}

int cortex::assign_diff_para(MTRand & mtrand,const run_parameter input)
{
    unsigned int no_spec = input.no_spec;
    if(no_spec != nuclear_array[mtrand.randInt(no_nucl-1)].get_no_spec())
    {
        cerr << "Parameter input.no_spec does not match no_spec in constructed nucleus!" << endl;
        exit(EXIT_FAILURE);
    }
    
    double * tmp_diff_para = new double[no_spec];
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            for(unsigned int spec = 0; spec < no_spec; spec++)
            {
                tmp_diff_para[spec] = input.diff_rate_const[spec+z_nucl*no_spec+r_nucl*no_z_nucl*no_spec];
            }
            nuclear_array[z_nucl+r_nucl*no_z_nucl].init_diff_para(tmp_diff_para);
        }
    }
    delete [] tmp_diff_para;

    return(EXIT_SUCCESS);
}

int cortex::assign_observables(MTRand & mtrand,const run_parameter input)
{
    unsigned int no_obs = input.no_obs;
    if(no_obs!= nuclear_array[mtrand.randInt(no_nucl-1)].get_no_obs())
    {
        cerr << "Parameter input.no_obs does not match no_obs in constructed nucleus!" << endl;
        exit(EXIT_FAILURE);
    }
    
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            nuclear_array[z_nucl+r_nucl*no_z_nucl].init_observables(input.obs);
        }
    }

    return(EXIT_SUCCESS);
}

int cortex::assign_configuration(MTRand & mtrand,const run_parameter input)
{
    unsigned int no_spec = input.no_spec;
    if(no_spec != nuclear_array[mtrand.randInt(no_nucl-1)].get_no_spec())
    {
        cerr << "Parameter input.no_spec does not match no_spec in constructed nucleus!" << endl;
        exit(EXIT_FAILURE);
    }
    
    double * tmp_config = new double[no_spec];
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            for(unsigned int spec = 0; spec < no_spec; spec++)
            {
                tmp_config[spec] = input.configuration[spec+z_nucl*no_spec+r_nucl*no_z_nucl*no_spec];
            }
            nuclear_array[z_nucl+r_nucl*no_z_nucl].init_config(tmp_config);
        }
    }
    delete [] tmp_config;

    return(EXIT_SUCCESS);
}

int cortex::init_nuclear_queue(void)
{
    double * tmp_reac_time = new double[no_nucl];

    for(unsigned int nucl = 0; nucl < no_nucl;nucl++)
    {
        tmp_reac_time[nucl] = nuclear_array[nucl].get_time();
    }
    nuclear_queue.init(no_nucl,tmp_reac_time);

    delete [] tmp_reac_time;

    return(EXIT_SUCCESS);
}

void cortex::gillespie_step(MTRand & mtrand)
{
    unsigned int neig_nucl = no_nucl, neig_spec;
    unsigned int this_nucl = nuclear_queue.get_index(0);
    double this_time = nuclear_queue.get_value(0);

    // Do the Gillespie step for the upmost nucleus in the tree
    nuclear_array[this_nucl].gillespie_step(mtrand,neig_nucl,neig_spec);

    // Get the updated nucleus' new next-event time
    unsigned int next_pos = 0;
    double new_time = nuclear_array[this_nucl].get_time();

    // Put the new time into the tree
    nuclear_queue.tree_updt(0,this_time,new_time);

    // Now accordingly update the neighbours if a diffusion event happened
    // (then neig_nucl is set to some value in [0,no_nucl), so the condition is true)
    if(neig_nucl < no_nucl)
    {
        nuclear_array[neig_nucl].add_molecule(neig_spec,this_time,mtrand);
        new_time = nuclear_array[neig_nucl].get_time();

        for(unsigned int nucl = 0; nucl < no_nucl; nucl++)
        {
            if(nuclear_queue.get_index(nucl) == neig_nucl)
            {
                next_pos = nucl;
                break;
            }
        }
        // Put the new time into the tree again
        nuclear_queue.tree_updt(next_pos,this_time,new_time);
    }
    
    // Rearrange the tree, putting the smallest next-time to the root
    nuclear_queue.tree_rsrt(next_pos);
    
}
