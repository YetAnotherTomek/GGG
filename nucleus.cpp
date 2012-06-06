//  ~~~~~~~~~~~~~~~~~~~~~
//  n u c l e u s . c p p
//  ~~~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~~~~
//    T.   Erdmann
//    T.R. Sokolowski    
//    2009 - 2012
// ~~~~~~~~~~~~~~~~~~~~~~

#include "nucleus.hpp"

using namespace std;


nucleus::nucleus()
{  
    // Create an empty dummy nucleus
    // nucleus.init() has to be run from upstream functions
    // to fully initialize the nuclear parameters.
    no_reac = 0;
    no_spec = 0;
    no_obs  = 0;
}

// Principal initialization: Assign number of reactions, species and observables
// for subsequent array creation
void nucleus::init(const int nr, const int ns, const int no)
{
    no_reac = nr;
    no_spec = ns;
    no_obs  = no;
}

// Initialize reaction network
void nucleus::init_reactions(const reaction_class* reaction_ext)
{
    reaction = new reaction_class[no_reac];
    for(unsigned int reac=0; reac<no_reac; reac++)   reaction[reac] = reaction_ext[reac];
    
}

// Initialize diffusion parameters
void nucleus::init_diff_para(const double* diff_para_ext)
{
    diff_para = new double[no_spec];
    for(unsigned int spec = 0; spec < no_spec; spec++)
    {
        diff_para[spec] = diff_para_ext[spec];
    }
}

// Initialize observables
void nucleus::init_observables(const observable* obs_ext)
{
    obs = new observable[no_obs];
    for(unsigned int s = 0; s < no_obs; s++)
    {
        obs[s] = obs_ext[s];
    }    
}

// Initialize configuration
void nucleus::init_config(const double* config_ext)
{
    configuration = new double[no_spec];
    for(unsigned int spec = 0; spec < no_spec; spec++)
    {
        configuration[spec] = config_ext[spec];
    }
    
    update_observables();
}

// Initialize reaction & diffusion propensities
void nucleus::init_propensities(MTRand &mtrand)
{
    // Propensities
    sum_reac_prop = 0.0;
    reac_prop = new double[no_reac];
    for(unsigned int reac = 0; reac < no_reac; reac++)
    {
        reac_prop[reac] = reaction[reac].prop(configuration);
        sum_reac_prop += reac_prop[reac];
    }
    sum_diff_prop = 0.0;
    diff_prop = new double[no_spec];
    for(unsigned int spec = 0; spec < no_spec; spec++)
    {
        diff_prop[spec] = configuration[spec]*diff_para[spec];
        sum_diff_prop += diff_prop[spec];
    }
    // Next reaction time
    if((sum_reac_prop+sum_diff_prop) > 0)
    {
        next_reac_time = mtrand.randExp(1.0/(sum_reac_prop+sum_diff_prop));
    }
    else
    {
        next_reac_time = ARTIFICIAL_INFINITY;
    }
}

// Set the neighbours of a nucleus
void nucleus::init_neighbours(const unsigned int no_neighbours)
{
    no_neig = no_neighbours;
    my_neig = new unsigned int[no_neig];
}

void nucleus::set_neighbour(const unsigned int index,const unsigned int position)
{
    if(index >= no_neig)
    {
        cerr << "Tried to assign more neighours than possible! " << endl;
        exit(EXIT_FAILURE);
    }
    my_neig[index] = position;
}

// Gillespie step
// This is the central routine called for subsequent Monte Carlo updates
int nucleus::gillespie_step(MTRand & mtrand,unsigned int & neig_nucl,unsigned int & neig_spec)
{
  
    // Check if at the last update everything went fine
    bool error = false;
    for(unsigned int spec = 0; spec < no_spec; spec++)
    {
        if(configuration[spec]<0)   error = true;
    }
    if(sum_reac_prop<=0.0)          error = true;
    
    if(error){
        cerr << "NEGATIVE COPY NUMBER AND/OR NON-POSITIVE PROPENSITY SUM! Something is wrong!" << endl;
        cerr << "sum_reac_prop = " << sum_reac_prop << endl;
        cerr << "sum_diff_prop = " << sum_diff_prop << endl;
        cerr << "last_reac = " << last_reac << endl;
        cerr << "Configuration:" << endl;
        for(unsigned int spec = 0; spec < no_spec; spec++)
        {
            cerr << configuration[spec] << " ";
        }
        cerr << endl;
        exit(EXIT_FAILURE);
    }
    
    if(mtrand.randDblExc()*(sum_reac_prop+sum_diff_prop) <= sum_reac_prop)
    {
        // Reaction step
        double rand_prop = mtrand.randDblExc()*sum_reac_prop;
        unsigned int this_reac;
        for(this_reac = 0; this_reac < no_reac; this_reac++)
        {
            rand_prop -= reac_prop[this_reac];
            if(rand_prop < 0) break;
        }
        // Fire the reaction = update configuration
        reaction[this_reac].fire(configuration);
        
        // Update propensities
        sum_reac_prop = 0.0;
        for(unsigned int reac = 0; reac < no_reac; reac++)
        {
            reac_prop[reac] = reaction[reac].prop(configuration);
            sum_reac_prop += reac_prop[reac];
        }
        sum_diff_prop = 0.0;
        for(unsigned int spec = 0; spec < no_spec; spec++)
        {
            diff_prop[spec] = configuration[spec]*diff_para[spec];
            sum_diff_prop += diff_prop[spec];
        }
        // Update next reaction time
        if((sum_reac_prop+sum_diff_prop) > 0)
        {
            next_reac_time += mtrand.randExp(1.0/(sum_reac_prop+sum_diff_prop));
        }
        else
        {
            next_reac_time = ARTIFICIAL_INFINITY;
        }
        // Remember which reaction fired
        last_reac = this_reac;
        // no neighbouring nucleus was bothered
    }
    else
    {	// Diffusion step
        unsigned int this_spec;
        double rand_prop = mtrand.randDblExc()*sum_diff_prop;
        for(this_spec = 0; this_spec < no_spec; this_spec++)
        {
            rand_prop -= diff_prop[this_spec];
            if(rand_prop < 0) break;
        }
        // Update configuration
        // ( only in this nucleus; the update in the neighbouring
        //   nucleus receiving the molecule is done by the upstream
        //   routine using nucleus::add_molecule )
        configuration[this_spec] -= 1;

        // Update propensities
        sum_reac_prop = 0.0;
        for(unsigned int reac = 0; reac < no_reac; reac++)
        {
            reac_prop[reac] = reaction[reac].prop(configuration);
            sum_reac_prop += reac_prop[reac];
        }
        sum_diff_prop = 0.0;
        diff_prop[this_spec] = configuration[this_spec] * diff_para[this_spec];
        for(unsigned int spec = 0; spec < no_spec; spec++)
        {
            sum_diff_prop += diff_prop[spec];
        }
        // Update next reaction time
        if((sum_reac_prop+sum_diff_prop) > 0)
        {
            next_reac_time += mtrand.randExp(1.0/(sum_reac_prop+sum_diff_prop));
        }
        else
        {
            next_reac_time = ARTIFICIAL_INFINITY;
        }
        // Identify receiving neighbour (chosen randomly from the neighbours)
        neig_nucl = my_neig[mtrand.randInt(no_neig-1)];
        neig_spec = this_spec;
    }
    
    update_observables();
    
    return(EXIT_SUCCESS);
}

int nucleus::add_molecule(const unsigned int neig_spec,const double this_time,MTRand &mtrand)
{
    // Update configuration and observables
    configuration[neig_spec] += 1;    
    update_observables();

    // Update propensities
    sum_reac_prop = 0.0;
    for(unsigned int reac = 0; reac < no_reac; reac++)
    {
        reac_prop[reac] = reaction[reac].prop(configuration);
        sum_reac_prop += reac_prop[reac];
    }
    sum_diff_prop = 0.0;
    diff_prop[neig_spec] = configuration[neig_spec]*diff_para[neig_spec];
    for(unsigned int spec = 0; spec < no_spec; spec++)
    {
        sum_diff_prop += diff_prop[spec];
    }
    // Update next reaction time
    if((sum_reac_prop+sum_diff_prop) > 0)
    {
        next_reac_time = this_time + mtrand.randExp(1.0/(sum_reac_prop+sum_diff_prop));
    }
    else
    {
        next_reac_time = ARTIFICIAL_INFINITY;
    }
    return(EXIT_SUCCESS);
}

void nucleus::update_observables(void)
{
    
    for(int s=0; s<no_obs; s++)
    {
        obs[s].value = 0.0;
        for(int c=0; c<obs[s].no_components; c++)   obs[s].value += 1.0*configuration[obs[s].component[c]];        
    }
}



