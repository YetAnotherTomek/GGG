//  ~~~~~~~~~~~~~~~~~~~~~~~~~~
//  p a r a m e t e rs . c p p
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~~~~~~~~~
//        T.R. Sokolowski    
//        T.   Erdmann
//        2009 - 2012
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include"parameters.hpp"

using namespace std ;

int run_parameter::read_filenames(int argc, char *argv[])
{
    ifstream in_file(argv[1]);
    if(!in_file)
    {
        cerr << "Failed to open in_file!" << endl;
        exit(EXIT_FAILURE);
    }
    in_file >> input_file;
    in_file >> new_input_file;
    in_file >> reactions_file;
    in_file >> observables_file;
    in_file >> temp_file;
    in_file >> traject_file;
    in_file >> average_file;
    in_file >> bnd_dst_file;
    in_file.close();

    return(EXIT_SUCCESS);
}

int run_parameter::read_parameter()
{
    double tmp;
    int i;
    char c;
    
    // Open general input file
    ifstream in_file(input_file);
    if(!in_file)
    {
        cerr << "Failure: INPUT_FILE was not opened!" << endl;
        exit(EXIT_FAILURE) ;
    }

    // Open reaction network file
    ifstream reac_file(reactions_file);
    if(!reac_file)
    {
        cerr << "Failure: REACTIONS_FILE file was not opened!" << endl;
        exit(EXIT_FAILURE) ;
    }
    while( '#' == reac_file.peek() )
    {
        reac_file.ignore(2048,'\n');    // discard initial comments
    }
          
    // Open file defining the observables
    ifstream obs_file(observables_file);
    if(!obs_file)
    {
        cerr << "Failure: OBSERVABLES file was not opened!" << endl;
        exit(EXIT_FAILURE) ;
    }
    while( '#' == obs_file.peek() )
    {
        obs_file.ignore(2048,'\n');    // discard initial comments
    }
    

    in_file >> what_to_do;
    in_file >> max_interrupts;
    in_file >> init_random_seed;
    random_seed = init_random_seed;

    in_file >> rlx_time;
    in_file >> main_time;
    in_file >> meas_time;    
    in_file >> rlx_steps;
    in_file >> main_steps;
    in_file >> meas_steps;
    in_file >> no_trajec;

    in_file >> subv_size;

    in_file >> no_r_nucl;
    in_file >> no_z_nucl;

    no_nucl = no_z_nucl*no_r_nucl;

    in_file >> no_reac;
    in_file >> no_para;
    in_file >> no_spec;
    in_file >> no_obs;

    in_file >> no_neig;

    threshold = new double[no_spec];
    for(unsigned int spec = 0; spec < no_spec; spec++)
    {
        in_file >> threshold[spec];
    }        

    reac_parameter = new double[no_r_nucl*no_z_nucl*no_para];
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            for(unsigned int para = 0; para < no_para; para++)
            {
                in_file >> reac_parameter[para+z_nucl*no_para+r_nucl*no_z_nucl*no_para];
            }
        }
    }
    
    diff_rate_const = new double[no_r_nucl*no_z_nucl*no_spec];
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            for(unsigned int spec = 0; spec < no_spec; spec++)
            {
                in_file >> diff_rate_const[spec+z_nucl*no_spec+r_nucl*no_z_nucl*no_spec];
            }
        }
    }
    
    configuration = new double[no_r_nucl*no_z_nucl*no_spec];
    for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
    {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            for(unsigned int spec = 0; spec < no_spec; spec++)
            {
                in_file >> configuration[spec+z_nucl*no_spec+r_nucl*no_z_nucl*no_spec];
            }
        }
    }
    
    in_file >> output_all;
    in_file >> instant_outputs;
    
    // Read reaction network
    reaction = new reaction_class[no_reac];
    for(unsigned int r = 0;  r < no_reac; r++)
    {
        reac_file >> reaction[r].number;
        reac_file >> c; // the reaction type code
        reac_file >> i; // the file contains the parameter no. only
        reac_file >> reaction[r].partner[0];
        reac_file >> reaction[r].partner[1];
        reac_file >> reaction[r].product[0];
        reac_file >> reaction[r].product[1];
        reac_file >> dump; // discard the user comment
        
        // Set the right reaction rate
        reaction[r].rate = reac_parameter[i];
        
        // Parse the reaction type code
        switch(c)
        {
            case 'p':   // PRODUCTION
              
              reaction[r].type = PRODUCTION;
              assert(   reaction[r].partner[0] >  -1 && reaction[r].partner[1] == -1 &&
                        reaction[r].product[0] >  -1 && reaction[r].product[1] >  -1    );
              break;
              
            case 'd':   // DECAY
              
              reaction[r].type = DECAY;
              assert(   reaction[r].partner[0] >  -1 && reaction[r].partner[1] == -1 &&
                        reaction[r].product[0] == -1 && reaction[r].product[1] == -1    );
              break;
           
            case 'b':   // BINDING
              
              reaction[r].type = BINDING;
              assert(   reaction[r].partner[0] >  -1 && reaction[r].partner[1] >  -1 &&
                        reaction[r].product[0] >  -1 && reaction[r].product[1] == -1    );
              break;
              
            case 'u':   // UNBINDING
              
              reaction[r].type = UNBINDING;
              assert(   reaction[r].partner[0] >  -1 && reaction[r].partner[1] == -1 &&
                        reaction[r].product[0] >  -1 && reaction[r].product[1] >  -1    );
              break;
              
            case 'D':   // DIMERIZATION
              
              reaction[r].type = DIMERIZATION;
              assert(   reaction[r].partner[0] >  -1 && reaction[r].partner[1] >  -1 &&
                        reaction[r].product[0] >  -1 && reaction[r].product[1] == -1    );
              break;              
        }
        
    }
    
    // Read observables set
    obs = new observable[no_obs];
    int no_comp;
    for(unsigned int s = 0; s < no_obs; s++)
    {
        obs_file >> obs[s].number;
        obs_file >> obs[s].name;
        obs_file >> obs[s].threshold;
        
        obs_file >> no_comp;
        if(no_comp > 100){
          
          // truncate to max. no. of components
          no_comp = 100;
          obs_file.ignore(1024,'\n');
          cerr << "Observable " << s << " has more than 100 components. Only reading first 100 entries." << endl;
          
        }
        obs[s].no_components = no_comp;
        
        for(unsigned int c = 0; c < no_comp; c++)
        {
            obs_file >> obs[s].component[c];
        }
        
        obs[s].value = 0.0; // will be calculated elsewhere      
    }
       
    // Finalize
    in_file.close();
    reac_file.close();
    obs_file.close();

    return(EXIT_SUCCESS);
}

int run_parameter::read_temp_file()
{
    double tmp;
    MTRand::uint32 ltmp;
    
    ifstream in_file(temp_file);
    if(!in_file)
    {
        cerr << "Failure: TEMP_FILE was not opened!" << endl;
        exit(EXIT_FAILURE) ;
    }

    // define local random variable arrays
    obs_config = new random_var*[no_obs];
    avg_obs_config = new random_var*[no_obs];
    for(unsigned int s=0; s<no_obs; s++)	// for all measured species
    {
        
        obs_config[s]		= new random_var[no_r_nucl*no_z_nucl];
        avg_obs_config[s]	= new random_var[no_r_nucl*no_z_nucl];
    }

    obs_boundary_pos = new histo1D[2*no_obs];
    obs_boundary_no = new histo1D[2*no_obs];
    for(unsigned int s=0; s<2*no_obs; s++)	// for all measured species
    {

        obs_boundary_pos[s].init(1,no_z_nucl-1);
        obs_boundary_no[s].init(1,no_z_nucl-1);
    }

    // Start reading from temporary file
    in_file >> interrupts;
    in_file >> random_seed;

    if(interrupts){

        // Load mersenne twister random number generator state
        in_file >> generator_state_len;
        generator_state = new MTRand::uint32[generator_state_len];
        for(unsigned int i=0; i<generator_state_len; i++)	in_file >> generator_state[i];

        // Next: load measurement variables / histograms
        for(unsigned int s=0; s<no_obs; s++)	// for all (measured) species
        {
          
        // Load random variables
        for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
        {
            for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
            {
            in_file >> tmp;	obs_config[s][z_nucl+r_nucl*no_z_nucl].events = (int)tmp;
            in_file >> tmp;	obs_config[s][z_nucl+r_nucl*no_z_nucl].weight = (double)tmp;
            in_file >> tmp;	obs_config[s][z_nucl+r_nucl*no_z_nucl].max_moments = (int)tmp;
            for(unsigned int j = 0; j<obs_config[s][z_nucl+r_nucl*no_z_nucl].max_moments; j++)
            {
                in_file >> tmp;	obs_config[s][z_nucl+r_nucl*no_z_nucl].moments[j] = (double)tmp;
            }
                

            in_file >> tmp;	avg_obs_config[s][z_nucl].events = (int)tmp;
            in_file >> tmp;	avg_obs_config[s][z_nucl].weight = (double)tmp;
            in_file >> tmp;	avg_obs_config[s][z_nucl].max_moments = (int)tmp;
            for(unsigned int j = 0; j<obs_config[s][z_nucl+r_nucl*no_z_nucl].max_moments; j++)
            {
                in_file >> tmp;	avg_obs_config[s][z_nucl].moments[j] = (double)tmp;
            }
            }
        }

        // Load histogram counts and weights
        for(unsigned int i=1;i<no_z_nucl;i++){

            // anterior boundaries
            in_file >> tmp;	obs_boundary_pos[s].events_to(i,(int)tmp);
            in_file >> tmp;	obs_boundary_pos[s].weight_to(i,(double)tmp);
            in_file >> tmp;	obs_boundary_no[s].events_to(i,(int)tmp);
            in_file >> tmp;	obs_boundary_no[s].weight_to(i,(double)tmp);
            // posterior boundaries
            in_file >> tmp;    obs_boundary_pos[no_obs+s].events_to(i,(int)tmp);
            in_file >> tmp;    obs_boundary_pos[no_obs+s].weight_to(i,(double)tmp);
            in_file >> tmp;    obs_boundary_no[no_obs+s].events_to(i,(int)tmp);
            in_file >> tmp;    obs_boundary_no[no_obs+s].weight_to(i,(double)tmp);
        }
        
        }	// species loop (s)

        // Load configurations
        configuration = new double[no_r_nucl*no_z_nucl*no_spec];
        event_time = new double[no_r_nucl*no_z_nucl];
        
        for(unsigned int r_nucl = 0; r_nucl < no_r_nucl; r_nucl++)
        {
        for(unsigned int z_nucl = 0; z_nucl < no_z_nucl; z_nucl++)
        {
            for(unsigned int spec = 0; spec < no_spec; spec++)	// for all species
                in_file >> configuration[spec+z_nucl*no_spec+r_nucl*no_z_nucl*no_spec];

            in_file >> event_time[z_nucl+r_nucl*no_z_nucl];
        }
        }

        // Load binary tree with next event times
        in_file >> ac_run_step;
        in_file >> ac_run_time;
        in_file >> last_time;
        in_file >> next_time;
        if(meas_time > 0.0 || meas_steps > 0){   // read this line only if averaging is on
          
            in_file >> ac_meas_step;
            in_file >> bgn_meas;
            in_file >> end_meas;
        }
        else{
          
            ac_meas_step = 0;
            bgn_meas = 0.0;
            end_meas = 0.0;
        }
        in_file >> tree_size;
        tree_value = new double[tree_size];
        tree_index = new int[tree_size];
        for(unsigned int n = 0; n < tree_size; n++)
        {
            in_file >> tree_value[n];
            in_file >> tree_index[n];
        }
  
    }	// end if(interrupts)

    in_file.close();

    return(EXIT_SUCCESS);
}

int run_parameter::open_output(ofstream & fout,const char * filename, bool header)
{
    fout.open(filename);
    if (!fout)
    {
        cerr << "Failure: OUTPUT_FILE " << filename << " was not opened!" << endl;
        exit(EXIT_FAILURE) ;
    }
    if(header)  write_header(fout,filename);
    fout.setf(ios::scientific,ios::floatfield);
    fout.precision(8);
    return(EXIT_SUCCESS);
}

int run_parameter::open_output(ofstream & fout,const char * filename, bool header, int prec)
// overloaded function to open outstream with specified precision
{
    fout.open(filename);
    if (!fout)
    {
        cerr << "Failure: OUTPUT_FILE " << filename << " was not opened!" << endl;
        exit(EXIT_FAILURE) ;
    }
    if(header)  write_header(fout,filename);
    fout.setf(ios::scientific,ios::floatfield);
    fout.precision(prec);
    return(EXIT_SUCCESS);
}

int run_parameter::open_output(ofstream & fout,const char * filename, ios_base::openmode mode, int prec)
// overloaded function for the case that a file shall not be opened in the default "overwrite" mode
// but e.g. "append" mode (trajectory files)
{
    fout.open(filename, mode);
    if (!fout)
    {
        cerr << "Failure: OUTPUT_FILE " << filename << " was not opened!" << endl;
        exit(EXIT_FAILURE) ;
    }
    fout.setf(ios::scientific,ios::floatfield);
    fout.precision(prec);
    return(EXIT_SUCCESS);
}

int run_parameter::write_header(ofstream & fout,const char * filename)
{
    time_t now;
    time(&now);
    struct tm * now_loc;
    now_loc = gmtime(&now);
    char * datum;
    datum = asctime(now_loc);
    fout << "# * * * * * * * * * * * * * * * * * * * * * * *"  << endl;
    fout << "# * MC simulations started on * * * * * * * * *"  << endl;
    fout << "# * * * * * * * * * "                  <<  datum  << endl;
    fout << "# * * * * * * * * * * * * * * * * * * * * * * *"  << endl << endl;
    fout << "# * * * * * * Files connected * * * * * * * * *"  << endl;
    fout << "# * Input..............: " << input_file          << endl;
    fout << "# * ...................: " << reactions_file      << endl;
    fout << "# * ...................: " << observables_file    << endl;
    fout << "# * ...................: " << temp_file           << endl;
    fout << "# * Output.............: " << temp_file           << endl;
    fout << "# * ...................: " << traject_file        << endl;
    fout << "# * ...................: " << bnd_dst_file        << endl;
    fout << "# * ...................: " << average_file        << endl;
    fout << "# * ...................: " << new_input_file      << endl;
    fout << "# * * * * * * * * * * * * * ** * * * * * * * * *" << endl << endl;
    fout << "# * * * * * * Parameters * * * * * * * * * * * *" << endl;
    fout << "# * what to do......: "<< what_to_do       << endl;
    fout << "# * no. trajectories: "<< no_trajec        << endl;
    fout << "# * max_interrupts..: "<< max_interrupts   << endl;
    fout << "# * instant_outputs.: "<< instant_outputs  << endl;
    fout << "# * output all?.....: "<< output_all       << endl;
    fout << "# * init_random_seed: "<< init_random_seed << endl;
    fout << "# * random_seed.....: "<< random_seed      << endl;
    fout << "# * relax. time.....: "<< rlx_time         << endl;
    fout << "# * avg. time.......: "<< meas_time        << endl;
    fout << "# * main time.......: "<< main_time        << endl;
    fout << "# * relax. steps....: "<< rlx_steps        << endl;
    fout << "# * avg. steps......: "<< meas_steps       << endl;
    fout << "# * main steps......: "<< main_steps       << endl;
    fout << "# * subv. size .....: "<< subv_size        << endl;
    fout << "# * no. subv........: "<< no_nucl          << endl;
    fout << "# * no. r-subv......: "<< no_r_nucl        << endl;
    fout << "# * no. z-subv......: "<< no_z_nucl        << endl;
    fout << "# * no. reac........: "<< no_reac          << endl;
    fout << "# * no. para........: "<< no_para          << endl;
    fout << "# * no. spec........: "<< no_spec          << endl;
    fout << "# * no. obs.........: "<< no_obs           << endl;
    fout << "# * no. neig........: "<< no_neig          << endl;
    fout << "# * thresholds......: ";
    for(int obs = 0; obs < no_obs; obs++) fout << threshold[obs] << " ";
    fout << endl;
    fout << "# * Reaction rate constants.: ";
    for(unsigned int i = 0; i < no_para; i++) fout << reac_parameter[i] << " " ;
    fout << endl;
    fout << "# * Diffusion rate constants.: ";
    for(int i = 0; i < no_spec; i++) fout << diff_rate_const[i] << " " ;
    fout << endl;
    fout << "# * Start configuration......: ";
    for(int i = 0; i < no_spec; i++) fout << configuration[i] << " " ;
    fout << endl;

    fout << "# * * * * * * * * * * * * * * * * * * *  * * * *" << endl << endl;
    fout.flush();

    return(EXIT_SUCCESS);
}

int run_parameter::close_output(ofstream & fout)
{
    time_t now;
    time(&now);
    struct tm * now_loc;
    now_loc = gmtime(&now);
    char * datum;
    datum = asctime(now_loc);
    fout << endl;
    fout << "# * Finished MC simulations at " << datum;
    fout.flush();
    fout.close();
    return(EXIT_SUCCESS);
}

