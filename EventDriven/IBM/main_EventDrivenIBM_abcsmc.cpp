#include "AbcSmc.h"
#include <vector>
#include "EventDriven_Sim.hpp"
#include <cstdlib>

const gsl_rng* RNG = gsl_rng_alloc(gsl_rng_taus2);



 vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp) {
    
    using namespace std;
    //const int numSims  = 1000;
    const int N        = 1000;
    const double DEATH = 0.012;//assumes max age of 85 y.o.
    const double BETA  = (double)args[0];//can't set to zero otherwise recovery time is inf -- do we want to change this?
    const double maxRunTime = 200; //units in years
    const int waningImmunity = 1;
    const double propVirusinWater = (double)args[1];
    

    EventDriven_MassAction_Sim sim(N,BETA,DEATH,maxRunTime,waningImmunity,propVirusinWater);
    sim.randomizePopulation(1);
    //100 year burn in time before metrics are gathered
    sim.runSimulation();
    double ageSum = 0.0;
    for(unsigned int k = 0; k < sim.printAvgAgeAtFirstInfect().size(); k++){
        ageSum+=sim.printAvgAgeAtFirstInfect()[k];
    }
    vector<double> metrics(1);
    
    metrics[0] = ageSum/sim.printAvgAgeAtFirstInfect().size();
    
    return metrics;

}

void usage(){
    cout<<"\n\t Usage: Input model parameters <beta> <propVirusinWater>";
    exit(-1);
}

int main(int argc, char* argv[]) {
    
    if (not(argc == 3 or argc == 5 or argc == 6)) usage();
    
    vector<double> initialValues(argc-1);
    for(int i = 1; i<argc; i++){
        initialValues[i-1] = atof(argv[i]);
    }
    
    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = -1;
    
    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) {
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            buffer_size = atoi(argv[++i]);
        } else {
            usage();
            exit(101);
        }
    }
    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (process_db) {
        gsl_rng_set(RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id
        abc->process_database(RNG);
    }
    
    if (simulate_db) {
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }
    
    return 0;
}