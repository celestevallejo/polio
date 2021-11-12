#define varname(s) #s

#include <iostream>
#include <iomanip>
#include <math.h>
#include <random>
#include <tuple>
#include <array>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <assert.h>
#include <sstream>
#include <map>
#include <sys/time.h>
#include <algorithm>
#include <cmath>
#include <cctype>


using namespace std;

std::string data = "Abc";
string tolower(const string input) {
    string output(input);
    transform(input.begin(), input.end(), output.begin(), [](unsigned char c){ return std::tolower(c); });
    return output;
}

//string output_dir = "/home/tjhladish/work/polio-small-pop/output/";
//string output_dir ="/Users/Celeste/Desktop/multipatch_model/sim_results/";
//string output_dir ="/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/new_data_after_review/";
string output_dir ="./";
//string output_dir = "/home/vallejo.26/";
string ext = "_test.csv";
const string SEP = ",";
const string parameterFileDatabase = "/home/vallejo.26/parameterDatabase.csv";
//const string parameterFileDatabase = "/Users/Celeste/Desktop/multipatch_model/multiPatch_repo/parameterDatabase.csv";

uniform_real_distribution<> runif(0.0, 1.0);

#define MAKE_ENUM(VAR) VAR,
#define MAKE_STRINGS(VAR) #VAR,


#define STATE_TYPE(VAR) \
    VAR(S) \
    VAR(I1) \
    VAR(R) \
    VAR(P) \
    VAR(IR) \
    VAR(NUM_OF_STATE_TYPES)
enum StateType{ STATE_TYPE(MAKE_ENUM) };
const char* const state_as_string[] = { STATE_TYPE(MAKE_STRINGS) };
inline std::ostream& operator<<(std::ostream& out, const StateType value){ return out << state_as_string[value]; }


#define EVENT_TYPE(VAR) \
    VAR(MOVE_EVENT) \
    VAR(FIRST_INFECTION_EVENT) \
    VAR(REINFECTION_EVENT) \
    VAR(RECOVERY_FROM_FIRST_INFECTION_EVENT) \
    VAR(RECOVERY_FROM_REINFECTION_EVENT) \
    VAR(WANING_EVENT) \
    VAR(DEATH_EVENT) \
    VAR(REINTRODUCTION_EVENT) \
    VAR(NUM_OF_EVENT_TYPES)
enum EventType{ EVENT_TYPE(MAKE_ENUM) };
const char* const event_as_string[] = { EVENT_TYPE(MAKE_STRINGS) };
inline std::ostream& operator<<(std::ostream& out, const EventType value){ return out << event_as_string[value]; }


#define QUANT_OUTPUT_TYPE(VAR) \
    VAR(TIME_OUT) \
    VAR(TIME_BETWEEN_PCASES_OUT) \
    VAR(TRANSMISSION_INTERVAL_OUT) \
    VAR(EXTINCTION_INTERVAL_OUT) \
    VAR(EXTINCTION_TIME_OUT) \
    VAR(NUM_OF_QUANTOUTPUT_TYPES)
enum QuantOutputType{ QUANT_OUTPUT_TYPE(MAKE_ENUM) };
const char* const quant_output_as_string[] = { QUANT_OUTPUT_TYPE(MAKE_STRINGS) };
inline std::ostream& operator<<(std::ostream& out, const QuantOutputType value){ return out << quant_output_as_string[value]; }


struct Params{
    double recovery;
    double beta;
    double birth;
    double death;
    double kappa;
    double rho;
    vector<int> Population;
};

struct VillageEvent{ //keeps track of which events are occurring to which village
    EventType event_type;
    int village;
    double time;
};

//fast waning parameters:
//kappa = 0.4179
//rho = 0.2

//intermediate waning parameters:
//kappa = 0.6383
//rho = 0.04

//slow waning parameters:
//kappa = 0.8434
//rho = 0.02

const double KAPPA                 = 0.4179; //waning depth parameter
const double RHO                   = 0.2; //waning speed parameter

//other parameters
const vector<int> village_pop      = {100000, 100000};
const int NUM_OF_VILLAGES          = village_pop.size(); //total number of villages under consideration
const int numDaysToRecover         = 28;
const double RECOVERY              = 365/numDaysToRecover;    //recovery rate (/year)
const double BETA                  = 135;   //contact rate (individuals/year)
const double lifespan              = 50;
const double BIRTH                 = 1/lifespan; //birth rate (per year)
const double DEATH                 = 1/lifespan; //death rate (per year)
const double PIR                   = 0.005;            //type 1 paralysis rate (naturally occurring cases)
const double DET_RATE              = 1.0;
const double expectedTimeUntilMove = 1.0/4.0; //years; was 0
const double MOVE_RATE             = expectedTimeUntilMove > 0 ? 1/expectedTimeUntilMove : 0;
const double REINTRODUCTION_RATE   = 1.0/10.0; // was 0
const double MIN_BURN_IN           = 10; //years
const double OBS_PERIOD            = 50; //years
const double seasonalAmp           = 0.0;
const int NUM_OF_SIMS              = 4;
const bool RETRY_SIMS              = false;

void cerr_state(const vector<vector<int>> &state_data, const int village) {
    for (size_t state = 0; state < NUM_OF_STATE_TYPES; ++state) {
        cerr << right << setw(8) << state_data[state][village] << " " << (StateType) state;
    }
    cerr << endl;
}

bool choose_event(double &ran, const double p) {
    if (ran < p) {
        return true;
    } else {
        ran -= p;
        return false;
    }
}

unsigned int rand_nonuniform_uint(const vector<int> weights, mt19937& RNG) {
    const double totalWeight = accumulate(weights.begin(), weights.end(), 0.0);
    double ran = totalWeight*runif(RNG);

    for (unsigned int idx = 0; idx < weights.size(); ++idx) {
        if (choose_event(ran, weights[idx])) {
            return idx;
            break;
        }
    }
    return weights.size(); // indicates failure to choose
}

vector<double> calculate_village_rates(const vector<vector<int>> &state_data, const int village, const double time) {

    const int S_  = state_data[S][village];
    const int I1_ = state_data[I1][village];
    const int R_  = state_data[R][village];
    const int P_  = state_data[P][village];
    const int IR_ = state_data[IR][village];
    double seasonalBeta = BETA*(1 + seasonalAmp*sin(time/(2*M_PI)));
    double foi = seasonalBeta*(I1_ + KAPPA*IR_)/village_pop[village];

    //vector<vector<double>> event_rates(NUM_OF_VILLAGES, vector<double>(NUM_OF_EVENT_TYPES, 0.0));
    vector<double> local_event_rates(NUM_OF_EVENT_TYPES, 0.0);
    local_event_rates[FIRST_INFECTION_EVENT]                = S_*foi;                          // first infection event
    local_event_rates[REINFECTION_EVENT]                    = KAPPA*P_*foi;                    // reinfection event
    local_event_rates[RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1_;                    // first infected revovery event
    local_event_rates[RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*IR_;            // reinfected recovery event
    local_event_rates[WANING_EVENT]                         = RHO*R_;                          // waning event
    local_event_rates[DEATH_EVENT]                          = DEATH*village_pop[village];      // natural death
    local_event_rates[MOVE_EVENT]                           = MOVE_RATE*village_pop[village];  // rate of movement from village

    if (time < MIN_BURN_IN) {
        local_event_rates[REINTRODUCTION_EVENT]             = REINTRODUCTION_RATE*village_pop[village]; // rate of reintroduction into system
    } else{
        local_event_rates[REINTRODUCTION_EVENT]             = 0;
    }
    return local_event_rates;
}

VillageEvent sample_event(mt19937& RNG, double& totalRate, const vector<vector<int>> &state_data, const double time) {
    vector<vector<double>> event_rates(NUM_OF_VILLAGES, vector<double>(NUM_OF_EVENT_TYPES));
    totalRate = 0.0;
    VillageEvent ve;
    for (int vil = 0; vil < NUM_OF_VILLAGES; vil++) {
        event_rates[vil] = calculate_village_rates(state_data, vil, time);
        totalRate += accumulate(event_rates[vil].begin(),event_rates[vil].end(),0.0);
    }
    exponential_distribution<> rexp(totalRate);
    ve.time = time + rexp(RNG);

    double ran = totalRate * runif(RNG);

    for (int event = 0; event < NUM_OF_EVENT_TYPES; ++event) {
        for (int vil = 0; vil < NUM_OF_VILLAGES; vil++) {
            if (choose_event(ran, event_rates[vil][event])) {
                ve.event_type = (EventType) event;
                ve.village = vil;
                return ve;
            }
        }
    }
    cerr << "Error: No event sampled" << endl;
    exit(-100);
}

bool is_infected_state(StateType state) { return state == I1 or state == IR; }

bool all_zeroes(const vector<int> &data) { return all_of(data.begin(), data.end(), [](int i) {return i==0;}); }

bool zero_infections(const vector<vector<int>> &state_data) { return all_zeroes(state_data[I1]) and all_zeroes(state_data[IR]); }

bool zero_local_infections(const vector<vector<int>> &state_data, const unsigned int village) { return state_data[I1][village] == 0 and state_data[IR][village] == 0; }

void move_from_A_to_B (vector<vector<int>> &state_data, const vector<unsigned int> &vil_pair, const vector<StateType> &state_pair, const double obs_time, vector<vector<double>> &villageExtinctionIntervals, vector<double> &villageExtinctionTimes) {
    const StateType state_A = state_pair[0];
    const StateType state_B = state_pair[1];
    const unsigned int vil_A = vil_pair[0];
    const unsigned int vil_B = vil_pair[1];

    --state_data[state_A][vil_A];

    if (obs_time > 0) {
        // does village A now have 0 infected, but lost net 1 infected person? == extinction
        // does village A now have 0 infected, but will gain net 1 infected person? == end of extinction
        if (zero_local_infections(state_data, vil_A)) {
            if (is_infected_state(state_A) and not is_infected_state(state_B)) {
                // Village A just lost its last infected person, and isn't getting one back
                villageExtinctionTimes[vil_A] = obs_time;
            }

            if (not is_infected_state(state_A) and is_infected_state(state_B)) {
                // Village A has had no infections, but is getting one
                if (villageExtinctionTimes[vil_A] != numeric_limits<double>::max()) {
                    villageExtinctionIntervals[vil_A].push_back(obs_time - villageExtinctionTimes[vil_A]);
                }
            }
        }
    }
    ++state_data[state_A][vil_B];
}


inline void process_death_event(vector<vector<int>> &state_data, const int village, mt19937& RNG, vector<double> &transmissionIntervals, const double obs_time, double &reinfectTime, vector<double> &extinctionIntervals, double &lastPCase, vector<double> &villageExtinctionTimes) {
    vector<int> rates(NUM_OF_STATE_TYPES, 0);
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) { rates[state] = state_data[state][village]; }

    StateType source_state = (StateType) rand_nonuniform_uint(rates, RNG);

    ++state_data[S][village];
    --state_data[source_state][village];

    if (is_infected_state(source_state) and obs_time > 0) {
        double extinctionTime = obs_time;
        if (zero_local_infections(state_data, village)) {
            villageExtinctionTimes[village] = extinctionTime;
        }
        if (zero_infections(state_data)) {
            if (reinfectTime != numeric_limits<double>::max()) {
                transmissionIntervals.push_back(extinctionTime - reinfectTime);
            }
            if (lastPCase != numeric_limits<double>::max()) {
                extinctionIntervals.push_back(extinctionTime - lastPCase);
            }
        }
    }
}

StateType sample_migrant(vector<vector<int>> &state_data, const int village, mt19937& RNG) {
    vector<int> weights(NUM_OF_STATE_TYPES, 0);
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) { weights[state] = state_data[state][village]; }
    return (StateType) rand_nonuniform_uint(weights, RNG);
}

inline void process_movement_event(vector<vector<int>> &state_data, vector<unsigned int> vil_pair, mt19937& RNG, const double obs_time, vector<vector<double>> &villageExtinctionIntervals, vector<double> &villageExtinctionTimes) {
    const unsigned int vil_A = vil_pair[0];
    const unsigned int vil_B = vil_pair[1];

    if (vil_A != vil_B) {
        vector<int> weights(NUM_OF_STATE_TYPES, 0);

        // Sample states of people to move from A to B and vice versa
        const StateType state_A = sample_migrant(state_data, vil_pair[0], RNG);
        const StateType state_B = sample_migrant(state_data, vil_pair[1], RNG);

        if (state_A != state_B) {
            move_from_A_to_B(state_data, {vil_A, vil_B}, {state_A, state_B}, obs_time, villageExtinctionIntervals, villageExtinctionTimes);
            move_from_A_to_B(state_data, {vil_B, vil_A}, {state_B, state_A}, obs_time, villageExtinctionIntervals, villageExtinctionTimes);
        }
    }
}

void process_reintroduction_helper(const StateType from, const StateType to, const double obs_time, vector<vector<int>> &state_data, const int village, double &reinfectTime, vector<vector<double>> &villageExtinctionIntervals, vector<double> &villageExtinctionTimes) {
    if (obs_time > 0) {
        if (zero_local_infections(state_data, village) and villageExtinctionTimes[village] != numeric_limits<double>::max() ) {
            villageExtinctionIntervals[village].push_back(obs_time - villageExtinctionTimes[village]);
        }
        if (zero_infections(state_data)) {
            reinfectTime = obs_time;
        }
    }
    --state_data[from][village]; ++state_data[to][village];
}

inline void process_reintroduction_event(vector<vector<int>> &state_data, const int village, mt19937& RNG, const double obs_time, vector<double> &timeBetweenPCases, double &reinfectTime, double &lastPCase, vector<vector<double>> &villageExtinctionIntervals,vector<double> &villageExtinctionTimes) {
    //treat reintroduction events as exposure events
    vector<int> weights(NUM_OF_STATE_TYPES, 0);

    // Sample state of person to which exposure event will occur
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) { weights[state] = state_data[state][village]; }
    StateType state = (StateType) rand_nonuniform_uint(weights, RNG);

    //exposure events only change state of S and P compartments
    if (state == S) {
        process_reintroduction_helper(S, I1, obs_time, state_data, village, reinfectTime, villageExtinctionIntervals, villageExtinctionTimes);

        if (obs_time > 0) {//start keeping track after burn in
            if (runif(RNG) < PIR * DET_RATE) {
                if (lastPCase != numeric_limits<double>::max()) {
                    timeBetweenPCases.push_back(obs_time - lastPCase);
                }
                lastPCase = obs_time;
            }
        }
    } else if (state == P) {
        process_reintroduction_helper(P, IR, obs_time, state_data, village, reinfectTime, villageExtinctionIntervals, villageExtinctionTimes);
    }

}

void output_results(vector<vector<stringstream>> &output_streams, vector<stringstream> &outputVillageExtinctionInterval_stream, vector<stringstream> &outputQuant_streams) {
    string numInEachVil = "";
    for (int i = 0; i < NUM_OF_VILLAGES; i++) {
        cout<<"num in each vil"<<village_pop[i]<<"\n";
        numInEachVil += to_string(int(village_pop[i]));
    }

    string base_filename = numInEachVil + "reintRate_"+to_string(REINTRODUCTION_RATE) + "migRate_"+to_string(MOVE_RATE)+ext;
    string parameters = +"beta_"+to_string(int(BETA))+"detect_rate_"+to_string(float(DET_RATE))+"rho_"+to_string(float(RHO))+ "NUM_OF_VILLAGES_"+to_string(NUM_OF_VILLAGES) + "migRate_"+to_string(float(MOVE_RATE))+"burnIn_"+to_string(MIN_BURN_IN)+"obsTime_"+to_string(OBS_PERIOD) + "seasonalAmp_"+to_string(seasonalAmp);
    //read parameters into database file
    fstream params;
    params.open(parameterFileDatabase, fstream::app);
    params<<base_filename<<SEP<<parameters<<endl;
    params.close();

    vector<vector<string>> output_filenames(NUM_OF_STATE_TYPES, vector<string>(NUM_OF_VILLAGES));
    vector<string> outputVillageExtinctionInterval_filenames(NUM_OF_VILLAGES);

    vector<string> outputQuant_filenames(NUM_OF_QUANTOUTPUT_TYPES);

    outputQuant_filenames[TIME_BETWEEN_PCASES_OUT ] = output_dir + "time_between_pcases_"+base_filename;
    outputQuant_filenames[TRANSMISSION_INTERVAL_OUT ] = output_dir + "transmission_interval_"+base_filename;
    outputQuant_filenames[EXTINCTION_INTERVAL_OUT ] = output_dir + "extinction_interval_"+base_filename;
    outputQuant_filenames[TIME_OUT ] = output_dir + "TIME_"+base_filename;
    outputQuant_filenames[EXTINCTION_TIME_OUT ] = output_dir + "extTime_"+base_filename;


    for (int vil = 0; vil < NUM_OF_VILLAGES; vil++) {
        int patch = vil + 1;
        output_filenames[S][vil]  = output_dir + "S"+to_string(patch)+"_"+base_filename;
        output_filenames[I1][vil] = output_dir + "I1"+to_string(patch)+"_"+base_filename;
        output_filenames[R][vil]  = output_dir + "R"+to_string(patch)+"_"+base_filename;
        output_filenames[P][vil]  = output_dir + "P"+to_string(patch)+"_"+base_filename;
        output_filenames[IR][vil] = output_dir + "IR"+to_string(patch)+"_"+base_filename;
        outputVillageExtinctionInterval_filenames[vil] = output_dir + "extInts_"+to_string(patch)+"_"+base_filename;
    }

    for (int ot_idx = 0; ot_idx < NUM_OF_QUANTOUTPUT_TYPES; ++ot_idx) {
        const QuantOutputType ot = (QuantOutputType) ot_idx;
        //const OutputType ot = CIRCULATION_INTERVAL_OUT;
        ofstream ofs;
        ofs.open(outputQuant_filenames[ot]);
        ofs << outputQuant_streams[ot].rdbuf();
        ofs.close();
    }
    for (int vil = 0; vil < NUM_OF_VILLAGES; vil++) {
        vector<ofstream> ofs(NUM_OF_STATE_TYPES);
        ofstream ofExt;

        for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) {
            ofs[state].open(output_filenames[state][vil]);
            ofs[state] << output_streams[state][vil].rdbuf();
            ofs[state].close();
        }

        ofExt.open(outputVillageExtinctionInterval_filenames[vil]);
        ofExt << outputVillageExtinctionInterval_stream[vil].rdbuf();
        ofExt.close();
    }
}

void log_recovery_event(const vector<vector<int>> &state_data, const unsigned int village, const double obs_time, vector<double> &villageExtinctionTimes, vector<double> &transmissionIntervals, vector<double> &extinctionIntervals, const double reinfectTime, const double lastPCase) {
    if (obs_time > 0) {
        const double extinctionTime = obs_time;
        if (state_data[I1][village] + state_data[IR][village] == 0) {
            villageExtinctionTimes[village] = extinctionTime;
        }
        if (zero_infections(state_data)) {
            if (reinfectTime!=numeric_limits<double>::max()) {
                transmissionIntervals.push_back(extinctionTime - reinfectTime);
            }
            if (lastPCase!= numeric_limits<double>::max()) {
                extinctionIntervals.push_back(extinctionTime - lastPCase);
            }
        }
    }
}

void append_output(vector<vector<vector<int>>> &output, const vector<vector<int>> &state_data) {
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) {
        for (int vil = 0; vil < NUM_OF_VILLAGES; vil++) {
            output[state][vil].push_back(state_data[state][vil]);
        }
    }
}

void event_handler(const VillageEvent &ve, vector<vector<int>> &state_data, vector<double> &timeBetweenPCases, vector<double> &transmissionIntervals, vector<double> &extinctionIntervals, const double obs_time, double &lastPCase, vector<double> &villageExtinctionTimes, vector<vector<double>> &villageExtinctionIntervals, double &reinfectTime, mt19937& RNG) {
    const unsigned int v = ve.village;
    switch(ve.event_type) {
        case FIRST_INFECTION_EVENT:{
            --state_data[S][v]; ++state_data[I1][v];
            //check to see if paralytic case occurs
            if (obs_time > 0) {//start keeping track after burn in
                if (runif(RNG) < PIR*DET_RATE) {
                    if (lastPCase != numeric_limits<double>::max()) {
                        timeBetweenPCases.push_back(obs_time - lastPCase);
                    }
                    lastPCase = obs_time;
                }
            }
        }
            break;
        case REINFECTION_EVENT:
            --state_data[P][v]; ++state_data[IR][v];
            break;
        case RECOVERY_FROM_FIRST_INFECTION_EVENT:
            --state_data[I1][v]; ++state_data[R][v];
            log_recovery_event(state_data, v, obs_time, villageExtinctionTimes, transmissionIntervals, extinctionIntervals, reinfectTime, lastPCase);
            break;
        case RECOVERY_FROM_REINFECTION_EVENT:
            --state_data[IR][v]; ++state_data[R][v];
            log_recovery_event(state_data, v, obs_time, villageExtinctionTimes, transmissionIntervals, extinctionIntervals, reinfectTime, lastPCase);
            break;
        case WANING_EVENT: --state_data[R][v]; ++state_data[P][v]; break;
        case DEATH_EVENT:
            process_death_event(state_data, v, RNG, transmissionIntervals, obs_time, reinfectTime, extinctionIntervals, lastPCase, villageExtinctionTimes);
            break;
        case MOVE_EVENT: {
            uint other = rand_nonuniform_uint(village_pop, RNG);
            process_movement_event(state_data, {v, other}, RNG, obs_time, villageExtinctionIntervals, villageExtinctionTimes); }
            break;
        case REINTRODUCTION_EVENT: {
            //treat reintroduction as exposure event that is density dependent
            uint village = rand_nonuniform_uint(village_pop,RNG); // TODO - this is a second declaration of a village.  what should this be called?  other?
            process_reintroduction_event(state_data, village, RNG, obs_time, timeBetweenPCases, reinfectTime, lastPCase, villageExtinctionIntervals, villageExtinctionTimes); }
            break;
        default:
            cerr << "ERROR: Unsupported event type" << endl;
            break;
    }
}


int main() {
    //vectors for holding information to be output
    vector<vector<vector<int>>> output(NUM_OF_STATE_TYPES, vector<vector<int>>(NUM_OF_VILLAGES));

    //keep track of all interparalytic case intervals for SC statistic calculation
    //interparalytic case interval for all villages
    vector<double> timeBetweenPCases;
    //keep track of all extinction intervals for SC statistic calculation
    //extinction interval for all villages (extinct in all villages at once)
    vector<double> extinctionIntervals;
    //keep track of time between reignition of infection and extinction
    vector<double> transmissionIntervals;
    //used to determine the length of a transmission interval
    double reinfectTime;
    //used to determine the length of an extinction interval
    //keeps track of time of last paralytic case
    double lastPCase;
    vector<double> timeVec;

    //keeps track of time interval of extinction for each village
    vector<vector<double>> villageExtinctionIntervals(NUM_OF_VILLAGES);
    //keep track of reinfection time in each village
    vector<double> villageExtinctionTimes(NUM_OF_VILLAGES);

    //vectors for outputing information
    vector<vector<stringstream>> output_streams;
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) {
        output_streams.push_back(vector<stringstream>(NUM_OF_VILLAGES));
    }
    //vector<vector<stringstream>> output_streams(NUM_OF_STATE_TYPES, vector<stringstream> (NUM_OF_VILLAGES)); // for some reason, compiler rejects this!
    vector<stringstream> outputQuant_stream(NUM_OF_QUANTOUTPUT_TYPES);

    vector<stringstream> outputVillageExtinctionInterval_stream(NUM_OF_VILLAGES);

    //initialize size of vector of vector of compartments
    //vector<vector<double>> compartments(NUM_OF_VILLAGES);

    //initialize size of vector of vector of initial values
    //vector<vector<int>> initialValues(NUM_OF_VILLAGES);

    //initialize size of vectors for individual compartments
    vector<vector<int>> state_data(NUM_OF_STATE_TYPES, vector<int>(NUM_OF_VILLAGES, 0.0));

    random_device rd;                       // generates a random real number for the seed
    //unsigned long int seed = rd();
    uint seed = 2186064846;
    cerr << "seed: " << seed << endl;
    mt19937 RNG(seed);                      // random number generator
    //mt19937 RNG(rd());                      // random number generator

    //find expected compartment size for each village
    /*for (int vil = 0; vil < NUM_OF_VILLAGES; vil++) {
        compartments[vil] = initialize_compartment(vil);
    }*/
    const int EPS_RES = 100; // resolution of endemic potential statistic, in divisions per year
    const int EPS_MAX = 50;  // circulation interval considered for EPS calculation, in years
    vector<int> eps_circ_ivls(EPS_MAX*EPS_RES, 0); // 50 yrs divided into 100 bins each
    vector<int> eps_intercase_ivls(EPS_MAX*EPS_RES, 0); // 50 yrs divided into 100 bins each

    //The Simulation
    for (int i = 0; i < NUM_OF_SIMS; i++) {
        double time         = 0.0;       // "absolute" time, relative to start of simulation
        double previousTime = 0.0;
        double burn_in      = -DBL_MAX; // don't know yet exactly how long burn-in will be.  burn_in ends with first event after MIN_BURN_IN
        double obs_time     = -DBL_MAX; // time measured since burn-in
        lastPCase = numeric_limits<double>::max(); //initialize to largest double
        reinfectTime = numeric_limits<double>::max(); //initialize to largest double

        for (int vil=0; vil<NUM_OF_VILLAGES; vil++) {
            villageExtinctionTimes[vil] = numeric_limits<double>::max(); //initialize to largest double
        }

        for (int vil = 0; vil < NUM_OF_VILLAGES; vil++) {
            //set initial values for each village using multinomial dist
            //initialValues[vil] = multinomial_Compartments(compartments[vil].size(),compartments[vil],vil,RNG());
            state_data[S][vil]        = (int) (0.99*village_pop[vil]);   //naive susceptible (no previous contact w/virus, moves into I1)
            state_data[I1][vil]       = village_pop[vil] - state_data[S][vil];  //first infected (only time paralytic case can occur, recovers into R)
            state_data[R][vil]        = 0;   //recovered (fully immune, wanes into P)
            state_data[P][vil]        = 0;   //partially susceptible (moves into IR)
            state_data[IR][vil]       = 0;  //reinfected (recovers into R)
        }

        long long int day_ct = -1;
        const double day = 1.0/365.0;

        for (int j = 0; j < 1e10; j++) {
            while (time/day > day_ct) {
                if (day_ct % 100 == 0) {
                    cerr << "step, time, day: " << right << setw(10) << j << ", " << setw(9) << setprecision(3) << time << ", " << setw(6) << day_ct << " | ";
                    cerr_state(state_data, 0);
                }
                ++day_ct;
            }
//            if (j % 10000 == 0) { cerr << "step, time, time/day, day_ct: " << j << ", " << time << ", " << time/day << ", " << day_ct << endl; }

            //keep track of previous time
            previousTime = time; // TH - why?
            double totalRate = 0.0;
            VillageEvent ve = sample_event(RNG, totalRate, state_data, time); // This is where 'time' gets updated
            event_handler(ve, state_data, timeBetweenPCases, transmissionIntervals, extinctionIntervals, obs_time, lastPCase, villageExtinctionTimes, villageExtinctionIntervals, reinfectTime, RNG);
            time = ve.time;

            //start collecting data after burn in
            if (time > MIN_BURN_IN) {
                if (burn_in < 0) {
                    burn_in = time;

                    //determine state of system at beginning of burn in
                    //if there are infected individuals start clock for transmission interval
                    //if not, wait until there are infected individuals
                    for (int vil = 0; vil < NUM_OF_VILLAGES; vil++) {
                        if (zero_local_infections(state_data, vil)) {
                            villageExtinctionTimes[vil] = 0; // t = 0 relative to observation period
                        }
                    }
                    if (not zero_infections(state_data)) {
                        reinfectTime = 0; // t = 0 relative to observation period
                    }
                }
                obs_time = time - burn_in; // here obs_time will always be >= 0.0

                double truncatePrevious = int(previousTime*10)/10.0; //truncate to 1 decimal place
                double truncateTime = int(time*10)/10.0;
                if (truncateTime > truncatePrevious) {
                    timeVec.push_back(obs_time);
                    append_output(output, state_data);
                }

                //stopping condition
                if (zero_infections(state_data) or (obs_time >= OBS_PERIOD)) {
                    //determine if there were no p cases in any villages

                    if (timeBetweenPCases.size() == 0) {
                        cerr << "Run failed, no cases observed\n\n";
                        if (RETRY_SIMS) { i--; } // CAN RESULT IN AN INFINITE LOOP!!!
                        //want at least two cases
                        extinctionIntervals.clear();
                    }

                    timeVec.push_back(obs_time);
                    append_output(output, state_data);
                    for (unsigned int vil = 0; vil < (unsigned) NUM_OF_VILLAGES; vil++) {
                        if (zero_local_infections(state_data, vil) and villageExtinctionTimes[vil] != numeric_limits<double>::max()) {
                            villageExtinctionIntervals[vil].push_back(obs_time - villageExtinctionTimes[vil]);
                        }
                    }
                    //output extinction time
                    outputQuant_stream[EXTINCTION_TIME_OUT] << obs_time <<endl;
                    //output time increments
                    for (unsigned int tLength = 0; tLength < timeVec.size(); tLength++) {
                        outputQuant_stream[TIME_OUT]<<timeVec[tLength]<<SEP;
                    }
                    outputQuant_stream[TIME_OUT]<<endl;

                    for (unsigned int vil = 0; vil < (unsigned) NUM_OF_VILLAGES; vil++) {
                        //output population counts
                        for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) {
                            for (unsigned int pop = 0; pop < output[state][vil].size(); pop++) { // What is pop? - tjh
                                output_streams[state][vil] << output[state][vil][pop] << SEP;
                            }
                            output_streams[state][vil] << endl;
                        }
                    }

                    //output extinction intervals for each village
                    for (unsigned int vil = 0; vil < (unsigned) NUM_OF_VILLAGES; vil++) {
                        unsigned int outputVecSize = villageExtinctionIntervals[vil].size();
                        for (unsigned int interval = 0; interval < outputVecSize; interval++) {
                            outputVillageExtinctionInterval_stream[vil]<<villageExtinctionIntervals[vil][interval]<<SEP;
                        }
                        outputVillageExtinctionInterval_stream[vil]<<endl;
                    }

                    //output intervals representing time between paralytic cases
                    for (unsigned int interval = 0; interval < timeBetweenPCases.size(); interval++) {
                        outputQuant_stream[TIME_BETWEEN_PCASES_OUT] << timeBetweenPCases[interval];
                        if (interval < timeBetweenPCases.size() - 1) {
                            outputQuant_stream[TIME_BETWEEN_PCASES_OUT]<< SEP;
                        }
                    }
                    outputQuant_stream[TIME_BETWEEN_PCASES_OUT] << endl;

                    //output transmission intervals
                    for (unsigned int interval = 0; interval < transmissionIntervals.size(); interval++) {
                        outputQuant_stream[TRANSMISSION_INTERVAL_OUT] << transmissionIntervals[interval];
                        if (interval < transmissionIntervals.size() - 1) {
                            outputQuant_stream[TRANSMISSION_INTERVAL_OUT] << SEP;
                        }
                    }
                    outputQuant_stream[TRANSMISSION_INTERVAL_OUT] << endl;

                    //output extinction intervals
                    for (unsigned int interval = 0; interval < extinctionIntervals.size(); interval++) {
                        assert(extinctionIntervals[interval] >= 0);
                        outputQuant_stream[EXTINCTION_INTERVAL_OUT] << extinctionIntervals[interval];
                        if (interval < extinctionIntervals.size() - 1) {
                            outputQuant_stream[EXTINCTION_INTERVAL_OUT]<< SEP;
                        }
                    }
                    outputQuant_stream[EXTINCTION_INTERVAL_OUT] << endl;

                    //clear vectors when done
                    timeVec.clear();
                    timeBetweenPCases.clear();
                    transmissionIntervals.clear();
                    extinctionIntervals.clear();
                    villageExtinctionTimes.clear();
                    output = vector<vector<vector<int>>>(NUM_OF_STATE_TYPES, vector<vector<int>>(NUM_OF_VILLAGES));
                    villageExtinctionIntervals = vector<vector<double>> (NUM_OF_VILLAGES);
                    break;
                }
            }
        }
    }
    assert(eps_intercase_ivls.size() == eps_circ_ivls.size());
    output_results(output_streams, outputVillageExtinctionInterval_stream, outputQuant_stream);
    return 0;
}
