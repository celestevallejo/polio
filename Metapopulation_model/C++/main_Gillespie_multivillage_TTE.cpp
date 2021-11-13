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

string tolower(const string input) {
    string output(input);
    transform(input.begin(), input.end(), output.begin(), [](unsigned char c){ return std::tolower(c); });
    return output;
}

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
    VAR(TIME) \
    VAR(EXTINCTION_TIME) \
    VAR(NUM_OF_QUANT_OUTPUT_TYPES)
enum QuantOutputType{ QUANT_OUTPUT_TYPE(MAKE_ENUM) };
const char* const quant_output_as_string[] = { QUANT_OUTPUT_TYPE(MAKE_STRINGS) };
inline std::ostream& operator<<(std::ostream& out, const QuantOutputType value){ return out << quant_output_as_string[value]; }


#define INTERVAL_TYPE(VAR) \
    VAR(PCASE_INTERVAL) \
    VAR(TRANSMISSION_INTERVAL) \
    VAR(EXTINCTION_INTERVAL) \
    VAR(NUM_OF_INTERVAL_TYPES)
enum IntervalType{ INTERVAL_TYPE(MAKE_ENUM) };
const char* const interval_as_string[] = { INTERVAL_TYPE(MAKE_STRINGS) };
inline std::ostream& operator<<(std::ostream& out, const IntervalType value){ return out << interval_as_string[value]; }


struct Params{
    double recovery;
    double beta;
    double birth;
    double death;
    double kappa;
    double rho;
    vector<int> Population;
};


struct VillageEvent{
    EventType event_type;
    int village;
    double time;
};

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

// fast waning parameters:
//kappa = 0.4179
//rho = 0.2

// intermediate waning parameters:
//kappa = 0.6383
//rho = 0.04

// slow waning parameters:
//kappa = 0.8434
//rho = 0.02

const double KAPPA                 = 0.4179;                // waning depth parameter
const double RHO                   = 0.2;                   // waning speed parameter

// other parameters
const vector<int> village_pop      = {100000, 100000};
const size_t NUM_OF_VILLAGES          = village_pop.size(); // total number of villages under consideration
const int numDaysToRecover         = 28;
const double RECOVERY              = 365/numDaysToRecover;  // recovery rate (/year)
const double BETA                  = 135;                   // contact rate (individuals/year)
const double lifespan              = 50;
const double BIRTH                 = 1/lifespan;            // birth rate (per year)
const double DEATH                 = 1/lifespan;            // death rate (per year)
const double PIR                   = 0.005;                 // type 1 paralysis rate (naturally occurring cases)
const double DET_RATE              = 1.0;
const double expectedTimeUntilMove = 1.0/4.0;               // years; was 0
const double MOVE_RATE             = expectedTimeUntilMove > 0 ? 1/expectedTimeUntilMove : 0;
const double REINTRODUCTION_RATE   = 1.0/10.0;              // was 0
const double MIN_BURN_IN           = 10;                    // years
const double OBS_PERIOD            = 50;                    // years
const double SEASONALITY           = 0.0;
const int NUM_OF_SIMS              = 4;
const bool RETRY_SIMS              = false;

bool extinction_observed(double val)     { return val != numeric_limits<double>::max(); }
bool reinfection_observed(double val)    { return val != numeric_limits<double>::max(); }
bool paralytic_case_observed(double val) { return val != numeric_limits<double>::max(); }
bool infection_observed(double val)      { return val != numeric_limits<double>::max(); }

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
    double seasonalBeta = BETA*(1 + SEASONALITY*sin(time/(2*M_PI)));
    double foi = seasonalBeta*(I1_ + KAPPA*IR_)/village_pop[village];

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
    for (size_t vil = 0; vil < NUM_OF_VILLAGES; ++vil) {
        event_rates[vil] = calculate_village_rates(state_data, vil, time);
        totalRate += accumulate(event_rates[vil].begin(),event_rates[vil].end(),0.0);
    }
    exponential_distribution<> rexp(totalRate);
    ve.time = time + rexp(RNG);

    double ran = totalRate * runif(RNG);

    for (size_t event = 0; event < NUM_OF_EVENT_TYPES; ++event) {
        for (size_t vil = 0; vil < NUM_OF_VILLAGES; ++vil) {
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
                if (extinction_observed(villageExtinctionTimes[vil_A])) {
                    villageExtinctionIntervals[vil_A].push_back(obs_time - villageExtinctionTimes[vil_A]);
                }
            }
        }
    }
    ++state_data[state_A][vil_B];
}

void log_event(const vector<vector<int>> &state_data, const unsigned int village, const double obs_time, vector<vector<double>> &intervals, vector<double> &villageExtinctionTimes, const double reinfectTime, const double lastPCase) {
    if (obs_time > 0) {
        const double extinctionTime = obs_time;
        if (zero_local_infections(state_data, village)) { villageExtinctionTimes[village] = extinctionTime; }

        if (zero_infections(state_data)) {
            if (reinfection_observed(reinfectTime)) { intervals[TRANSMISSION_INTERVAL].push_back(extinctionTime - reinfectTime); }
            if (paralytic_case_observed(lastPCase)) { intervals[EXTINCTION_INTERVAL].push_back(extinctionTime - lastPCase); }
        }
    }
}

inline void process_death_event(vector<vector<int>> &state_data, const int village, const double obs_time, vector<vector<double>> &intervals, double &reinfectTime, double &lastPCase, vector<double> &villageExtinctionTimes, mt19937& RNG) {
    vector<int> rates(NUM_OF_STATE_TYPES, 0);
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) { rates[state] = state_data[state][village]; }

    StateType source_state = (StateType) rand_nonuniform_uint(rates, RNG);

    ++state_data[S][village];
    --state_data[source_state][village];

    if (is_infected_state(source_state)) {
        log_event(state_data, village, obs_time, intervals, villageExtinctionTimes, reinfectTime, lastPCase);
    }
}

StateType sample_migrant(vector<vector<int>> &state_data, const int village, mt19937& RNG) {
    vector<int> weights(NUM_OF_STATE_TYPES, 0);
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) { weights[state] = state_data[state][village]; }
    return (StateType) rand_nonuniform_uint(weights, RNG);
}


inline void process_movement_event(vector<vector<int>> &state_data, vector<unsigned int> vil_pair, const double obs_time, vector<double> &villageExtinctionTimes, vector<vector<double>> &villageExtinctionIntervals, mt19937& RNG) {
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
        if (zero_local_infections(state_data, village) and extinction_observed(villageExtinctionTimes[village])) {
            villageExtinctionIntervals[village].push_back(obs_time - villageExtinctionTimes[village]);
        }
        if (zero_infections(state_data)) {
            reinfectTime = obs_time;
        }
    }
    --state_data[from][village]; ++state_data[to][village];
}

inline void process_reintroduction_event(vector<vector<int>> &state_data, const int village, const double obs_time, vector<double> &pcaseInterval, double &reinfectTime, double &lastPCase, vector<double> &villageExtinctionTimes, vector<vector<double>> &villageExtinctionIntervals, mt19937& RNG) {
    // treat reintroduction events as exposure events
    vector<int> weights(NUM_OF_STATE_TYPES, 0);

    // Sample state of person to which exposure event will occur
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) { weights[state] = state_data[state][village]; }
    StateType state = (StateType) rand_nonuniform_uint(weights, RNG);

    // exposure events only change state of S and P compartments
    if (state == S) {
        process_reintroduction_helper(S, I1, obs_time, state_data, village, reinfectTime, villageExtinctionIntervals, villageExtinctionTimes);

        if (obs_time > 0) {// start keeping track after burn in
            if (runif(RNG) < PIR * DET_RATE) {
                if (paralytic_case_observed(lastPCase)) {
                    pcaseInterval.push_back(obs_time - lastPCase);
                }
                lastPCase = obs_time;
            }
        }
    } else if (state == P) {
        process_reintroduction_helper(P, IR, obs_time, state_data, village, reinfectTime, villageExtinctionIntervals, villageExtinctionTimes);
    }

}

void output_results(vector<vector<stringstream>> &output_streams, vector<stringstream> &outputVillageExtinctionInterval_stream, vector<stringstream> &outputQuant_streams, vector<stringstream> &output_interval_streams) {
    string numInEachVil = "";
    for (size_t i = 0; i < NUM_OF_VILLAGES; ++i) {
        cout << "num in each vil" << village_pop[i] << "\n";
        numInEachVil += to_string(int(village_pop[i]));
    }

    string base_filename = numInEachVil + "reintRate_" + to_string(REINTRODUCTION_RATE) + "migRate_" + to_string(MOVE_RATE)+ext;
    string parameters = "beta_" + to_string(int(BETA)) + "detect_rate_" + to_string(float(DET_RATE)) + "rho_" + to_string(float(RHO)) + 
                        "NUM_OF_VILLAGES_" + to_string(NUM_OF_VILLAGES) + "migRate_" + to_string(float(MOVE_RATE)) + "burnIn_" + to_string(MIN_BURN_IN) +
                        "obsTime_" + to_string(OBS_PERIOD) + "seasonality_" + to_string(SEASONALITY);
    // read parameters into database file
    fstream params;
    params.open(parameterFileDatabase, fstream::app);
    params<<base_filename<<SEP<<parameters<<endl;
    params.close();

    vector<vector<string>> output_filenames(NUM_OF_STATE_TYPES, vector<string>(NUM_OF_VILLAGES));
    vector<string> outputVillageExtinctionInterval_filenames(NUM_OF_VILLAGES);

    vector<string> outputQuant_filenames(NUM_OF_QUANT_OUTPUT_TYPES);
    for (size_t output_type = 0; output_type < outputQuant_filenames.size(); ++output_type) {
        outputQuant_filenames[output_type] = output_dir + quant_output_as_string[output_type] + base_filename;
    }

    vector<string> output_interval_filenames(NUM_OF_INTERVAL_TYPES);
    for (size_t output_type = 0; output_type < output_interval_filenames.size(); ++output_type) {
        output_interval_filenames[output_type] = output_dir + interval_as_string[output_type] + base_filename;
    }

    for (size_t vil = 0; vil < NUM_OF_VILLAGES; ++vil) {
        for (size_t state_type = 0; state_type < NUM_OF_STATE_TYPES; ++state_type) {
            output_filenames[state_type][vil] = output_dir + state_as_string[state_type] + to_string(vil + 1) + "_" + base_filename;
        }
        outputVillageExtinctionInterval_filenames[vil] = output_dir + "extInts_" + to_string(vil + 1) + "_" + base_filename;
    }

    ofstream ofs;
    for (int ot = 0; ot < NUM_OF_QUANT_OUTPUT_TYPES; ++ot) {
        ofs.open(outputQuant_filenames[ot]);
        ofs << outputQuant_streams[ot].rdbuf();
        ofs.close();
    }

    for (int ot = 0; ot < NUM_OF_INTERVAL_TYPES; ++ot) {
        ofs.open(output_interval_filenames[ot]);
        ofs << output_interval_streams[ot].rdbuf();
        ofs.close();
    }

    for (size_t vil = 0; vil < NUM_OF_VILLAGES; ++vil) {
        for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) {
            ofs.open(output_filenames[state][vil]);
            ofs << output_streams[state][vil].rdbuf();
            ofs.close();
        }

        ofs.open(outputVillageExtinctionInterval_filenames[vil]);
        ofs << outputVillageExtinctionInterval_stream[vil].rdbuf();
        ofs.close();
    }
}

void append_output(vector<vector<vector<int>>> &output, const vector<vector<int>> &state_data) {
    for (size_t state = 0; state < NUM_OF_STATE_TYPES; ++state) {
        for (size_t vil = 0; vil < NUM_OF_VILLAGES; ++vil) {
            output[state][vil].push_back(state_data[state][vil]);
        }
    }
}

void event_handler(const VillageEvent &ve, vector<vector<int>> &state_data, vector<vector<double>> &intervals, const double obs_time, double &lastPCase, vector<double> &villageExtinctionTimes, vector<vector<double>> &villageExtinctionIntervals, double &reinfectTime, mt19937& RNG) {
    const unsigned int v = ve.village;
    switch(ve.event_type) {
        case FIRST_INFECTION_EVENT:{
            --state_data[S][v]; ++state_data[I1][v];
            // check to see if paralytic case occurs
            if (obs_time > 0) {// start keeping track after burn in
                if (runif(RNG) < PIR*DET_RATE) {
                    if (paralytic_case_observed(lastPCase)) {
                        intervals[PCASE_INTERVAL].push_back(obs_time - lastPCase);
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
            log_event(state_data, v, obs_time, intervals, villageExtinctionTimes, reinfectTime, lastPCase);
            break;
        case RECOVERY_FROM_REINFECTION_EVENT:
            --state_data[IR][v]; ++state_data[R][v];
            log_event(state_data, v, obs_time, intervals, villageExtinctionTimes, reinfectTime, lastPCase);
            break;
        case WANING_EVENT: --state_data[R][v]; ++state_data[P][v]; break;
        case DEATH_EVENT:
            process_death_event(state_data, v, obs_time, intervals, reinfectTime, lastPCase, villageExtinctionTimes, RNG);
            break;
        case MOVE_EVENT: {
            uint other = rand_nonuniform_uint(village_pop, RNG);
            process_movement_event(state_data, {v, other}, obs_time, villageExtinctionTimes, villageExtinctionIntervals, RNG); }
            break;
        case REINTRODUCTION_EVENT: {
            // treat reintroduction as exposure event that is density dependent
            uint village = rand_nonuniform_uint(village_pop, RNG); // TODO - this is a second declaration of a village.  what should this be called?  other?
            process_reintroduction_event(state_data, village, obs_time, intervals[PCASE_INTERVAL], reinfectTime, lastPCase, villageExtinctionTimes, villageExtinctionIntervals, RNG); }
            break;
        default:
            cerr << "ERROR: Unsupported event type" << endl;
            break;
    }
}


int main() {
    // TODO -- review what data is being used and for what.  we may not need everything we're outputting
    // vectors for holding information to be output
    vector<vector<vector<int>>> output(NUM_OF_STATE_TYPES, vector<vector<int>>(NUM_OF_VILLAGES));

    // interparalytic case intervals,  transmission intervals, and extinction intervals for overall population; for SC statistic calculation
    vector<vector<double>> intervals(NUM_OF_INTERVAL_TYPES);

    double reinfectTime;    // time between reintroduction and extinction, used to determine the length of a transmission interval
    double lastPCase;       // used to determine the length of an extinction interval
    vector<double> timeVec;

    // keeps track of time interval of extinction for each village
    vector<vector<double>> villageExtinctionIntervals(NUM_OF_VILLAGES);
    // keep track of reinfection time in each village
    vector<double> villageExtinctionTimes(NUM_OF_VILLAGES);

    // vectors for outputing information
    // vector<vector<stringstream>> output_streams(NUM_OF_STATE_TYPES, vector<stringstream> (NUM_OF_VILLAGES)); // for some reason, compiler rejects this!
    vector<vector<stringstream>> output_streams;
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) { output_streams.push_back(vector<stringstream>(NUM_OF_VILLAGES)); }

    vector<stringstream> outputQuant_stream(NUM_OF_QUANT_OUTPUT_TYPES);
    vector<stringstream> output_interval_streams(NUM_OF_INTERVAL_TYPES);
    vector<stringstream> outputVillageExtinctionInterval_stream(NUM_OF_VILLAGES);

    // initialize size of vectors for individual compartments
    vector<vector<int>> state_data(NUM_OF_STATE_TYPES, vector<int>(NUM_OF_VILLAGES, 0.0));

    random_device rd;                                   // generates a random real number for the seed
    // unsigned long int seed = rd();
    uint seed = 2186064846;
    cerr << "seed: " << seed << endl;
    mt19937 RNG(seed);                                  // random number generator

    const int EPS_RES = 100;                            // resolution of endemic potential statistic, in divisions per year
    const int EPS_MAX = 50;                             // circulation interval considered for EPS calculation, in years
    vector<int> eps_circ_ivls(EPS_MAX*EPS_RES, 0);      // 50 yrs divided into 100 bins each
    vector<int> eps_intercase_ivls(EPS_MAX*EPS_RES, 0); // 50 yrs divided into 100 bins each

    // The Simulation
    for (int i = 0; i < NUM_OF_SIMS; ++i) {
        double time         = 0.0;                      // "absolute" time, relative to start of simulation
        double previousTime = 0.0;
        double burn_in      = -DBL_MAX;                 // don't know yet exactly how long burn-in will be.  burn_in ends with first event after MIN_BURN_IN
        double obs_time     = -DBL_MAX;                 // time measured since burn-in
        lastPCase           = DBL_MAX;
        reinfectTime        = DBL_MAX;
        villageExtinctionTimes = vector<double>(NUM_OF_VILLAGES, DBL_MAX);
        
        for (size_t vil = 0; vil < NUM_OF_VILLAGES; ++vil) {
            // set initial values for each village using multinomial dist
            //initialValues[vil] = multinomial_Compartments(compartments[vil].size(),compartments[vil],vil,RNG());
            state_data[S][vil]        = (int) (0.99*village_pop[vil]);         // naive susceptible (no previous contact w/virus, moves into I1)
            state_data[I1][vil]       = village_pop[vil] - state_data[S][vil]; // first infected (only time paralytic case can occur, recovers into R)
            state_data[R][vil]        = 0;                                     // recovered (fully immune, wanes into P)
            state_data[P][vil]        = 0;                                     // partially susceptible (moves into IR)
            state_data[IR][vil]       = 0;                                     // reinfected (recovers into R)
        }

        long long int day_ct = -1;
        const double day = 1.0/365.0;

        for (int j = 0; j < 1e10; ++j) {
            while (time/day > day_ct) {
                if (day_ct % 100 == 0) {
                    cerr << "step, time, day: " << right << setw(10) << j << ", " << setw(9) << setprecision(3) << time << ", " << setw(6) << day_ct << " | ";
                    cerr_state(state_data, 0);
                }
                ++day_ct;
            }
//            if (j % 10000 == 0) { cerr << "step, time, time/day, day_ct: " << j << ", " << time << ", " << time/day << ", " << day_ct << endl; }

            // keep track of previous time
            previousTime = time; // TH - why?
            double totalRate = 0.0;
            VillageEvent ve = sample_event(RNG, totalRate, state_data, time); // This is where 'time' gets updated
            event_handler(ve, state_data, intervals, obs_time, lastPCase, villageExtinctionTimes, villageExtinctionIntervals, reinfectTime, RNG);
            time = ve.time;

            // start collecting data after burn in
            if (time > MIN_BURN_IN) {
                if (burn_in < 0) {
                    burn_in = time;

                    // determine state of system at beginning of burn in
                    // if there are infected individuals start clock for transmission interval
                    // if not, wait until there are infected individuals
                    for (size_t vil = 0; vil < NUM_OF_VILLAGES; ++vil) {
                        if (zero_local_infections(state_data, vil)) {
                            villageExtinctionTimes[vil] = 0; // t = 0 relative to observation period
                        }
                    }
                    if (not zero_infections(state_data)) {
                        reinfectTime = 0; // t = 0 relative to observation period
                    }
                }
                obs_time = time - burn_in; // here obs_time will always be >= 0.0

                double truncatePrevious = int(previousTime*10)/10.0; // truncate to 1 decimal place
                double truncateTime = int(time*10)/10.0;
                if (truncateTime > truncatePrevious) {
                    timeVec.push_back(obs_time);
                    append_output(output, state_data);
                }

                // stopping condition
                if (zero_infections(state_data) or (obs_time >= OBS_PERIOD)) {
                    // determine if there were no p cases in any villages

                    if (intervals[PCASE_INTERVAL].size() == 0) {
                        cerr << "Run failed, no cases observed\n\n";
                        if (RETRY_SIMS) { i--; } // CAN RESULT IN AN INFINITE LOOP!!!
                        // want at least two cases
                        intervals[EXTINCTION_INTERVAL].clear();
                    }

                    timeVec.push_back(obs_time);
                    append_output(output, state_data);
                    for (unsigned int vil = 0; vil < (unsigned) NUM_OF_VILLAGES; ++vil) {
                        if (zero_local_infections(state_data, vil) and extinction_observed(villageExtinctionTimes[vil])) {
                            villageExtinctionIntervals[vil].push_back(obs_time - villageExtinctionTimes[vil]);
                        }
                    }

                    for (unsigned int vil = 0; vil < (unsigned) NUM_OF_VILLAGES; ++vil) {
                        // output population counts
                        for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) {
                            for (size_t compartment_size: output[state][vil]) {
                                output_streams[state][vil] << compartment_size << SEP;
                            }
                            output_streams[state][vil] << endl;
                        }

                        // output extinction intervals for each village
                        for (size_t i = 0; i < villageExtinctionIntervals[vil].size(); ++i) {
                            outputVillageExtinctionInterval_stream[vil] << villageExtinctionIntervals[vil][i] << SEP;
                        }
                        outputVillageExtinctionInterval_stream[vil] << endl;

                    }

                    // output extinction time
                    outputQuant_stream[EXTINCTION_TIME] << obs_time <<endl;

                    // output time increments
                    for (unsigned int tLength = 0; tLength < timeVec.size(); ++tLength) { outputQuant_stream[TIME] << timeVec[tLength] << SEP; }
                    outputQuant_stream[TIME] << endl;

                    // output interval data
                    for (size_t interval_type = 0; interval_type < NUM_OF_INTERVAL_TYPES; ++interval_type) {
                        for (size_t i = 0; i < intervals[interval_type].size(); ++i) {
                            output_interval_streams[interval_type] << intervals[interval_type][i] << SEP;
                        }
                        output_interval_streams[interval_type] << endl;
                    }

                    // clear vectors when done
                    timeVec.clear();
                    intervals = vector<vector<double>>(NUM_OF_INTERVAL_TYPES);
                    villageExtinctionTimes.clear();
                    output = vector<vector<vector<int>>>(NUM_OF_STATE_TYPES, vector<vector<int>>(NUM_OF_VILLAGES));
                    villageExtinctionIntervals = vector<vector<double>> (NUM_OF_VILLAGES);
                    break;
                }
            }
        }
    }
    assert(eps_intercase_ivls.size() == eps_circ_ivls.size());
    output_results(output_streams, outputVillageExtinctionInterval_stream, outputQuant_stream, output_interval_streams);
    return 0;
}
