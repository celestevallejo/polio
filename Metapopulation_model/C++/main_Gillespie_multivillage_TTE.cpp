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

const bool RETRY_SIMS              = false;

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


#define REPORTABLE_EVENT_TYPE(VAR) \
    VAR(FIRST_INFECTION) \
    VAR(REINFECTION)\
    VAR(EXTINCTION) \
    VAR(NUM_OF_REPORTABLE_EVENTS)
enum ReportableEventType{ REPORTABLE_EVENT_TYPE(MAKE_ENUM) };
const char* const reportable_event_as_string[] = { REPORTABLE_EVENT_TYPE(MAKE_STRINGS) };
inline std::ostream& operator<<(std::ostream& out, const ReportableEventType value){ return out << reportable_event_as_string[value]; }


#define DETECTION_TYPE(VAR) \
    VAR(AFP_DETECTION) \
    VAR(ES_DETECTION) \
    VAR(EXTINCTION_DETECTION) \
    VAR(NUM_OF_DETECTION_TYPES)
enum DetectionType{ DETECTION_TYPE(MAKE_ENUM) };
const char* const detection_as_string[] = { DETECTION_TYPE(MAKE_STRINGS) };
inline std::ostream& operator<<(std::ostream& out, const DetectionType value){ return out << detection_as_string[value]; }


const map<EventType, pair<StateType, StateType>> TRANSITIONS = {{FIRST_INFECTION_EVENT,               {S, I1}},
                                                                {REINFECTION_EVENT,                   {P, IR}},
                                                                {RECOVERY_FROM_FIRST_INFECTION_EVENT, {I1, R}},
                                                                {RECOVERY_FROM_REINFECTION_EVENT,     {IR, R}},
                                                                {WANING_EVENT,                        {R, P}}};

class Parameters {
  public:
    Parameters(double _kappa, double _rho, double _numDaysToRecover, double _beta, double _deathRate, double _PIR, double _AFP_det, double _reintroRate,
               double _minBurnIn, double _obsPeriod, double _seasonalAmp, size_t _numSims, size_t _movModel, double _moveRate, double _ES_det, double _vacRate, vector<int> _village_pop) {
        kappa = _kappa;
        rho = _rho;
        numDaysToRecover = _numDaysToRecover;
        beta = _beta;
        deathRate = _deathRate;
        PIR = _PIR;
        AFP_det = _AFP_det;
        reintroRate = _reintroRate;
        minBurnIn = _minBurnIn;
        obsPeriod = _obsPeriod;
        seasonalAmp = _seasonalAmp;
        numSims = _numSims;
        movModel = _movModel;
        moveRate = _moveRate;
        ES_det = _ES_det;
        vacRate = _vacRate;
        village_pop = _village_pop;
        numVillages = village_pop.size();
        totalPop = accumulate(village_pop.begin(), village_pop.end(), 0.0);
    }

    double kappa, rho, PIR, reintroRate, seasonalAmp, ES_det, vacRate, beta, deathRate, numDaysToRecover, AFP_det, minBurnIn, obsPeriod, moveRate;
    size_t numSims, movModel, numVillages;
    double totalPop;
    vector<int> village_pop;
};

class DailyDetectedEvents {
  public:
    DailyDetectedEvents() {
        day = 0;
        detection_events = vector<int>(NUM_OF_DETECTION_TYPES, 0);
    };

    DailyDetectedEvents(int d) {
        day = d;
        detection_events = vector<int>(NUM_OF_DETECTION_TYPES, 0);
    };

    int day;
    vector<int> detection_events;
};

vector<int> vectorize_village_pop(string village_pop_str) {
    // removes all curly braces from the string
    // expected input format eg: {32000,32000}
    village_pop_str.erase(remove(village_pop_str.begin(), village_pop_str.end(), '{'), village_pop_str.end());
    village_pop_str.erase(remove(village_pop_str.begin(), village_pop_str.end(), '}'), village_pop_str.end());

    // splits string by commas and saves as a vector of ints
    vector<int> village_pop;
    stringstream ss(village_pop_str);
    string token;
    while (getline(ss, token, ',')) {
        village_pop.push_back(stoi(token));
    }

    return village_pop;
}

Parameters* parse_params(string filename, size_t par_line) {
    double kappa, rho, PIR, reintroRate, seasonalAmp, ES_det, vacRate, numDaysToRecover, beta, deathRate, AFP_det, minBurnIn, obsPeriod, moveRate;
    size_t numSims, movModel;
    string village_pop_str;

    ifstream iss(filename);
    if (!iss) {
        cerr << "ERROR: " << filename << " not found." << endl;
        exit(-101);
    }

    string buffer;
    istringstream line;
    int line_ct = -1;
    Parameters* par = nullptr;
    while (getline(iss, buffer)) {
        ++line_ct;
        if (line_ct == (int) par_line) {
            line.clear();
            line.str(buffer);

            if (line >> kappa >> rho >> numDaysToRecover >> beta >> deathRate >> PIR >> AFP_det >> reintroRate >> minBurnIn >> obsPeriod >> seasonalAmp
                     >> numSims >> movModel >> moveRate >> ES_det >> vacRate >> village_pop_str) {
                vector<int> village_pop = vectorize_village_pop(village_pop_str);
                par = new Parameters(kappa, rho, numDaysToRecover, beta, deathRate, PIR, AFP_det, reintroRate, minBurnIn, obsPeriod, seasonalAmp, numSims, movModel, moveRate, ES_det, vacRate, village_pop);
                break;
            } else {
                cerr << "ERROR: did not find valid parameter combination.  Found: " << buffer << endl;
                exit(-102);
            }
        }
    }
    return par;
}

struct VillageEvent{
    EventType event_type;
    int village;
    double time;
};

//string output_dir = "/home/tjhladish/work/polio-small-pop/output/";
//string output_dir ="/Users/Celeste/Desktop/multipatch_model/sim_results/";
//string output_dir ="/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/new_data_after_review/";
string output_dir ="./";
string ext = "_test.csv";
const string SEP = ",";

uniform_real_distribution<> runif(0.0, 1.0);

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

size_t rand_nonuniform_uint(const vector<int> weights, mt19937& RNG) {
    const double totalWeight = accumulate(weights.begin(), weights.end(), 0.0);
    double ran = totalWeight*runif(RNG);

    for (size_t idx = 0; idx < weights.size(); ++idx) {
        if (choose_event(ran, weights[idx])) {
            return idx;
            break;
        }
    }
    return weights.size(); // indicates failure to choose
}

vector<double> calculate_village_rates(const Parameters* par, const vector<vector<int>> &state_data, const int village, const double time) {
    // fast waning parameters:
    //kappa = 0.4179
    //rho = 0.2

    // intermediate waning parameters:
    //kappa = 0.6383
    //rho = 0.04

    // slow waning parameters:
    //kappa = 0.8434
    //rho = 0.02

//        vacRate = _vacRate; // TODO -- implement this

    const double KAPPA                 = par->kappa;            //0.4179;                // waning depth parameter
    const double RHO                   = par->rho;              //0.2;                   // waning speed parameter
    const double local_pop             = par->village_pop[village];
    const double RECOVERY              = 365.0/par->numDaysToRecover;// recovery rate (/year)
    const double BETA                  = par->beta;             //135.0;                   // contact rate (individuals/year)
    const double DEATH                 = par->deathRate;        //1.0/lifespan;          // death rate (per year)
    const double MOVE_RATE             = par->moveRate;         //expectedTimeUntilMove > 0 ? 1/expectedTimeUntilMove : 0;
    const double REINTRODUCTION_RATE   = par->reintroRate;      //1.0/10.0;              // was 0
//    const double VAC_RATE              = par->vacRate;
//    const double VAC_FRAC              = 0.20;                  // percentage of population vaccinated at birth
    const bool STRICT_TRAVEL           = (bool) par->movModel;  //true;                  // true: movement means moving to another patch; false: movement (implicitly) can mean moving within a patch
                                                            // original results were produced with the latter (STRICT_TRAVEL==false) interpretation

    const int S_  = state_data[S][village];
    const int I1_ = state_data[I1][village];
    const int R_  = state_data[R][village];
    const int P_  = state_data[P][village];
    const int IR_ = state_data[IR][village];
    double seasonalBeta = BETA*(1 + par->seasonalAmp*sin(time/(2*M_PI)));
    double foi = seasonalBeta*(I1_ + KAPPA*IR_)/local_pop;

    vector<double> local_event_rates(NUM_OF_EVENT_TYPES, 0.0);
    local_event_rates[FIRST_INFECTION_EVENT]                = S_*foi;               // first infection event
    local_event_rates[REINFECTION_EVENT]                    = KAPPA*P_*foi;         // reinfection event
    local_event_rates[RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1_;         // first infected revovery event
    local_event_rates[RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*IR_; // reinfected recovery event
    local_event_rates[WANING_EVENT]                         = RHO*R_;               // waning event
    local_event_rates[DEATH_EVENT]                          = DEATH*local_pop;      // natural death
    local_event_rates[MOVE_EVENT]                           = MOVE_RATE*local_pop;  // rate of movement from village
    if (not STRICT_TRAVEL) {
        local_event_rates[MOVE_EVENT] *= (par->totalPop - local_pop) / par->totalPop;   // if people can "move" to their own village
    }

    if (time < par->minBurnIn) {
        local_event_rates[REINTRODUCTION_EVENT]             = REINTRODUCTION_RATE*local_pop;// rate of reintroduction into system
    } else{
        local_event_rates[REINTRODUCTION_EVENT]             = 0;
    }
    return local_event_rates;
}

VillageEvent sample_event(const Parameters* par, const vector<vector<int>> &state_data, const double time, mt19937& RNG) {
    vector<vector<double>> event_rates(par->numVillages, vector<double>(NUM_OF_EVENT_TYPES));
    double totalRate = 0.0;
    VillageEvent ve;
    for (size_t vil = 0; vil < par->numVillages; ++vil) {
        event_rates[vil] = calculate_village_rates(par, state_data, vil, time);
        totalRate += accumulate(event_rates[vil].begin(), event_rates[vil].end(), 0.0);
    }
    exponential_distribution<> rexp(totalRate);
    ve.time = time + rexp(RNG);

    double ran = totalRate * runif(RNG);

    for (size_t event = 0; event < NUM_OF_EVENT_TYPES; ++event) {
        for (size_t vil = 0; vil < par->numVillages; ++vil) {
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

void move_from_A_to_B (vector<vector<int>> &state_data, const vector<size_t> &vil_pair, StateType state) {
    const size_t vil_A = vil_pair[0];
    const size_t vil_B = vil_pair[1];

    --state_data[state][vil_A];
    ++state_data[state][vil_B];
}


StateType choose_villager(const vector<vector<int>> &state_data, const size_t village, mt19937 &RNG) {
    vector<int> weights(NUM_OF_STATE_TYPES, 0);
    for (size_t state = 0; state < NUM_OF_STATE_TYPES; ++state) {
        weights[state] = state_data[state][village];
    }

    return (StateType) rand_nonuniform_uint(weights, RNG);
}


void process_movement_event(const Parameters* par, vector<vector<int>> &state_data, const size_t vil_A, mt19937& RNG) {
    vector<int> weights(par->village_pop);
    weights[vil_A] = 0;
    const size_t vil_B = rand_nonuniform_uint(weights, RNG);

    // Sample states of people to move from A to B and vice versa
    const StateType state_A = choose_villager(state_data, vil_A, RNG);
    const StateType state_B = choose_villager(state_data, vil_B, RNG);

    if (state_A != state_B) {
        move_from_A_to_B(state_data, {vil_A, vil_B}, state_A);
        move_from_A_to_B(state_data, {vil_B, vil_A}, state_B);
    }
}


void process_reintroduction_event(vector<vector<int>> &state_data, const size_t village, mt19937 &RNG) {
    // treat reintroduction as exposure event that is density dependent
    // Sample state of person to whom exposure event will occur
    StateType state = choose_villager(state_data, village, RNG);

    // exposure events only change state of S and P compartments
    EventType event = state == S ? FIRST_INFECTION_EVENT :
                      state == P ? REINFECTION_EVENT :
                      NUM_OF_EVENT_TYPES;

    if (TRANSITIONS.count(event) > 0) {
        --state_data[TRANSITIONS.at(event).first][village];
        ++state_data[TRANSITIONS.at(event).second][village];
    }
}


void append_output(const Parameters* par, vector<vector<vector<int>>> &output, const vector<vector<int>> &state_data) {
    for (size_t state = 0; state < NUM_OF_STATE_TYPES; ++state) {
        for (size_t vil = 0; vil < par->numVillages; ++vil) {
            output[state][vil].push_back(state_data[state][vil]);
        }
    }
}


void event_handler(const Parameters* par, const VillageEvent &ve, vector<vector<int>> &state_data, const double obs_time, vector<int> &reportable_event_ct, mt19937& RNG) {
    const size_t village  = ve.village;
    const EventType event = ve.event_type;
    StateType prev_state  = NUM_OF_STATE_TYPES;
    StateType new_state   = NUM_OF_STATE_TYPES;

    // transitions that are known exactly at this point are handled here
    if (TRANSITIONS.count(event) > 0) {
        prev_state = TRANSITIONS.at(event).first;
        new_state  = TRANSITIONS.at(event).second;
        --state_data[prev_state][village];
        ++state_data[new_state][village];
    }

    switch(ve.event_type) {
        case FIRST_INFECTION_EVENT:{
            if (obs_time > 0) {
                reportable_event_ct[FIRST_INFECTION]++; // will be determined later whether paralytic and detected
            }
        }
            break;
        case RECOVERY_FROM_FIRST_INFECTION_EVENT: //[[fallthrough]]
        case RECOVERY_FROM_REINFECTION_EVENT:
            if (obs_time > 0 and zero_infections(state_data)) {
                reportable_event_ct[EXTINCTION]++;
            }
            break;
        case DEATH_EVENT:
            {
                StateType deceased = choose_villager(state_data, village, RNG);
                --state_data[deceased][village];
                ++state_data[S][village];
                if (is_infected_state(deceased) and obs_time > 0 and zero_infections(state_data)) {
                    reportable_event_ct[EXTINCTION]++;
                }
            }
            break;
        case MOVE_EVENT:
            process_movement_event(par, state_data, village, RNG);
            break;
        case REINTRODUCTION_EVENT:
            if (obs_time > 0) {
                process_reintroduction_event(state_data, village, RNG);
            }
            break;
        case REINFECTION_EVENT:
            if (obs_time > 0)  {
                reportable_event_ct[REINFECTION]++; // will be determined later whether environmental surveillance occurred
            }
            break;
        case WANING_EVENT: break;
        default:
            cerr << "ERROR: Unsupported event type" << endl;
            break;
    }
}


void log_output(const vector<DailyDetectedEvents> /*&detected_event_ct_ts*/, const vector<vector<int>> /*&initial_states*/) {
//serial par_combo replicate village_id village_size S I1 R P IR

//serial par_combo replicate day_of_event event_type event_type_ct
// reportable events include: population-wide extinction (I1 + IR -> 0 during observation period); S->I1 and P->IR
// possibly will want patch ID as a output dimension

//for (size_t state = 0; state < initial_states.size(); ++state) {
//   for(size_t village = 0; village < initial_states[state].size(); ++village) {
//   }
//}

}


int main(int argc, char** argv) {
    assert(argc == 3);
    string par_file = argv[1];
    size_t par_idx  = atoi(argv[2]); // which line in the par file to use
    Parameters* par = parse_params(par_file, par_idx);

//cerr << par->kappa << ' ' << par->rho << ' ' << par->numDaysToRecover << ' ' << par->beta << ' ' << par->PIR << ' ' << par->AFP_det << ' ' << par->reintroRate << ' ' << par->minBurnIn << ' ' << par->obsPeriod
//     << ' ' << par->seasonalAmp << ' ' << par->numSims << ' ' << par->movModel << ' ' << par->moveRate << ' ' << par->ES_det << ' ' << par->vacRate << ' ';
//for (int vp : par->village_pop) {
//    cerr << vp << ' ';
//}
//cerr << endl;
//exit(-565);

    // TODO -- review what data is being used and for what.  we may not need everything we're outputting
    // vectors for holding information to be output

    //random_device rd;                                   // generates a random real number for the seed
    // unsigned long int seed = rd();
    const size_t initial_seed = 8675309;
    //cerr << "seed: " << seed << endl;

    // The Simulation
    for (size_t i = 0; i < par->numSims; ++i) {
        // initialize size of vectors for individual compartments
        vector<vector<int>> state_data(NUM_OF_STATE_TYPES, vector<int>(par->numVillages, 0.0));

        mt19937 RNG(initial_seed + i);                  // random number generator
        double time         = 0.0;                      // "absolute" time, relative to start of simulation
        double burn_in      = -DBL_MAX;                 // don't know yet exactly how long burn-in will be.  burn_in ends with first event after MIN_BURN_IN
        double obs_time     = -DBL_MAX;                 // time measured since burn-in

        for (size_t vil = 0; vil < par->numVillages; ++vil) {
            // set initial values for each village using multinomial dist
            //initialValues[vil] = multinomial_Compartments(compartments[vil].size(),compartments[vil],vil,RNG());
            const double local_pop    = par->village_pop[vil];
            state_data[S][vil]        = (int) (0.99*local_pop);         // naive susceptible (no previous contact w/virus, moves into I1)
            state_data[I1][vil]       = local_pop - state_data[S][vil]; // first infected (only time paralytic case can occur, recovers into R)
            state_data[R][vil]        = 0;                              // recovered (fully immune, wanes into P)
            state_data[P][vil]        = 0;                              // partially susceptible (moves into IR)
            state_data[IR][vil]       = 0;                              // reinfected (recovers into R)
        }

        //long long int day_ct = -1;
        const double day_length = 1.0/365.0;

        // Two new output files w columns:
        // initial conditions file: serial par_combo replicate village_id village_size S I1 R P IR
        // event data:              serial par_combo replicate day_of_event event_type event_type_ct
        vector<vector<int>> initial_states;                                                 // for each village, number of people in each compartment just after burn-in
        vector<int> reportable_event_ct         = vector<int>(NUM_OF_REPORTABLE_EVENTS, 0); // things that *might* be detected
        vector<DailyDetectedEvents> detected_event_ct_ts;                                   // time series of things that *were* detected

        int prev_day                            = 0;
        int num_days_with_detections            = 0;

        while (true) {
            const double day = time/day_length;
            if (obs_time > 0 and (int) day > prev_day) { // 'tis a new day!  tally what happened yesterday
                DailyDetectedEvents dde(prev_day);

                binomial_distribution<int> AFP_det_binom(reportable_event_ct[FIRST_INFECTION], par->PIR * par->AFP_det);
                binomial_distribution<int> ES_det_binom(reportable_event_ct[FIRST_INFECTION] + reportable_event_ct[REINFECTION], par->ES_det);

                dde.detection_events[AFP_DETECTION]        = AFP_det_binom(RNG);
                dde.detection_events[ES_DETECTION]         = ES_det_binom(RNG);
                dde.detection_events[EXTINCTION_DETECTION] = reportable_event_ct[EXTINCTION];

                bool had_detections = dde.detection_events[AFP_DETECTION] + dde.detection_events[ES_DETECTION] > 0;
                num_days_with_detections += had_detections;

                detected_event_ct_ts.push_back(dde);
                reportable_event_ct = vector<int>(NUM_OF_REPORTABLE_EVENTS, 0); // things that *might* be detected
            }

            /* while (time/day_length > day_ct) { // increment day counter as needed
                // TODO -- if after burn-in and isNewDay, report extinctions and *detected* infections as appropriate, then clear counts
                if (day_ct % 100 == 0) {
                    cerr << "step, time, day: " << right << setw(9) << setprecision(3) << time << ", " << setw(6) << day_ct << " | ";
                    cerr_state(state_data, 0);
                }
                ++day_ct;
            } */

            if (time > par->minBurnIn) {
                if (burn_in < 0) { // this is the first event after the burn-in has been completed
                    burn_in = time;
                    initial_states = state_data; // capture state of system at end of burn-in
                }
                obs_time = time - burn_in; // here obs_time will always be >= 0.0

                // stopping condition
                if (zero_infections(state_data) or (obs_time >= par->obsPeriod)) {
                    if (num_days_with_detections < 2) { // no intercase intervals to report
                        cerr << "Run failed, <2 days with detected infections\n\n";
                        if (RETRY_SIMS) { i--; } // CAN RESULT IN AN INFINITE LOOP!!!
                    }
                    break;
                }
            }

            VillageEvent ve = sample_event(par, state_data, time, RNG); // This is where 'time' gets updated
            event_handler(par, ve, state_data, obs_time, reportable_event_ct, RNG);
            time = ve.time;
        }
        log_output(detected_event_ct_ts, initial_states);
    }
//    output_results();
    return 0;
}
