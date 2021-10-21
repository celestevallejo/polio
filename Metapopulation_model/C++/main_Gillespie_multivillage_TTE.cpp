#include <iostream>
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


using namespace std;

//string output_dir = "/home/tjhladish/work/polio-small-pop/output/";
//string output_dir ="/Users/Celeste/Desktop/multipatch_model/sim_results/";
string output_dir ="/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/new_data_after_review/";
//string output_dir = "/home/vallejo.26/";
string ext = "_test.csv";
const string SEP = ",";
const string parameterFileDatabase = "/home/vallejo.26/parameterDatabase.csv";
//const string parameterFileDatabase = "/Users/Celeste/Desktop/multipatch_model/multiPatch_repo/parameterDatabase.csv";

uniform_real_distribution<> unifdis(0.0, 1.0);

enum StateType {S,
                I1,
                R,
                P,
                IR,
                NUM_OF_STATE_TYPES};

enum EventType {MOVE_EVENT,
                FIRST_INFECTION_EVENT,
                REINFECTION_EVENT,
                RECOVERY_FROM_FIRST_INFECTION_EVENT,
                RECOVERY_FROM_REINFECTION_EVENT,
                WANING_EVENT,
                DEATH_EVENT,
                REINTRODUCTION_EVENT,
                NUM_OF_EVENT_TYPES};

enum QuantOutputType {TIME_OUT,
                      TIME_BETWEEN_PCASES_OUT,
                      TRANSMISSION_INTERVAL_OUT,
                      EXTINCTION_INTERVAL_OUT,
                      EXTINCTION_TIME_OUT,
                      NUM_OF_QUANTOUTPUT_TYPES};


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
    EventType event;
    int village;
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
const vector<int> village_pop      = {1000, 1000};
const int numVillages              = village_pop.size(); //total number of villages under consideration
const int numDaysToRecover         = 28;
const double RECOVERY              = 365/numDaysToRecover;    //recovery rate (/year)
const double BETA                  = 135;   //contact rate (individuals/year)
const double lifespan              = 50;
const double BIRTH                 = 1/lifespan; //birth rate (per year)
const double DEATH                 = 1/lifespan; //death rate (per year)
const double PIR                   = 0.005;            //type 1 paralysis rate (naturally occurring cases)
const double DET_RATE              = 1.0;
const double expectedTimeUntilMove = 0; //years
const double MOVE_RATE             = expectedTimeUntilMove > 0 ? 1/expectedTimeUntilMove : 0;
const double REINTRODUCTION_RATE      = 0;
const double burnInTime            = 0; //years
const double obsTime               = 0; //years
const double seasonalAmp           = 0.05;
const int numSims                  = 1;

vector<vector<double>> event_rates(numVillages, vector<double>(NUM_OF_EVENT_TYPES, 0.0));

bool choose_event(double &ran, const double p) {
    if (ran < p) {
        return true;
    } else {
        ran -= p;
        return false;
    }
}

unsigned int rand_nonuniform_uint(const vector<int> weights, mt19937& gen) {
    const double totalWeight = accumulate(weights.begin(), weights.end(), 0.0);
    double ran = totalWeight*unifdis(gen);

    for (unsigned int idx = 0; idx < weights.size(); ++idx) {
        if (choose_event(ran, weights[idx])) {
            return idx;
            break;
        }
    }
    return weights.size(); // indicates failure to choose
}

void calculate_rates(const vector<vector<int>> &state_data, const int i, const double time) {
    const int S_  = state_data[S][i];
    const int I1_ = state_data[I1][i];
    const int R_  = state_data[R][i];
    const int P_  = state_data[P][i];
    const int IR_ = state_data[IR][i];
    //double seasonalBeta = BETA*(1+seasonalAmp*sin(time/(2*M_PI)));
    //double foi = seasonalBeta*(I1[i]+ KAPPA*IR[i])/village_pop[i];
    double foi = BETA*(I1_+ KAPPA*IR_)/village_pop[i];
    
    event_rates[i][FIRST_INFECTION_EVENT]                = S_*foi;                          // first infection event
    event_rates[i][REINFECTION_EVENT]                    = KAPPA*P_*foi;                    // reinfection event
    event_rates[i][RECOVERY_FROM_FIRST_INFECTION_EVENT]  = RECOVERY*I1_;                    // first infected revovery event
    event_rates[i][RECOVERY_FROM_REINFECTION_EVENT]      = (RECOVERY/KAPPA)*IR_;            // reinfected recovery event
    event_rates[i][WANING_EVENT]                         = RHO*R_;                          // waning event
    event_rates[i][DEATH_EVENT]                          = DEATH*village_pop[i];            // natural death
    event_rates[i][MOVE_EVENT]                           = MOVE_RATE*village_pop[i];        // rate of movement from village i

    if(time < burnInTime){
        event_rates[i][REINTRODUCTION_EVENT]             = REINTRODUCTION_RATE*village_pop[i]; // rate of reintroduction into system
    } else{
        event_rates[i][REINTRODUCTION_EVENT]             = 0;
    }
}

VillageEvent sample_event(mt19937& gen, double& totalRate) {
    totalRate = 0.0;
    VillageEvent ve;
    for(int vil = 0; vil < numVillages; vil++){
        totalRate += accumulate(event_rates[vil].begin(),event_rates[vil].end(),0.0);
    }
    double ran = totalRate*unifdis(gen);

    for (int event = 0; event < NUM_OF_EVENT_TYPES; ++event) {
        for(int vil = 0; vil < numVillages; vil++){
            if (choose_event(ran, event_rates[vil][event])) {
                ve.event = (EventType) event;
                ve.village = vil;
                return ve;
            }
        }
    }
    cerr << "Error: No event sampled" << endl;
    exit(-100);
}

//void update_compartments(vector<int> &S, vector<int> &I1, vector<int> &R, vector<int> &P, vector<int> &IR, const uint A, const uint B, const StateType state_from_A, const StateType state_from_B, const double time, const double burnInTime, const double timeAfterBurnIn, vector<vector<double>> &villageExtinctionIntervals,vector<double> &villageExtinctionTime) {
void update_compartments(vector<vector<int>> &state_data, vector<unsigned int> vil_pair, const vector<StateType> &state_pair, const double time, const double burnInTime, const double timeAfterBurnIn, vector<vector<double>> &villageExtinctionIntervals,vector<double> &villageExtinctionTime) {

    for (int i = 0; i < 2; ++i) { // foreach village (vil_pair reversed at end of first loop)
        unsigned int A = vil_pair[0];
        unsigned int B = vil_pair[1];
        switch (state_pair[i]) {
            case S:  --state_data[S][A];  ++state_data[S][B];  break;
            case I1:{
                if(time > burnInTime){
                    if((state_data[I1][B]+state_data[IR][B]==0) and villageExtinctionTime[B]!=numeric_limits<double>::max()){
                        villageExtinctionIntervals[B].push_back((time-timeAfterBurnIn) - villageExtinctionTime[B]);
                    }
                }
                --state_data[I1][A]; ++state_data[I1][B];
                double extinctionTime = time - timeAfterBurnIn;
                if(state_data[I1][A]+state_data[IR][A]==0 and time > burnInTime){
                    villageExtinctionTime[A] = extinctionTime;
                }
            }
                break;
            case R:  --state_data[R][A];  ++state_data[R][B];  break;
            case P:  --state_data[P][A];  ++state_data[P][B];  break;
            case IR:{
                if(time > burnInTime){
                    if((state_data[I1][B]+state_data[IR][B]==0) and villageExtinctionTime[B]!=numeric_limits<double>::max()){
                        villageExtinctionIntervals[B].push_back((time-timeAfterBurnIn) - villageExtinctionTime[B]);
                    }
                }
                --state_data[IR][A]; ++state_data[IR][B];
                double extinctionTime = time - timeAfterBurnIn;
                if(state_data[I1][A]+state_data[IR][A]==0 and time > burnInTime){
                    villageExtinctionTime[A] = extinctionTime;
                }
            }
                break;
            default: break;
        }

        vil_pair = {vil_pair[1], vil_pair[0]};
    }

    calculate_rates(state_data, vil_pair[0], time);
    calculate_rates(state_data, vil_pair[1], time);
}

bool all_zeroes(const vector<int> &data) { return all_of(data.begin(), data.end(), [](int i){return i==0;}); }

bool no_infections(const vector<vector<int>> &state_data) { return all_zeroes(state_data[I1]) and all_zeroes(state_data[IR]); }

bool no_local_infections(const vector<vector<int>> &state_data, const unsigned int village) { return state_data[I1][village] == 0 and state_data[IR][village] == 0; }

inline void process_death_event(vector<vector<int>> &state_data, const int chosenVillage, mt19937& gen, vector<double> &transmissionIntervals, const double time, const double timeAfterBurnIn, double &reinfectTime, const double burnInTime, vector<double> &extinctionIntervals, double &lastPCase, vector<double> &villageExtinctionTime) {
    const int j = chosenVillage;
    vector<int> rates(NUM_OF_STATE_TYPES, 0);
    rates[S]  = state_data[S][j];
    rates[I1] = state_data[I1][j];
    rates[R]  = state_data[R][j];
    rates[P]  = state_data[P][j];
    rates[IR] = state_data[IR][j];

    StateType source_state = (StateType) rand_nonuniform_uint(rates, gen);

    if (source_state == S) return; // important to bail now, since nothing happens in this case

    ++state_data[S][j];
    switch(source_state) {
        case I1:{
            --state_data[I1][j];
            if(time > burnInTime){
                double extinctionTime = time - timeAfterBurnIn;
                if(state_data[I1][j]+state_data[IR][j]==0){
                    villageExtinctionTime[j] = extinctionTime;
                }

                if(no_infections(state_data)){
                    if(reinfectTime!=numeric_limits<double>::max()){
                        transmissionIntervals.push_back(extinctionTime - reinfectTime);
                    }
                    if(lastPCase!= numeric_limits<double>::max()){
                        extinctionIntervals.push_back(extinctionTime - lastPCase);
                    }
                }
            }
        }
            break;
        case R:
            --state_data[R][j];
            break;
        case P:
            --state_data[P][j];
            break;
        case IR:{
            --state_data[IR][j];
            if(time > burnInTime){
                double extinctionTime = time - timeAfterBurnIn;
                if(state_data[I1][j]+state_data[IR][j]==0){
                    villageExtinctionTime[j] = extinctionTime;
                }
                if(no_infections(state_data)){
                    if(reinfectTime!=numeric_limits<double>::max()){
                        transmissionIntervals.push_back(extinctionTime - reinfectTime);
                    }
                    if(lastPCase!= numeric_limits<double>::max()){
                        extinctionIntervals.push_back(extinctionTime - lastPCase);
                    }
                }
            }
        }
            break;
        default:
            break;
    }
}

StateType sample_migrant(vector<vector<int>> &state_data, const int village, mt19937& gen) {
    vector<int> weights(NUM_OF_STATE_TYPES, 0);
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) { weights[state] = state_data[state][village]; } 
    return (StateType) rand_nonuniform_uint(weights, gen);
}

inline void process_movement_event(vector<vector<int>> &state_data, vector<unsigned int> vil_pair, mt19937& gen, const double time, const double burnInTime, const double timeAfterBurnIn, vector<vector<double>> &villageExtinctionIntervals, vector<double> &villageExtinctionTime){
    vector<int> weights(NUM_OF_STATE_TYPES, 0);
    vector<StateType> state_pair;

    // Sample states of people to move from A to B and vice versa
    state_pair.push_back(sample_migrant(state_data, vil_pair[0], gen));
    state_pair.push_back(sample_migrant(state_data, vil_pair[1], gen));

    if (state_pair[0] == state_pair[1]) {
        return; // nothing actually happens
    } else {
        update_compartments(state_data, vil_pair, state_pair, time, burnInTime, timeAfterBurnIn, villageExtinctionIntervals, villageExtinctionTime);
    }
}

void process_reintroduction_helper(const StateType from, const StateType to, const double time, vector<vector<int>> &state_data, const int village, double &timeAfterBurnIn, double &reinfectTime, vector<vector<double>> &villageExtinctionIntervals, vector<double> &villageExtinctionTime) {
    if (time > burnInTime){
        if ( no_local_infections(state_data, village) and villageExtinctionTime[village] != numeric_limits<double>::max() ) {
            villageExtinctionIntervals[village].push_back((time - timeAfterBurnIn) - villageExtinctionTime[village]);
        }
        if(no_infections(state_data)){
            reinfectTime = time - timeAfterBurnIn;
        }
    }
    --state_data[from][village]; ++state_data[to][village];
    calculate_rates(state_data, village, time);
}

inline void process_reintroduction_event(vector<vector<int>> &state_data, const int village, mt19937& gen, vector<double> &timeBetweenPCases, double time, double &timeAfterBurnIn, double &reinfectTime, double &lastPCase, vector<vector<double>> &villageExtinctionIntervals,vector<double> &villageExtinctionTime){
    //treat reintroduction events as exposure events
    vector<int> weights(NUM_OF_STATE_TYPES, 0);

    // Sample state of person to which exposure event will occur
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) { weights[state] = state_data[state][village]; }
    StateType state = (StateType) rand_nonuniform_uint(weights, gen);

    //exposure events only change state of S and P compartments
    if (state == S){
        process_reintroduction_helper(S, I1, time, state_data, village, timeAfterBurnIn, reinfectTime, villageExtinctionIntervals, villageExtinctionTime);

        double rr = unifdis(gen);
        if(time > burnInTime){//start keeping track after burn in
            if(rr<(PIR*DET_RATE)){
                assert(timeAfterBurnIn > burnInTime);
                if(lastPCase != numeric_limits<double>::max()){
                    timeBetweenPCases.push_back(time - timeAfterBurnIn - lastPCase);
                }
                lastPCase = time - timeAfterBurnIn;
            }
        }
    } else if (state == P) {
        process_reintroduction_helper(P, IR, time, state_data, village, timeAfterBurnIn, reinfectTime, villageExtinctionIntervals, villageExtinctionTime);
    }
    
}

void output_results(vector<vector<stringstream>> &output_streams, vector<stringstream> &outputVillageExtinctionInterval_stream, vector<stringstream> &outputQuant_streams) {
    string numInEachVil = "";
    for(int i = 0; i < numVillages; i++){
        cout<<"num in each vil"<<village_pop[i]<<"\n";
        numInEachVil += to_string(int(village_pop[i]));
    }

    string base_filename = numInEachVil + "reintRate_"+to_string(REINTRODUCTION_RATE) + "migRate_"+to_string(MOVE_RATE)+ext;
    string parameters = +"beta_"+to_string(int(BETA))+"detect_rate_"+to_string(float(DET_RATE))+"rho_"+to_string(float(RHO))+ "numVillages_"+to_string(numVillages) + "migRate_"+to_string(float(MOVE_RATE))+"burnIn_"+to_string(burnInTime)+"obsTime_"+to_string(obsTime) + "seasonalAmp_"+to_string(seasonalAmp);
    //read parameters into database file
    fstream params;
    params.open(parameterFileDatabase, fstream::app);
    params<<base_filename<<SEP<<parameters<<endl;
    params.close();
    
    vector<vector<string>> output_filenames(NUM_OF_STATE_TYPES, vector<string>(numVillages));
//    vector<string> outputS_filenames(numVillages);
//    vector<string> outputI1_filenames(numVillages);
//    vector<string> outputR_filenames(numVillages);
//    vector<string> outputP_filenames(numVillages);
//    vector<string> outputIR_filenames(numVillages);
    vector<string> outputVillageExtinctionInterval_filenames(numVillages);
    
    vector<string> outputQuant_filenames(NUM_OF_QUANTOUTPUT_TYPES);

    outputQuant_filenames[TIME_BETWEEN_PCASES_OUT ] = output_dir + "time_between_pcases_"+base_filename;
    outputQuant_filenames[TRANSMISSION_INTERVAL_OUT ] = output_dir + "transmission_interval_"+base_filename;
    outputQuant_filenames[EXTINCTION_INTERVAL_OUT ] = output_dir + "extinction_interval_"+base_filename;
    outputQuant_filenames[TIME_OUT ] = output_dir + "TIME_"+base_filename;
    outputQuant_filenames[EXTINCTION_TIME_OUT ] = output_dir + "extTime_"+base_filename;
    

    for(int vil = 0; vil < numVillages; vil++){
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
    for(int vil = 0; vil < numVillages; vil++){
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

void process_recovery_event(const vector<vector<int>> &state_data, const unsigned int village, const double time, const double timeAfterBurnIn, vector<double> &villageExtinctionTime, vector<double> &transmissionIntervals, vector<double> &extinctionIntervals, const double reinfectTime, const double lastPCase) {
    calculate_rates(state_data, village, time);
    if(time > burnInTime){
        const double extinctionTime = time - timeAfterBurnIn;
        if(state_data[I1][village] + state_data[IR][village] == 0){
            villageExtinctionTime[village] = extinctionTime;
        }
        if(no_infections(state_data)){
            if(reinfectTime!=numeric_limits<double>::max()){
                transmissionIntervals.push_back(extinctionTime - reinfectTime);
            }
            if(lastPCase!= numeric_limits<double>::max()){
                extinctionIntervals.push_back(extinctionTime - lastPCase);
            }
        }
    }
}

void append_output(vector<vector<vector<int>>> &output, const vector<vector<int>> &state_data) {
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) {
        for(int vil = 0; vil < numVillages; vil++){
            output[state][vil].push_back(state_data[state][vil]);
        }
    }
}


int main(){
    
    //vectors for holding information to be output
    vector<vector<vector<int>>> output(NUM_OF_STATE_TYPES, vector<vector<int>>(numVillages));
    //vector<vector<int>> outputS(numVillages);
    //vector<vector<int>> outputI1(numVillages);
    //vector<vector<int>> outputR(numVillages);
    //vector<vector<int>> outputP(numVillages);
    //vector<vector<int>> outputIR(numVillages);
    
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
    vector<vector<double>> villageExtinctionIntervals(numVillages);
    //keep track of reinfection time in each village
    vector<double> villageExtinctionTime(numVillages);
    
    //vectors for outputing information
    vector<vector<stringstream>> output_streams;
    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) {
        output_streams.push_back(vector<stringstream>(numVillages));
    }
    //vector<vector<stringstream>> output_streams(NUM_OF_STATE_TYPES, vector<stringstream> (numVillages)); // for some reason, compiler rejects this!
    //vector<stringstream> outputS_stream(numVillages);
    //vector<stringstream> outputI1_stream(numVillages);
    //vector<stringstream> outputR_stream(numVillages);
    //vector<stringstream> outputP_stream(numVillages);
    //vector<stringstream> outputIR_stream(numVillages);
    vector<stringstream> outputQuant_stream(NUM_OF_QUANTOUTPUT_TYPES);
    
    vector<stringstream> outputVillageExtinctionInterval_stream(numVillages);

    //initialize size of vector of vector of compartments
    //vector<vector<double>> compartments(numVillages);

    //initialize size of vector of vector of initial values
    //vector<vector<int>> initialValues(numVillages);

    //initialize size of vectors for individual compartments
    vector<vector<int>> state_data(NUM_OF_STATE_TYPES, vector<int>(numVillages, 0.0));

    random_device rd;                       // generates a random real number for the seed
    unsigned long int seed = rd();
    //uint seed = 2186064846;
    cerr << "seed: " << seed << endl;
    mt19937 gen(seed);                      // random number generator
    //mt19937 gen(rd());                      // random number generator

    //find expected compartment size for each village
    /*for(int vil = 0; vil < numVillages; vil++){
        compartments[vil] = initialize_compartment(vil);
    }*/
    const int EPS_RES = 100; // resolution of endemic potential statistic, in divisions per year
    const int EPS_MAX = 50;  // circulation interval considered for EPS calculation, in years
    vector<int> eps_circ_ivls(EPS_MAX*EPS_RES, 0); // 50 yrs divided into 100 bins each
    vector<int> eps_intercase_ivls(EPS_MAX*EPS_RES, 0); // 50 yrs divided into 100 bins each

    //The Simulation
    for(int i = 0; i < numSims; i++){
        double time = 0.0;
        double previousTime = 0.0;
        bool firstBurnInLoop = true; //keeps track of first time entering burn in loop
        double timeAfterBurnIn; //keeps track of time after burn in over
        lastPCase = numeric_limits<double>::max(); //initialize to largest double
        reinfectTime = numeric_limits<double>::max(); //initialize to largest double
        
        for(int vil=0; vil<numVillages; vil++){
            villageExtinctionTime[vil] = numeric_limits<double>::max(); //initialize to largest double
        }

        for(int vil = 0; vil < numVillages; vil++){
            //set initial values for each village using multinomial dist
            //initialValues[vil] = multinomial_Compartments(compartments[vil].size(),compartments[vil],vil,gen());
            state_data[S][vil]        = (int)0.99*village_pop[vil];   //naive susceptible (no previous contact w/virus, moves into I1)
            state_data[I1][vil]       = village_pop[vil] - state_data[S][vil];  //first infected (only time paralytic case can occur, recovers into R)
            state_data[R][vil]        = 0;   //recovered (fully immune, wanes into P)
            state_data[P][vil]        = 0;   //partially susceptible (moves into IR)
            state_data[IR][vil]       = 0;  //reinfected (recovers into R)
            calculate_rates(state_data, vil, time);
        }

        for(int j = 0; j < 1e10; j++){
            double totalRate = 0.0;
            VillageEvent ve = sample_event(gen, totalRate);
            EventType event_type = ve.event;
            const uint A = ve.village;

            switch(event_type){
                case FIRST_INFECTION_EVENT:{
                    --state_data[S][A]; ++state_data[I1][A];
                    //check to see if paralytic case occurs
                    double rr = unifdis(gen);
                    if(time > burnInTime){//start keeping track after burn in
                        if(rr<(PIR*DET_RATE)){
                            assert(timeAfterBurnIn > burnInTime);
                            if(lastPCase != numeric_limits<double>::max()){
                                timeBetweenPCases.push_back(time - timeAfterBurnIn - lastPCase);
                            }
                            lastPCase = time - timeAfterBurnIn;
                        }
                    }
                    
                    calculate_rates(state_data, A, time);
                }
                    break;
                case REINFECTION_EVENT:{
                    --state_data[P][A]; ++state_data[IR][A];
                    calculate_rates(state_data, A, time);
                }
                    break;
                case RECOVERY_FROM_FIRST_INFECTION_EVENT:{
                    --state_data[I1][A]; ++state_data[R][A];
                    process_recovery_event(state_data, A, time, timeAfterBurnIn, villageExtinctionTime, transmissionIntervals, extinctionIntervals, reinfectTime, lastPCase);
                }
                    break;
                case RECOVERY_FROM_REINFECTION_EVENT:{
                    --state_data[IR][A]; ++state_data[R][A];
                    process_recovery_event(state_data, A, time, timeAfterBurnIn, villageExtinctionTime, transmissionIntervals, extinctionIntervals, reinfectTime, lastPCase);
                }
                    break;
                case WANING_EVENT:{
                    --state_data[R][A]; ++state_data[P][A];
                    calculate_rates(state_data, A, time);
                }
                    break;
                case DEATH_EVENT:{
                    process_death_event(state_data, A, gen, transmissionIntervals, time, timeAfterBurnIn, reinfectTime, burnInTime, extinctionIntervals, lastPCase, villageExtinctionTime);
                    calculate_rates(state_data, A, time);
                }
                    break;
                case MOVE_EVENT:{
                    uint B = rand_nonuniform_uint(village_pop,gen);
                    //if no movement occurs then no rates need to be updated
                    if(A != B){
                        process_movement_event(state_data, {A, B}, gen, time, burnInTime, timeAfterBurnIn, villageExtinctionIntervals, villageExtinctionTime);
                        //rates are calculated after movement
                    }
                }
                    break;
                case REINTRODUCTION_EVENT:{
                    //treat reintroduction as exposure event that is density dependent
                    uint village = rand_nonuniform_uint(village_pop,gen);
                    process_reintroduction_event(state_data, village, gen, timeBetweenPCases, time, timeAfterBurnIn, reinfectTime, lastPCase, villageExtinctionIntervals, villageExtinctionTime);
                    //rates are recalculated after event is processed
                }
                    break;
                default:
                    cerr << "ERROR: Unsupported event type" << endl;
                    break;
            }

            //keep track of previous time
            previousTime = time;
            //generate the time at which the event occurs
            exponential_distribution<>rng(totalRate);
            time+=rng(gen);

            //start collecting data after burn in
            if(time > burnInTime){
                if(firstBurnInLoop){
                    timeAfterBurnIn = time;
                    firstBurnInLoop = false; //to keep timeAfterBurnIn constant

                    //determine state of system at beginning of burn in
                    //if there are infected individuals start clock for transmission interval
                    //if not, wait until there are infected individuals
                    for(int vil = 0; vil < numVillages; vil++){
                        if (no_local_infections(state_data, vil)) {
                            villageExtinctionTime[vil] = time - timeAfterBurnIn;
                        }
                    }
                    if(not no_infections(state_data)){
                        reinfectTime = time - timeAfterBurnIn;
                    }
                }

                double truncatePrevious = int(previousTime*10)/10.0; //truncate to 1 decimal place
                double truncateTime = int(time*10)/10.0;
                if(truncateTime > truncatePrevious){
                    timeVec.push_back(time - timeAfterBurnIn);
                    append_output(output, state_data);
                }
            }

            //stopping condition
            if((no_infections(state_data) and time > burnInTime) or ((time - timeAfterBurnIn) >= obsTime)){
                //determine if there were no p cases in any villages
                
                if(timeBetweenPCases.size() == 0){
                    i--;
                    //want at least two cases
                    extinctionIntervals.clear();
                }

                timeVec.push_back(time - timeAfterBurnIn);
                append_output(output, state_data);
                for(unsigned int vil = 0; vil < (unsigned) numVillages; vil++){
                    if (no_local_infections(state_data, vil) and villageExtinctionTime[vil] != numeric_limits<double>::max()){
                        villageExtinctionIntervals[vil].push_back((time-timeAfterBurnIn) - villageExtinctionTime[vil]);
                    }
                }
                //output extinction time
                outputQuant_stream[EXTINCTION_TIME_OUT]<<(time - timeAfterBurnIn)<<endl;
                //output time increments
                for(unsigned int tLength = 0; tLength < timeVec.size(); tLength++){
                    outputQuant_stream[TIME_OUT]<<timeVec[tLength]<<SEP;
                }
                outputQuant_stream[TIME_OUT]<<endl;

                for(unsigned int vil = 0; vil < (unsigned) numVillages; vil++){
                    //output population counts
                    for (unsigned int state = 0; state < NUM_OF_STATE_TYPES; ++state) {
                        for(unsigned int pop = 0; pop < output[state][vil].size(); pop++){ // What is pop? - tjh
                            output_streams[state][vil] << output[state][vil][pop] << SEP;
                        }
                        output_streams[state][vil] << endl;
                    }
                }
                
                //output extinction intervals for each village
                for(unsigned int vil = 0; vil < (unsigned) numVillages; vil++){
                    unsigned int outputVecSize = villageExtinctionIntervals[vil].size();
                    for(unsigned int interval = 0; interval < outputVecSize; interval++){
                        outputVillageExtinctionInterval_stream[vil]<<villageExtinctionIntervals[vil][interval]<<SEP;
                    }
                    outputVillageExtinctionInterval_stream[vil]<<endl;
                }
                    
                //output intervals representing time between paralytic cases
                for(unsigned int interval = 0; interval < timeBetweenPCases.size(); interval++) {
                    outputQuant_stream[TIME_BETWEEN_PCASES_OUT] << timeBetweenPCases[interval];
                    if(interval < timeBetweenPCases.size() - 1) {
                        outputQuant_stream[TIME_BETWEEN_PCASES_OUT]<< SEP;
                    }
                }
                outputQuant_stream[TIME_BETWEEN_PCASES_OUT] << endl;
                
                //output transmission intervals
                for(unsigned int interval = 0; interval < transmissionIntervals.size(); interval++) {
                    outputQuant_stream[TRANSMISSION_INTERVAL_OUT] << transmissionIntervals[interval];
                    if(interval < transmissionIntervals.size() - 1) {
                        outputQuant_stream[TRANSMISSION_INTERVAL_OUT] << SEP;
                    }
                }
                outputQuant_stream[TRANSMISSION_INTERVAL_OUT] << endl;

                //output extinction intervals
                for(unsigned int interval = 0; interval < extinctionIntervals.size(); interval++) {
                    assert(extinctionIntervals[interval] >= 0);
                    outputQuant_stream[EXTINCTION_INTERVAL_OUT] << extinctionIntervals[interval];
                    if(interval < extinctionIntervals.size() - 1){
                        outputQuant_stream[EXTINCTION_INTERVAL_OUT]<< SEP;
                    }
                }
                outputQuant_stream[EXTINCTION_INTERVAL_OUT] << endl;
                
                //clear vectors when done
                timeVec.clear();
                timeBetweenPCases.clear();
                transmissionIntervals.clear();
                extinctionIntervals.clear();
                villageExtinctionTime.clear();
                output = vector<vector<vector<int>>>(NUM_OF_STATE_TYPES, vector<vector<int>>(numVillages));
                villageExtinctionIntervals = vector<vector<double>> (numVillages);
                break;
            }
        }
    }
    assert(eps_intercase_ivls.size() == eps_circ_ivls.size());
    output_results(output_streams, outputVillageExtinctionInterval_stream, outputQuant_stream);
    return 0;
}
