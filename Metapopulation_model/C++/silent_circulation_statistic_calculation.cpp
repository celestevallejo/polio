#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <assert.h>
#include <sstream>
#include <random>

using namespace std;

int main(int argc, char** argv) {
    assert(argc == 3);
    string circ_ivl_filename = argv[1]; //intercase interval file first
    string ext_ivl_filename = argv[2]; //extinction interval file second
    vector<double> intercase_intervals;
    vector<double> extinction_intervals;
    int required_pcases = 1;
    //string dir ="/Users/Celeste/Desktop/polio-small-pop/polio-small-pop/new_data_after_review/";
    string dir = "/Users/Celeste/Desktop/multipatch_model/SC_statistic/1pcase_";
    string extension = "_paper_testing_may_delete.csv";
    string pop ="16x4000_migRate_0.2";
    string sc_output_filename = dir + pop + extension;
    string sc_output_filename_numExt = dir + pop + "_numExtIntLessThanDeltaT" + extension;
    string sc_output_filename_numIntercase = dir + pop + "_numInterCaseIntGretThanDeltaT" + extension;
    string sc_output_filename_totalNumExt = dir + pop + "_totalNumExtInts" + extension;

    //string sc_output_filename = dir + "sc_1pcase_N_15000_numVil_2_det_1_beta_135_fast_multipatch_burnIn_50_reint_0.001.out";
    random_device rd;                       // generates a random real number for the seed
    mt19937 gen(rd());                      // random number generator

    //first read in intercase intervals
    string intercase_line;
    ifstream circ_fh(circ_ivl_filename);
    double val;
    if (circ_fh) {
        while (getline(circ_fh, intercase_line)) {
            vector<double> intercase_times;
            stringstream ss(intercase_line);
            while (ss >> val) {
                intercase_times.push_back(val);
                if (ss.peek() == ','){
                    ss.ignore();
                }
            }
            if (intercase_times.size() > 0) {
                intercase_intervals.insert(intercase_intervals.end(), intercase_times.begin(), intercase_times.end());
            }
        }
    }
    circ_fh.close();
    
    //second read in extinction intervals
    string ext_line;
    ifstream ext_fh(ext_ivl_filename);
    double ext_val;
    if (ext_fh) {
        while (getline(ext_fh, ext_line)) {
            vector<double> ext_times;
            stringstream ss(ext_line);
            while (ss >> ext_val) {
                ext_times.push_back(ext_val);
                if (ss.peek() == ','){
                    ss.ignore();
                }
            }
            if(ext_times.size() > 0){
                extinction_intervals.insert(extinction_intervals.end(),ext_times.begin(),ext_times.end());
            }
        }
            
    }
    ext_fh.close();
    
    vector<double> totalVec = extinction_intervals;
    totalVec.insert(totalVec.end(),intercase_intervals.begin(),intercase_intervals.end());
    
    sort(extinction_intervals.begin(),extinction_intervals.end());
    sort(intercase_intervals.begin(),intercase_intervals.end());
    sort(totalVec.begin(),totalVec.end());//put all times together to get intervals
    

    const int num_ei = extinction_intervals.size();//number of extinction intervals (used in denom)
    const int num_ii = intercase_intervals.size();
    const int num_tv = totalVec.size();

    int last_ecase_idx = 0;
    unsigned int ecase_idx = 0;
    int last_icase_idx = 0;
    unsigned int icase_idx = 0;


    ofstream fo(sc_output_filename);
    ofstream fo_numExt(sc_output_filename_numExt);
    ofstream fo_numIntercase(sc_output_filename_numIntercase);
    ofstream fo_totalExt(sc_output_filename_totalNumExt);
    //calculate SC statistic using 1 - (# of ext<=delta t /(# intercase >= delta t + # ext))
    for (unsigned int interval = 0; interval < num_tv; ++interval) {
        double interval_duration = totalVec[interval];
        //find all extinction intervals <= interval_duration (for numerator)
        for (ecase_idx = last_ecase_idx; ecase_idx < extinction_intervals.size(); ++ecase_idx) {
            if (extinction_intervals[ecase_idx] > interval_duration) {
                break;
            }
        }
        last_ecase_idx = ecase_idx;
        
        //find all intercase intervals < interval_duration (for denominator)
        for (icase_idx = last_icase_idx; icase_idx < intercase_intervals.size(); ++icase_idx) {
            if (intercase_intervals[icase_idx] >= interval_duration) {
                break;
            }
        }
        last_icase_idx = icase_idx;
        double numerator = ecase_idx;
        double denominator = num_ei + num_ii - icase_idx;
        fo << interval_duration << " " << (1 - (numerator/denominator)) << endl;
        fo_numExt << interval_duration<< " "<<numerator<<endl;
        fo_totalExt << interval_duration<<" "<<num_ei<<endl;
        fo_numIntercase<<interval_duration<<" "<<num_ii - icase_idx<<endl;
    }
    fo.close();
    fo_numExt.close();
    fo_totalExt.close();
    fo_numIntercase.close();
}
