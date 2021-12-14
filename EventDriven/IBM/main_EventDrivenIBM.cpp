#include "EventDriven_Sim.hpp"
#include <cstdlib>
#include <string>



int main(int argc, char** argv) {
    if (argc != 3) { cerr << "\t\nUsage: ./polio <pop_size> <waning_function>\n \t waning_function: Enter <FAMULARE> for Famulare waning and <TEUNIS> for Teunis waning\n\n"; exit(-1);}
    
    using namespace std;
    string output_dir = "/Users/Celeste/Desktop/C++_Polio_IBM/";
    
    //string output_dir = "../polio_data/";
    //string output_dir = "./";
    ofstream myfile4;
    ofstream myfile5;
    ofstream myfile6;
    ofstream myfile7;
    ofstream myfile8;

    const int numSims  = 1;
    const int N        = atoi(argv[1]);
    const double DEATH = 0.012;//assumes max age of 85 y.o.
    //adding seasonality to transmission term (highest transmission in summer months, lowest in winter months)
    //vector<double> beta(12) = {1000};
    const double BETA  = 1000;//can't set to zero otherwise recovery time is inf -- do we want to change this?
    const double maxRunTime = 5; //units in years
    const int waningImmunity = convertCLArg(argv[2]);
    assert(waningImmunity!=NUM_OF_WANING_IMMUNITY);
    const double propVirusinWater = 1e-9;
    string ext = "N_" + to_string(N)+"_beta_"+to_string(BETA)+"_wan_immune_"+to_string(waningImmunity)+ "_test.csv";
    
    myfile4.open (output_dir + "I1Vec_" + ext);
    myfile6.open (output_dir + "IEnvVec_" + ext);
    myfile5.open (output_dir + "Svec_" + ext);
    myfile7.open (output_dir + "time_" + ext);
    myfile8.open (output_dir + "end_age_dist_valid3_" + ext);

    //for counts at end of sim
    int births=0;
    int deaths=0;
    int I1inf=0;
    int Irinf=0;
    int nums=0;
    int numi1=0;
    int numr=0;
    //int seed = 2;
    
    for(int i=0; i<numSims; i++) {
        EventDriven_MassAction_Sim sim(N,BETA,DEATH,maxRunTime,waningImmunity,propVirusinWater);
        sim.randomizePopulation(2);//this input is number of initially infecteds...make this a command line input variable?
        sim.runSimulation();
        for(unsigned int k=0;k<sim.printVectorIDC().size();k++){
        myfile4 << sim.printVectorIDC()[k] << " , ";
        myfile5 << sim.printVectorNonInf()[k] << " , ";
        myfile6 << sim.printVectorIE()[k]  << " , ";
        myfile7 << sim.printTimeVector()[k]<<" , ";
        }
        for(unsigned int k = 0; k < sim.printAgeDist().size();k++){
            myfile8<<sim.printAgeDist()[k]<<" , ";
        }
        myfile5 << "\n";
        myfile4 << "\n";
        myfile7 << "\n";
        myfile6 << "\n";
        myfile8 << "\n";
        /*myfile9 << "\n";
        myfile10 << "\n";
        myfile11 << "\n";*/
        /*for(int k=0;k<sim.printVectorS().size();k++){
            myfile4 << sim.printTimeVector()[k] << " , ";
            myfile5 << sim.printTimeVector()[k] << " , ";
            myfile6 << sim.printTimeVector()[k]  << " , ";
            myfile7 << sim.printTimeVector()[k]  << " , ";
        }*/
        births+=sim.NumBirths();
        deaths+=sim.NumDeaths();
        I1inf+=sim.NumDCInf();
        Irinf+=sim.NumEInf();
        nums+=sim.numNonInf();
        numi1+=sim.numInfDC();
        numr+=sim.numInfE();
    }
    myfile4.close();
    myfile5.close();
    myfile6.close();
    myfile7.close();
    myfile8.close();
    
    cout<<"num births "<<births/(double)numSims<<"\n";
    cout<<"num deaths "<<deaths/(double)numSims<<"\n";
    cout<<"num DC Inf "<<I1inf/(double)numSims<<"\n";
    cout<<"num E Inf "<<Irinf/(double)numSims<<"\n";

    cout<<"distribution of pop\n";
    cout<<"num Non Inf "<<nums/(double)numSims<<"\n";
    cout<<"num Inf DC "<<numi1/(double)numSims<<"\n";
    cout<<"num Inf E "<<numr/(double)numSims<<"\n";


    return 0;
}
