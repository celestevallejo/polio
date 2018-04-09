#include "EventDriven_Sim.hpp"
#include <cstdlib>



int main(int argc, char** argv) {
    if (argc != 3) { cerr << "\t\nUsage: ./polio <pop_size> <waning_function>\n \t waning_function: Enter 1 for Famulare waning and 2 for Teunis waning\n\n"; exit(-1);}
    
    using namespace std;
    string output_dir = "/Users/Celeste/Desktop/C++_Polio_IBM/";
    //string output_dir = "../polio_data/";
    //string output_dir = "./";
    ofstream myfile4;
    ofstream myfile5;
    ofstream myfile6;
    ofstream myfile7;
    ofstream myfile8;

    myfile4.open (output_dir + "N=1000_I1vec_waningfn_nonexprec_ateq_minTiter_fixedsheddingfn_test.csv");
    myfile6.open (output_dir + "N=1000_IEnvVec_waningfn_nonexprec_ateq_minTiter_fixedsheddingfn_test.csv");
    myfile5.open (output_dir + "N=1000_Svec_waningfn_nonexprec_ateq_minTiter_fixedsheddingfn_test.csv");
    myfile7.open (output_dir + "N=1000_time_waningfn_nonexprec_ateq_minTiter_fixedsheddingfn_test.csv");
    myfile8.open (output_dir + "end_age_dist.csv");
    const int numSims  = 1000;
    const int N        = atoi(argv[1]);
    const double DEATH = 0.012;//assumes max age of 85 y.o.
    const double BETA  = 1000;//can't set to zero otherwise recovery time is inf -- do we want to change this?
    const double maxRunTime = 100; //units in years
    const int waningImmunity = atoi(argv[2]);

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
        EventDriven_MassAction_Sim sim(N,BETA,DEATH,maxRunTime,waningImmunity);
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