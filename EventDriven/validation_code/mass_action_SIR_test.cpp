#include "EventDriven_Sim_exp_recovery.hpp"


int main() {
    const int numSims = 1000;
    const int N        = 10000;
    const double GAMMA = 13;
    
    const double KAPPA = .4179;
    const double RHO = .2;
    const double BIRTH = 0;
    const double DEATH = 0;
    const double BETA  = 270;

    std::vector<double> TTE;
    double vecSumm=0.0;
    double firstInfecSum=0.0;
    vector<double> R0;
    double R0sum = 0.0;
    int seed = 0;
    std::array<double,1000> AVG_I1vec={0};
    std::array<double,1000> AVG_Irvec={0};
    int ii=0;
    
    for(int i=0; i<numSims; i++) {
        EventDriven_MassAction_Sim sim(N,BETA,BIRTH,DEATH,GAMMA,RHO,KAPPA);
        sim.randomizePopulation();//randomize with # of infecteds
        sim.runSimulation();
        TTE.push_back(sim.FinalTime());
        vecSumm+=sim.FinalTime();
        sim.printVector(ii);
        ii++;
        R0sum+= sim.NumInfected();
        for(int k = 0;k<sim.printVectorI1().size();k++){
            AVG_I1vec[k]= (AVG_I1vec[k]*i+sim.printVectorI1()[k])/(double)(i+1);
        }
        for(int k = 0;k<sim.printVectorIr().size();k++){
            AVG_Irvec[k]= (AVG_Irvec[k]*i+sim.printVectorIr()[k])/(double)(i+1);
        }
    }
    
    std::ofstream myfile7;
    myfile7.open ("/Users/Celeste/Desktop/C++IBM/N=10000,beta=135,fast_I1vec_valid_time=6_test.csv");
    for(int i=0;i<AVG_I1vec.size();++i){
        myfile7<<AVG_I1vec[i]<<"\n";
    }
    myfile7.close();
    std::ofstream myfile8;
    myfile8.open ("/Users/Celeste/Desktop/C++IBM/N=10000,beta=135,fast_Irvec_valid_time=6_test.csv");
    for(int i=0;i<AVG_Irvec.size();++i){
        myfile8<<AVG_Irvec[i]<<"\n";
    }
    myfile8.close();

    cout<<"TTE "<<vecSumm/numSims<<"\n";

    return 0;
}
