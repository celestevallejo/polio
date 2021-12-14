#include "EventDriven_Sim_nonexp_recovery_waningfn_ENV.hpp"



int main() {
    //using namespace std;
    string output_dir = "/Users/Celeste/Desktop/C++_Polio_IBM/";
    //string output_dir = "./";
    std::ofstream myfile4;
    std::ofstream myfile5;
    std::ofstream myfile6;
    std::ofstream myfile7;
    std::ofstream myfile8;
    /*std::ofstream myfile9;
    std::ofstream myfile10;
    std::ofstream myfile11;*/
    myfile4.open (output_dir + "N=1000_I1vec_waningfn_nonexprec_ateq_minTiter_fixedsheddingfn_test.csv");
    myfile5.open (output_dir + "N=1000_IRvec_waningfn_nonexprec_ateq_minTiter_fixedsheddingfn_test.csv");
    myfile6.open (output_dir + "N=1000_Svec_waningfn_nonexprec_ateq_minTiter_fixedsheddingfn_test.csv");
    myfile7.open (output_dir + "N=1000_Rvec_waningfn_nonexprec_ateq_minTiter_fixedsheddingfn_test.csv");
    myfile8.open (output_dir + "N=1000_time_waningfn_nonexprec_ateq_minTiter_fixedsheddingfn_test.csv");
    const int numSims = 1;
    const int N        = 1000;
    const double DEATH = 0.02;
    const double BETA  =135*2;//multiply by 2
    const double maxRunTime = .2;
    int seed = 1;
    
    //initial conditions
    const int Seq=1;
    const int I1eq=2;
    const int Req=369;
    const int Ireq=628;
    const int I1env=0;
    const int Irenv=0;

    //for counts at end of sim
    int births=0;
    int deaths=0;
    int I1inf=0;
    int Irinf=0;
    int recI1=0;
    int recIr=0;
    int nums=0;
    int numi1=0;
    int numr=0;
    int nump=0;
    int numir=0;
    int numi1e=0;
    int numire=0;
    
    for(int i=0; i<numSims; i++) {
        EventDriven_MassAction_Sim sim(N,BETA,DEATH,Seq,I1eq,Req,Ireq,I1env,Irenv,maxRunTime,seed);
        sim.randomizePopulation();
        sim.runSimulation();
        for(int k=0;k<sim.printVectorI1().size();k++){
        myfile4 << sim.printVectorI1()[k] << " , ";
        myfile5 << sim.printVectorIr()[k] << " , ";
        myfile6 << sim.printVectorS()[k]  << " , ";
        myfile7 << sim.printVectorR()[k]  << " , ";
        myfile8 << sim.printTimeVector()[k]<<" , ";
        /*myfile9 << sim.printAntTit50()[k]<<" , ";
        myfile10 << sim.printAntTit100()[k]<<" , ";
        myfile11 << sim.printAntTit2048()[k]<<" , ";*/
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
        I1inf+=sim.NumI1Inf();
        Irinf+=sim.NumIRInf();
        recI1+=sim.NumRecI1();
        recIr+=sim.NumRecIr();
        nums+=sim.numS();
        numi1+=sim.numI1();
        numr+=sim.numR();
        nump+=sim.numP();
        numir+=sim.numIr();
        numi1e+=sim.numI1E();
        numire+=sim.numIrE();
    }
    myfile4.close();
    myfile5.close();
    myfile6.close();
    myfile7.close();
    myfile8.close();
    /*myfile9.close();
    myfile10.close();
    myfile11.close();*/
    
    cout<<"num births "<<births/(double)numSims<<"\n";
    cout<<"num deaths "<<deaths/(double)numSims<<"\n";
    cout<<"num I1 inf "<<I1inf/(double)numSims<<"\n";
    cout<<"num Ir inf "<<Irinf/(double)numSims<<"\n";
    cout<<"num rec I1 "<<recI1/(double)numSims<<"\n";
    cout<<"num rec Ir "<<recIr/(double)numSims<<"\n";

    cout<<"distribution of pop\n";
    cout<<"num S "<<nums/(double)numSims<<"\n";
    cout<<"num I1 "<<numi1/(double)numSims<<"\n";
    cout<<"num R "<<numr/(double)numSims<<"\n";
    cout<<"num Ir "<<numir/(double)numSims<<"\n";
    cout<<"num I1E "<<numi1e/(double)numSims<<"\n";
    cout<<"num IRE "<<numire/(double)numSims<<"\n";


    return 0;
}
