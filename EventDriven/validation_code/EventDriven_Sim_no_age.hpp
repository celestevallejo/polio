//
//  EventDriven_Sim.hpp
//  RtoC++practice
//
//  Created by Celeste on 5/3/17.
//  Copyright © 2017 Celeste. All rights reserved.
//

#ifndef EDMASSIM_H
#define EDMASSIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <queue>
#include <random>
#include <assert.h>
#include <array>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
#include <limits>
#include "EventDriven_parameters.hpp"




using namespace std;

enum EventType {INFECTIOUS_CONTACT, INFECTION_RECOVERY, DEATH, BIRTH, BEGIN_SHEDDING,NUM_OF_EVENT_TYPES};

enum InfectionStatus{I,V,IE,NA,D}; //NA means never infected, D means dead

class Person{
    friend class Event;
private:
    //attributes
    double m_age;
    double m_titerLevel;
    double m_timeAtInfection;//infection means either OPV or WPV
    int m_numInfections;
    int m_numVaccinations;
    int m_numInfectionsEnv;
    InfectionStatus m_infectionStatus; //I, V, or I_E (infection by environment) - cause of most recent infection
    const int m_index;
    double m_timeToRec;
    double m_timeToShed;
    
public:
    //default constructor
    Person(int idx, double age=0, double titerLevel=1.0, double timeAtInfection=numeric_limits<double>::max(), InfectionStatus infStat = NA, int numInf=0,int numVacc=0,int numInfEnv=0):m_age(age),m_titerLevel(titerLevel),m_timeAtInfection(timeAtInfection), m_infectionStatus(infStat),m_numInfections(numInf),m_numVaccinations(numVacc),m_numInfectionsEnv(numInfEnv),m_index(idx){
        
    }
    int getIndex () const {
        return m_index;
    }
    
    void setAge(double age){
        m_age = age;
    }
    void updateAge(double age){
        m_age +=age;
    }
    double deathTime(){
        return maxAge-m_age;
    }
    void reset(){
        m_age = 0.0;
        m_titerLevel = minTiter;
        m_infectionStatus=NA;
        m_timeAtInfection=numeric_limits<double>::max();
        m_numInfections=0.0;
        m_numVaccinations=0.0;
        m_numInfectionsEnv=0.0;
        m_timeToRec=0.0;
        m_timeToShed=0.0;
    }
    double getTimeToRec(){
        return m_timeToRec;
    }
    double getTimeToShed(){
        return m_timeToShed;
    }
    void setTimeToRec(double time){
        m_timeToRec = time;
    }
    void setTimeToShed(double time){
        m_timeToShed = time;
    }
    int getNumInfections(){
        return m_numInfections;
    }
    int getNumInfectionsEnv(){
        return m_numInfectionsEnv;
    }
    int getNumVaccinations(){
        return m_numVaccinations;
    }
    void setTiterLevel(double titerLevel){
        m_titerLevel = std::min(maxTiter,titerLevel);
    }
    void setTimeAtInfection(double t){
        m_timeAtInfection = t;
    }
    void setInfectionStatus(InfectionStatus infectionStatus){
        m_infectionStatus = infectionStatus;
    }
    void setNumInfections(){
        m_numInfections++;
    }
    void setNumVaccinations(){
        m_numVaccinations++;
    }
    void setNumInfectionsEnv(){
        m_numInfectionsEnv++;
    }
    double sheddingThreshold(InfectionStatus infectionStatus){
        if(infectionStatus==I or infectionStatus==IE){//need to change I_E if incorrect
            return WPVrecThresh;
        }
        else{
            return OPVrecThresh;
        }
    }
    double getAge(){
        return m_age;
    }
    double getTiterLevel(){
        return m_titerLevel;
    }
    double getTimeAtInfection(){
        return m_timeAtInfection;
    }
    InfectionStatus getInfectionStatus(){
        return m_infectionStatus;
    }
    
    double probInfGivenDose(InfectionStatus infstat){//dose will vary based on if it is WPV or OPV, used mean of all shape parameters
        if(m_infectionStatus!=D){
            if(infstat==I){
                return 1-pow((1.0+(infDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
            }
            else if(infstat == V){
                return 1-pow((1.0+(vaccDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
            }
            else{//must be infstat=="I_E"
                return 1-pow((1.0+(envDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
            }
        }
        else{
            return 0.0;
        }
    }
    void waning(double t){
        const double tnew = convertToMonths(t-m_timeAtInfection); //time needs to be in months post infection
        if(tnew>=1){//only wanes after one month post infection
            m_titerLevel= max(minTiter, m_titerLevel*pow((tnew),-waningLambda));
        }
        return;
    }
    double peakShedding(){//age needs to be converted to months
        if(m_age<newBorn){
            return Smax;
        }
        else{
            return ((Smax-Smin)*exp((newBorn-(convertToMonths(m_age)))/tau)+Smin);
        }
    }
    
    double probShedding(double t){
        if(m_infectionStatus!=D){
            const double tnew = convertToDays(t - m_timeAtInfection);//***t needs to be time since infection (days)
            if(m_infectionStatus==I or m_infectionStatus==IE){
                return .5*erfc((log(tnew)-(log(muWPV)-log(deltaShedding)*log(m_titerLevel)))/sqrt(2.0)*log(sigmaWPV));
            }
            else{//m_infectionStatus== V
                return .5*erfc((log(tnew)-(log(muOPV)-log(deltaShedding)*log(m_titerLevel)))/sqrt(2.0)*log(sigmaOPV));
            }
        }
        else{
            return 0.0;
        }
    }
    
    double stoolViralLoad(double t){
        if(m_infectionStatus!=D){
            const double tnew = convertToDays(t-m_timeAtInfection);//***t needs to be time since infection (days)
            return max(pow(10.0,2.6),pow(10,((1-k*log2(m_titerLevel))*log(peakShedding())))*(exp(eta-(pow(nu,2)/(double)2)-(pow(log(tnew)-eta,2)/(double)2*pow(nu+xsi*log(tnew),2)))/tnew));
        //**units are in TCID50/g
        }
        else{
            return 0.0;
        }
    }
    double convertToDays(double t){
        return t*365;
    }
    double convertToMonths(double t){
        return t*12;
    }
    

};

class Event {
    friend class Person;
public:
    double time;
    EventType type;
    Person* person;
    
    
    // Event(const Event& o) {time = o.time; type=o.type;}//*** this is a copy constructor
    Event(double t, EventType e, Person* p):time(t),type(e),person(p){} // *** this is a constructor
    // Event& operator=(const Event& o) {
    //     if(this != &o){
    //         time = o.time;
    //         type=o.type;
    //     }
    //     return *this;
    // } //*** this is copy assignment operator
    
};

class compTime {
public:
    //two in case we use pointers at some point
    bool operator() (const Event* lhs, const Event* rhs) const {
        
        return (lhs->time>rhs->time); //(->) is a reference operator for pointers
    }
    
    bool operator() (const Event& lhs, const Event& rhs) const {
        return (lhs.time>rhs.time);
    }
};


class EventDriven_MassAction_Sim {
public:
    // constructor
    EventDriven_MassAction_Sim(const int n, const double beta, const double birth, const double death, int seed = (random_device())()): rng(seed), N(n), BETA(beta), BIRTH_RATE(birth), DEATH_RATE(death), unif_real(0.0,1.0),unif_int(0,n-2),people_counter(0){
        people = vector<Person*>(n);
        for (Person* &p: people) p = new Person(people_counter++);
        Now=0.0;
        ii=0;
        k=0;
        IIvec = {0};
        I1vec = {0};
        queueCount = 1;
        event_counter = vector<int>(NUM_OF_EVENT_TYPES, 0);
        numBirths=0;
        numDeaths=0;
        numI1inf=0;
        numIrinf=0;
        numRec=0;
        numWane=0;
        avgInfI1rate=0;
        numRemainingInfectives=0;

    }
    
    ~EventDriven_MassAction_Sim() {
        for (Person* &p: people) delete p;
    }
    
    
    const double BETA;
    const double BIRTH_RATE;
    const double DEATH_RATE;
    
    exponential_distribution<double> exp_beta;
    exponential_distribution<double> exp_death;
    exponential_distribution<double> exp_vacc;
    uniform_real_distribution<double> unif_real;
    uniform_int_distribution<int> unif_int;
    mt19937 rng;
    
    //containers to keep track of various pieces of information
    vector<Person*> people;
    priority_queue<Event, vector<Event>, compTime > EventQ;
    array<double,1>previousTime;
    array<double,1>finalTime;
    vector<double> timeOfParalyticCase;
    vector<double> avgAgeOfFirstInfect;
    unordered_set<Person*> sheddingPeople;
    
    
    int people_counter;
    double Now;
    int numInfected;
    double counter=0.0;
    int ii;
    vector<int> event_counter;
    const int N;
    array<double,100000> IIvec;
    array<double,100000> I1vec;
    int k;
    int queueCount;
    int numBirths;
    int numDeaths;
    int numI1inf;
    int numIrinf;
    int numRec;
    int numWane;
    double avgInfI1rate;
    unsigned long numRemainingInfectives;


    
    
    
    void runSimulation(){
        IIvec[0] = 0;
        while(nextEvent()>=0){
            if(nextEvent()==0){
                numRemainingInfectives = sheddingPeople.size();
                TTE=Now;
                break;
            }
            continue;
        }
    }
    int NumInfected(){
        return numInfected;
    }
    
    void randomizePopulation(int infected){
        //change initial conditions depending on how we want population to start
        sheddingPeople.clear();
        for(Person* p: people) {

            death(p); //set when each individual will die

            p->setInfectionStatus(NA); //sets most recent cause of infection

            p->setTiterLevel(maxTiter);

            p->setTimeAtInfection(numeric_limits<double>::max());//set time at last infection..if infection status is NA then time at infection is largest int
        }
        for(int i=0;i<infected;i++){
            infect(people[i]);//infects first k people in vector
        }
    }
    
    void printPopultion(){
        for(Person* p:people){
            cout<<"Person attributes\n";
            cout<<"infection status "<<p->getInfectionStatus()<<"\n";
            cout<<"titer level "<<p->getTiterLevel()<<"\n";
        }
    }
    
    double FinalTime(){
        return TTE;
    }
    
    array<double,100000> printVectorII(){
        return IIvec;
    }
    array<double,100000> printVectorI1(){
        return I1vec;
    }
    int NumBirths(){
        return numBirths;
    }
    int NumDeaths(){
        return numDeaths;
    }
    int NumI1Inf(){
        return numI1inf;
    }
    int NumRec(){
        return numRec;
    }
    unsigned long NumRemainingInfectives(){
        return numRemainingInfectives;
    }
    
    
    /*void printVector(int ii){
        if(ii==0){
        std::ofstream myfile7;
        myfile7.open ("/Users/Celeste/Desktop/C++IBM/N=3500,beta=135,fast_I1vec_valid_time=6_test_insideh_first.csv");
        for(int i=0;i<I1vec.size();++i){
            myfile7<<I1vec[i]<<"\n";
        }
        myfile7.close();
        }
        else{
            std::ofstream myfile7;
            myfile7.open ("/Users/Celeste/Desktop/C++IBM/N=3500,beta=135,fast_I1vec_valid_time=6_test_insideh_second.csv");
            for(int i=0;i<I1vec.size();++i){
                myfile7<<I1vec[i]<<"\n";
            }
            myfile7.close();
        }
        
    }*/
    
    
    void infect(Person* p) {
        assert(p->getInfectionStatus()!=D);
        assert(sheddingPeople.count(p)==0);
        cout<<"in Infect loop\n";
        cout<<"number of infections before update "<<numI1inf<<"\n";
        numI1inf++;
        cout<<"number of infections after update "<<numI1inf<<"\n";
        p->setNumInfections();
        p->setInfectionStatus(I);
        p->setTimeAtInfection(Now);
        p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold

        //time to shedding (able to cause infections)
        double sheddingTime = Now + (1/(double)365);
        cout<<"index of shedder "<<p->getIndex()<<"\n";
        cout<<"shedding time "<<sheddingTime<<"\n";
        EventQ.emplace(sheddingTime,BEGIN_SHEDDING,p);//assume shedding begins a day later--can't be instantaneous
        event_counter[BEGIN_SHEDDING]++;
        p->setTimeToShed(sheddingTime);

        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;

        while (p->probShedding(Tc)>WPVrecThresh) {     // does contact occur before recovery?
            if(Tc > sheddingTime){//only make infectious contacts after shedding begins
                EventQ.emplace(Tc,INFECTIOUS_CONTACT,p);
                event_counter[INFECTIOUS_CONTACT]++;
            }
            Tc += exp_beta(rng);
        }
        
        //time to recovery
        double Tr = Tc;//once the contact time gets late enough then it is a recovery time (trying to make it not look like a bug)
        EventQ.emplace(Tr,INFECTION_RECOVERY,p);
        event_counter[INFECTION_RECOVERY]++;
        p->setTimeToRec(Tr);
        return;
    }

    void death(Person* p){
        assert(p->getInfectionStatus()!=D);
        exponential_distribution<double> exp_death(DEATH_RATE);
        double Td = exp_death(rng) + Now;
        EventQ.emplace(Td,DEATH,p);
        event_counter[DEATH]++;
        return;
    }
    void birth(Person* p){
        assert(p->getInfectionStatus()==D);
        exponential_distribution<double> exp_birth(BIRTH_RATE);
        double Tb = exp_birth(rng) + Now;
        EventQ.emplace(Tb,BIRTH,p);
        event_counter[BIRTH]++;
        return;
    }
    
    
    int nextEvent() {
        if(((sheddingPeople.size()==0) or Now > 5) and Now > (1/(double)365)){
            return 0;
        }
        cout<<"number of infecteds "<<sheddingPeople.size()<<"\n";
        Event event = EventQ.top();
        Now = event.time;
        cout<<"event "<<event.type<<"\n";
        if(k%30==0){
            int firstInfCount=0;
            for(Person* p: people){
                if(p->getNumInfections()==1){
                    firstInfCount++;
                }
            }
            I1vec[queueCount] = firstInfCount;
            IIvec[queueCount] = sheddingPeople.size();
            queueCount++;
        }
        k++;
         //if(Now>counter){
             cout<<"Loop "<< ii<<"\n";
             cout<<"Now "<<Now<<"\n";
             cout<<" queue size "<<EventQ.size()<<"\n";
             cout << "\tINFECTIOUS_CONTACT: " << event_counter[INFECTIOUS_CONTACT] << endl;
             cout << "\tBEGIN_SHEDDING: " << event_counter[BEGIN_SHEDDING] << endl;
             cout << "\tINFECTION_RECOVERY: " << event_counter[INFECTION_RECOVERY] << endl;
             cout << "\tDEATH: " << event_counter[DEATH] << endl;
             cout << "\tBIRTH: " << event_counter[BIRTH] << endl;
             //counter+=.1;
             ii++;
         //}
        Person* individual = event.person;
        if(event.type==INFECTIOUS_CONTACT){
            if(individual->getInfectionStatus()==I){
                int contact_idx = unif_int(rng);
                contact_idx = contact_idx >= individual->getIndex() ? contact_idx++ : contact_idx;
                Person* contact = people[contact_idx];
                cout<<"index of contact "<<contact->getIndex()<<"\n";
                if(contact->getInfectionStatus()!=D and contact->getTimeToShed() < Now and sheddingPeople.count(contact)==0){
                    cout<<"time to shed "<<contact->getTimeToShed()<<"\n";
                    assert(contact->getInfectionStatus()!=D);
                    assert(sheddingPeople.count(contact)==0);
                    //first wane the immunity of the person that is contacted
                    if(contact->getTimeAtInfection()!=numeric_limits<double>::max()){
                        contact->waning(Now);
                    }
                    double r1 = unif_real(rng);
                    double r2 = unif_real(rng);
                    if(r1 < 0.5){//Assume 50% probability a contact is considered a dose?
                        if(r2 < contact->probInfGivenDose(I)){
                            cout<<"infection\n";
                            assert(sheddingPeople.count(contact)==0);
                            infect(contact);
                        }
                    }
                }
            }
        }
        else if(event.type==BEGIN_SHEDDING){
            if(individual->getInfectionStatus()==I){
                cout<<"index of shedder "<<individual->getIndex()<<"\n";
                cout<<"time of shedding "<<individual->getTimeToShed()<<"\n";
                cout<<"Now "<<Now<<"\n";
                assert(sheddingPeople.count(individual)==0);
                sheddingPeople.insert(individual);
                cout<<"Now "<<Now<<"\n";
                cout<<"time to shed "<<individual->getTimeToShed()<<"\n";
            }
        }
        else if (event.type == INFECTION_RECOVERY) {
            if(individual->getInfectionStatus()==I){
                assert(individual->getInfectionStatus()!=D);
                TTE = Now;
                sheddingPeople.erase(individual);
                numRec++;
                //numI1inf--;
                cout<<"Now "<<Now <<"\n";
                cout<<"time to recover "<<individual->getTimeToRec()<<"\n";
            }
        }
        else if(event.type==DEATH){
            numDeaths++;
            assert(individual->getInfectionStatus()!=D);
            if(sheddingPeople.count(individual)>0){
                numRec++;
                //numI1inf--;
                cout<<"infected died!!!!!!!!!!!!!!\n";
                cout<<"!!!!!!!!!!!!!!!!!\n";
                cout<<"!!!!!!!!!!!!!!!!!!!\n";
            }
            sheddingPeople.erase(individual);//removes if they are contained in set
            individual->setInfectionStatus(D);
            birth(individual);
        }
        else if(event.type==BIRTH){
            assert(individual->getInfectionStatus()==D);
            individual->reset();
            numBirths++;
            death(individual);

        }
        cout<<"num I1 inf variable "<<numI1inf<<"\n";
        cout<<"num rec variable "<<numRec<<"\n";
        EventQ.pop();
        event_counter[event.type]--;
        cout<<"num I inf "<<sheddingPeople.size()<<"\n";
        cout<<"num recoveries in queue "<<event_counter[INFECTION_RECOVERY]<<"\n";
        return 1;
    }
    
};
#endif

