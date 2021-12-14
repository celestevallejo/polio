//
//  EventDriven_Sim_exp_recovery_waningfn.hpp
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

enum EventType {INFECTIOUS_CONTACT, INFECTION_RECOVERY, DEATH, BIRTH, BEGIN_SHEDDING, CHECK_ENVIRONMENT, ENVIRONMENT_CONTACT, NUM_OF_EVENT_TYPES};

enum InfectionStatus{S, I1, R, IR, I1_E, IR_E};

class Person{
    friend class Event;
private:
    //attributes
    double m_initialTiterLevel;//need initial titer level after infection for waning fn
    double m_titerLevel;
    double m_timeAtInfection;
    int m_numInfections; //number of person-person infections
    int m_numInfectionsEnv; //number of infections due to env
    InfectionStatus m_infectionStatus; //designates compartment -- eventually will relax to be last cause of infection
    double m_timeToRec;
    double m_timeToShed;
    const int m_index;
    
public:
    //default constructor
    Person(int idx, double initialTiterLevel = 1.0, double titerLevel=1.0, double timeAtInfection=numeric_limits<double>::max(), InfectionStatus infStat = S, int numInf=0,int numInfEnv=0, double timeToShed = numeric_limits<double>::max(), double timeToRec = numeric_limits<double>::max()):m_initialTiterLevel(initialTiterLevel), m_titerLevel(titerLevel),m_timeAtInfection(timeAtInfection), m_infectionStatus(infStat),m_numInfections(numInf),m_numInfectionsEnv(numInfEnv),m_index(idx),m_timeToRec(timeToRec),m_timeToShed(timeToShed){
        
    }
    int getIndex () const {
        return m_index;
    }
    double getTimeToShed(){
        return m_timeToShed;
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
    void setTiterLevel(double titerLevel){
        m_titerLevel = std::min(maxTiter,titerLevel);
    }
    void setInitialTiterLevel(double titerLevel){
        m_initialTiterLevel = std::min(maxTiter, titerLevel);
        m_titerLevel = m_initialTiterLevel; //set equal for the purposes of further calculation
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
    void setNumInfectionsEnv(){
        m_numInfectionsEnv++;
    }
    double getTiterLevel(){
        return m_titerLevel;
    }
    double getInitialTiterLevel(){
        return m_initialTiterLevel;
    }
    double getTimeAtInfection(){
        return m_timeAtInfection;
    }
    InfectionStatus getInfectionStatus(){
        return m_infectionStatus;
    }
    
    void waning(double t){
        const double tnew = convertToMonths(abs(t-m_timeAtInfection)); //time needs to be in months post infection

        if(tnew>=1){//only wanes after one month post infection
            m_titerLevel= max(minTiter, m_initialTiterLevel*pow((tnew),-waningLambda));
        }
        return;
    }
    double probInfGivenDose(InfectionStatus infstat){//dose will vary based on if it is WPV direct contact or env, used mean of all shape parameters
        if(infstat==I1 or IR){
            return 1-pow((1.0+(infDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
        }
        else{//must be infstat=="I1_E" or "IR_E"
            return 1-pow((1.0+(envDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
        }
        
    }
    double probShedding(double t){
        const double tnew = convertToDays(t - m_timeAtInfection);//***t needs to be time since infection (days)
        return .5*erfc((log(tnew)-(log(muWPV)-log(deltaShedding)*log(m_titerLevel)))/(sqrt(2.0)*log(sigmaWPV)));
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
    
    Event(double t, EventType e, Person* p):time(t),type(e),person(p){}
    //Event keeps track of time of event, type of event, and person to whom event occurs
    
};

class compTime {
public:
    //two in case we use pointers at some point
    bool operator() (const Event* lhs, const Event* rhs) const {
        
        return (lhs->time>rhs->time);
    }
    
    bool operator() (const Event& lhs, const Event& rhs) const {
        return (lhs.time>rhs.time);
    }
};


class EventDriven_MassAction_Sim {
public:
    // constructor
    EventDriven_MassAction_Sim(const int n, const double beta, const double death, const int seq, const int i1eq, const int req, const int ireq, const int i1env, const int irenv, int seed = (random_device())()): rng(seed), N(n), BETA(beta), DEATH_RATE(death), Ssum(seq), I1sum(i1eq), Rsum(req), Irsum(ireq), I1Esum(i1env), IrEsum(irenv), unif_real(0.0,1.0),unif_int(0,n-2),people_counter(0){
        people = vector<Person*>(n);
        for (Person* &p: people) p = new Person(people_counter++);
        Now=0.0;//simulation time
        event_counter = vector<int>(NUM_OF_EVENT_TYPES, 0);
        I1vec = {0};
        Irvec = {0};
        Svec = {0};
        Rvec = {0};
        I1Evec={0};
        IrEvec={0};
        timevec={0};
        antTit50={0};
        antTit100={0};
        antTit2048={0};
        queueCount = 1;
        numBirths=0;
        numDeaths=0;
        numI1inf=0;
        numIrinf=0;
        numRecI1=0;
        numRecIr=0;
        I1vec[0] = I1sum;
        Irvec[0] = Irsum;
        Svec[0] = Ssum;
        Rvec[0] = Rsum;
        I1Evec[0]=I1Esum;
        IrEvec[0]=IrEsum;
        delta=0;

    }
    
    ~EventDriven_MassAction_Sim() {
        for (Person* &p: people) delete p;
    }
    
    
    const double BETA;
    const double DEATH_RATE;
    double virusCon;
    double environment;

    
    exponential_distribution<double> exp_beta;
    uniform_real_distribution<double> unif_real;
    uniform_int_distribution<int> unif_int;
    mt19937 rng;
    
    //containers to keep track of various pieces of information
    vector<Person*> people;
    priority_queue<Event, vector<Event>, compTime > EventQ;
    unordered_set<Person*> sheddingPeople;
    
    
    int people_counter;
    double Now;
    double previousTime; //need for time step calculation
    double counter=0.0;
    int Ssum;
    int I1sum;
    int Rsum;
    int Psum;
    int Irsum;
    int I1Esum;
    int IrEsum;
    const int N;
    vector<int> event_counter;
    vector<double> I1vec;
    vector<double> Irvec;
    vector<double> Svec;
    vector<double> Rvec;
    vector<double> I1Evec;
    vector<double> IrEvec;
    vector<double>timevec;
    vector<double> antTit50;
    vector<double> antTit100;
    vector<double> antTit2048;
    int queueCount;
    int numBirths;
    int numDeaths;
    int numI1inf;
    int numIrinf;
    int numRecI1;
    int numRecIr;
    double delta;


    
    
    
    void runSimulation(){
        double maxRunTime = 1.0;//units in years
        cout<<"population at start: "<<"\n";
        cout<<" num S "<<Ssum<<"\n";
        cout<<"num I1 "<<I1sum<<"\n";
        cout<<" num R "<<Rsum<<"\n";
        cout<<" num Ir "<<Irsum<<"\n";
        while(nextEvent() and EventQ.top().time < maxRunTime) {
        }
    }
    
    void randomizePopulation(){
        int j=0;
        //discretize antibody titer level
        int TiterLevel50=0; //level 1-50
        int TiterLevel100=0;//level 51-100
        int TiterLevel2048=0;//level 101+
        for(Person* p: people) {
            if(j<Ssum){
                //add invididual to compartment
                p->setInfectionStatus(S);
                
                //set initial titer level
                //uniform_real_distribution<double>titer(minTiter,maxTiter);
                p->setInitialTiterLevel(minTiter);
                
                //set time at infection (used to calculate time since infection)
                p->setTimeAtInfection(numeric_limits<double>::max());
                j++;
            }
            else if(j<(Ssum + I1sum)){
                //set initial titer level
                //uniform_real_distribution<double>titer(minTiter,maxTiter);
                p->setInitialTiterLevel(minTiter);
                
                //add individual to compartment
                infectByDirectContact(p);
                j++;
            }
            else if(j<(Ssum + I1sum + Rsum)){
                //add individual to compartment
                p->setInfectionStatus(R);
                
                //set initial titer level
                //uniform_real_distribution<double>titer(minTiter,maxTiter);
                p->setInitialTiterLevel(minTiter);
                
                //set time at infection (used to calculate time since infection)
                p->setTimeAtInfection(numeric_limits<double>::max());
                j++;
            }
            else{
                //since these are reinfected individuals need to set numinfections>1
                p->setNumInfections();
                
                //set initial titer level
                //uniform_real_distribution<double>titer(minTiter,maxTiter);
                p->setInitialTiterLevel(minTiter);
                infectByDirectContact(p);
                j++;
            }
            death(p); //set when each individual will die
            
            //set environment contact--occurs daily
            //random so that at initialization everyone isn't contacting at exact same time
            double envCon = unif_real(rng)/(double)365;
            //EventQ.emplace(envCon,ENVIRONMENT_CONTACT,p);
            //event_counter[ENVIRONMENT_CONTACT]++;
            
            //antibody titer level distribution
            if(p->getTiterLevel()<=50){
                TiterLevel50++;
            }
            else if(p->getTiterLevel()<=100){
                TiterLevel100++;
            }
            else{
                TiterLevel2048++;
            }
        }
        
        //Environmental Surveillance occurs monthly
        double chkEnv = unif_real(rng)/(double)12;
        //EventQ.emplace(chkEnv,CHECK_ENVIRONMENT,nullptr);
        //event_counter[CHECK_ENVIRONMENT]++;
        
        antTit50[0]=TiterLevel50;
        antTit100[0]=TiterLevel100;
        antTit2048[0]=TiterLevel2048;
    }
    
    //these functions are used to get data at end of sim
    vector<double> printVectorI1(){
        if(I1vec.size()<100){
            I1vec.resize(100,I1vec.back());
        }
        return I1vec;
    }
    vector<double> printVectorS(){
        if(Svec.size()<100){
            Svec.resize(100,Svec.back());
        }
        return Svec;
    }
    vector<double> printVectorR(){
        if(Rvec.size()<100){
            Rvec.resize(100,Rvec.back());
        }
        return Rvec;
    }
    vector<double> printVectorIr(){
        if(Irvec.size()<100){
            Irvec.resize(100,Irvec.back());
        }
        return Irvec;
    }
    vector<double> printTimeVector(){
        if(timevec.size()<100){
            timevec.resize(100,timevec.back());
        }
        return timevec;
    }
    vector<double> printAntTit50(){
        return antTit50;
    }
    vector<double> printAntTit100(){
        return antTit100;
    }
    vector<double> printAntTit2048(){
        return antTit2048;
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
    int NumIRInf(){
        return numIrinf;
    }
    int NumRecI1(){
        return numRecI1;
    }
    int NumRecIr(){
        return numRecIr;
    }
    int numS(){
        return Ssum;
    }
    int numI1(){
        return I1sum;
    }
    int numR(){
        return Rsum;
    }
    int numP(){
        return Psum;
    }
    int numIr(){
        return Irsum;
    }
    int numI1E(){
        return I1Esum;
    }
    int numIRE(){
        return IrEsum;
    }
    
    void infectByDirectContact(Person* p) {
        p->setNumInfections();
        if(p->getNumInfections()==1){
            p->setInfectionStatus(I1);
        }
        else{
            p->setInfectionStatus(IR);
        }
        p->setTimeAtInfection(Now);
        p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
        
        //time to shedding (able to cause infections)
        /*double sheddingTime = Now + (1/(double)365);
        EventQ.emplace(sheddingTime, BEGIN_SHEDDING,p); //assume shedding begins a day later -- can't be instantaneous
        p->setTimeToShed(sheddingTime);*/
        
        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;

        while (p->probShedding(Tc)>WPVrecThresh) {
            //if(Tc > sheddingTime){//only make infectious contacts after shedding begins
                EventQ.emplace(Tc,INFECTIOUS_CONTACT,p);
                //event_counter[INFECTIOUS_CONTACT_I1]++;
            //}
            Tc += exp_beta(rng);
        }
        
        //time to recovery
        double Tr = Tc;//once the contact time is late enough such that the probability of shedding is below WPVrecThresh then it is a recovery time (trying to make it not look like a bug)
        EventQ.emplace(Tr,INFECTION_RECOVERY,p);
        //event_counter[INFECTION_RECOVERY_I1]++;
        return;
    }
    void infectByEnvironment(Person* p) {
        p->setNumInfections();
        if(p->getNumInfections()==1){
            p->setInfectionStatus(I1_E);
        }
        else{
            p->setInfectionStatus(IR_E);
        }
        p->setTimeAtInfection(Now);
        p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
        
        //time to shedding (able to cause infections)
        double sheddingTime = Now + (1/(double)365);
        EventQ.emplace(sheddingTime, BEGIN_SHEDDING,p); //assume shedding begins a day later -- can't be instantaneous
        p->setTimeToShed(sheddingTime);
        
        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;
        
        while (p->probShedding(Tc)>WPVrecThresh) {
            if(Tc > sheddingTime){//only make infectious contacts after shedding begins
                EventQ.emplace(Tc,INFECTIOUS_CONTACT,p);
                //event_counter[INFECTIOUS_CONTACT_I1]++;
            }
            Tc += exp_beta(rng);
        }
        
        //time to recovery
        double Tr = Tc;//once the contact time is late enough such that the probability of shedding is below WPVrecThresh then it is a recovery time (trying to make it not look like a bug)
        EventQ.emplace(Tr,INFECTION_RECOVERY,p);
        //event_counter[INFECTION_RECOVERY_I1]++;
        return;
    }

    void death(Person* p){
        exponential_distribution<double> exp_death(DEATH_RATE);
        double Td = exp_death(rng) + Now;
        EventQ.emplace(Td,DEATH,p);
        event_counter[DEATH]++;
        return;
    }
    
    void environmentalSurveillance(){
        if(virusCon>detectionRate){
            cout<<"detected pathogen in water\n";
        }
        else{
            cout<<"did not detect pathogen in water\n";
        }
        //after each ES event generate time to next one--monthly
        EventQ.emplace(Now+(1/(double)12),CHECK_ENVIRONMENT,nullptr);
        event_counter[CHECK_ENVIRONMENT]++;
        return;
    }
    
    
    int nextEvent() {
        Event event = EventQ.top();
        EventQ.pop();
        Now = event.time;
        if(Now-delta>.01){
            //count population every 0.01 time steps
            I1vec.push_back(I1sum);
            Irvec.push_back(Irsum);
            Svec.push_back(Ssum);
            Rvec.push_back(Rsum);
            timevec.push_back(Now);
            int TiterLevel50 = 0;
            int TiterLevel100 = 0;
            int TiterLevel2048 = 0;
            for(Person* p: people){
                if(p->getTiterLevel()<=50){
                    TiterLevel50++;
                }
                else if(p->getTiterLevel()<=100){
                    TiterLevel100++;
                }
                else{
                    TiterLevel2048++;
                }
            }
            antTit50.push_back(TiterLevel50);
            antTit100.push_back(TiterLevel100);
            antTit2048.push_back(TiterLevel2048);
            delta=Now;
        }
        //left this in case we want to look at number of events in queue
         /*if(Now>counter){
             cout<<"Loop "<< ii<<"\n";
             cout<<"Now "<<Now<<"\n";
             cout<<" queue size "<<EventQ.size()<<"\n";
             cout << "\tINFECTIOUS_CONTACT: " << event_counter[INFECTIOUS_CONTACT_I1] << endl;
             cout << "\tINFECTIOUS_CONTACT_IR: " << event_counter[INFECTIOUS_CONTACT_IR] << endl;
             cout << "\tINFECTION_RECOVERY: " << event_counter[INFECTION_RECOVERY_I1] << endl;
             cout << "\tINFECTION_RECOVERY_IR: " << event_counter[INFECTION_RECOVERY_IR] << endl;
             cout << "\tDEATH: " << event_counter[DEATH] << endl;
             counter+=.1;
             ii++;
         }*/
        
        //calculate concentration of virus particles in environment
        double timeStep = Now - previousTime;
        if(timeStep > 0){
            environment+=exp(-inactivationRate*timeStep*365);
        }
        
        Person* individual = event.person;
        if(event.type==INFECTIOUS_CONTACT){
            if(individual->getInfectionStatus()==I1 or individual->getInfectionStatus()==IR){
            //if(sheddingPeople.count(individual)>0){//check this event still applicable
                int contact_idx = unif_int(rng);
                contact_idx = contact_idx >= individual->getIndex() ? contact_idx++ : contact_idx;
                Person* contact = people[contact_idx];
                if(contact->getInfectionStatus()==S or contact->getInfectionStatus()==R){
                //if(contact->getTimeToShed()< Now and sheddingPeople.count(contact)==0){//no simultaneous infections
                    
                    //wane immunity
                    if(contact->getTimeAtInfection()!=numeric_limits<double>::max()){//can't wane if never infected according to Famulare fn defn
                        contact->waning(Now);
                    }
                    
                    //check susceptibility
                    double r2 = unif_real(rng);
                    if(r2<contact->probInfGivenDose(I1)){
                        if(contact->getInfectionStatus()==S){
                            Ssum--;
                            I1sum++;
                            infectByDirectContact(contact);
                        }
                        else{
                            Rsum--;
                            Irsum++;
                            infectByDirectContact(contact);
                        }
                        if(individual->getInfectionStatus()==I1){
                            numI1inf++;
                        }
                        else{
                            numIrinf++;
                        }
                    
                    }
                }
            }
        }
        else if(event.type==ENVIRONMENT_CONTACT){
            //define a dose of virus from environment
            virusCon = environment; //assumes all virus particles shed end up in water source
            envDose = 2*virusCon; //2 is num L water drank per day
            
            //check susceptibility
            double r2 = unif_real(rng);
            if(individual->probInfGivenDose(I1_E)){
                infectByEnvironment(individual);
            }
            
            //update environment with additions
            if(sheddingPeople.count(individual)>0){
                environment+=gramsFeces*100;//just a place holder since real function to determine viral load depends on age
            }
            
            //set next time to contact environment
            EventQ.emplace(Now + (1/(double)365),ENVIRONMENT_CONTACT,individual);
            event_counter[ENVIRONMENT_CONTACT]++;
        }
        else if(event.type==BEGIN_SHEDDING){
            if(individual->getInfectionStatus()==I1 or individual->getInfectionStatus()==IR or individual->getInfectionStatus()==I1_E or individual->getInfectionStatus()==IR_E){//check event still applicable
                sheddingPeople.insert(individual);
            }
        }
        else if (event.type==INFECTION_RECOVERY) {
            if(individual->getInfectionStatus()==IR){
                Irsum--;
                Rsum++;
                individual->setInfectionStatus(R);
                numRecIr++;
            }
            else if(individual->getInfectionStatus()==I1){
                I1sum--;
                Rsum++;
                individual->setInfectionStatus(R);
                numRecI1++;
            }
            else if(individual->getInfectionStatus()==I1_E){
                I1Esum--;
                Rsum++;
                individual->setInfectionStatus(R);
            }
            else if(individual->getInfectionStatus()==IR_E){
                IrEsum--;
                Rsum++;
                individual->setInfectionStatus(R);
            }
        }
        else if(event.type==DEATH){
            //death event generates simultaneous birth event
            //all individuals born as naive susceptible (S)
            numDeaths++;
            numBirths++;
            switch (individual->getInfectionStatus()) {
                case S:
                {
                    death(individual);
                    individual->setInitialTiterLevel(minTiter);
                    break;
                }
                case I1:
                {
                    I1sum--;
                    Ssum++;
                    individual->setInfectionStatus(S);
                    death(individual);
                    individual->setInitialTiterLevel(minTiter);
                    break;
                }
                case R:
                {
                    Rsum--;
                    Ssum++;
                    individual->setInfectionStatus(S);
                    death(individual);
                    individual->setInitialTiterLevel(minTiter);
                    break;
                }
                case IR:
                {
                    Irsum--;
                    Ssum++;
                    individual->setInfectionStatus(S);
                    death(individual);
                    individual->setInitialTiterLevel(minTiter);
                    break;
                }
                case I1_E:
                {
                    I1Esum--;
                    Ssum++;
                    individual->setInfectionStatus(S);
                    death(individual);
                    individual->setInitialTiterLevel(minTiter);
                    break;
                }
                case IR_E:
                {
                    IrEsum--;
                    Ssum++;
                    individual->setInfectionStatus(S);
                    death(individual);
                    individual->setInitialTiterLevel(minTiter);
                    break;
                }
            }
        }
        if(I1sum<0 or Irsum<0 or Ssum<0 or Rsum<0){
            cout<<"last event "<<event.type<<"\n";
            cout<<"num S "<<Ssum<<"\n";
            cout<<"num I1 "<<I1sum<<"\n";
            cout<<"num R "<<Rsum<<"\n";
            cout<<"num Ir "<<Irsum<<"\n";
            assert(I1sum>=0);
            assert(Irsum>=0);
            assert(Ssum>=0);
            assert(Rsum>=0);
        }
        event_counter[event.type]--;
        
        return 1;
    }
    
};
#endif

