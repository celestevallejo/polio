//
//  EventDriven_Sim_nonexp_recovery_waningfn_ENV.hpp
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

enum EventType {INFECTIOUS_CONTACT, INFECTION_RECOVERY, DEATH, BEGIN_SHEDDING, CHECK_ENVIRONMENT, ENVIRONMENT_CONTACT, AGING, NUM_OF_EVENT_TYPES};

enum InfectionStatus{IDC, IE, NA};
//infection status means the most recent cause of infection
//NA means never infected
//DC suffix -> infection by direct contact
//E suffix -> infection by contacting environment

enum AgeClass{AGE5,AGE15,AGE100}; //discretize age class
//Age5 = 0-5
//Age15 = 5-15
//Age100 = 15+

class Person{
    friend class Event;
private:
    //attributes
    double m_initialTiterLevel;//need initial titer level after infection for waning fn
    double m_titerLevel;
    double m_timeAtInfection;
    int m_numInfectionsDC; //DC means direct contact
    int m_numInfectionsEnv;
    InfectionStatus m_infectionStatus;
    double m_timeToShed;
    int m_age;
    AgeClass m_ageClass;
    const int m_index;
    
public:
    //default constructor
    Person(int idx, double initialTiterLevel = 1.0, double titerLevel=1.0, double timeAtInfection=numeric_limits<double>::max(), InfectionStatus infStat = S, int numInf=0, int numInfEnv=0, double timeToShed = 0.0):m_initialTiterLevel(initialTiterLevel), m_titerLevel(titerLevel),m_timeAtInfection(timeAtInfection), m_infectionStatus(infStat),m_numInfectionsDC(numInf),m_numInfectionsEnv(numInfEnv),m_index(idx), m_timeToShed(timeToShed){
        
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

    int getNumInfectionsDC(){ //DC means direct contact
        return m_numInfectionsDC;
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
    void setNumInfectionsDC(int numInf){
        m_numInfectionsDC = numInf;
    }
    void setNumInfectionsEnv(int numInfEnv){
        m_numInfectionsEnv = numInfEnv;
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
    void reset(){
        m_infectionStatus = S;
        m_initialTiterLevel = minTiter;
        m_numInfectionsDC = 0;
        m_numInfectionsEnv = 0;
    }
    
    
    void waning(double t){
        const double tnew = convertToMonths(abs(t-m_timeAtInfection)); //time needs to be in months post infection
                if(tnew>=1){//only wanes after one month post infection
                    m_titerLevel= max(minTiter, m_initialTiterLevel*pow((tnew),-waningLambda));
                }
        return;
    }
    double probInfGivenDose(InfectionStatus infstat){//dose will vary based on cause of infection, used mean of all shape parameters
        if(infstat == I1DC or infstat == IRDC){
            return 1-pow((1.0+(infDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
        }
        else if(infstat == I1E or infstat == IRE){
            return 1-pow((1.0+(envDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
        }
        else{
            return 0;
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
    EventDriven_MassAction_Sim(const int n, const double beta, const double death, const int seq, const int i1eq, const int req, const int ireq, const int i1env, const int irenv, const double maxRunTime, const int wanIm, int seed = (random_device())()): rng(seed), N(n), BETA(beta), DEATH_RATE(death),Ssum(seq), I1sum(i1eq), Rsum(req), Irsum(ireq), I1Esum(i1env), IrEsum(irenv), maxRunTime(maxRunTime), waningImmunityScenario(wanIm),unif_real(0.0,1.0),unif_int(0,n-2),people_counter(0){
        people = vector<Person*>(n);
        for (Person* &p: people) p = new Person(people_counter++);
        Now=0.0;
        ii=0;//counter for displaying events in queue
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
        event_counter = vector<int>(NUM_OF_EVENT_TYPES, 0);
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
        timevec[0]=0;
        delta=0;
    }
    
    ~EventDriven_MassAction_Sim() {
        for (Person* &p: people) delete p;
    }
    
    
    const double BETA;
    const double DEATH_RATE;
    double virusCon;
    double environment;
    const int waningImmunityScenario;
    
    exponential_distribution<double> exp_beta;
    exponential_distribution<double> exp_death;
    uniform_real_distribution<double> unif_real;
    uniform_int_distribution<int> unif_int;
    mt19937 rng;
    
    //containers to keep track of various pieces of information
    vector<Person*> people;
    priority_queue<Event, vector<Event>, compTime > EventQ;
    unordered_set<Person*> sheddingPeople;
    
    
    int people_counter;
    const double maxRunTime;
    const double reportingThreshold = 0.01;//used to determine how often to count compartments
    double Now; //simulation time
    double previousTime; //used for time step calculation
    
    //used for displaying events in queue
    double counter=0.0;
    int ii;
    vector<int> event_counter;
    
    //used to keep track of number in each compartment -- to be relaxed in next code iteration
    int Ssum;
    int I1sum;
    int Rsum;
    int Psum;
    int Irsum;
    int I1Esum;
    int IrEsum;
    
    const int N; //total population size
    
    //vectors used for export of compartment counts to csv at end of sim
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
    
    
    //used to keep track of number of events that occur
    int numBirths;
    int numDeaths;
    int numI1inf;
    int numIrinf;
    int numRecI1;
    int numRecIr;
    
    //used to determine how often to add compartment counts to above vectors
    double delta;

    
    void runSimulation(){
        while(nextEvent() and EventQ.top().time < maxRunTime) {
        }
    }
    
    void randomizePopulation(){
        int j=0;
        int TiterLevel50=0;
        int TiterLevel100=0;
        int TiterLevel2048=0;
        for(Person* p: people) {
            if(j<Ssum){
                //add individual to compartment
                p->setInfectionStatus(S);
                
                //set initial titer level
                //uniform_real_distribution<double>titer(minTiter,maxTiter);
                p->setInitialTiterLevel(maxTiter);
                
                //set time at infection (used to calculate time since infection)
                p->setTimeAtInfection(numeric_limits<double>::max());
                j++;
            }
            else if(j<(Ssum+I1sum)){
                //set initial titer level
                //uniform_real_distribution<double>titer(minTiter,maxTiter);
                p->setInitialTiterLevel(maxTiter);
                
                //add individual to compartment
                infectByDirectContact(p);
                j++;
            }
            else if(j<(Ssum+I1sum+Rsum)){
                //add individual to compartment
                p->setInfectionStatus(R);
                
                //set initial number of infections > 0
                p->setNumInfectionsDC(p->getNumInfectionsDC()+1);

                //set initial titer level
                //uniform_real_distribution<double>titer(minTiter,maxTiter);
                p->setInitialTiterLevel(maxTiter);

                //set time at infection (used to calculate time since infection)
                p->setTimeAtInfection(numeric_limits<double>::max());
                j++;
            }
            else{
                //since these are reinfected individuals need to set numInfections > 0
                p->setNumInfectionsDC(p->getNumInfectionsDC()+1);
                
                //set initial titer level
                //uniform_real_distribution<double>titer(minTiter,maxTiter);
                p->setInitialTiterLevel(maxTiter);
                
                //add individual to compartment
                infectByDirectContact(p);
                j++;
            }
            death(p); //set when each individual will die
            
            //set environment contact--occurs daily
            //random so that at initialization everyone isn't contacting at exact same time
            double envCon = unif_real(rng)/(double)365;
            EventQ.emplace(envCon,ENVIRONMENT_CONTACT,p);
            event_counter[ENVIRONMENT_CONTACT]++;
            
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
        EventQ.emplace(chkEnv,CHECK_ENVIRONMENT,nullptr);
        event_counter[CHECK_ENVIRONMENT]++;
        
        antTit50[0]=TiterLevel50;
        antTit100[0]=TiterLevel100;
        antTit2048[0]=TiterLevel2048;
    }
    
    //use these functions to retrieve end of simulation data
    const double maxVectorLength = maxRunTime/reportingThreshold;
    vector<double> printVectorI1(){
        if(I1vec.size() < maxVectorLength){
            I1vec.resize(maxVectorLength,I1vec.back());
        }
        return I1vec;
    }
    vector<double> printVectorS(){
        if(Svec.size() < maxVectorLength){
            Svec.resize(maxVectorLength,Svec.back());
        }
        return Svec;
    }
    vector<double> printVectorR(){
        if(Rvec.size() < maxVectorLength){
            Rvec.resize(maxVectorLength,Rvec.back());
        }
        return Rvec;
    }
    vector<double> printVectorIr(){
        if(Irvec.size() < maxVectorLength){
            Irvec.resize(maxVectorLength,Irvec.back());
        }
        return Irvec;
    }
    vector<double> printTimeVector(){
        if(timevec.size() < maxVectorLength){
            timevec.resize(maxVectorLength,timevec.back());
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
    int numIrE(){
        return IrEsum;
    }

    void infectByDirectContact(Person* p) {
        p->setNumInfectionsDC(p->getNumInfectionsDC()+1);
        if((p->getNumInfectionsDC() + p->getNumInfectionsEnv()) == 1){
            p->setInfectionStatus(I1DC);
            numI1inf++;
        }
        else{
            p->setInfectionStatus(IRDC);
            numIrinf++;
        }
        p->setTimeAtInfection(Now);
        p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
        
        //time to shedding (able to cause infections)
        double sheddingTime = Now + (1/(double)365);
        EventQ.emplace(sheddingTime, BEGIN_SHEDDING,p);//assume shedding begins a day later--can't be instantaneous
        event_counter[BEGIN_SHEDDING]++;
        p->setTimeToShed(sheddingTime);
        
        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;

        while (p->probShedding(Tc)>WPVrecThresh) {
            if(Tc > sheddingTime){//only make infectious contacts after shedding begins
                EventQ.emplace(Tc,INFECTIOUS_CONTACT,p);
                event_counter[INFECTIOUS_CONTACT]++;
            }
            Tc += exp_beta(rng);
        }
        //time to recovery
        double Tr = Tc;//once the contact time is late enough such that the probability of shedding is below WPVrecThresh then it is a recovery time (trying to make it not look like a bug)
        EventQ.emplace(Tr,INFECTION_RECOVERY,p);
        event_counter[INFECTION_RECOVERY]++;
        return;
    }
    void infectByEnvironment(Person* p) {
        p->setNumInfectionsEnv(p->getNumInfectionsEnv()+1);
        if((p->getNumInfectionsDC() + p->getNumInfectionsEnv()) == 1){
            p->setInfectionStatus(I1E);
            numI1inf++;
        }
        else{
            p->setInfectionStatus(IRE);
            numIrinf++;
        }
        p->setTimeAtInfection(Now);
        p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
        
        //time to shedding (able to cause infections)
        double sheddingTime = Now + (1/(double)365);
        EventQ.emplace(sheddingTime, BEGIN_SHEDDING,p);//assume shedding begins a day later--can't be instantaneous
        event_counter[BEGIN_SHEDDING]++;
        p->setTimeToShed(sheddingTime);
        
        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;
        
        while (p->probShedding(Tc)>WPVrecThresh) {
            if(Tc > sheddingTime){//only make infectious contacts after shedding begins
                EventQ.emplace(Tc,INFECTIOUS_CONTACT,p);
                event_counter[INFECTIOUS_CONTACT]++;
            }
            Tc += exp_beta(rng);
        }
        //time to recovery
        double Tr = Tc;//once the contact time is late enough such that the probability of shedding is below WPVrecThresh then it is a recovery time (trying to make it not look like a bug)
        EventQ.emplace(Tr,INFECTION_RECOVERY,p);
        event_counter[INFECTION_RECOVERY]++;
        return;
    }
    
    void environmentalSurveillance(){
/*        if(virusCon>detectionRate){
            cout<<"detected pathogen in water\n";
        }
        else{
            cout<<"did not detect pathogen in water\n";
        }*/
        
        //after each ES event generate time to next one--monthly
        EventQ.emplace(Now+(1/(double)12),CHECK_ENVIRONMENT,nullptr);
        event_counter[CHECK_ENVIRONMENT]++;
        return;
    }

    void death(Person* p){
        exponential_distribution<double> exp_death(DEATH_RATE);
        double Td = exp_death(rng) + Now;
        EventQ.emplace(Td,DEATH,p);
        event_counter[DEATH]++;
        return;
    }
    
    
    int nextEvent() {
        switch (waningImmunityScenario) {
            case 1:{//Famulare waning function
                Event event = EventQ.top();
                EventQ.pop();
                Now = event.time;
                
                //take compartment counts for output
                if(Now-delta>reportingThreshold){
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
                //display events in queue
                if(Now>counter){
                    /*             cout<<"Loop "<< ii<<"\n";
                     cout<<"Now "<<Now<<"\n";
                     cout<<" queue size "<<EventQ.size()<<"\n";
                     cout << "\tINFECTIOUS_CONTACT: " << event_counter[INFECTIOUS_CONTACT] << endl;
                     cout << "\tINFECTION_RECOVERY: " << event_counter[INFECTION_RECOVERY] << endl;
                     cout<< "\tBEGIN_SHEDDING: " << event_counter[BEGIN_SHEDDING] << endl;
                     cout<<"\tENVIRONMENT_CONTACT " << event_counter[ENVIRONMENT_CONTACT]<< endl;
                     cout<<"\tCHECK_ENVIRONMENT " << event_counter[CHECK_ENVIRONMENT]<< endl;
                     cout << "\tDEATH: " << event_counter[DEATH] << endl;*/
                    counter+=.1;
                    ii++;
                }
                //update virus particles in environment
                double timeStep = Now - previousTime;
                if(timeStep > 0 and environment > 0){
                    environment+=exp(-inactivationRate*timeStep*365);
                    if(environment < 0){
                        environment = 0;
                    }
                }
                
                Person* individual = event.person;
                if(event.type==INFECTIOUS_CONTACT){
                    if(sheddingPeople.count(individual)>0){//check this event is still applicable
                        
                        //choose person to contact
                        int contact_idx = unif_int(rng);
                        contact_idx = contact_idx >= individual->getIndex() ? contact_idx++ : contact_idx;
                        Person* contact = people[contact_idx];
                        
                        //no simultaneous infections
                        if(contact->getTimeToShed() < Now and sheddingPeople.count(contact)==0){
                            
                            //wane immunity if applicable
                            if(contact->getTimeAtInfection()!=numeric_limits<double>::max()){
                                contact->waning(Now);
                            }
                            
                            //check susceptibility
                            double r2 = unif_real(rng);
                            if(r2<contact->probInfGivenDose(I1DC)){//just need suffix, prefix doesn't change dose
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
                            }
                        }
                        
                    }
                }
                else if (event.type==INFECTION_RECOVERY) {
                    if(individual->getInfectionStatus()!=S or individual->getInfectionStatus()!=R){//check that this event is still applicable
                        if(individual->getInfectionStatus()==I1DC){
                            I1sum--;
                            Rsum++;
                            individual->setInfectionStatus(R);
                            sheddingPeople.erase(individual);
                            numRecI1++;
                        }
                        else if(individual->getInfectionStatus()==I1E){
                            I1Esum--;
                            Rsum++;
                            individual->setInfectionStatus(R);
                            sheddingPeople.erase(individual);
                            numRecI1++;
                        }
                        else if(individual->getInfectionStatus()==IRE){
                            IrEsum--;
                            Rsum++;
                            individual->setInfectionStatus(R);
                            sheddingPeople.erase(individual);
                            numRecIr++;
                        }
                        else if(individual->getInfectionStatus()==IRDC){
                            Irsum--;
                            Rsum++;
                            individual->setInfectionStatus(R);
                            sheddingPeople.erase(individual);
                            numRecIr++;
                        }
                    }
                }
                else if(event.type == BEGIN_SHEDDING){
                    if(individual->getInfectionStatus()!=S or individual->getInfectionStatus()!=R){//check event still applicable
                        sheddingPeople.insert(individual);
                    }
                }
                else if(event.type == ENVIRONMENT_CONTACT){
                    //define a dose of virus from environment
                    virusCon = environment; //assumes all virus particles shed end up in water source -- will relax this assumption in another code iteration
                    envDose = 2*virusCon; //2 is num L water drank per day
                    
                    //no simultaneous infections
                    if(individual->getTimeToShed() < Now and sheddingPeople.count(individual)==0){
                        
                        //check susceptibility
                        double r2 = unif_real(rng);
                        if(r2<individual->probInfGivenDose(I1E)){//just need suffix, prefix doesn't change dose
                            if(individual->getInfectionStatus()==S){
                                Ssum--;
                                I1Esum++;
                                infectByEnvironment(individual);
                            }
                            else{
                                Rsum--;
                                IrEsum++;
                                infectByEnvironment(individual);
                            }
                        }
                    }
                    
                    //update environment with additions
                    if(sheddingPeople.count(individual) > 0){
                        environment+=gramsFeces*100;//100 is just a place holder since real function to determine viral load depends on age
                    }
                    
                    //set next time to contact environment
                    EventQ.emplace(Now + (1/(double)365),ENVIRONMENT_CONTACT,individual);
                    event_counter[ENVIRONMENT_CONTACT]++;
                }
                else if(event.type==CHECK_ENVIRONMENT){
                    environmentalSurveillance();
                }
                else if(event.type == DEATH){
                    numDeaths++;
                    numBirths++;
                    switch (individual->getInfectionStatus()) {
                        case S:
                        {
                            //don't need to update compartment counts
                            break;
                        }
                        case I1DC:
                        {
                            I1sum--;
                            Ssum++;
                            break;
                        }
                        case R:
                        {
                            Rsum--;
                            Ssum++;
                            break;
                        }
                        case IRDC:
                        {
                            Irsum--;
                            Ssum++;
                            break;
                        }
                        case I1E:
                        {
                            I1Esum--;
                            Ssum++;
                            break;
                        }
                        case IRE:
                        {
                            IrEsum--;
                            Ssum++;
                            break;
                        }
                    }
                    individual->reset();
                    
                    //set time to death
                    death(individual);
                    
                }
                event_counter[event.type]--;
                previousTime = Now;
                

                break;
            }
            case 2:{
                Event event = EventQ.top();
                EventQ.pop();
                Now = event.time;
                
                //take compartment counts for output
                if(Now-delta>reportingThreshold){
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
                //display events in queue
                if(Now>counter){
                    /*             cout<<"Loop "<< ii<<"\n";
                     cout<<"Now "<<Now<<"\n";
                     cout<<" queue size "<<EventQ.size()<<"\n";
                     cout << "\tINFECTIOUS_CONTACT: " << event_counter[INFECTIOUS_CONTACT] << endl;
                     cout << "\tINFECTION_RECOVERY: " << event_counter[INFECTION_RECOVERY] << endl;
                     cout<< "\tBEGIN_SHEDDING: " << event_counter[BEGIN_SHEDDING] << endl;
                     cout<<"\tENVIRONMENT_CONTACT " << event_counter[ENVIRONMENT_CONTACT]<< endl;
                     cout<<"\tCHECK_ENVIRONMENT " << event_counter[CHECK_ENVIRONMENT]<< endl;
                     cout << "\tDEATH: " << event_counter[DEATH] << endl;*/
                    counter+=.1;
                    ii++;
                }
                //update virus particles in environment
                double timeStep = Now - previousTime;
                if(timeStep > 0 and environment > 0){
                    environment+=exp(-inactivationRate*timeStep*365);
                    if(environment < 0){
                        environment = 0;
                    }
                }
                
                Person* individual = event.person;
                if(event.type==INFECTIOUS_CONTACT){
                    if(sheddingPeople.count(individual)>0){//check this event is still applicable
                        
                        //choose person to contact
                        int contact_idx = unif_int(rng);
                        contact_idx = contact_idx >= individual->getIndex() ? contact_idx++ : contact_idx;
                        Person* contact = people[contact_idx];
                        
                        //no simultaneous infections
                        if(contact->getTimeToShed() < Now and sheddingPeople.count(contact)==0){
                            
                            //wane immunity if applicable
                            if(contact->getTimeAtInfection()!=numeric_limits<double>::max()){
                                contact->waning(Now);
                            }
                            
                            //check susceptibility
                            double r2 = unif_real(rng);
                            if(r2<contact->probInfGivenDose(I1DC)){//just need suffix, prefix doesn't change dose
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
                            }
                        }
                        
                    }
                }
                else if (event.type==INFECTION_RECOVERY) {
                    if(individual->getInfectionStatus()!=S or individual->getInfectionStatus()!=R){//check that this event is still applicable
                        if(individual->getInfectionStatus()==I1DC){
                            I1sum--;
                            Rsum++;
                            individual->setInfectionStatus(R);
                            sheddingPeople.erase(individual);
                            numRecI1++;
                        }
                        else if(individual->getInfectionStatus()==I1E){
                            I1Esum--;
                            Rsum++;
                            individual->setInfectionStatus(R);
                            sheddingPeople.erase(individual);
                            numRecI1++;
                        }
                        else if(individual->getInfectionStatus()==IRE){
                            IrEsum--;
                            Rsum++;
                            individual->setInfectionStatus(R);
                            sheddingPeople.erase(individual);
                            numRecIr++;
                        }
                        else if(individual->getInfectionStatus()==IRDC){
                            Irsum--;
                            Rsum++;
                            individual->setInfectionStatus(R);
                            sheddingPeople.erase(individual);
                            numRecIr++;
                        }
                    }
                }
                else if(event.type == BEGIN_SHEDDING){
                    if(individual->getInfectionStatus()!=S or individual->getInfectionStatus()!=R){//check event still applicable
                        sheddingPeople.insert(individual);
                    }
                }
                else if(event.type == ENVIRONMENT_CONTACT){
                    //define a dose of virus from environment
                    virusCon = environment; //assumes all virus particles shed end up in water source -- will relax this assumption in another code iteration
                    envDose = 2*virusCon; //2 is num L water drank per day
                    
                    //no simultaneous infections
                    if(individual->getTimeToShed() < Now and sheddingPeople.count(individual)==0){
                        
                        //check susceptibility
                        double r2 = unif_real(rng);
                        if(r2<individual->probInfGivenDose(I1E)){//just need suffix, prefix doesn't change dose
                            if(individual->getInfectionStatus()==S){
                                Ssum--;
                                I1Esum++;
                                infectByEnvironment(individual);
                            }
                            else{
                                Rsum--;
                                IrEsum++;
                                infectByEnvironment(individual);
                            }
                        }
                    }
                    
                    //update environment with additions
                    if(sheddingPeople.count(individual) > 0){
                        environment+=gramsFeces*100;//100 is just a place holder since real function to determine viral load depends on age
                    }
                    
                    //set next time to contact environment
                    EventQ.emplace(Now + (1/(double)365),ENVIRONMENT_CONTACT,individual);
                    event_counter[ENVIRONMENT_CONTACT]++;
                }
                else if(event.type==CHECK_ENVIRONMENT){
                    environmentalSurveillance();
                }
                else if(event.type == DEATH){
                    numDeaths++;
                    numBirths++;
                    switch (individual->getInfectionStatus()) {
                        case S:
                        {
                            //don't need to update compartment counts
                            break;
                        }
                        case I1DC:
                        {
                            I1sum--;
                            Ssum++;
                            break;
                        }
                        case R:
                        {
                            Rsum--;
                            Ssum++;
                            break;
                        }
                        case IRDC:
                        {
                            Irsum--;
                            Ssum++;
                            break;
                        }
                        case I1E:
                        {
                            I1Esum--;
                            Ssum++;
                            break;
                        }
                        case IRE:
                        {
                            IrEsum--;
                            Ssum++;
                            break;
                        }
                    }
                    individual->reset();
                    
                    //set time to death
                    death(individual);
                    
                }
                event_counter[event.type]--;
                previousTime = Now;
                break;
            }
            default:
                cout<<"incorrect input";
                exit(-1);
                break;
        }
        return 1;

    }
    
};
#endif

