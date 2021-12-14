//
//  EventDriven_Sim_exp_recovery.hpp
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

enum EventType {INFECTIOUS_CONTACT_I1, INFECTION_RECOVERY_I1, DEATH, INFECTIOUS_CONTACT_IR,INFECTION_RECOVERY_IR, BIRTH, WANE, NUM_OF_EVENT_TYPES};

enum InfectionStatus{S,I1, R, P, IR, NA};

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
    
public:
    //default constructor
    Person(int idx, double age=0, double titerLevel=1.0, double timeAtInfection=numeric_limits<double>::max(), InfectionStatus infStat = S, int numInf=0,int numVacc=0,int numInfEnv=0):m_age(age),m_titerLevel(titerLevel),m_timeAtInfection(timeAtInfection), m_infectionStatus(infStat),m_numInfections(numInf),m_numVaccinations(numVacc),m_numInfectionsEnv(numInfEnv),m_index(idx){
        
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
    EventDriven_MassAction_Sim(const int n, const double beta, const double birth, const double death, const double gamma, const double rho, const double kappa, int seed = (random_device())()): rng(seed), N(n), BETA(beta), BIRTH_RATE(birth), DEATH_RATE(death), GAMMA(gamma), RHO(rho), KAPPA(kappa), unif_real(0.0,1.0),unif_int(0,n-2),people_counter(0){
        people = vector<Person*>(n);
        for (Person* &p: people) p = new Person(people_counter++);
        Now=0.0;
        ii=0;
        k=0;
        I1vec = {0};
        Irvec = {0};
        Svec = {0};
        Rvec = {0};
        Pvec = {0};
        timevec={0};
        Ssum=Seq;
        I1sum=I1eq;
        Rsum=Req;
        Psum =Peq;
        Irsum=Ireq;
        queueCount = 1;
        event_counter = vector<int>(NUM_OF_EVENT_TYPES, 0);
        numBirths=0;
        numDeaths=0;
        numI1inf=0;
        numIrinf=0;
        numRecI1=0;
        numRecIr=0;
        numWane=0;
        avgInfI1rate=0;
        delta = 0;
        I1vec[0] = I1eq;
        Irvec[0] = Ireq;
        Svec[0] = Seq;
        Rvec[0] = Req + Peq;

    }
    
    ~EventDriven_MassAction_Sim() {
        for (Person* &p: people) delete p;
    }
    
    
    const double BETA;
    const double BIRTH_RATE;
    const double DEATH_RATE;
    const double GAMMA;
    const double KAPPA;
    const double RHO;
    
    exponential_distribution<double> exp_beta;
    exponential_distribution<double> exp_rho;
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
    const int Seq=2;
    const int I1eq=15;
    const int Req=2731;
    const int Peq=962;
    const int Ireq=6290;
    int Ssum;
    int I1sum;
    int Rsum;
    int Psum;
    int Irsum;
    const int N;
    vector<double> I1vec;
    vector<double> Irvec;
    vector<double> Svec;
    vector<double> Rvec;
    vector<double> Pvec;
    vector<double> timevec;
    int k;
    int queueCount;
    unordered_set<Person*> susIndivS;
    unordered_set<Person*> infIndivI1;
    unordered_set<Person*> recIndivR;
    unordered_set<Person*> susIndivP;
    unordered_set<Person*> infIndivIr;
    int numBirths;
    int numDeaths;
    int numI1inf;
    int numIrinf;
    int numRecI1;
    int numRecIr;
    int numWane;
    double avgInfI1rate;
    double delta;


    
    
    
    void runSimulation(){
        while(nextEvent()>=0){
            if(nextEvent()==0){
                TTE=Now;
                break;
            }
            continue;
        }
    }
    int NumInfected(){
        return numInfected;
    }
    
    void randomizePopulation(){
        int j=0;
        for(Person* p: people) {
            if(j<Seq){
                p->setInfectionStatus(S);
                j++;
                susIndivS.emplace(p);
            }
            else if(j<(Seq+I1eq)){
                p->setInfectionStatus(I1);
                infectI1(p);
                infIndivI1.insert(p);
                j++;
            }
            else if(j<(Seq+I1eq+Req)){
                p->setInfectionStatus(R);
                wane(p);
                recIndivR.insert(p);
                j++;
            }
            else if(j<(Seq+I1eq+Req+Peq)){
                p->setInfectionStatus(P);
                j++;
                susIndivP.emplace(p);
            }
            else{
                p->setInfectionStatus(IR);
                infectIr(p);
                infIndivIr.insert(p);
                j++;
            }
            death(p); //set when each individual will die
            //birth(p);
        }
    }
    
    
    
    double FinalTime(){
        return TTE;
    }
    
    vector<double> printVectorI1(){
        return I1vec;
    }
    vector<double> printVectorS(){
        return Svec;
    }
    vector<double> printVectorR(){
        return Rvec;
    }
    vector<double> printVectorP(){
        return Pvec;
    }
    vector<double> printVectorIr(){
        return Irvec;
    }
    vector<double> printTimeVector(){
        return timevec;
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
    int NumWane(){
        return numWane;
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

    
    
    
    void infectI1(Person* p) {
        int numI1contacts=0;
        numI1inf++;
        p->setInfectionStatus(I1);
        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;
        //time to recovery
        exponential_distribution<double> exp_gamma(GAMMA);
        double Tr = exp_gamma(rng) + Now;
        while (Tc<Tr) {     // does contact occur before recovery?
            EventQ.emplace(Tc,INFECTIOUS_CONTACT_I1,p);
            event_counter[INFECTIOUS_CONTACT_I1]++;
            Tc += exp_beta(rng);
            numI1contacts++;
        }
        //cout<<"num I1 contacts "<<numI1contacts<<"\n";
        EventQ.emplace(Tr,INFECTION_RECOVERY_I1,p);
        event_counter[INFECTION_RECOVERY_I1]++;
        return;
    }
    void infectIr(Person* p) {
        int numIrcontacts=0;
        numIrinf++;
        p->setInfectionStatus(IR);//IR means reinfected
        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;
        //time to recovery
        exponential_distribution<double> exp_gamma((GAMMA/KAPPA)); //IR individuals recover sooner than I1
        double Tr = exp_gamma(rng) + Now;
        while (Tc<Tr) {     // does contact occur before recovery?
            EventQ.emplace(Tc,INFECTIOUS_CONTACT_IR,p);
            event_counter[INFECTIOUS_CONTACT_IR]++;
            Tc += exp_beta(rng);
            numIrcontacts++;
        }
        //cout<<"num Ir contacts "<<numIrcontacts<<"\n";
        EventQ.emplace(Tr,INFECTION_RECOVERY_IR,p);
        event_counter[INFECTION_RECOVERY_IR]++;
        return;
    }

    void death(Person* p){
        exponential_distribution<double> exp_death(DEATH_RATE);
        double Td = exp_death(rng) + Now;
        EventQ.emplace(Td,DEATH,p);
        event_counter[DEATH]++;
        return;
    }
    void birth(Person* p){
        //if((Ssum+I1sum+Rsum+Psum+Irsum)<N){
            exponential_distribution<double> exp_birth(BIRTH_RATE*(Ssum+I1sum+Rsum+Psum+Irsum));
            double Tb = exp_birth(rng) + Now;
            EventQ.emplace(Tb,BIRTH,p);
            event_counter[BIRTH]++;
        //}
        return;
        
    }
    void wane(Person* p){
        exponential_distribution<double> exp_rho(RHO);
        double Tw = exp_rho(rng) + Now;
        EventQ.emplace(Tw,WANE,p);
        event_counter[WANE]++;
    }
    
    
    int nextEvent() {
        if(Now > 1){
            return 0;
        }
        Event event = EventQ.top();
        EventQ.pop();
        Now = event.time;
        //cout<<"event "<<event.type<<"\n";
        if(Now-delta>.01){
            /*I1vec[queueCount] = I1sum;
             Irvec[queueCount] = Irsum;
             Svec[queueCount] = Ssum;
             Rvec[queueCount] = Rsum;
             timevec[queueCount]=Now;
             ++queueCount;*/
            I1vec.push_back(I1sum);
            Irvec.push_back(Irsum);
            Svec.push_back(Ssum);
            Rvec.push_back(Rsum);
            timevec.push_back(Now);
            delta=Now;
        }
         /*if(Now>counter){
             cout<<"Loop "<< ii<<"\n";
             cout<<"Now "<<Now<<"\n";
             cout<<" queue size "<<EventQ.size()<<"\n";
             cout << "\tINFECTIOUS_CONTACT: " << event_counter[INFECTIOUS_CONTACT_I1] << endl;
             cout << "\tINFECTIOUS_CONTACT_IR: " << event_counter[INFECTIOUS_CONTACT_IR] << endl;
             cout << "\tINFECTION_RECOVERY: " << event_counter[INFECTION_RECOVERY_I1] << endl;
             cout << "\tINFECTION_RECOVERY_IR: " << event_counter[INFECTION_RECOVERY_IR] << endl;
             cout << "\tDEATH: " << event_counter[DEATH] << endl;
             cout << "\tBIRTH: " << event_counter[BIRTH] << endl;
             cout << "\tWANE: " << event_counter[WANE] << endl;
             counter+=.1;
             ii++;
         }*/
        Person* individual = event.person;
        if(event.type==INFECTIOUS_CONTACT_I1){
            if((susIndivS.size()+susIndivP.size())!=0 and individual->getInfectionStatus()==I1){
                int contact_idx = unif_int(rng);
                contact_idx = contact_idx >= individual->getIndex() ? contact_idx++ : contact_idx;
                Person* contact = people[contact_idx];
                double r2 = unif_real(rng);
                if(contact->getInfectionStatus()==S){
                    cout<<"contact has S infection status\n";
                    cout<<"num s "<<Ssum<<"\n";
                    cout<<"num I1 "<<I1sum<<"\n";
                    cout<<"num P "<<Psum<<"\n";
                    if(r2<(.5)){
                        cout<<"contact infection status: "<<contact->getInfectionStatus()<<"\n";
                        assert(contact->getInfectionStatus()==S);
                        Ssum--;
                        I1sum++;
                        Person* q = *susIndivS.begin();
                        infectI1(q);
                        susIndivS.erase(q);
                        infIndivI1.insert(q);
                    }
                }
                else if(contact->getInfectionStatus()==P){
                    if(r2<(.5*KAPPA)){
                        cout<<"contact infection status: "<<contact->getInfectionStatus()<<"\n";
                        assert(contact->getInfectionStatus()==P);
                        Psum--;
                        Irsum++;
                        Person* q = *susIndivP.begin();
                        infectIr(q);
                        susIndivP.erase(q);
                        infIndivIr.insert(q);
                    }
                }
                
            }
        }
        if(event.type==INFECTIOUS_CONTACT_IR){
            if((susIndivS.size()+susIndivP.size())!=0 and individual->getInfectionStatus()==IR){
                int contact_idx = unif_int(rng);
                contact_idx = contact_idx >= individual->getIndex() ? contact_idx++ : contact_idx;
                Person* contact = people[contact_idx];
                double r2= unif_real(rng);
                if(contact->getInfectionStatus()==S){
                    if(r2<(.5*KAPPA)){
                        Ssum--;
                        I1sum++;
                        Person* q = *susIndivS.begin();
                        infectI1(q);
                        susIndivS.erase(q);
                        infIndivI1.insert(q);
                    }
                }
                else if(contact->getInfectionStatus()==P){
                    if(r2<(.5*pow(KAPPA,2))){
                        Psum--;
                        Irsum++;
                        Person* q = *susIndivP.begin();
                        infectIr(q);
                        susIndivP.erase(q);
                        infIndivIr.insert(q);
                    }
                }
            }
        }
        else if (event.type==INFECTION_RECOVERY_IR) {
            if(individual->getInfectionStatus()==IR){
                Irsum--;
                Rsum++;
                individual->setInfectionStatus(R);
                wane(individual);
                infIndivIr.erase(individual);
                recIndivR.insert(individual);
                numRecIr++;
            }
        }
        else if (event.type == INFECTION_RECOVERY_I1) {
            if(individual->getInfectionStatus()==I1){
                I1sum--;
                Rsum++;
                individual->setInfectionStatus(R);
                wane(individual);
                infIndivI1.erase(individual);
                recIndivR.insert(individual);
                numRecI1++;
            }
        }
        else if(event.type==WANE){
            if(individual->getInfectionStatus()==R){
                Rsum--;
                Psum++;
                individual->setInfectionStatus(P);
                recIndivR.erase(individual);
                susIndivP.insert(individual);
                numWane++;
            }
        }
        else if(event.type==DEATH){
            numDeaths++;
            numBirths++;
            switch (individual->getInfectionStatus()) {
                case S:
                {
                    //Ssum--;
                    //susIndivS.erase(individual);
                    //individual->setInfectionStatus(NA);
                    //birth(individual);
                    death(individual);
                    break;
                }
                case I1:
                {
                    I1sum--;
                    Ssum++;
                    infIndivI1.erase(individual);
                    //individual->setInfectionStatus(NA);
                    //birth(individual);
                    individual->setInfectionStatus(S);
                    susIndivS.insert(individual);
                    death(individual);
                    break;
                }
                case R:
                {
                    Rsum--;
                    Ssum++;
                    recIndivR.erase(individual);
                    //individual->setInfectionStatus(NA);
                    //birth(individual);
                    individual->setInfectionStatus(S);
                    susIndivS.insert(individual);
                    death(individual);
                    break;
                }
                case P:
                {
                    Psum--;
                    Ssum++;
                    susIndivP.erase(individual);
                    //individual->setInfectionStatus(NA);
                    //birth(individual);
                    individual->setInfectionStatus(S);
                    susIndivS.insert(individual);
                    death(individual);
                    break;
                }
                case IR:
                {
                    Irsum--;
                    Ssum++;
                    infIndivIr.erase(individual);
                    //individual->setInfectionStatus(NA);
                    //birth(individual);
                    individual->setInfectionStatus(S);
                    susIndivS.insert(individual);
                    death(individual);
                    break;
                }
                case NA:
                    break;
            }
            //individual->setInfectionStatus(S);
            //susIndivS.insert(individual);
            //numBirths++;
            //Ssum++;
            //death(individual);
        }
        else if(event.type==BIRTH){
            if(individual->getInfectionStatus()==NA){
                if((Ssum+I1sum+Rsum+Psum+Irsum)<N){
                    individual->setInfectionStatus(S);
                    Ssum++;
                    numBirths++;
                    susIndivS.insert(individual);
                    death(individual);
                }
                else{
                    birth(individual);
                }
            }

        }
        if(I1sum<0 or Irsum<0 or Ssum<0 or Rsum<0){
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

