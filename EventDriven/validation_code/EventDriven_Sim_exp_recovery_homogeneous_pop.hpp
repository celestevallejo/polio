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
extern std:: ofstream &myfile4;
extern std:: ofstream &myfile5;
extern std:: ofstream &myfile6;
extern std:: ofstream &myfile7;
extern std:: ofstream &myfile8;



using namespace std;

enum EventType {INFECTIOUS_CONTACT_I1, INFECTION_RECOVERY_I1, DEATH, INFECTIOUS_CONTACT_IR,INFECTION_RECOVERY_IR, BIRTH, WANE, NUM_OF_EVENT_TYPES};

enum InfectionStatus{S, I1, R, P, IR, NA};

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
        numRec=0;
        numWane=0;
        avgInfI1rate=0;
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
    const int Req=2825;
    const int Peq=649;
    const int Ireq=6509;
    int Ssum;
    int I1sum;
    int Rsum;
    int Psum;
    int Irsum;
    const int N;
    array<double,5000> I1vec;
    array<double,5000> Irvec;
    array<double,5000> Svec;
    array<double,5000> Rvec;
    array<double,5000> Pvec;
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
    int numRec;
    int numWane;
    double avgInfI1rate;

    
    
    
    void runSimulation(){
        I1vec[0] = I1eq;
        Irvec[0] = Ireq;
        Svec[0] = Seq;
        Rvec[0] = Req;
        Pvec[0] = Peq;
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
                infIndivI1.emplace(p);
                j++;
            }
            else if(j<(Seq+I1eq+Req)){
                p->setInfectionStatus(R);
                recIndivR.emplace(p);
                j++;
            }
            else if(j<(Seq+I1eq+Req+Peq)){
                p->setInfectionStatus(P);
                j++;
                susIndivP.emplace(p);
            }
            else{
                p->setInfectionStatus(IR);
                infIndivIr.emplace(p);
                j++;
            }
        }
    }
    
    
    
    double FinalTime(){
        return TTE;
    }
    
    array<double,5000> printVectorI1(){
        return I1vec;
    }
    array<double,5000> printVectorS(){
        return Svec;
    }
    array<double,5000> printVectorR(){
        return Rvec;
    }
    array<double,5000> printVectorP(){
        return Pvec;
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
    int NumRec(){
        return numRec;
    }
    int NumWane(){
        return numWane;
    }
    
    void printVector(int ii){
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
        
    }
    array<double,5000> printVectorIr(){
        return Irvec;
    }
    
    int nextEvent() {
        if(((infIndivI1.size()+infIndivI1.size())==0) or Now > 0.1){
            return 0;
        }
        double birth;
        if(((susIndivS.size()+susIndivP.size()+infIndivI1.size()+infIndivIr.size()+recIndivR.size()))< N){
            birth = BIRTH_RATE*(susIndivS.size()+susIndivP.size()+infIndivI1.size()+infIndivIr.size()+recIndivR.size());
        }
        else{
            birth= 0;
        }
        double infect1 = ((BETA*susIndivS.size())/(double)N)*(infIndivI1.size()+KAPPA*infIndivIr.size());
        double recover1 = GAMMA*infIndivI1.size();
        double wane = RHO*recIndivR.size();
        double infectr = ((KAPPA*BETA*susIndivP.size())/(double)N)*(infIndivI1.size()+KAPPA*infIndivIr.size());
        double recover2 = (GAMMA/KAPPA)*infIndivIr.size();
        double deathS = DEATH_RATE*susIndivS.size();
        double deathI1 = DEATH_RATE*infIndivI1.size();
        double deathR = DEATH_RATE*recIndivR.size();
        double deathP = DEATH_RATE*susIndivP.size();
        double deathIr = DEATH_RATE*infIndivIr.size();
        
        double totalRate = (birth+infect1+recover1+wane+infectr+recover2+deathS+deathI1+deathR+deathP+deathIr);
        
        exponential_distribution<double>timeStep(totalRate);
        double t1 = timeStep(rng);
        Now+=t1;
        
        double r = unif_real(rng);
        
        if(r<(birth/totalRate)){
            Ssum++;
            numBirths++;
            for(Person* &p:people){
                if(p->getInfectionStatus()==NA){
                    susIndivS.insert(p);
                    p->setInfectionStatus(S);
                    
                    break;
                }
            }
        }
        else if(r<(birth+infect1)/totalRate){
            Ssum--;
            I1sum++;
            Person* q = *susIndivS.begin();
            susIndivS.erase(q);
            infIndivI1.insert(q);
            q->setInfectionStatus(I1);
            numI1inf++;
        }
        else if(r<((birth+infect1+recover1)/totalRate)){
            I1sum--;
            Rsum++;
            Person* q = *infIndivI1.begin();
            infIndivI1.erase(q);
            recIndivR.insert(q);
            q->setInfectionStatus(R);
            numRec++;
        }
        else if(r<((birth+infect1+recover1+wane)/totalRate)){
            Rsum--;
            Psum++;
            Person* q = *recIndivR.begin();
            recIndivR.erase(q);
            susIndivP.insert(q);
            q->setInfectionStatus(P);
            numWane++;
        }
        else if(r<(birth+infect1+recover1+wane+infectr)/totalRate){
            Psum--;
            Irsum++;
            Person* q = *susIndivP.begin();
            susIndivP.erase(q);
            infIndivIr.insert(q);
            q->setInfectionStatus(IR);
            numIrinf++;
        }
        else if(r<((birth+infect1+recover1+wane+infectr+recover2)/totalRate)){
            Irsum--;
            Rsum++;
            Person* q = *infIndivIr.begin();
            infIndivIr.erase(q);
            recIndivR.insert(q);
            q->setInfectionStatus(R);
            numRec++;
        }
        else if(r<(birth+infect1+recover1+wane+infectr+recover2+deathS)/totalRate){
            Ssum--;
            Person* q = *susIndivS.begin();
            susIndivS.erase(q);
            q->setInfectionStatus(NA);
            numDeaths++;
        }
        else if(r<(birth+infect1+recover1+wane+infectr+recover2+deathS+deathI1)/totalRate){
            I1sum--;
            Person* q = *infIndivI1.begin();
            infIndivI1.erase(q);
            q->setInfectionStatus(NA);
            numDeaths++;
        }
        else if(r<(birth+infect1+recover1+wane+infectr+recover2+deathS+deathI1+deathR)/totalRate){
            Rsum--;
            Person* q = *recIndivR.begin();
            recIndivR.erase(q);
            q->setInfectionStatus(NA);
            numDeaths++;
        }
        else if(r<(birth+infect1+recover1+wane+infectr+recover2+deathS+deathI1+deathR+deathP)/totalRate){
            Psum--;
            Person* q = *susIndivP.begin();
            susIndivP.erase(q);
            q->setInfectionStatus(NA);
            numDeaths++;
        }
        else{
            Irsum--;
            Person* q = *infIndivIr.begin();
            infIndivIr.erase(q);
            q->setInfectionStatus(NA);
            numDeaths++;
        }
        if(k%150==0){
            I1vec[queueCount] = I1sum;
            Irvec[queueCount] = Irsum;
            Svec[queueCount] = Ssum;
            Pvec[queueCount] = Psum;
            Rvec[queueCount] = Rsum;
            queueCount++;
        }
        k++;
        return 1;
    }
    
};
#endif

