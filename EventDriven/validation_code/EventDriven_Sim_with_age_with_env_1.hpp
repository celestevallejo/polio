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

enum EventType {INFECTIOUS_CONTACT, INFECTION_RECOVERY, DEATH, BIRTH, BEGIN_SHEDDING,AGING,CHECK_ENVIRONMENT,ENVIRONMENT_CONTACT,NUM_OF_EVENT_TYPES};

enum InfectionStatus{I,V,IE,NA,D}; //NA means never infected, D means dead

enum AgeClass{AGE14,AGE24,AGE54,AGE64,AGE100};

class Person{
    friend class Event;
private:
    //attributes
    double m_age;
    AgeClass m_ageclass;
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
    Person(int idx, double age=0, double titerLevel=1.0, double timeAtInfection=numeric_limits<double>::max(), InfectionStatus infStat = NA, int numInf=0,int numVacc=0,int numInfEnv=0,AgeClass ac = AGE14):m_age(age),m_titerLevel(titerLevel),m_timeAtInfection(timeAtInfection), m_infectionStatus(infStat),m_numInfections(numInf),m_numVaccinations(numVacc),m_numInfectionsEnv(numInfEnv), m_ageclass(ac),m_index(idx){
        
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
        m_ageclass=AGE14;
        m_titerLevel = minTiter;
        m_infectionStatus=NA;
        m_timeAtInfection=numeric_limits<double>::max();
        m_numInfections=0.0;
        m_numVaccinations=0.0;
        m_numInfectionsEnv=0.0;
        m_timeToRec=0.0;
        m_timeToShed=0.0;
    }
    void setAgeClass(AgeClass ac){
        m_ageclass = ac;
    }
    AgeClass getAgeClass(){
        return m_ageclass;
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
        Ivec = {0};
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
        Environment=0.0;
        virusCon=0.0;
        numIEInf=0;

    }
    
    ~EventDriven_MassAction_Sim() {
        for (Person* &p: people) delete p;
    }
    
    
    const double BETA;
    const double BIRTH_RATE;
    const double DEATH_RATE;
    double wellVolume;
    double Environment;
    double virusCon;
    
    exponential_distribution<double> exp_beta;
    exponential_distribution<double> exp_death;
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
    array<double,100000> Ivec;
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
    int numIEInf;


    
    
    
    void runSimulation(){
        wellVolume = M_PI*.00043*people.size();//assume well volume stays constant relative to pop size
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
    int age14(){
        int Age14 = 0;
        for(Person* p: people){
            if(p->getAge()<=14){
                Age14++;
            }
        }
        return Age14;
    }
    int age24(){
        int Age24 = 0;
        for(Person* p: people){
            if(p->getAge()<=24 and p->getAge()>=15){
                Age24++;
            }
        }
        return Age24;
    }
    int age54(){
        int Age54 = 0;
        for(Person* p: people){
            if(p->getAge()<=54 and p->getAge()>=25){
                Age54++;
            }
        }
        return Age54;
    }
    int age64(){
        int Age64 = 0;
        for(Person* p: people){
            if(p->getAge()<=64 and p->getAge()>=55){
                Age64++;
            }
        }
        return Age64;
    }
    int age100(){
        int Age100 = 0;
        for(Person* p: people){
            if(p->getAge()<=100 and p->getAge()>=65){
                Age100++;
            }
        }
        return Age100;
    }
    void randomizePopulation(int infected){
        //change initial conditions depending on how we want population to start
        sheddingPeople.clear();
        int age14=0;
        int age24=0;
        int age54=0;
        int age64=0;
        int age100=0;
        for(Person* p: people) {
            discrete_distribution<int> age {0.3199,0.2131,0.3687,0.0543,0.044}; //Pakistan age dist (indexmundi.com)
            int personAge = age(rng);
            switch (personAge) {
                case 0:
                {
                    uniform_int_distribution<int> AGE0(0,14);
                    p->setAge(AGE0(rng));
                    p->setAgeClass(AGE14);
                    age14++;
                    break;
                }
                case 1:
                {
                    uniform_int_distribution<int> AGE1(15,24);
                    p->setAge(AGE1(rng));
                    p->setAgeClass(AGE24);
                    age24++;
                    break;
                }
                case 2:
                {
                    uniform_int_distribution<int> AGE2(25,54);
                    p->setAge(AGE2(rng));
                    p->setAgeClass(AGE54);
                    age54++;
                    break;
                }
                case 3:
                {
                    uniform_int_distribution<int> AGE3(55,64);
                    p->setAge(AGE3(rng));
                    p->setAgeClass(AGE64);
                    age64++;
                    break;
                }
                case 4:
                {
                    uniform_int_distribution<int> AGE4(65,85);
                    p->setAge(AGE4(rng));
                    p->setAgeClass(AGE100);
                    age100++;
                    break;
                }
                    
            }
            aging(p);//sets time to age
            death(p); //set when each individual will die
            p->setInfectionStatus(NA); //sets most recent cause of infection
            p->setTiterLevel(maxTiter);
            p->setTimeAtInfection(numeric_limits<double>::max());//set time at last infection..if infection status is NA then time at infection is largest int
            double envcon = unif_real(rng)/(double)365;
            cout<<"env contact from initialization "<<envcon<<"\n";
            EventQ.emplace(envcon,ENVIRONMENT_CONTACT,p);//set environmental contact--daily
            event_counter[ENVIRONMENT_CONTACT]++;
        }
        cout<<" 0 - 14 "<< age14<<"\n";
        cout<<"15-24 " <<age24<<"\n";
        cout<<"25-54 "<<age54<<"\n";
        cout<<"55-64 "<<age64<<"\n";
        cout<<"65+ "<<age100<<"\n";
        for(int i=0;i<infected;i++){
            infect(people[i]);//infects first k people in vector
        }
        //ES occurs monthly
        double chkenv = unif_real(rng)/(double) 12;
        cout<<"check env time "<<chkenv<<"\n";
        EventQ.emplace(chkenv,CHECK_ENVIRONMENT,nullptr);//nullptr since checking env does not involve a person
        event_counter[CHECK_ENVIRONMENT]++;
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
    
    array<double,100000> printVectorI(){
        return Ivec;
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
    int NumIEInf(){
        return numIEInf;
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
    
    int findOldest(AgeClass ac){
        int oldestAGE =0;
        int oldestIndex =0;
        for(Person* p: people){
            if(p->getAgeClass()==ac and p->getInfectionStatus()!=D){
                if(p->getAge()> oldestAGE){
                    oldestAGE= p->getAge();
                    oldestIndex = p->getIndex();
                }
            }
        }
        return oldestIndex;
        
    }
    void incrementAge(int age){
        for(Person* p: people){
            if(p->getInfectionStatus()!=D and p->getAge() < age){
                p->setAge(p->getAge()+1);
                 if (p->getAge()>=14 and p->getAge()<=24) {
                 p->setAgeClass(AGE24);
                 }
                 else if(p->getAge()>=24 and p->getAge()<=54){
                 p->setAgeClass(AGE54);
                 }
                 else if(p->getAge()>=55 and p->getAge()<=64){
                 p->setAgeClass(AGE64);
                 }
                 else if(p->getAge()>=65){
                 p->setAgeClass(AGE100);
                 }
            }
        }
    }
    
    double ageLifeExpectancy(double age){
        return (67.36697778 - 0.7881155556*age);//formula fit using data from knoema.com
    }
    
    
    void infect(Person* p) {
        assert(p->getInfectionStatus()!=D);
        assert(sheddingPeople.count(p)==0);
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
                cout<<"infection contact time "<<Tc<<"\n";
            }
            Tc += exp_beta(rng);
        }
        
        //time to recovery
        double Tr = Tc;//once the contact time gets late enough then it is a recovery time (trying to make it not look like a bug)
        EventQ.emplace(Tr,INFECTION_RECOVERY,p);
        event_counter[INFECTION_RECOVERY]++;
        cout<<"infection recovery time "<<Tr<<"\n";
        p->setTimeToRec(Tr);
        return;
    }
    
    void infectByEnvironment(Person* p){
        assert(p->getInfectionStatus()!=D);
        assert(sheddingPeople.count(p)==0);
        p->setNumInfections();
        p->setInfectionStatus(IE);
        p->setTimeAtInfection(Now);
        p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
        
        //time to shedding (able to cause infections)
        double sheddingTime = Now + (1/(double)365);
        cout<<"IE shed time "<<sheddingTime<<"\n";
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
                cout<<"IE contact time "<<Tc<<"\n";
            }
            Tc += exp_beta(rng);
        }
        
        //time to recovery
        double Tr = Tc;//once the contact time gets late enough then it is a recovery time (trying to make it not look like a bug)
        cout<<"IE recovery time "<<Tr<<"\n";
        EventQ.emplace(Tr,INFECTION_RECOVERY,p);
        event_counter[INFECTION_RECOVERY]++;
        p->setTimeToRec(Tr);
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
        cout<<"check env time "<<Now + (1/(double)12)<<"\n";
        event_counter[CHECK_ENVIRONMENT]++;
        return;
    }

    void death(Person* p){
        assert(p->getInfectionStatus()!=D);
        exponential_distribution<double> exp_death(1/ageLifeExpectancy(p->getAge()));
        double randDeath = exp_death(rng);
        cout<<"individual age "<<p->getAge()<<"\n";
        cout<<"age-specific life expectancy "<<randDeath<<"\n";
        double Td = randDeath + Now;
        cout<<"time of death "<<Td<<"\n";
        if((p->getAge()+Td)> maxAge){//if the individual is going to die after the set maxAge
            Td = p->deathTime() + Now;
            cout<<"new death time "<<Td<<"\n";
        }
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
    
    void aging(Person* p){
        double Ta = Now + 1;//age every year
        EventQ.emplace(Ta,AGING,p);
        event_counter[AGING]++;
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
        assert(Now >=0);
        if(k%10==0){
            int firstInfCount=0;
            for(Person* p: people){
                if(p->getNumInfections()==1){
                    firstInfCount++;
                }
            }
            I1vec[queueCount] = firstInfCount;
            Ivec[queueCount] = sheddingPeople.size();
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
             cout << "\tAGING: " << event_counter[AGING] << endl;
             cout << "\tENVIRONMENT_CONTACT: " << event_counter[ENVIRONMENT_CONTACT] << endl;
            // counter+=.1;
             ii++;
         //}
        double timeStep = Now - previousTime[0];
        if(timeStep > 0){
            Environment*=exp(-inactivationRate*timeStep*365);
        }
        Person* individual = event.person;
        if(event.type==INFECTIOUS_CONTACT){
            if(individual->getInfectionStatus()==I or individual->getInfectionStatus()==IE){
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
            if(individual->getInfectionStatus()==I or individual->getInfectionStatus()==IE){
                cout<<"index of shedder "<<individual->getIndex()<<"\n";
                cout<<"time of shedding "<<individual->getTimeToShed()<<"\n";
                cout<<"Now "<<Now<<"\n";
                assert(sheddingPeople.count(individual)==0);
                sheddingPeople.insert(individual);
                cout<<"Now "<<Now<<"\n";
                cout<<"time to shed "<<individual->getTimeToShed()<<"\n";
                if(individual->getInfectionStatus()==I){
                    numI1inf++;
                }
                else{
                    numIEInf++;
                }
            }
        }
        else if (event.type == INFECTION_RECOVERY) {
            if(individual->getInfectionStatus()==I or individual->getInfectionStatus()==IE){
                assert(individual->getInfectionStatus()!=D);
                TTE = Now;
                sheddingPeople.erase(individual);
                numRec++;
            }
        }
        else if(event.type==ENVIRONMENT_CONTACT){
            if(individual->getInfectionStatus()!=D){
                //next update virus concentration in water source to exclude inactivated virus
                virusCon = (propVirusinWater*Environment/wellVolume);
                envDose = 2*virusCon;//2 is num Liters of water drank per day
                cout<<"envDose "<<envDose<<"\n";
                double rand1 = unif_real(rng);
                double rand2 = unif_real(rng);
                if(rand1 < 0.5 and sheddingPeople.count(individual)==0 and individual->getTimeToShed() < Now){//assume 50% probability contact
                    if(rand2 < individual->probInfGivenDose(IE)){
                        infectByEnvironment(individual);
                        cout<<"infection by environment!!!!!!!\n";
                    }
                }
                
                //update Environment with additions
                if(sheddingPeople.count(individual)>0){
                    Environment+=gramsFeces*individual->stoolViralLoad(Now);
                }
                cout<<"new additions to environment "<<Environment<<"\n";
                //set next time to contact environment
                EventQ.emplace(Now + (1/(double)365),ENVIRONMENT_CONTACT,individual);
                event_counter[ENVIRONMENT_CONTACT]++;
            }
        }
        else if(event.type==DEATH){
            numDeaths++;
            assert(individual->getInfectionStatus()!=D);
            if(sheddingPeople.count(individual)>0){
                numRec++;
            }
            sheddingPeople.erase(individual);//removes if they are contained in set
            individual->setInfectionStatus(D);
            birth(individual);
        }
        else if(event.type==BIRTH){
            incrementAge(individual->getAge());
            switch (individual->getAgeClass()) {
                case AGE14:{
                    break;
                }
                case AGE24:{
                    int agingIndex = findOldest(AGE14);
                    Person* aging = people[agingIndex];
                    aging->setAgeClass(AGE24);
                    aging->setAge(15);//lowest age in AGE24
                    break;
                }
                case AGE54:{
                    int agingIndex = findOldest(AGE14);
                    Person* aging = people[agingIndex];
                    aging->setAgeClass(AGE24);
                    aging->setAge(15);//lowest age in AGE24
                    int agingIndex1 = findOldest(AGE24);
                    Person* aging1 = people[agingIndex1];
                    aging1->setAgeClass(AGE54);
                    aging1->setAge(25); //lowest age in AGE54
                    break;
                }
                case AGE64:{
                    int agingIndex = findOldest(AGE14);
                    Person* aging = people[agingIndex];
                    aging->setAgeClass(AGE24);
                    aging->setAge(15);//lowest age in AGE24
                    int agingIndex1 = findOldest(AGE24);
                    Person* aging1 = people[agingIndex1];
                    aging1->setAgeClass(AGE54);
                    aging1->setAge(25); //lowest age in AGE54
                    int agingIndex2 = findOldest(AGE54);
                    Person* aging2 = people[agingIndex2];
                    aging2->setAgeClass(AGE64);
                    aging2->setAge(55); //lowest age in AGE64
                    break;
                }
                case AGE100:{
                    int agingIndex = findOldest(AGE14);
                    Person* aging = people[agingIndex];
                    aging->setAgeClass(AGE24);
                    aging->setAge(15);//lowest age in AGE24
                    int agingIndex1 = findOldest(AGE24);
                    Person* aging1 = people[agingIndex1];
                    aging1->setAgeClass(AGE54);
                    aging1->setAge(25); //lowest age in AGE54
                    int agingIndex2 = findOldest(AGE54);
                    Person* aging2 = people[agingIndex2];
                    aging2->setAgeClass(AGE64);
                    aging2->setAge(55); //lowest age in AGE64
                    int agingIndex3 = findOldest(AGE64);
                    Person* aging3 = people[agingIndex3];
                    aging3->setAgeClass(AGE100);
                    aging3->setAge(65); //lowest age in AGE100
                    break;
                }
            }
            numBirths++;
            individual->reset();
            death(individual);
            aging(individual);
            
            //set time to contact environment
            EventQ.emplace(Now + (1/(double)365),ENVIRONMENT_CONTACT,individual);
            event_counter[ENVIRONMENT_CONTACT]++;

        }
        else if(event.type==AGING){
            if(individual->getInfectionStatus()!=D){
                individual->setAge(individual->getAge()+1);
                if (individual->getAge()>=14 and individual->getAge()<=24) {
                    individual->setAgeClass(AGE24);
                }
                else if(individual->getAge()>=24 and individual->getAge()<=54){
                    individual->setAgeClass(AGE54);
                }
                else if(individual->getAge()>=55 and individual->getAge()<=64){
                    individual->setAgeClass(AGE64);
                }
                else if(individual->getAge()>=65){
                    individual->setAgeClass(AGE100);
                }
                aging(individual);//set next time to age
            }
        }
        else if(event.type==CHECK_ENVIRONMENT){
            environmentalSurveillance();
        }
        cout<<"num I1 inf variable "<<numI1inf<<"\n";
        cout<<"num IE inf variable "<<numIEInf<<"\n";
        cout<<"num rec variable "<<numRec<<"\n";
        
        previousTime[0]=event.time;
        EventQ.pop();
        event_counter[event.type]--;
        cout<<"num I inf "<<sheddingPeople.size()<<"\n";
        cout<<"num recoveries in queue "<<event_counter[INFECTION_RECOVERY]<<"\n";
        return 1;
    }
    
};
#endif

