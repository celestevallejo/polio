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
#include <numeric>
#include "EventDriven_parameters.hpp"


using namespace std;

enum EventType {INFECTIOUS_CONTACT, INFECTION_RECOVERY, DEATH, BEGIN_SHEDDING, CHECK_ENVIRONMENT, ENVIRONMENT_CONTACT, SHED_INTO_ENVIRONMENT,NUM_OF_EVENT_TYPES};

enum InfectionStatus{IDC, IE, NA, NUM_OF_INFECTION_STATUS};
//infection status means the most recent cause of infection
//NA means never infected
//DC suffix -> infection by direct contact
//E suffix -> infection by contacting environment

enum AgeClass{AGE5, AGE15, AGE100, NUM_OF_AGE_CLASSES}; //discretize age class
//Age5 = 0-5
//Age15 = 5-15
//Age100 = 15+

enum WaningImmunity{FAMULARE, TEUNIS, NUM_OF_WANING_IMMUNITY};

WaningImmunity convertCLArg(const string wi){
    if(wi == "FAMULARE") return FAMULARE;
    else if(wi == "TEUNIS") return TEUNIS;
    else return NUM_OF_WANING_IMMUNITY; //I don't think this will get caught until infectious contact
}

class Person{
    friend class Event;
private:
    //Famulare specific attributes
    double m_initialTiterLevel;//need initial titer level after infection for waning fn
    double m_titerLevel;


    //General attributes
    double m_timeAtInfection;
    InfectionStatus m_infectionStatus;
    int m_numInfectionsDC; //DC means direct contact
    int m_numInfectionsEnv;
    const int m_index;
    double m_timeToShed;
    int m_age;
    AgeClass m_ageClass;

    //Teunis specific attributes
    double m_durationInfection;
    double m_initialPathogen;
    double m_initialImmunityLevel;
    double m_peakImmunityLevel;
    double m_immunityLevel;
    double m_pathogenLevel;

public:
    //default constructor
    Person(int idx, double initialTiterLevel = 1.0, double titerLevel=1.0, double timeAtInfection=numeric_limits<double>::max(), InfectionStatus infStat = NA, int numInf=0, int numInfEnv=0, double timeToShed = numeric_limits<double>::max(), int age =0.0, AgeClass a = AGE5, double durInf=0, double intPat = 0.0, double intAnt=1.0, double pkAnt=1.0, double antLvl = 1.0, double patLvl = 0.0):m_initialTiterLevel(initialTiterLevel), m_titerLevel(titerLevel),m_timeAtInfection(timeAtInfection), m_infectionStatus(infStat),m_numInfectionsDC(numInf),m_numInfectionsEnv(numInfEnv),m_index(idx), m_timeToShed(timeToShed), m_age(age), m_ageClass(a), m_durationInfection(durInf), m_initialPathogen(intPat), m_initialImmunityLevel(intAnt), m_peakImmunityLevel(pkAnt), m_immunityLevel(antLvl), m_pathogenLevel(patLvl){

    }
    //General Functions
    int getIndex () const {
        return m_index;
    }
    int getAge(){
        return m_age;
    }
    void setAge(int a){
        m_age = a;
    }
    int getAgeClass(){
        return m_ageClass;
    }
    void setAgeClass(AgeClass a){
        m_ageClass = a;
    }
    int getNumInfectionsDC(){
        return m_numInfectionsDC;
    }
    int getNumInfectionsEnv(){
        return m_numInfectionsEnv;
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
    InfectionStatus getInfectionStatus(){
        return m_infectionStatus;
    }
    double getTimeToShed(){
        return m_timeToShed;
    }
    void setTimeToShed(double time){
        m_timeToShed = time;
    }
    void reset(int waningImmunityScenario){
        setNumInfectionsDC(0);
        setNumInfectionsEnv(0);
        setInfectionStatus(NA);
        setAge(0);
        setAgeClass(AGE5);
        setTimeToShed(numeric_limits<double>::max());
        setTimeAtInfection(numeric_limits<double>::max());
        if(waningImmunityScenario==FAMULARE){
            setInitialTiterLevel(minTiter);
        }
        else if(waningImmunityScenario==TEUNIS){
            setInitialImmunityLevel(1.0);
            setInitialPathogen(0.0);
            setDurationInfection();//sets peak antibody level
        }
        
    }
    double convertToDays(double t){
        return t*365;
    }
    double convertToMonths(double t){
        return t*12;
    }


    //Famulare specific functions
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
    double getTiterLevel(){
        return m_titerLevel;
    }
    double getInitialTiterLevel(){
        return m_initialTiterLevel;
    }
    double getTimeAtInfection(){
        return m_timeAtInfection;
    }
    void waningFamulare(double t){
        const double tnew = convertToMonths(t-m_timeAtInfection); //time needs to be in months post infection
        if(tnew>=1){//only wanes after one month post infection
            m_titerLevel= max(minTiter, m_initialTiterLevel*pow((tnew),-waningLambda));
        }
        return;
    }
    double probInfGivenDoseFamulare(InfectionStatus infstat){//dose will vary based on cause of infection, used mean of all shape parameters
        if(infstat == IDC){
            return 1-pow((1.0+(infDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
        }
        else if(infstat == IE){
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
    double peakShedding(int age){//age needs to be converted to months
        double ageInMonths = convertToMonths(age);
        if(ageInMonths <newBorn){
            return Smax;
        }
        else{
            return ((Smax-Smin)*exp((newBorn-ageInMonths)/tau)+Smin);
        }
    }

    double stoolViralLoad(double t, int age){
        const double tnew = convertToDays(t-m_timeAtInfection);
        return max(pow(10.0,2.6),pow(10.0,((1.0-k*log2(m_titerLevel))*log(peakShedding(age))))*(exp(eta-(pow(nu,2.0)/2.0)-(pow(log(tnew)-eta,2.0)/(2.0*pow(nu+xsi*log(tnew),2.0))))/tnew));
            //**units are in TCID50/g
    }


    //Teunis specific functions
    double getDurationInfection(){
        return m_durationInfection;
    }
    void setDurationInfection(){
        //sets duration of infection, peak antibody level
        m_durationInfection = (1/(antibodyGrowth - pathogenGrowth))*log(1+((antibodyGrowth - pathogenGrowth)*m_initialPathogen)/(double)(pathogenClearance*m_initialImmunityLevel));
        m_peakImmunityLevel = m_initialImmunityLevel*exp(antibodyGrowth*m_durationInfection);
    }
    double getPeakImmunityLevel(){
        return m_peakImmunityLevel;
    }
    void setInitialPathogen(double pat){
        m_initialPathogen = pat;
    }
    double timePeakPathogen(){
        return (1/(antibodyGrowth - pathogenGrowth))*log((pathogenGrowth/antibodyGrowth)*(1+((antibodyGrowth-pathogenGrowth)*m_initialPathogen)/(pathogenClearance*m_initialImmunityLevel)));
    }
    double getInitialImmunityLevel(){
        return m_initialImmunityLevel;
    }
    void setInitialImmunityLevel(double ant){
        m_initialImmunityLevel = ant;
    }
    double getInitialPathogen(){
        return m_initialPathogen;
    }
    double getImmunityLevel(double t){
        const double tnew = convertToDays(t - m_timeAtInfection);
        const double tnow = convertToDays(t);
        if(m_pathogenLevel > 0){
            m_immunityLevel = m_initialImmunityLevel*exp(antibodyGrowth*tnew);
        }
        else{
            m_immunityLevel = m_peakImmunityLevel*pow((1+(r - 1)*pow(m_peakImmunityLevel,(r-1))*antibodyDecay*(tnow-m_durationInfection)),(-1/(r-1)));
        }
        return m_immunityLevel;
    }
    double getPathogenLevel(double t){
        const double tnew = convertToDays(t);
        if(t <= (m_timeAtInfection + m_durationInfection)){
            m_pathogenLevel = m_initialPathogen*exp(pathogenGrowth*tnew) - (pathogenClearance*m_initialImmunityLevel*(exp(antibodyGrowth*tnew)-exp(pathogenGrowth*tnew))/(antibodyGrowth-pathogenGrowth));
        }
        else{
            m_pathogenLevel = 0.0;
        }
        return m_pathogenLevel;
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
    EventDriven_MassAction_Sim(const int n, const double beta, const double death, const double maxRunTime, const int wanIm, const double virIntWat,int seed=(random_device())()): rng(seed), N(n), BETA(beta), DEATH_RATE(death), maxRunTime(maxRunTime), waningImmunityScenario(wanIm), propVirusinWater(virIntWat),unif_real(0.0,1.0),unif_int(0,n-2),people_counter(0){
        people = vector<Person*>(n);
        deathTime = vector<double>(n);
        birthTime = vector<double>(n);
        for (Person* &p: people) p = new Person(people_counter++);
        Now=0.0;
        ii=0;//counter for displaying events in queue
        IDCvec={0};
        IEvec={0};
        NonInfvec={0};
        timevec={0};
        antTit50={0};
        antTit100={0};
        antTit2048={0};
        ageDist = {0};
        event_counter = vector<int>(NUM_OF_EVENT_TYPES, 0);
        numBirths=0;
        numDeaths=0;
        numDCInf=0;
        numEInf=0;
        delta=0;
        NonInfsum=n;
        IDCsum=0;
        IEsum=0;
        environment = 0;
        previousTime = 0;
        sheddingPeople.clear();
    }

    ~EventDriven_MassAction_Sim() {
        for (Person* &p: people) delete p;
    }


    mt19937 rng;
    const int N; //total population size
    const double BETA;
    const double DEATH_RATE;
    double virusCon;
    double environment;
    const double maxRunTime;
    const int waningImmunityScenario;
    const double propVirusinWater;
    

    exponential_distribution<double>  exp_beta;
    exponential_distribution<double>  exp_death;
    uniform_real_distribution<double> unif_real;
    uniform_int_distribution<int>     unif_int;

    //containers to keep track of various pieces of information
    vector<Person*> people;
    vector<double> deathTime;
    vector<double> birthTime;
    priority_queue<Event, vector<Event>, compTime> EventQ;
    unordered_set<Person*> sheddingPeople;
    vector<int> ageAtFirstInfect;


    int people_counter;
    const double reportingThreshold = 0.01;//used to determine how often to count compartments
    double Now; //simulation time
    double previousTime; //used for time step calculation

    //used for displaying events in queue
    double counter=0.0;
    int ii;
    vector<int> event_counter;

    //used to keep track of number in each compartment -- to be relaxed in next code iteration
    int NonInfsum; //all individuals not infected
    int IDCsum;
    int IEsum;

    //vectors used for export of compartment counts to csv at end of sim
    vector<double> IDCvec;
    vector<double> IEvec;
    vector<double> NonInfvec;
    vector<double> timevec;
    vector<double> antTit50;
    vector<double> antTit100;
    vector<double> antTit2048;
    vector<int>    ageDist;


    //used to keep track of number of events that occur
    int numBirths;
    int numDeaths;
    int numDCInf;
    int numEInf;

    //used to determine how often to add compartment counts to above vectors
    double delta;



    void runSimulation(){
        IDCvec[0] = IDCsum;
        IEvec[0] = IEsum;
        NonInfvec[0] = NonInfsum;
        timevec[0]=0;

        while(nextEvent() and EventQ.top().time < maxRunTime) {
        }
        if(EventQ.top().time >= maxRunTime){
            for(Person* p: people){
                assert((int)((calculateCurrentAge(p))/lengthAgeBuckets) < 20);
                ageDist[(int)((calculateCurrentAge(p))/lengthAgeBuckets)]++;
            }
        }
    }

    void randomizePopulation(int infected){
        int TiterLevel50=0;
        int TiterLevel100=0;
        int TiterLevel2048=0;
        ageDist.resize(numAgeGroups);

        for(Person* p: people) {

            //age distribution based on Pakistan data (pbs.gov.pk - last updated 4/18/18)
            //5 year age buckets
            //age distribution within age buckets taken from fitting age distribution to exponential function
            discrete_distribution<int> age {18.33,15.16,12.54,10.38,8.58,7.10,5.88,4.86,4.02,3.33,2.75,2.28,1.88,1.56,1.29,1.07,0.88,0.73,0.60,0.50};
            int personAge = age(rng);
            p->setAge(chooseAge(personAge,"age"));
            
            if(calculateCurrentAge(p)<=5){
                p->setAgeClass(AGE5);
            }
            else if(calculateCurrentAge(p)<=15){
                p->setAgeClass(AGE15);
            }
            else{
                p->setAgeClass(AGE100);
            }

            p->setInfectionStatus(NA);
            deathTimeEvent(p);
            birthTime[p->getIndex()] = Now;
            p->setTimeAtInfection(numeric_limits<double>::max());
            
            if(waningImmunityScenario == FAMULARE){
                p->setInitialTiterLevel(maxTiter);
            }
            else if(waningImmunityScenario == TEUNIS){
                p->setInitialImmunityLevel(1.0);
            }

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
        
        for(int i = 0; i < infected; i++){
            if(waningImmunityScenario == FAMULARE){
                NonInfsum--;
                IDCsum++;
                infectByDirectContactFamulare(people[i]);
            }
            else if(waningImmunityScenario == TEUNIS){
                NonInfsum--;
                IDCsum++;
                people[i]->setInitialPathogen(1000);//input variable to change
                infectByDirectContactTeunis(people[i]);//infect first "infected" number of people in vector
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
    vector<double> printVectorIDC(){
        if(IDCvec.size() < maxVectorLength){
            IDCvec.resize(maxVectorLength,IDCvec.back());
        }
        return IDCvec;
    }
    vector<double> printVectorNonInf(){
        if(NonInfvec.size() < maxVectorLength){
            NonInfvec.resize(maxVectorLength,NonInfvec.back());
        }
        return NonInfvec;
    }
    vector<double> printVectorIE(){
        if(IEvec.size() < maxVectorLength){
            IEvec.resize(maxVectorLength,IEvec.back());
        }
        return IEvec;
    }
    vector<double> printTimeVector(){
        if(timevec.size() < maxVectorLength){
            timevec.resize(maxVectorLength,timevec.back());
        }
        return timevec;
    }
    vector<int> printAgeDist(){
        return ageDist;
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
    vector<int> printAvgAgeAtFirstInfect(){
        return ageAtFirstInfect;
    }
    
    int NumBirths(){
        return numBirths;
    }
    int NumDeaths(){
        return numDeaths;
    }
    int NumDCInf(){
        return numDCInf;
    }
    int NumEInf(){
        return numEInf;
    }
    int numNonInf(){
        return NonInfsum;
    }
    int numInfDC(){
        return IDCsum;
    }
    int numInfE(){
        return IEsum;
    }
    
    int chooseAge(int personAge,string event){
        int ageGroupIndex = numeric_limits<int>::max();
        //chooses age at death
        if(event=="death"){
            discrete_distribution<int> ageVec(age_Death[personAge].begin(),age_Death[personAge].end());
            ageGroupIndex = ageVec(rng);
        }
        //chooses age at initialization
        else if(event=="age"){
            discrete_distribution<int> ageVec(age_Aging[personAge].begin(),age_Aging[personAge].end());
            ageGroupIndex = ageVec(rng);
        }
        assert(ageGroupIndex != numeric_limits<int>::max());
        int age = personAge*lengthAgeBuckets + ageGroupIndex;
        return age;
    }
    
    int calculateCurrentAge(Person* p){
        return (int)(p->getAge() + Now - birthTime[p->getIndex()]);
    }
    
    void infectByDirectContactFamulare(Person* p) {

        //update number of infections
        p->setNumInfectionsDC(p->getNumInfectionsDC()+1);
        numDCInf++;
        //find average age at first infection (starts counting after half time to let dynamics settle)
        if((p->getNumInfectionsDC() + p->getNumInfectionsEnv())==1 && Now >= maxRunTime/2){
            ageAtFirstInfect.push_back(calculateCurrentAge(p));
        }
        p->setInfectionStatus(IDC);
        p->setTimeAtInfection(Now);
        p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold

        //set time to shedding (able to cause infections) if it occurs before death
        double sheddingTime = Now + (1/(double)365);//assume shedding begins a day later--can't be instantaneous
        p->setTimeToShed(sheddingTime);
        if(p->getTimeToShed()< deathTime[p->getIndex()]){
            EventQ.emplace(sheddingTime, BEGIN_SHEDDING,p);
            event_counter[BEGIN_SHEDDING]++;
        }
        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;
        while (p->probShedding(Tc)>WPVrecThresh) {
            if((Tc > sheddingTime) and (Tc < deathTime[p->getIndex()])){//only make infectious contacts after shedding begins and before death
                EventQ.emplace(Tc,INFECTIOUS_CONTACT,p);
                event_counter[INFECTIOUS_CONTACT]++;
            }
            Tc += exp_beta(rng);
        }
        //set time to recovery if it occurs before death
        double Tr = Tc;//once the contact time is late enough such that the probability of shedding is below WPVrecThresh then it is a recovery time (trying to make it not look like a bug)
        if(Tr < deathTime[p->getIndex()]){
            EventQ.emplace(Tr,INFECTION_RECOVERY,p);
            event_counter[INFECTION_RECOVERY]++;
        }
        return;
    }
    void infectByDirectContactTeunis(Person* p){
        
        //set initial immunity level based on previous immunity level -- setDurationInfection uses to determine boosting
        if((p->getNumInfectionsDC() + p->getNumInfectionsEnv()) > 1){
            p->setInitialImmunityLevel(p->getImmunityLevel(Now));
        }

        //update number of infections
        p->setNumInfectionsDC(p->getNumInfectionsDC()+1);
        numDCInf++;
        p->setTimeAtInfection(Now);
        
        //find average age at first infection (starts counting after half time to let dynamics settle)
        if((p->getNumInfectionsDC() + p->getNumInfectionsEnv())==1 && Now >= maxRunTime/2){
            ageAtFirstInfect.push_back(calculateCurrentAge(p));
        }
        p->setInfectionStatus(IDC);
        p->setDurationInfection();

        //set time to shedding (able to cause infections) if it occurs before death
        double sheddingTime = Now + (1/(double)365);
        p->setTimeToShed(sheddingTime);
        if(p->getTimeToShed() < deathTime[p->getIndex()]){
            EventQ.emplace(sheddingTime, BEGIN_SHEDDING,p);//assume shedding begins a day later
            event_counter[BEGIN_SHEDDING]++;
        }

        //set time to recovery if it occurs before death
        double Tr = p->getDurationInfection() + Now;
        if(Tr < deathTime[p->getIndex()]){
            EventQ.emplace(Tr,INFECTION_RECOVERY,p);
            event_counter[INFECTION_RECOVERY]++;
        }

        //time to next human-human contact
        exponential_distribution<double> exp_beta(BETA);
        double Tc = exp_beta(rng) + Now;

        while ((Tr>Tc) and (Tc < deathTime[p->getIndex()])) {//if contact occurs before recovery and death
            EventQ.emplace(Tc, INFECTIOUS_CONTACT,p);
            event_counter[INFECTIOUS_CONTACT]++;
            Tc+=exp_beta(rng);
        }

    }

    void infectByEnvironment(Person* p) {

        //update number of infections
        p->setNumInfectionsEnv(p->getNumInfectionsEnv()+1);
        numEInf++;
        //find average age at first infection
        if((p->getNumInfectionsDC() + p->getNumInfectionsEnv())==1 && Now >= maxRunTime/2){
            ageAtFirstInfect.push_back(calculateCurrentAge(p));
        }
        p->setInfectionStatus(IE);
        p->setTimeAtInfection(Now);
        

        //set time to shedding (able to cause infections) if it occurs before death
        double sheddingTime = Now + (1/(double)365);//assume shedding begins a day later--can't be instantaneous
        p->setTimeToShed(sheddingTime);
        if(p->getTimeToShed() < deathTime[p->getIndex()]){
            EventQ.emplace(sheddingTime, BEGIN_SHEDDING,p);
            event_counter[BEGIN_SHEDDING]++;
        }

        // time to next human-human contact
        exponential_distribution<double> exp_beta(BETA); //BETA is contact rate/individual/year
        double Tc = exp_beta(rng) + Now;
        
        //initialize recovery time to be set depending on waning immunity scenario
        double Tr = numeric_limits<double>::max();
        
        if(waningImmunityScenario == FAMULARE){
            //boost immunity 10 fold if Famulare waning model
            p->setTiterLevel(11.0*p->getTiterLevel());
            
            //determine infectious contacts (only can occur after shedding begins and before death)
            while (p->probShedding(Tc)>WPVrecThresh) {
                if((Tc > p->getTimeToShed()) and (Tc < deathTime[p->getIndex()])){
                    EventQ.emplace(Tc,INFECTIOUS_CONTACT,p);
                    event_counter[INFECTIOUS_CONTACT]++;
                }
                Tc += exp_beta(rng);
            }
            //set time to recovery if it occurs before death
            Tr = Tc;//once the contact time is late enough such that the probability of shedding is below WPVrecThresh then it is a recovery time (trying to make it not look like a bug)
        }
        else if(waningImmunityScenario == TEUNIS){
            //set initial immunity level to determine boosting amount
            if((p->getNumInfectionsDC() + p->getNumInfectionsEnv()) > 1){
                p->setInitialImmunityLevel(p->getImmunityLevel(Now));
            }
            p->setDurationInfection();
            //set time to recovery
            Tr = p->getDurationInfection() + Now;
            
            //determine infectious contacts
            while ((Tr>Tc) and (Tc < deathTime[p->getIndex()])) {//if contact occurs before recovery and death
                EventQ.emplace(Tc, INFECTIOUS_CONTACT,p);
                event_counter[INFECTIOUS_CONTACT]++;
                Tc+=exp_beta(rng);
            }
        }

        assert(Tr != numeric_limits<double>::max());
        if(Tr < deathTime[p->getIndex()]){
            EventQ.emplace(Tr,INFECTION_RECOVERY,p);
            event_counter[INFECTION_RECOVERY]++;
        }
        return;
    }

    void environmentalSurveillance(){
        /*if(virusCon>detectionRate){
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

    void deathTimeEvent(Person* p){
        
        double rand = unif_real(rng);
        int deathAgeGroup = numeric_limits<int>::max();//do i want to put a default value just in case?
        for(unsigned int i = 0; i < deathCDF.size(); i++){
            if(rand < deathCDF[i]){
                deathAgeGroup = i;
                break;
            }
        }
        assert(deathAgeGroup != numeric_limits<int>::max());
        int deathAge = chooseAge(deathAgeGroup,"death");
        //calculate age to set age at death
        int timeToDeath = deathAge - calculateCurrentAge(p);
        assert((timeToDeath + calculateCurrentAge(p)) <=99);
        
        double Td;
        if(timeToDeath < 0){
            Td = Now;
        }
        else{
            Td = timeToDeath + Now;
        }
        deathTime[p->getIndex()] = Td;
        EventQ.emplace(Td,DEATH,p);
        event_counter[DEATH]++;
        return;
    }
    
    void infectiousContact(Person* p){
        switch (waningImmunityScenario) {
            case FAMULARE:
            {
                //choose person to contact
                int contact_idx = unif_int(rng);
                if(contact_idx >= p->getIndex()) contact_idx++;
                Person* contact = people[contact_idx];
                //no simultaneous infections
                if((contact->getTimeToShed() < Now or contact->getTimeToShed()== numeric_limits<double>::max()) and sheddingPeople.count(contact)==0){
                    
                    //wane immunity if applicable
                    if(contact->getTimeAtInfection()!=numeric_limits<double>::max()){
                        contact->waningFamulare(Now);
                    }
                    
                    //check susceptibility
                    double r2 = unif_real(rng);
                    if(r2<contact->probInfGivenDoseFamulare(IDC)){
                        NonInfsum--;
                        IDCsum++;
                        infectByDirectContactFamulare(contact);
                    }
                }
                break;
                
            }
            case TEUNIS:
            {
                
                //choose person to contact
                int contact_idx = unif_int(rng);
                if(contact_idx >= p->getIndex()) contact_idx++;
                Person* contact = people[contact_idx];
                
                //no simultaneous infections
                //if(contact->getTimeToShed() < Now and sheddingPeople.count(contact)==0){
                if((contact->getTimeToShed() < Now or contact->getTimeToShed()== numeric_limits<double>::max()) and sheddingPeople.count(contact)==0){
                    
                    
                    double dose = p->getPathogenLevel(Now);//assumes all of individual's pathogen is transmitted
                    
                    //check susceptibility (is pathogen concentration greather than immunity level?)
                    if(((dose*exp(pathogenGrowth*Now)- (pathogenClearance/(antibodyGrowth-pathogenGrowth))*(exp(antibodyGrowth*Now)-exp(pathogenGrowth*Now)))- p->getImmunityLevel(Now)*exp(antibodyGrowth*Now)) > 0){
                        NonInfsum--;
                        IDCsum++;
                        infectByDirectContactTeunis(contact);
                    }
                }
                break;
            }
            default:
            {
                cout<<"incorrect input\n";
                exit(-1);
                break;
            }
        }
    }
    
    void infectionRecovery(Person* p){
        //no else statement to incorporate vacc later
        if(p->getInfectionStatus()==IDC){
            IDCsum--;
        }
        else if(p->getInfectionStatus()==IE){
            IEsum--;
        }
        NonInfsum++;
        sheddingPeople.erase(p);
        p->setInfectionStatus(NA);
    }
    
    void environmentContact(Person* p){
        //define a dose of virus from environment
        virusCon = propVirusinWater*environment;
        
        //assumes all virus particles shed end up in water source -- will relax this assumption in another code iteration
        envDose = numLitersDrinkWater*virusCon;
        
        //no simultaneous infections
        if((p->getTimeToShed() < Now or p->getTimeToShed()== numeric_limits<double>::max()) and sheddingPeople.count(p)==0){
            
            //check susceptibility
            if(waningImmunityScenario == FAMULARE){
                //wane immunity if applicable
                if(p->getTimeAtInfection()!=numeric_limits<double>::max()){
                    p->waningFamulare(Now);
                }
                double r2 = unif_real(rng);
                if(r2<p->probInfGivenDoseFamulare(IE)){
                    NonInfsum--;
                    IEsum++;
                    infectByEnvironment(p);
                }
            }
            else if(waningImmunityScenario == TEUNIS){
                //is pathogen concentration greater than immunity level?
                if(((envDose*exp(pathogenGrowth*Now)- (pathogenClearance/(antibodyGrowth-pathogenGrowth))*(exp(antibodyGrowth*Now)-exp(pathogenGrowth*Now)))- p->getImmunityLevel(Now)*exp(antibodyGrowth*Now)) > 0){
                    NonInfsum--;
                    IEsum++;
                    infectByEnvironment(p);
                }
            }
        }
        
        //update environment with additions
        if(sheddingPeople.count(p) > 0){
            if(waningImmunityScenario == FAMULARE){
                environment+=gramsFeces*p->stoolViralLoad(Now,calculateCurrentAge(p));
            }
            else if(waningImmunityScenario == TEUNIS){
                //assumes equal pathogen level shed is equal to pathogen level in individual -- probably should change this assumption
                environment += p->getPathogenLevel(Now);
            }
        }
        
        //set next time to contact environment if it occurs before death
        double nextEnvCon = Now + 1/(double)365;
        if(nextEnvCon < deathTime[p->getIndex()]){
            EventQ.emplace(nextEnvCon,ENVIRONMENT_CONTACT,p);
            event_counter[ENVIRONMENT_CONTACT]++;
        }

    }
    
    void deathEvent(Person* p){
        numDeaths++;
        numBirths++;
        switch (p->getInfectionStatus()) {
            case NA:
            {
                //don't need to update compartment counts
                break;
            }
            case IDC:
            {
                IDCsum--;
                NonInfsum++;
                break;
            }
            case IE:
            {
                IEsum--;
                NonInfsum++;
                break;
            }
            default:
            {
                cout<<"Invalid input\n";
                exit(-1);
                break;
            }
        }
        sheddingPeople.erase(p);
        p->reset(waningImmunityScenario);
        deathTimeEvent(p);
        birthTime[p->getIndex()] = Now;
    }
    
    void updateEnvironment(double timeStep){
        if(timeStep > 0 and environment > 0){
            environment-=exp(inactivationRate*timeStep*365);
            if(environment < 0){
                environment = 0;
            }
        }
    }


    int nextEvent() {
        Event event = EventQ.top();
        EventQ.pop();
        Now = event.time;

        //take compartment counts for output
        if(Now-delta>reportingThreshold){
            IDCvec.push_back(IDCsum);
            IEvec.push_back(IEsum);
            NonInfvec.push_back(NonInfsum);
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
        /*if(Now>counter){
          cout<<"Loop "<< ii<<"\n";
          cout<<"Now "<<Now<<"\n";
          cout<<" queue size "<<EventQ.size()<<"\n";
          cout << "\tINFECTIOUS_CONTACT: " << event_counter[INFECTIOUS_CONTACT] << endl;
          cout << "\tINFECTION_RECOVERY: " << event_counter[INFECTION_RECOVERY] << endl;
          cout<< "\tBEGIN_SHEDDING: " << event_counter[BEGIN_SHEDDING] << endl;
          cout<<"\tENVIRONMENT_CONTACT " << event_counter[ENVIRONMENT_CONTACT]<< endl;
          cout<<"\tCHECK_ENVIRONMENT " << event_counter[CHECK_ENVIRONMENT]<< endl;
          cout << "\tDEATH: " << event_counter[DEATH] << endl;
          counter+=.1;
          ii++;
          }*/
        //update virus particles in environment
        double timeStep = Now - previousTime;
        
        updateEnvironment(timeStep);

        Person* individual = event.person;
        if(event.type==INFECTIOUS_CONTACT){
            infectiousContact(individual);
        }
        else if(event.type==INFECTION_RECOVERY){
            infectionRecovery(individual);
        }
        else if(event.type==BEGIN_SHEDDING){
            sheddingPeople.insert(individual);
        }
        else if(event.type == ENVIRONMENT_CONTACT){
            environmentContact(individual);
        }
        else if(event.type==CHECK_ENVIRONMENT){
            environmentalSurveillance();
        }
        else if(event.type == DEATH){
            deathEvent(individual);
        }
        event_counter[event.type]--;
        previousTime = Now;
        return 1;
    }
};
#endif

