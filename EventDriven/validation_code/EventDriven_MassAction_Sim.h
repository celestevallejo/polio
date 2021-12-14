#ifndef EDMASSIM_H
#define EDMASSIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
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

enum EventType { VACCINATED_CONTACT, INFECTIOUS_CONTACT, ENVIRONMENT_CONTACT, INFECTION_RECOVERY, VACCINE_RECOVERY, CHECK_ENVIRONMENT, VACCINATE, DEATH, NUM_OF_EVENT_TYPES};

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
    string m_infectionStatus; //I, V, or I_E (infection by environment) - cause of most recent infection
    const int m_index;

    public:
    //default constructor
    Person(int idx, double age=0, double titerLevel=1.0, double timeAtInfection=numeric_limits<double>::max(), string infStat = " ", int numInf=0,int numVacc=0,int numInfEnv=0):m_age(age),m_titerLevel(titerLevel),m_timeAtInfection(timeAtInfection), m_infectionStatus(infStat),m_numInfections(numInf),m_numVaccinations(numVacc),m_numInfectionsEnv(numInfEnv),m_index(idx){

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
        m_titerLevel = 1.0;
        m_infectionStatus=" ";
        m_timeAtInfection=numeric_limits<double>::max();
        m_numInfections=0.0;
        m_numVaccinations=0.0;
        m_numInfectionsEnv=0.0;
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
    void setInfectionStatus(string infectionStatus){
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
    double sheddingThreshold(string infectionStatus){
        if(infectionStatus=="I" or infectionStatus=="I_E"){//need to change I_E if incorrect
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
    string getInfectionStatus(){
        return m_infectionStatus;
    }

    double probInfGivenDose(string infstat){//dose will vary based on if it is WPV or OPV, used mean of all shape parameters
        if(infstat=="I"){
            return 1-pow((1.0+(infDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
        }
        else if(infstat == "V"){
            return 1-pow((1.0+(vaccDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
        }
        else{//must be infstat=="I_E"
            return 1-pow((1.0+(envDose/betaDose)),-alphaDose*pow(m_titerLevel,-gammaDose));
        }
    }
    void waning(double t){
        const double tnew = (t-m_timeAtInfection)*12; //time needs to be in months post infection
        if(tnew>=1){//only wanes after one month post infection
            m_titerLevel= max(1.0, m_titerLevel*pow((tnew),-waningLambda));
        }
        return;
    }
    double peakShedding(){//age needs to be converted to months
        if(m_age<(7/(double)12)){
            return Smax;
        }
        else{
            return ((Smax-Smin)*exp((7-(m_age*12))/tau)+Smin);
        }
    }

    double probShedding(double t){
        const double tnew = (t - m_timeAtInfection)*365;//***t needs to be time since infection (days)
        if(m_infectionStatus=="I" or m_infectionStatus=="I_E"){
            return .5*erfc((log(tnew)-(log(muWPV)-log(deltaShedding)*log(m_titerLevel)))/sqrt(2.0)*log(sigmaWPV));
        }
        else{//m_infectionStatus=="V"
            return .5*erfc((log(tnew)-(log(muOPV)-log(deltaShedding)*log(m_titerLevel)))/sqrt(2.0)*log(sigmaOPV));
        }
    }

    double stoolViralLoad(double t){
        const double tnew = (t - m_timeAtInfection)*365;//***t needs to be time since infection (days)
        if(tnew<30){
            return 0;
        }
        else{
        return max(pow(10.0,2.6),pow(10,((1-k*log2(m_titerLevel))*log(peakShedding())))*(exp(eta-(pow(nu,2)/(double)2)-(pow(log(tnew)-eta,2)/(double)2*pow(nu+xsi*log(tnew),2)))/tnew));
        //**units are in TCID50/g
        }
    }
    double gramsDailyFeces(double t){
        const double tnew = (t - m_timeAtInfection)*365;
        return gramsFeces*tnew;
    }

    void print(int i){
        std::cout<<"Attributes for person: "<<(i+1)<<" : "<<m_age<<" , "<<m_titerLevel<<" ,"<<m_infectionStatus<<" , "<<m_timeAtInfection<<"\n";
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
        //EventDriven_MassAction_Sim(const int n, const double beta, const double birth, const double death, const double betaenv): rng((random_device())()), BETA(beta), BIRTH_RATE(birth), DEATH_RATE(death), BETAENV(betaenv),unif_real(0.0,1.0),unif_int(0,n-2),people_counter(0){
        EventDriven_MassAction_Sim(const int n, const double beta, const double birth, const double death, const double betaenv, int seed = (random_device())()): rng(seed), BETA(beta), BIRTH_RATE(birth), DEATH_RATE(death), BETAENV(betaenv),unif_real(0.0,1.0),unif_int(0,n-2),people_counter(0){
            previousTime = {0};
            people = vector<Person*>(n);
            for (Person* &p: people) p = new Person(people_counter++);
            Now=0.0;
            Environment=0.0;
            virusCon=0.0;
            numInfected=0;//keeps track of number of "active" infected and vacc individuals (active means prob shedding >.15,.24)
            //wellVolume = M_PI*.00043*n;//well is right cylinder with radius 1 m and depth is 21.5 m/50,000 people (Quality of Well Water in Owo, SW Nigeria paper)
            event_counter = vector<int>(NUM_OF_EVENT_TYPES, 0);
        }

        ~EventDriven_MassAction_Sim() {
            for (Person* &p: people) delete p;
        }


        const double BETA;
        const double BIRTH_RATE; // currently isn't used
        const double DEATH_RATE;
        const double BETAENV;
        double wellVolume;
        double Environment;
        double virusCon;
        //double evapRate; //**units in L/day

        exponential_distribution<double> exp_beta;
        exponential_distribution<double> exp_death;
        exponential_distribution<double> exp_betaEnvironment;
        exponential_distribution<double> exp_checkEnvironment;
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


        int people_counter;
        double Now;
        int numInfected;
    double counter=0.0;
    int ii = 1;
    vector<int> event_counter;



        void runSimulation(){
            wellVolume = M_PI*.00043*people.size();//assume well volume stays constant
            //assume evaporation occurs at a rate of 12 mm/day
            //evapRate = .12*people.size(); //**units in L/day
            while(nextEvent()>=0){
                if(nextEvent()==0){
                    finalTime[0]=Now;
                    break;
                }
                continue;
            }
        }

        void randomizePopulation(int k){
        //temporary initial conditions: these can be changed at a later date
            exponential_distribution<double> age (1/meanAge);
            //vector<int> age5indx;
            for(Person* p: people) {
                p->setAge(min(100.0,age(rng)));
                p->setInfectionStatus(" ");
                p->setTiterLevel(1.0);
                p->setTimeAtInfection(numeric_limits<double>::max());
               /* if(p->getAge()<5){
                    age5indx.push_back(p->getIndex());
                }*/
                death(p); //set when each individual will die
            }
           /* double currentIndexCounter=age5indx.size();
            for(auto i=age5indx.rbegin();i!=age5indx.rend();++i,--currentIndexCounter){
                uniform_int_distribution<> dis(0,currentIndexCounter-1);
                const int randomIndex = dis(rng);

                if(*i !=age5indx.at(randomIndex)){
                    swap(age5indx.at(randomIndex),*i);
                }
            }
            for(int i=0;i<rint(propToVacc*age5indx.size());i++){
                vaccinate(people[age5indx[i]]);
            }*/
            for(int i=0;i<k;i++){
                infect(people[i]);//infects first person in vector
            }
          /*  exponential_distribution<double> checkEnvironment(chkEnvRate);
            double Te = checkEnvironment(rng) + Now;
            EventQ.push(Event(Te,CHECK_ENVIRONMENT,nullptr));//nullptr since checking env does not involve a person
            event_counter[CHECK_ENVIRONMENT]++;
           */
        }

        void printPeople(){
            for(unsigned int i = 0; i<people.size(); ++i) {
                people[i]->print(i);
            }
        }


        double FinalTime(){
            return TTE;
        }

        double avgFirstInfec(){
            double vecSum=0.0;
            for(int i=0;i<avgAgeOfFirstInfect.size();i++){
                vecSum+=avgAgeOfFirstInfect[i];
            }
            return vecSum/avgAgeOfFirstInfect.size();
        }
        double R0(){
            return numInfected;
        }


        void infect(Person* p) {
            //cout<<"infection\n";
            numInfected++;
            p->setNumInfections();
            //cout<<"person "<<p->getIndex()<<"\n";
            //cout<<"num of infections "<<p->getNumInfections()<<"\n";
            p->setInfectionStatus("I");
            p->setTimeAtInfection(Now);
            p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
            if(p->getNumInfections()==1){
                //cout<<"in first infection loop: age "<<p->getAge()<<"\n";
                avgAgeOfFirstInfect.push_back(p->getAge());
            }
            //is this a paralytic case?
            double r1 = unif_real(rng);
            if(r1<PIR){
                timeOfParalyticCase.push_back(Now);
            }
            // time to next human-human contact
            int numContacts=0;
            exponential_distribution<double> exp_beta(BETA);
            double Tc = exp_beta(rng) + Now;
            while (p->probShedding(Tc)>WPVrecThresh) {     // does contact occur before recovery?
                EventQ.push(Event(Tc,INFECTIOUS_CONTACT,p));
                event_counter[INFECTIOUS_CONTACT]++;
                Tc += exp_beta(rng);
                numContacts++;
            }
            cout<<"numContacts "<<numContacts<<"\n";
            // tjh - not crazy about pushing Tc as a recovery event.  it *looks* like a bug even if it's not
            EventQ.push(Event(Tc,INFECTION_RECOVERY,p));//recovery used for decrementing # of infecteds
            event_counter[INFECTION_RECOVERY]++;
            cout<<"recovery time "<<Tc<<"\n";
            return;
        }

        void vaccinate(Person* p){
           // cout<<"in vaccinate\n";
            numInfected++;
            p->setInfectionStatus("V");
            p->setNumVaccinations();
            //cout<<"person "<<p->getIndex()<<"\n";
            //cout<<"num vaccinations "<<p->getNumVaccinations()<<"\n";
            p->setTimeAtInfection(Now);
            p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
            // time to next contact
            exponential_distribution<double> exp_beta(BETA);//**temporary contact rate
            double Tc = exp_beta(rng) + Now;
            while (p->probShedding(Tc)>OPVrecThresh) {     // does contact occur before recovery?
                EventQ.push(Event(Tc,VACCINATED_CONTACT,p));
                event_counter[VACCINATED_CONTACT]++;
                Tc += exp_beta(rng);
            }
            EventQ.push(Event(Tc,VACCINE_RECOVERY,p));
            event_counter[VACCINE_RECOVERY]++;
            return;

        }

        void infectByEnvironment(Person* p){
            numInfected++;
            p->setNumInfectionsEnv();
            //cout<<"person "<<p->getIndex()<<"\n";
            //cout<<"num environmental infect "<<p->getNumInfectionsEnv()<<"\n";
            p->setInfectionStatus("I_E");
            p->setTimeAtInfection(Now);
            p->setTiterLevel(11.0*p->getTiterLevel());//boost 10 fold
            if(p->getNumInfectionsEnv()==1){
                cout<<"in first infection loop: age "<<p->getAge()<<"\n";
                avgAgeOfFirstInfect.push_back(p->getAge());
            }
            //is this a paralytic case?
            double r1 = unif_real(rng);
            if(r1<PIR){
                timeOfParalyticCase.push_back(Now);
                cout<<"paralytic case\n";
            }
            exponential_distribution<double> exp_beta(BETA);//contact rate
            double Tc = exp_beta(rng) + Now;
            while (p->probShedding(Tc)>WPVrecThresh) {     // does contact occur before recovery?
                EventQ.push(Event(Tc,INFECTIOUS_CONTACT,p));
                event_counter[INFECTIOUS_CONTACT]++;
                Tc += exp_beta(rng);
            }
            EventQ.push(Event(Tc,INFECTION_RECOVERY,p));
            event_counter[INFECTION_RECOVERY]++;
            cout<<"recovery\n";
            return;
        }
        void death(Person* p){
            exponential_distribution<double> exp_death(DEATH_RATE);
            double deathAge = exp_death(rng);
            double Td = deathAge + Now;
            if((p->getAge()+deathAge)>100.0){
                Td = p->deathTime() + Now;
            }
            EventQ.push(Event(Td,DEATH,p));
            event_counter[DEATH]++;
            return;
        }

        void environmentalSurveillance(){

            if(virusCon>detectionRate){//Environment means virus particles
                cout<<"detected pathogen in water\n";
            }
            else{
                cout<<"did not detect pathogen in water\n";
            }
            //after each ES event generate time to next one
            exponential_distribution<double> checkEnvironment(chkEnvRate);//**temp check rate--needs to change
            double Te = checkEnvironment(rng) + Now;
            EventQ.push(Event(Te,CHECK_ENVIRONMENT,nullptr));
            event_counter[CHECK_ENVIRONMENT]++;
            return;
        }

        int nextEvent() {
            if(numInfected==0 or EventQ.empty()) return 0;
            Event event = EventQ.top();
            Now = event.time;
            //cout<<"Now "<<Now<<"\n";
            //cout<<"event "<<event.type<<"\n";
            if(Now>counter){
                cout<<"Loop "<< ii<<"\n";
                cout<<"Now "<<Now<<"\n";
                cout<<" queue size "<<EventQ.size()<<"\n";
                cout << "\tVACCINATED_CONTACT: " << event_counter[VACCINATED_CONTACT] << endl;
                cout << "\tINFECTIOUS_CONTACT: " << event_counter[INFECTIOUS_CONTACT] << endl;
                cout << "\tENVIRONMENT_CONTACT: " << event_counter[ENVIRONMENT_CONTACT] << endl;
                cout << "\tINFECTION_RECOVERY: " << event_counter[INFECTION_RECOVERY] << endl;
                cout << "\tVACCINE_RECOVERY: " << event_counter[VACCINE_RECOVERY] << endl;
                cout << "\tCHECK_ENVIRONMENT: " << event_counter[CHECK_ENVIRONMENT] << endl;
                cout << "\tVACCINATE: " << event_counter[VACCINATE] << endl;
                cout << "\tDEATH: " << event_counter[DEATH] << endl;
                counter+=.1;
                ii++;
            }
            double timeStep = Now - previousTime[0];
            if(timeStep!=0){
                exponential_distribution<double> exp_virusDeath(inactivationRate*timeStep*365);
                minusEnvironment += exp_virusDeath(rng);; //keeps track of decay rate until it is needed (i.e. used for an event) then is reset
            }
            for(Person* p: people) {
               /* if(p->getAge()<5){//only care about age up to 5 for vacc
                    p->updateAge(timeStep);
                }*/
                // does this person contact the environment?
                exponential_distribution<double> exp_betaEnvironment(BETAENV);//**temporary contact rate
                double rand2 = unif_real(rng);
                double TceStep = exp_betaEnvironment(rng);
                double Tce = TceStep+Now;
                if(rand2 < TceStep){
                    EventQ.push(Event(Tce,ENVIRONMENT_CONTACT,p));
                    event_counter[ENVIRONMENT_CONTACT]++;
                }
            }
            Person* individual = event.person;
            if(event.type==INFECTIOUS_CONTACT){//includes contact with infected and infected by environment
                //cout<<"in contact\n";
                int contact_idx = unif_int(rng);
                contact_idx = contact_idx >= individual->getIndex() ? contact_idx++ : contact_idx;
                Person* contact = people[contact_idx];
                double r1 = unif_real(rng);
                //first wane the immunity of the person that is contacted
                if(contact->getTimeAtInfection()!=numeric_limits<double>::max()){
                    contact->waning(Now);
                }
                //contacted person is infected if probability of infection given appropriate dose is sufficiently low
                if(r1 < contact->probInfGivenDose("I")){
                    infect(contact);
                }

            }
            else if(event.type==VACCINATED_CONTACT){//contact with vaccinated individual
                int contact_idx = unif_int(rng);
                contact_idx = contact_idx >= individual->getIndex() ? contact_idx++ : contact_idx;
                Person* contact = people[contact_idx];
                double r1 = unif_real(rng);
                //first wane the immunity of the person that is contacted
                if(contact->getTimeAtInfection()!=numeric_limits<double>::max()){
                    contact->waning(Now);
                }
                //contacted person is infected if probability of infection given appropriate dose is sufficiently low
                if(r1 < contact->probInfGivenDose("I")){//contact w/vacc individual can induce infection
                    infect(contact);
                }
            }
            else if (event.type == INFECTION_RECOVERY) {//recovery from vacc and WPV may be different
                TTE=Now;
                numInfected--;
                if(individual->getIndex()==0){
                    return 0;
                }
            }
            else if(event.type== VACCINE_RECOVERY){
                TTE=Now;
                numInfected--;
            }
            else if(event.type==DEATH){
                //cout<<"death\n";
                if(individual->getIndex()==0){
                    return 0;
                }
                //updates numInfected count if individual that dies is still infected
                if(individual->probShedding(Now)<individual->sheddingThreshold(individual->getInfectionStatus())){
                    numInfected--;
                }
                individual->reset();//keeps population constant
                death(individual);//give each new person a death
            }
            else if(event.type==ENVIRONMENT_CONTACT){
                //first update Environment with additions
                double r3 = unif_real(rng);
                for(Person* p: people){
                    if(p->getInfectionStatus()!= " " and r3 < p->probShedding(Now)){
                        Environment+=p->gramsDailyFeces(Now)*p->stoolViralLoad(Now);//1/1000th of virus particles end up in well
                    }
                }
                //next update virus concentration in water source to exclude inactivated virus
                virusCon = ((1/(double)100000)*Environment/wellVolume)*minusEnvironment;
                minusEnvironment=0;
                double rand2 = unif_real(rng);
                if(rand2 < individual->probInfGivenDose("I_E")){
                    envDose = .5*virusCon;//.5 is num Liters of water drank in a single sitting
                    infectByEnvironment(individual);
                }
            }
            else if(event.type==VACCINATE){
                vector<int> age5indx;
                for(Person* p: people) {
                    if(p->getAge()<5){
                        age5indx.push_back(p->getIndex());
                    }
                }
                double currentIndexCounter=age5indx.size();
                for(auto i=age5indx.rbegin();i!=age5indx.rend();++i,--currentIndexCounter){
                    uniform_int_distribution<> dis(0,currentIndexCounter-1);
                    const int randomIndex = dis(rng);

                    if(*i !=age5indx.at(randomIndex)){
                        swap(age5indx.at(randomIndex),*i);
                    }
                }
                for(int i=0;i<rint(propToVacc*age5indx.size());i++){
                    vaccinate(people[age5indx[i]]);
                }
                //set next time of vaccination
                exponential_distribution<double> exp_vacc(vaccRate);
                double Tv = exp_vacc(rng) + Now;
                EventQ.push(Event(Tv,VACCINATE,nullptr));
                event_counter[VACCINATE]++;
              }
            else if(event.type==CHECK_ENVIRONMENT){
                //first update Environment with additions
                double r3 = unif_real(rng);
                for(Person* p: people){
                    if(p->getInfectionStatus()!=" " and r3 < p->probShedding(Now)){
                        cout<<"stool viral load "<< p->stoolViralLoad(Now)<<"\n";
                        Environment+=p->gramsDailyFeces(Now)*p->stoolViralLoad(Now);//1/1000th of virus particles end up in well
                    }
                }
                //next update virus concentration in water source to exclude inactivated virus
                virusCon = ((1/(double)100000)*Environment/wellVolume)*minusEnvironment;
                minusEnvironment=0;
                environmentalSurveillance();
            }
            previousTime[0]=event.time;
            EventQ.pop();
            event_counter[event.type]--;

            return 1;
        }

};
#endif
