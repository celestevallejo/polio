//
//  main.cpp
//  Practice
//
//  Created by Celeste on 2/16/17.
//  Copyright © 2017 Celeste. All rights reserved.
//

#include <iostream>
#include <array>
#include <random>
#include <algorithm>
#include <math.h>

//create a class for each person in the simulation
class Person{
private:
    //attributes:
    int m_age;
    double m_timeSinceInfection;//infection means any contact with virus
    double m_titerLevel; //does this need to be a double?
    char m_infectionStatus; //'S','I','V','R'
public:
    //default constructor
    Person(int age=0, double timeSinceInfection=0.0,double titerLevel=1.0,char infectionStatus='S'):m_age(age),m_timeSinceInfection(timeSinceInfection),m_titerLevel(titerLevel),m_infectionStatus(infectionStatus){
        //defaults to age 0, 0 time since infection, 1 titer level, and naive susceptible
        
    }
    //set values for attributes
    void setAge(int age){
        m_age = age;
    }
    void setTimeSinceInfection(double timeSinceInfection){
        m_timeSinceInfection += timeSinceInfection;
    }
    void setTiterLevel(double titerLevel){
        m_titerLevel = titerLevel;//max titer level is 2048
    }
    void setInfectionStatus(char infectionStatus){
        m_infectionStatus = infectionStatus;
    }
    int getAge(){
        return m_age;
    }
    double getTimeSinceInfection(){
        return m_timeSinceInfection;
    }
    double getTiterLevel(){
        return m_titerLevel;
    }
    char getInfectionStatus(){
        return m_infectionStatus;
    }
    double probInfGivenDose(double dose){//dose will vary based on if it is WPV or OPV
        if(m_infectionStatus=='S'){
            return 1-pow((1.0+(dose/14.0)),-.44*pow(m_titerLevel,0.55));//put in mean and shape parameters
        }
        else{
            return 0.0;
        }
    }
    double waning(){
        if(m_infectionStatus=='R'){
            return std::max(1.0, 3.0*pow(m_timeSinceInfection,-.75));///what is baseline line immmunity one month post infection??--replace 3.0 with this
        }
        else{
            return 0.0;
        }
    }
    
    void print(int i){
        std::cout<<"Attributes for person: "<<(i+1)<<" : "<<m_age<<" , "<<m_timeSinceInfection<<" , "<<m_titerLevel<<" ,"<<m_infectionStatus<<"\n";
    }
    
};



int main(){
  
    //define number of individuals
    const int numPeople = 5;
    //define number of times to run sim
    const int numSims = 5;
    
    //create fixed array of number of class Person
    std::array<Person,numPeople> individual;
    
    //generate a random real number
    std::random_device rd;//this is the seed
    std::mt19937 gen(rd());//this generates the rn with the above seed
    std::uniform_real_distribution<> unifdis(0.0, 1.0);//generates unif rn
    std::uniform_real_distribution<> titerdis(1.0, 2048.0);//generates unif rn for titer level
    
    //randomly initialize population
    for(int i=0;i<numPeople;++i){
        //first sets age
        double rr=unifdis(gen);
        if(rr<.8){
            individual[i].setAge(1);
        }
        else if(rr<.9){
            individual[i].setAge(2);
        }
        else{
            individual[i].setAge(3);
        }
        //second sets time since last infection
        individual[i].setTimeSinceInfection(unifdis(gen));
        
        //third sets titer level
        individual[i].setTiterLevel(titerdis(gen));

        //fourth sets infection status
    //    double r2 = dis(gen);
    //    if(r2<.1){
    //        individual[i].setInfectionStatus('I');
    //    }
    //    else if(r2<.7){
    //        individual[i].setInfectionStatus('R');
    //    }
    }
    //set one infected individual in the population
    individual[1].setInfectionStatus('I');
  
    
    //parameters
    //const int beta = 200;
    const int recoverInf = 13; //recover from WPV
    const int recoverVacc= 10; //recover from OPV
   // const double wane = .2;
    //const double age = .0025;
   // const double lambda = .75;
   // const double phi = .3; //vaccination rate
    double infDose = 5; //number of doses from WPV infection
    double vaccDose = 3; //number of doses from OPV/vacc
    
    
    //events
    for (int l=0;l<numSims;++l){
        //resets sums after each simulation
        int sSum =0;
        int iSum=0;
        int vSum=0;
        int rSum=0;
        double time =0.0;
        double ageSum=0.0;
        double waneSum=0.0;
        double infectionSum=0.0;
        double vaccinationSum=0.0;
        double chosenEventSum=0.0;
        
        for(int k=0;k<5;++k){//k gives number of steps in each sim
            for(int i=0;i<numPeople;++i){
                //first get total num in each class
                switch(individual[i].getInfectionStatus()){
                    case 'S':
                        sSum+=(1/individual[i].getTiterLevel()); //lower titer level=more susc
                        break;
                    case 'I':
                        iSum+=1;
                        break;
                    case 'V':
                        vSum+=1;
                        break;
                    case 'R':
                        rSum+=1;
                        break;
                }
                //get total aging rate
                ageSum += 1.0/individual[i].getAge();
                //get total waning rate
                waneSum += individual[i].waning();
                //get total infection rate
                infectionSum += individual[i].probInfGivenDose(infDose);
                //get total vaccination rate
                if(individual[i].getAge()==1){
                    vaccinationSum+=individual[i].probInfGivenDose(vaccDose);
                }
                
            }
            
            // get total WPV recovery rate
            double recoverInfSum = recoverInf*iSum;
            //get total OPV recovery rate
            double recoverVaccSum = recoverVacc*vSum;
            
            //sum all rates and choose the event that will occur
            double eventSum = ageSum + infectionSum + vaccinationSum + recoverInfSum + recoverVaccSum + waneSum;
            double r4 = unifdis(gen);
            
            //chooses WPV infection to occur
            if(r4 < (infectionSum/eventSum)){
                //now choose the individual to which the infection will occur
                std::cout<<"infection\n";
                double r5 = unifdis(gen);
                double sum = 0.0;
                for(int i = 0;i<numPeople;++i){
                    if(r5<(sum+(individual[i].probInfGivenDose(infDose)/infectionSum))){
                        individual[i].setInfectionStatus('I');
                        individual[i].setTimeSinceInfection(0.0);
                        individual[i].setTiterLevel(std::min(2048.0,11.0*individual[i].getTiterLevel()));//boost 10 fold 
                        chosenEventSum = infectionSum;
                        break;
                    }
                    sum += individual[i].probInfGivenDose(infDose)/infectionSum;
                }
            }
            //chooses recovery from WPV to occur
            else if(r4<((recoverInfSum+infectionSum)/eventSum)){
                //now choose the individual that will recover
                std::cout<<"recover\n";
                for(int i=0;i<numPeople;++i){//.11 is 40 days
                    if(individual[i].getTimeSinceInfection()>0.11 and individual[i].getInfectionStatus()=='I'){
                        individual[i].setInfectionStatus('R');
                        chosenEventSum = recoverInfSum;
                        break;
                    }
                    
                    
                }
            }
            //chooses waning to occur
            else if(r4<((waneSum+recoverInfSum+infectionSum)/eventSum)){
                //now choose the individual whose immunity wanes
                std::cout<<"wane\n";
                double r6 = unifdis(gen);
                double sum = 0.0;
                for(int i=0;i<numPeople;++i){
                    if(r6<(sum+(individual[i].waning()/waneSum))){
                        individual[i].setTiterLevel(individual[i].waning());
                        individual[i].setInfectionStatus('S');
                        chosenEventSum = waneSum;
                        break;
                    }
                    sum += individual[i].waning()/waneSum;
                }
            }
            //chooses vaccination to occur
            else if(r4<((vaccinationSum+waneSum+recoverInfSum+infectionSum)/eventSum)){
                std::cout<<"vacc\n";
                //now choose the individual to which the vaccination occurs
                double r8 = unifdis(gen);
                double sum = 0.0;
                for(int i=0;i<numPeople;++i){
                    if(r8<(sum+(individual[i].probInfGivenDose(vaccDose)/vaccinationSum))){
                        individual[i].setInfectionStatus('V');
                        individual[i].setTimeSinceInfection(0.0);
                        individual[i].setTiterLevel(std::min(2048.0,11.0*individual[i].getTiterLevel()));//boost 10 fold increase
                        chosenEventSum = vaccinationSum;
                        break;
                    }
                    sum += individual[i].probInfGivenDose(vaccDose)/vaccinationSum;
                }
            }
            //chooses recover from vaccination to occur
            else if(r4<((recoverVaccSum+vaccinationSum+waneSum+recoverInfSum+infectionSum)/eventSum)){
                std::cout<<"recover from vacc\n";
                //now choose individual that recovers from OPV
                //double r9 = dis(gen);
                for(int i=0;i<numPeople;++i){//.11 is 40 days
                    if(individual[i].getTimeSinceInfection()>.11&&individual[i].getInfectionStatus()=='V'){
                        individual[i].setInfectionStatus('R');
                        chosenEventSum = recoverVaccSum;
                        break;
                    }
                    
                }
            }
            //chooses aging to occur
            else{
                //now choose which individual will age
                std::cout<<"age\n";
                double r7 = unifdis(gen);
                double sum = 0.0;
                for(int i =0;i<numPeople;++i){
                    if(r7<(sum+(individual[i].getAge()/ageSum))){
                        individual[i].setAge((individual[i].getAge()+1));
                        chosenEventSum = ageSum;
                        break;
                    }
                    sum += individual[i].getAge()/ageSum;
                }
            }
            //generate exp rn
            std::exponential_distribution<>expdis((eventSum*chosenEventSum));
            double t1 = expdis(gen);
            
            //update all individual's time since infection
            for(int i=0;i<numPeople;++i){
                individual[i].setTimeSinceInfection(t1);
            }
            //update simulation time
            time += t1;

            for(int i=0;i<numPeople;++i){
                individual[i].print(i);
            }
        }

    }
    return 0;
}

