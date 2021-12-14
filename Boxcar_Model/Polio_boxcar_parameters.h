//
//  Polio_boxcar_parameters.h
//  
//
//  Created by Celeste on 12/11/17.
//
//

#ifndef Polio_boxcar_parameters_h
#define Polio_boxcar_parameters_h

#include <math.h>
#include <vector>
#include <string>

using namespace std;


//set number of compartments
int nbins = 13; //total compartments
int sbins = 1;
int ibins = 3;
int rbins = nbins - (sbins + ibins);
int subrbins = 2;

const double numContacts = 100; //per year
const double birthRate = 0.02;

//parameter for Famulare probability of infection model
const double beta_fam = 14;    //mean from Famulare
const double alpha_fam = .44; //mean from Famulare
const double gamma_fam = .55;//mean from Famulare

//parameters for use in Teunis model

//lb -> lower bound, ub-> upper bound
//dpi -> days past infection
const double lb_dpi_fast = 0.0;//for fastWaningRate
const double ub_dpi_fast = 50.0;
const double lb_dpi_slow = 50.0;//for slowWaningRate
const double ub_dpi_slow = 1000.0;

const double antGrowth = 0.2; //antibody growth rate
const double waningShape = 1.2;//r in Teunis model
const double waningRate = 1.5; //nu in Teunis model

vector<double> intAnt(sbins + rbins);//initial antibody concentrations for S and R compartments

void setIntAnt(){
    for(unsigned int i = 0; i < intAnt.size();i++){
        intAnt[i] = pow(10,i);//need to change for real simulations
    }
}

vector<double> peakAnt(sbins + rbins); //peak antibody concentration for S and R compartments
void setPeakAnt(){
    for(unsigned int i = 0; i < peakAnt.size(); i++){
        peakAnt[i] = intAnt[i]*11.0;//assume 10 fold boosting
    }
}

vector<double> infectPeriod(ibins); //infectious period of I compartments
void setInfectPeriod(){
    for(unsigned int i = 0; i< infectPeriod.size();i++){
        int recIndex = i + 1;//I compartment has same antibody level as associated R compartment (i.e. I_1 and R_1)
        infectPeriod[i] = (1.0/antGrowth)*log(peakAnt[recIndex]/intAnt[recIndex]);
    }
}

double getInfectPeriod(int index){
    return infectPeriod[index];
}


//waning functions from Teunis model
double fastWaningRate(int i){//assumes fast waning period occurs between upper and lower bound days past infection for fast waning
    //takes average of upper and lower bound waning rate
    int recIndex = i + 1;//increment to account for S class in peakAnt vector
    return (peakAnt[recIndex]*pow((1.0+(waningShape-1.0)*pow(peakAnt[recIndex],(waningShape-1.0))*waningRate*ub_dpi_fast),(-1.0/(waningShape-1.0)))-peakAnt[recIndex]*pow((1.0+(waningShape-1.0)*pow(peakAnt[recIndex],(waningShape-1.0))*waningRate*lb_dpi_fast),(-1.0/(waningShape-1.0))))/(ub_dpi_fast-lb_dpi_fast);
}
double slowWaningRate(int i){//slow waning time window between upper and lower bound days past infection for slow waning
    //takes average of upper and lower bound waning rate
    int recIndex = i + 1;//increment to account for S class in peakAnt vector
    return (peakAnt[recIndex]*pow((1.0+(waningShape-1.0)*pow(peakAnt[recIndex],(waningShape-1.0))*waningRate*ub_dpi_slow),(-1.0/(waningShape-1.0)))-peakAnt[recIndex]*pow((1.0+(waningShape-1.0)*pow(peakAnt[recIndex],(waningShape-1.0))*waningRate*lb_dpi_slow),(-1.0/(waningShape-1.0))))/(ub_dpi_slow-lb_dpi_slow);
}

//probability of infection model from Famulare
double probInfection(double antibody){
    return 1-pow(1+(1000.0/beta_fam),-alpha_fam*pow(antibody,-gamma_fam));
}

/*array<double,numComp> probInfection(){//dose response model from Famulare with dose fixed at 1000 TCID50
    array<double,numComp> pInfection;
    for(unsigned int i = 0; i < intAnt.size(); i++){
    cout<<"antibody "<<antibody<<"\n";
    cout<<1-pow(1+(1000.0/beta_fam),(-alpha_fam*pow(intAnt[i],2)))<<"\n";
        pInfection[i] = 1-pow(1+(1000.0/beta_fam),-alpha_fam*pow(intAnt[i],-gamma_fam));
    //cout<<"prob infection "<<1-pow(1+(1000.0/beta_fam),-alpha_fam*pow(antibody,-gamma_fam))<<"\n";
    }
    return pInfection;
}*/



#endif /* Polio_boxcar_parameters_h */
