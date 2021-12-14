//
//  adjustable_boxcar_model.h
//  
//
//  Created by Celeste on 1/11/18.
//
//

#ifndef adjustable_boxcar_model_h
#define adjustable_boxcar_model_h

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include "DIFFEQ_SIM.h"
#include <math.h>
#include <vector>
#include <array>
#include <cassert>

using namespace std;

class adjustable_boxcar_model : public DiffEq_Sim {
    
private:
    //birth/death rate
    const double mu;
    //effective contact rates
    vector<double> b;
    //recovery rates
    vector<double> g;
    //waning rates
    vector<double> r;
    vector<double> fastR;
    vector<double> slowR;
    //number infected compartments
    unsigned int nibins;
    //number recovered compartments
    unsigned int nrbins;
    //number of sub-recovered compartments
    unsigned int nsubrbins;
    
    
public:
    adjustable_boxcar_model(): mu(0.0){
        nbins = 0;
        nibins = 0;
        nrbins = 0;
        nsubrbins = 0;
        b = {0};
        g = {0};
        fastR = {0};
        slowR = {0};
    }
    adjustable_boxcar_model(double bd, vector<double> infect, vector<double> recover, vector<double> wane, vector<double> fastWane, vector<double> slowWane, int totalBins, int ibins, int rbins, int subrbins): DiffEq_Sim(totalBins), mu(bd), b(infect), g(recover), r(wane), fastR(fastWane), slowR(slowWane), nibins(ibins), nrbins(rbins), nsubrbins(subrbins){}
    ~adjustable_boxcar_model() {};
    
    void initialize(vector<double>initialValues) {
        x = new double[nbins];
        for(unsigned int i = 0; i < nbins; i++){
            if(i >= (initialValues.size()-1)){
                x[i] = initialValues.back()/(double)nrbins;//last element in initial values vector is sum of all R compartments
            }
            else{
                x[i] = initialValues[i];
            }
        }
    }
    
    void derivative(double const x[], double dxdt[]){
        const int S_idx = 0;
        double totalPop = 0.0;
        double infectPop = 0.0;

        for(unsigned int i = 0; i < nbins; i++){
            assert(x[i]>=0);
            totalPop += x[i];
            if(i > 0 and i < (nibins+1)){
                infectPop += x[i];
            }
        }
        assert(totalPop!= 0);
        //suscpetible compartment is index 0
        //equation for naive susceptible (S) compartment
        dxdt[S_idx] = mu*totalPop - (b[S_idx]*infectPop/totalPop + mu)*x[S_idx];

        //infected compartments start at index 1
        //equations for infected compartments
        const int num_Cols = nibins;
        const int startRow = 2;//start at row 2 since this is the sub R compartments
        for(unsigned int i = 0; i < nibins; i++){
            const int I_idx = i + 1; // offset to account for S compartment
            double subRsum = 0.0;
            if(I_idx == 1){
                for(unsigned int j = 0; j < nsubrbins; j++){
                    subRsum += b[num_Cols*(1+j)+(I_idx-1)]*x[I_idx + startRow*num_Cols + num_Cols*j];
                }
                dxdt[I_idx] = (b[S_idx]*x[S_idx] + subRsum)*infectPop/totalPop - (g[S_idx] + mu)*x[I_idx];
            }
            else{
                for(unsigned int j = 0; j < nsubrbins; j++){
                    subRsum += b[num_Cols*(1+j)+(I_idx-1)]*x[I_idx + startRow*num_Cols + num_Cols*j];
                }
                dxdt[I_idx] = (b[I_idx-1]*x[num_Cols+(I_idx-1)] + subRsum)*infectPop/totalPop - (g[I_idx-1] + mu)*x[I_idx];
                
            }
        }
        //equations for recovered compartments
        int beta_counter = 1;
        for(unsigned int i = 0; i < nrbins; i++){
            const int R_idx = i + nibins + 1;//offset to account for S and I compartments
            const int gamma_counter = i;
            const int fast_wane_counter = i;
            const int slow_wane_counter = i - num_Cols;
            if(R_idx < startRow*num_Cols){//first row of R compartments
                dxdt[R_idx] = g[gamma_counter]*x[R_idx-num_Cols] - (b[beta_counter]*infectPop/totalPop+fastR[fast_wane_counter] + mu)*x[R_idx];
                beta_counter++;
            }
            else if(R_idx == startRow*num_Cols){//only R compartment with complete immunity
                dxdt[R_idx] = g[gamma_counter]*x[R_idx-num_Cols] - (fastR[fast_wane_counter] + mu)*x[R_idx];
            }
            else if(R_idx <= num_Cols*(startRow + 1)){//second row of R compartments
                dxdt[R_idx] = fastR[fast_wane_counter-num_Cols]*x[R_idx-num_Cols] - (b[beta_counter]*infectPop/totalPop+slowR[slow_wane_counter] + mu)*x[R_idx];//change slowR_counter
                beta_counter++;
            }
            else{//last row of R compartments
                dxdt[R_idx] = slowR[slow_wane_counter - num_Cols]*x[R_idx-num_Cols] - (b[beta_counter]*infectPop/totalPop + mu)*x[R_idx];
                beta_counter++;
            }
        }
    }
    void printX(){
        for(unsigned int i=0; i < nbins; i++){
            cout<<x[i]<<" ";
        }
        cout<<endl;
        delete [] x;
        x=NULL;
    }
    
};


#endif /* adjustable_boxcar_model_h */
