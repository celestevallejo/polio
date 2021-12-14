//
//  SIR_H.h
//  
//
//  Created by Celeste on 10/1/17.
//
//

#ifndef Polio_boxcar_model_h
#define Polio_boxcar_model_h

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include "DIFFEQ_SIM.h"
#include <math.h>
#include <vector>

using namespace std;

class Polio_boxcar_model : public DiffEq_Sim {
    
private:
    //birth/death rate
    const double mu;
    
    const int nbins;
    //effective contact rates
    vector<double> b;
    //recovery rates
    vector<double> g;
    //waning rates
    vector<double> r;
    
public:
    //Polio_boxcar_model() : mu(0.0), b1(0.0), b2(0.0),b3(0.0),b4(0.0), g1(0.0), g2(0.0), g3(0.0), g4(0.0), r1(0.0), r2(0.0), r3(0.0) {nbins=9;}
    Polio_boxcar_model(): mu(0.0), nbins(0){
        b = {0};
        g = {0};
        r = {0};
    }
    //Polio_boxcar_model(double bd, double infect1, double infect2, double infect3, double infect4, double rec1, double rec2, double rec3, double rec4, double wane1, double wane2, double wane3): mu(bd), b1(infect1), b2(infect2), b3(infect3), b4(infect4), g1(rec1), g2(rec2), g3(rec3), g4(rec4), r1(wane1), r2(wane2), r3(wane3) {nbins=9;}
    Polio_boxcar_model(double bd, vector<double> beta, vector<double> gamma, vector<double> recovery, int totalBins): mu(bd), nbins(totalBins){
        b.resize(beta.size());
        g.resize(gamma.size());
        r.resize(recovery.size());
        for(int i = 0; i < beta.size(); i++){
            b[i] = beta[i];
        }
        for(int i = 0; i < gamma.size(); i++){
            g[i] = gamma[i];
        }
        for(int i = 0; i < recovery.size(); i++){
            r[i] = recovery[i];
        }
    }
    ~Polio_boxcar_model() {};
    
    void initialize(vector<double> initialValues, double rbins) {
        x = new double[nbins];
        for(int i = 0; i < nbins; i++){
            if(i >= (initialValues.size()-1)){
                x[i] = initialValues.back()/(double)rbins;//last element in initial values vector is sum of all R compartments
            }
            else{
                x[i] = initialValues[i];
            }
        }
    }
    
    void derivative(double const x[], double dxdt[]) {
        //key:
        //x[0] -> S         x[5] -> R_1
        //x[1] -> I_1       x[6] -> R_2
        //x[2] -> I_2       x[7] -> R_3
        //x[3] -> I_3       x[8] -> R_4
        //x[4] -> I_4
        cout<<"next run\n";
        cout<<"dxdt[0] "<<dxdt[0]<<"\n";
        cout<<"dxdt[1] "<<dxdt[1]<<"\n";
        cout<<"dxdt[2] "<<dxdt[2]<<"\n";
        cout<<"dxdt[3] "<<dxdt[3]<<"\n";
        cout<<"dxdt[4] "<<dxdt[4]<<"\n";
        cout<<"dxdt[5] "<<dxdt[5]<<"\n";
        cout<<"dxdt[6] "<<dxdt[6]<<"\n";
        cout<<"dxdt[7] "<<dxdt[7]<<"\n";
        cout<<"dxdt[8] "<<dxdt[8]<<"\n";
        double totalPop = 0;
        double infectPop = 0;
        for(int i = 0; i < nbins; i++){
            totalPop += x[i];
            //cout<<"x total pop "<<x[i]<<"\n";
            if(i >0 and i < 5){
                infectPop += x[i];
                //cout<<"x infect pop "<<x[i]<<"\n";
            }
        }
        dxdt[0] = mu*totalPop - (b[0]*infectPop/totalPop + mu)*x[0];
        dxdt[1] = b[0]*x[0]*infectPop/totalPop - (g[0] + mu)*x[1];
        dxdt[2] = b[1]*x[5]*infectPop/totalPop - (g[1] + mu)*x[2];
        dxdt[3] = b[2]*x[6]*infectPop/totalPop - (g[2] + mu)*x[3];
        dxdt[4] = b[3]*x[7]*infectPop/totalPop - (g[3] + mu)*x[4];
        dxdt[5] = g[0]*x[1] - (b[1]*infectPop/totalPop + mu)*x[5] + r[0]*x[6];
        dxdt[6] = g[1]*x[2] - (b[2]*infectPop/totalPop + r[0] + mu)*x[6] + r[1]*x[7];
        dxdt[7] = g[2]*x[3] - (b[3]*infectPop/totalPop + r[1] + mu)*x[7] + r[2]*x[8];
        dxdt[8] = g[3]*x[4] - (r[2] + mu)*x[8];
        cout<<"dxdt[0] "<<dxdt[0]<<"\n";
        cout<<"dxdt[1] "<<dxdt[1]<<"\n";
        cout<<"dxdt[2] "<<dxdt[2]<<"\n";
        cout<<"dxdt[3] "<<dxdt[3]<<"\n";
        cout<<"dxdt[4] "<<dxdt[4]<<"\n";
        cout<<"dxdt[5] "<<dxdt[5]<<"\n";
        cout<<"dxdt[6] "<<dxdt[6]<<"\n";
        cout<<"dxdt[7] "<<dxdt[7]<<"\n";
        cout<<"dxdt[8] "<<dxdt[8]<<"\n";
        
    }
    void printX(){for(int i=0; i < nbins; i++) { cout << x[i] << " ";} cout << endl; }
    
};


#endif /* Polio_boxcar_model_h */
