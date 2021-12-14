//
//  EventDriven_parameters.hpp
//  EventDrivenModel
//
//  Created by Celeste on 4/1/17.
//  Copyright © 2017 Celeste. All rights reserved.
//

#ifndef EventDriven_parameters_hpp
#define EventDriven_parameters_hpp

#include <stdio.h>
#include <math.h>
#include <vector>



const double PIR                  = 0.001; //type 3 paralysis incidence rate
const double maxTiter             = 2048.0;
const double minTiter             = 1.0;
const double maxAge               = 70.0;// (max age in simulation so that mean of dist not zero (in death function)
const double meanAge              = 18;//Nigeria mean age
double TTE                        = 0;
const double vaccRate             = 0.01;
//thresholds used to determine when an individual recovers -- needed to end contact events
const double WPVrecThresh         = 0.17;//used mean shedding duration (for WPV) in prob shedding at t function .17
const double OPVrecThresh         = 0.24;//used mean shedding duration (for OPV)in prob shedding at t function .24
const double propToVacc           = 0.25;
const double numLitersDrinkWater = 2;

//Environment contact parameters
const double inactivationRate     = 2*pow(10,0.223161);//TCID50 per day assumes surface water at 22 deg C (converted from PFU to TCID5)(Hurst,1989)
//const double feces = .128; //Liters/day
//const double urine = 1.4; //Liters/day
const double gramsFeces           = 128.0;//avg grams of feces produced per day
const double chkEnvRate           = 12.0;//once per month
//const double fiftyPerInf=pow(10,7);//units TCID50/L used saturating incidence function -- concentration of PV at which infection rate is 50% of max
//const double propVirusinWater     = 1e-9;

//probability of shedding parameters (these are all means from Famulare paper) - units in days
const double muOPV                = 30.3;
const double sigmaOPV             = 1.86;
const double muWPV                = 43.0;
const double sigmaWPV             = 1.69;
const double deltaShedding        = 1.16;

//peak shedding concentration (these are all means from Famulare paper)
const double Smax                 = 6.7; //units TCID50/gram
const double Smin                 = 4.3; //units TCID50/gram
const double tau                  = 10.0; //units months
const double newBorn              = 7.0/12.0;//age in months

//stool viral load parameters (these are all means from Famulare paper)
const double k                    = 0.056; //unitless
const double eta                  = 1.65; //unitless
const double nu                   = 0.17; //unitless
const double xsi                  = 0.32; //unitless

//probability of infection given dose parameters (these are all means from Famulare paper)
const double betaDose             = 14.0; //TCID50
const double alphaDose            = 0.44; // unitless shape parameter
const double gammaDose            = 0.55; //unitless immunity-dependent shape parameter
const double infDose              = 1000.0; //units TCID50 per infection (upper bound on number of virus particles need to initiate infection)
const double vaccDose             = 1e5; //units TCID50 per vaccine
double envDose;//changes based on concentration of virus in environment (water source)

//waning parameters (these are means from Famulare paper)
const double waningLambda         = 0.75;//unitless

//environmental detection parameters (4/3/2017 - these are just guesses)
const double detectionRate        = 1.0/5.0; //can detect 1 virus particle in 5 liters of water?

//death rate parameters (indexmundi.com, Pakistan as of 2016)
const double infantMortalityRate  = 53.9/1000.0; //0-1 year olds
const double childMortalityRate   = 87.0/1000.0; // <5 year olds
const double adultMortalityRate   = 6.4/1000.0; //>5

//Teunis boosting/waning model
const double pathogenClearance    = 5.0;//needs to be fit
const double pathogenGrowth       = 4.0;//needs to be fit
const double antibodyGrowth       = 5.0;//needs to be fit

//shape parameter for waning fn
const double r                    = 2.19; //needs to be fit
const double antibodyDecay        = 1.49; //needs to be fit

//age determination parameters
int lengthAgeBuckets = 5;
int numAgeGroups = 20;
//CDF for death age based on Pakistan lifetable
std::vector<double> deathCDF {0.178,0.324,0.446,0.546,0.630,0.698,0.755,0.802,0.841,0.874,0.900,0.922,0.941,0.956,0.968,0.979,0.987,0.994,1.000};
//weights for each age based on Pakistan life table - fit to exponential function
std::vector<std::vector<double> >age_Aging = {{18.33,17.64,16.99,16.36,15.75},{15.16,14.60,14.05,13.53,13.03},{12.54,12.08,11.63,11.19,10.78},{10.38,9.99,9.62,9.26,8.92},{8.58,8.27,7.96,7.66,7.38},{7.10,6.84,6.58,6.34,6.10},{5.88,5.66,5.45,5.24,5.05},{4.86,4.68,4.51,4.34,4.18},{4.02,3.87,3.73,3.59,3.46},{3.33,3.20,3.08,2.97,2.86},{2.75,2.65,2.55,2.46,2.36},{2.28,2.19,2.11,2.03,1.96},{1.88,1.81,1.75,1.68,1.62},{1.56,1.50,1.44,1.39,1.34},{1.29,1.24,1.20,1.15,1.11},{1.07,1.03,0.99,0.95,0.92},{0.88,0.85,0.82,0.79,0.76},{0.73,0.70,0.68,0.65,0.63},{0.60,0.58,0.56,0.54,0.52},{0.50,0.48,0.46,0.45,0.43}};
//cdf for age distribution
std::vector<std::vector<double> >age_Death = {{18.33,35.97,52.96,69.31,85.06},{100.22,114.82,128.87,142.40,155.42},{167.97,180.05,191.68,202.87,213.65},{224.03,234.02,243.63,252.90,261.81},{270.40,278.66,286.62,294.28,301.66},{308.76,315.60,322.18,328.52,334.62},{340.50,346.15,351.60,356.84,361.89},{366.75,371.43,375.94,380.28,384.45},{388.47,392.35,396.07,399.66,403.12},{406.44,409.65,412.73,415.70,418.56},{421.31,423.96,426.51,428.97,431.33},{433.61,435.80,437.91,439.94,441.90},{443.79,445.60,447.34,449.03,450.64},{452.20,453.70,455.15,456.54,457.88},{459.17,460.41,461.60,462.75,463.86},{464.93,465.96,466.94,467.90,468.81},{469.69,470.54,471.36,472.15,472.91},{473.64,474.34,475.02,475.67,476.30},{476.90,477.48,478.04,478.58,479.10},{479.60,480.08,480.54,480.99,481.42}};
#endif /* EventDriven_parameters_hpp */
