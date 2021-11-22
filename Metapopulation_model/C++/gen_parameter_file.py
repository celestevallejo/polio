#!/usr/bin/python3

# STATIC PARAMETERS
kappa            = 0.4179
rho              = 0.2
numDaysToRecover = 28
beta             = 135
pir              = 0.005
detRate          = 1
reintroduceRate  = 0.001
burnInTime       = 50
obsTime          = 100
seasonalAmp      = 0.05
numSims          = 10000

# DYNAMIC PARAMETERS
movementModel             = [0, 1]
expectedTimeUntilMove     = [0, 1, 5, 10, 20]
environmentalSurveillance = [0, 0.001, 0.01, 0.1]
vaccinationRate           = [0, 0.05, 0.2, 0.5]
villagePop                = [[64000],
                             [32000,32000],
                             [16000,16000,16000,16000],
                             [8000,8000,8000,8000,8000,8000,8000,8000],
                             [4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000],
                             [2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000]]

# OUTPUT FILE HANDLING
outputFilename = 'parameters.txt'

with open(outputFilename, 'w') as outFile:
    outFile.write('kappa rho numDaysToRecover beta pir detRate reintroduceRate burnInTime obsTime seasonalAmp numSims movModel expTimeUntilMov envSurv vacRate vilPop\n')
    for mm in movementModel:
        for etum in expectedTimeUntilMove:
            for es in environmentalSurveillance:
                for vr in vaccinationRate:
                    for vp in villagePop:
                        vilPopStr = '{'
                        for v in range(len(vp)):
                            vilPopStr += str(vp[v])
                            if v < (len(vp) - 1):
                                vilPopStr += ','
                    
                        vilPopStr += '}'
                        outStr = f'{kappa} {rho} {numDaysToRecover} {beta} {pir} {detRate} {reintroduceRate} {burnInTime} {obsTime} {seasonalAmp} {numSims} {mm} {etum} {es} {vr} {vilPopStr}\n'
                        outFile.write(outStr)
