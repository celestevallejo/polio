#!/usr/bin/python3

import zlib

# STATIC PARAMETERS
kappa            = 0.4179
rho              = 0.2
numDaysToRecover = 28
beta             = 135
pir              = 0.005
AFP_Detection    = 1
reintroRate      = 0.001
minBurnIn        = 50
obsPeriod        = 100
seasonalAmp      = 0.05
numSims          = 10000
lifespan         = 50.0
deathRate        = 1.0/lifespan

# DYNAMIC PARAMETERS
movementModel     = [0, 1]
movementRate      = [0, 0.05, 0.1, 0.2, 1.0]
ES_Detection      = [0, 0.001, 0.01, 0.1]
vaccinationRate   = [0, 0.05, 0.2, 0.5]
villagePop        = [[64000],
                     [32000,32000],
                     [16000,16000,16000,16000],
                     [8000,8000,8000,8000,8000,8000,8000,8000],
                     [4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000],
                     [2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000],
                     [1024000,64000],
                     [1024000,32000,32000],
                     [1024000,16000,16000,16000,16000],
                     [1024000,8000,8000,8000,8000,8000,8000,8000,8000],
                     [1024000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000,4000],
                     [1024000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000],
                     [32000,8000,8000,8000,8000],
                     [32000,4000,4000,4000,4000,4000,4000,4000,4000]]

# OUTPUT FILE HANDLING
outputFilename = 'parameters.txt'

with open(outputFilename, 'w') as outFile:
    outFile.write('cksum kappa rho numDaysToRecover beta deathRate pir AFP_Detection reintroRate minBurnIn obsPeriod seasonalAmp numSims movModel moveRate ES_Detection vacRate vilPop\n')
    for mm in movementModel:
        for mr in movementRate:
            for es in ES_Detection:
                for vr in vaccinationRate:
                    for vp in villagePop:
                        vilPopStr = '{'
                        for v in range(len(vp)):
                            vilPopStr += str(vp[v])
                            if v < (len(vp) - 1):
                                vilPopStr += ','
                    
                        vilPopStr += '}'
                        paramStr = f'{kappa} {rho} {numDaysToRecover} {beta} {deathRate} {pir} {AFP_Detection} {reintroRate} {minBurnIn} {obsPeriod} {seasonalAmp} {numSims} {mm} {mr} {es} {vr} {vilPopStr}'

                        cksum = hex(zlib.crc32(paramStr.encode()))[2:]
                        outStr = f'{cksum} {paramStr}\n'
                        outFile.write(outStr)
