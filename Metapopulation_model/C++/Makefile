INCLUDE = -I.
ifdef HPC_GSL_INC
INCLUDE += -I$$HPC_GSL_INC
endif

LIBS = -lm -lgsl -lgslcblas

default: main_Gillespie_multivillage_TTE.cpp
	g++ -O2 -std=c++17 -Wall -Wextra -Wno-deprecated-declarations --pedantic main_Gillespie_multivillage_TTE.cpp $(INCLUDE)
