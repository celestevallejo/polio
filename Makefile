FLAGS = -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
output_dir=.

ABCDIR = $(HOME)/work/AbcSmc
GSL_PATH = $(HOME)/work/AbcSmc/gsl_local
SQLDIR = $(ABCDIR)/sqdb
ABC_LIB = -L$(ABCDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -lpthread -ldl

INCLUDE = -I$(ABCDIR) -I$(GSL_PATH)/include/ 

polio_ed: EventDriven/IBM/main_EventDrivenIBM.cpp EventDriven/IBM/EventDriven_Sim.hpp EventDriven/IBM/EventDriven_parameters.hpp | $(output_dir)/polio_data
	g++ $(FLAGS) EventDriven/IBM/main_EventDrivenIBM.cpp -o polio_ed

polio_abc: EventDriven_parameters.hpp EventDriven_Sim.hpp main_EventDrivenIBM_abcsmc.cpp
	g++ $(FLAGS) $(INCLUDE) -I$(SQLDIR) main_EventDrivenIBM_abcsmc.cpp -o polio_abc $(ABC_LIB) $(GSL_LIB)

debug: main_Gillespie.cpp
	g++ $(FLAGS) main_Gillespie.cpp -o debug

polio: main.cpp
	g++ $(FLAGS) main.cpp -o polio

edma: EventDriven_MassAction_Sim.h mass_action_SIR.cpp
	g++ $(FLAGS) mass_action_SIR.cpp -o edma

pde: PDE_Simulator.cpp
	g++ $(FLAGS) PDE_Simulator.cpp -o pde

pde_abc: libabc PDE_Simulator_abc.cpp
	g++ $(FLAGS) $(INCLUDE) -I$(SQLDIR) PDE_Simulator_abc.cpp -o pde_abc $(ABC_LIB) $(GSL_LIB)

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

polio_boxcar: Boxcar_Model/Polio_boxcar_simulator.cpp Boxcar_Model/Polio_boxcar_model_extended_waning.h Boxcar_Model/DIFFEQ_SIM.h Boxcar_Model/Polio_boxcar_parameters.h
	g++ $(FLAGS) $(INCLUDE) Boxcar_Model/Polio_boxcar_simulator.cpp -o polio_boxcar $(GSL_LIB)

polio_boxcar_cel: Boxcar_Model/Polio_boxcar_simulator.cpp Boxcar_Model/Polio_boxcar_model_extended_waning.h Boxcar_Model/DIFFEQ_SIM.h Boxcar_Model/Polio_boxcar_parameters.h Boxcar_Model/adjustable_boxcar_model.h
	g++ $(FLAGS) Boxcar_Model/Polio_boxcar_simulator.cpp -o polio_boxcar_cel -lgsl

boxcar_abc: libabc Boxcar_Model/Polio_boxcar_simulator.cpp
	g++ $(FLAGS) $(INCLUDE) -I$(SQLDIR) Boxcar_Model/Polio_boxcar_simulator.cpp -o boxcar_abc $(ABC_LIB) $(GSL_LIB)


$(output_dir)/polio_data:
	mkdir $@
