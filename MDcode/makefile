#compile = /apps/openmpi/bin/mpicc -O3 -std=c99
compile = mpicc -O3 -std=c99  
#compile = gcc -O3
#link= -lm -std=c99
link= -lm  -lgsl -lgslcblas 

info: 		#before colon is the command we need to give with make
	@echo "type make run to create executables , make dataclean to clear all data directories and clean to clear all object files" 

object= extern.o main.o cell_list+forces.o initialise.o time_evolve.o print.o Parallel.o active.o

run: $(object)                            #first line checks for version
	$(compile) $(object) $(link) -o run
run1: $(object)                            #first line checks for version
	$(compile) $(object) $(link) -o run1
run2: $(object)                            #first line checks for version
	$(compile) $(object) $(link) -o run2
run3: $(object)                            #first line checks for version
	$(compile) $(object) $(link) -o run3

main.o: main.c global.h subroutine.h
	$(compile) -c main.c -o main.o

print.o: print.c global.h subroutine.h
	$(compile) -c print.c -o print.o

	
cell_list+forces.o: cell_list+forces.c global.h subroutine.h
	$(compile) -c cell_list+forces.c -o cell_list+forces.o


extern.o: extern.c global.h subroutine.h
	$(compile) -c extern.c -oextern.o

time_evolve.o: time_evolve.c global.h subroutine.h
	$(compile) -c time_evolve.c -o time_evolve.o

initialise.o: initialise.c global.h subroutine.h
	$(compile) -c initialise.c -o initialise.o

Parallel.o: Parallel.c global.h subroutine.h
	$(compile) -c Parallel.c -o Parallel.o
 
active.o: active.c global.h subroutine.h
	$(compile) -c active.c -o active.o
clean:
	@ rm -f $(object) run*
	@ echo "All object files are now deleted"
dataclean:
	@ rm -rf data_0* $(object) run*
	@ rm *.dat 
	@ echo "All object files and data directories are now deleted"
