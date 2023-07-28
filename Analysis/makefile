#compile = /hpcapps/mpi/openmpi/gcc/1.6.5/bin/mpicc -w -std=c99 -O3 
compile = mpicc	-std=c99 -O3
#compile = gcc -std=c99 -O3
#link= -lm -std=c99
link= -lm

info: 		#before colon is the command we need to give with make
	@echo "type make analyse to create executables , make dataclean to clear all data directories and clean to clear all object files" 

object= main.o avg_and_print.o guu.o read_all_data.o msd_pos.o overlap.o S4Qt.o Fsqt.o
analyse: $(object)                            #first line checks for version
	$(compile) $(object) $(link) -o analyse
analyse1: $(object)                            #first line checks for version
	$(compile) $(object) $(link) -o analyse1

analysis.o: analysis.c global.h subroutine.h
	$(compile) -c analysis.c -o analysis.o

main.o: main.c global.h subroutine.h
	$(compile) -c main.c -o main.o
avg_and_print.o: avg_and_print.c global.h subroutine.h
	$(compile) -c avg_and_print.c -o avg_and_print.o

guu.o: guu.c global.h subroutine.h
	$(compile) -c guu.c -o guu.o

S4Qt.o: S4Qt.c global.h subroutine.h
	$(compile) -c S4Qt.c -o S4Qt.o
Fsqt.o: Fsqt.c global.h subroutine.h
	$(compile) -c Fsqt.c -o Fsqt.o
	
read_all_data.o: read_all_data.c global.h subroutine.h
	$(compile) -c read_all_data.c -o read_all_data.o


msd_pos.o: msd_pos.c global.h subroutine.h
	$(compile) -c msd_pos.c -o msd_pos.o

overlap.o: overlap.c global.h subroutine.h
	$(compile) -c overlap.c -o overlap.o

clean:
	@ rm -f $(object) analyse* *.o
	@ echo "All object files are now deleted"
dataclean:
	@ rm -rf data_0* $(object) analyse
	@ rm *.dat 
	@ echo "All object files and data directories are now deleted"
