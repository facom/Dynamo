CC=g++
CFLAGS=-c -I. -Iutil -g
LFLAGS=-lm -lgsl -lgslcblas
PROGRAMS=\
bin2sim.out transform.out updatesim.out slice.out filter.out sort.out \
stats.out fit.out scalarmap.out map2dist.out vectormap.out derintegral.out

all:dynamo.h
	make $(PROGRAMS)

%.out:dynamo.o functions.o %.o
	$(CC) $^ $(LFLAGS) -o $@

%.o:%.cpp dynamo.h functions.h
	$(CC) $< $(CFLAGS) -o $@

%.h:%.cpp %.hpp
	./.hpp2h $(@:.h=)

cleandata:
	@echo -n "Cleaning data..."
	@rm -rf *.dst *.slc *.fil *.map *.sim *.trn *.fit
	@echo "Done."

clean:
	@echo -n "Cleaning..."
	@rm -rf *.o *.out *~ 
	@echo "Done."

edit:
	emacs -nw *.cpp *.hpp run .hpp2h

updaterepo:
	@echo -n "Updating repo..."
	git commit -am "Commit..."
	git push origin master
