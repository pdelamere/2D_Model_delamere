%COMP=mpiifort
%FLAGS=-O4
COMP=mpif90
FLAGS=-O2 -mmacosx-version-min=10.6

OBJECTS=functions.o debug.o reactions.o varTypes.o global.o cfit.o rrfit.o dielectronic.o \
  timeStep.o readEmis.o inputs.o ftmix.o para.o output.o diffusion.o dimensions.o
MODS=functions.mod debug.mod reactions.mod vartypes.mod global.mod fitc.mod fitr.mod dielectronic.mod \
  timeStep.mod readEmis.mod inputs.mod ftmix.mod para.mod output.mod diffusion.mod dimensions.mod
EXE= torus

all: onebox

onebox: $(OBJECTS)
	$(COMP) $(FLAGS) -o $(EXE) onebox.f90 $(OBJECTS)

inputs.o:
	$(COMP) $(FLAGS) -c inputs.f90

timeStep.o: functions.o
	$(COMP) $(FLAGS) -c timeStep.f90

ftmix.o:
	$(COMP) $(FLAGS) -c ftmix.f90

readEmis.o: varTypes.o
	$(COMP) $(FLAGS) -c readEmis.f90

functions.o: reactions.o varTypes.o global.o cfit.o rrfit.o dielectronic.o debug.o para.o inputs.o
	$(COMP) $(FLAGS) -c functions.f90

para.o: varTypes.o
	$(COMP) $(FLAGS) -o para.o -c ParallelVars.f90

debug.o: varTypes.o
	$(COMP) $(FLAGS) -c debug.f90

varTypes.o: dimensions.o
	$(COMP) $(FLAGS) -c varTypes.f90

reactions.o:
	$(COMP) $(FLAGS) -c reactions.f90

global.o:
	$(COMP) $(FLAGS) -c global.f90

cfit.o:
	$(COMP) $(FLAGS) -c cfit.f90

rrfit.o:
	$(COMP) $(FLAGS) -c rrfit.f90

dielectronic.o: para.o
	$(COMP) $(FLAGS) -c dielectronic.f90

output.o:
	$(COMP) $(FLAGS) -c output.f90

diffusion.o:
	$(COMP) $(FLAGS) -c diffusion.f90

dimensions.o:
	$(COMP) $(FLAGS) -c dimensions.f90

clean:
	\rm *.o *.mod $(EXE)

try: all
	\rm *.o *.mod $(EXE)
