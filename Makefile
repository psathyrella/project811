CC=g++ -O2  -g
CCFLAGS = `root-config --cflags` -I /afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms4/include
LDFLAGS = `root-config --libs` -L /afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms4/lib  -l RooFit -l RooFitCore -L $(PWD) -lProject
# note: RooClassFactory::makePdf("EtaPdf","theta","","1/(2*3.1415926535*sin(theta))");

LIB = libProject.a
EXE = test.exe
OBJ = MitStyleRemix.o EtaPdf_cc.so Simulator.o Detector.o hi.o

all : lib exe
exe : $(EXE)
obj : $(OBJ)

lib : $(LIB)
$(LIB) : $(OBJ)
	ar cr $@ $^

MitStyleRemix.o : MitStyleRemix.cc MitStyleRemix.h
	$(CC) $(CCFLAGS) -c MitStyleRemix.cc -o MitStyleRemix.o $(LDFLAGS)

Simulator.o : Simulator.cc Simulator.h
	$(CC) $(CCFLAGS) -c Simulator.cc -o Simulator.o $(LDFLAGS)

Detector.o : Detector.cc Detector.h
	$(CC) $(CCFLAGS) -c Detector.cc -o Detector.o $(LDFLAGS)

hi.o : hi.cc hi.h
	$(CC) $(CCFLAGS) -c hi.cc -o hi.o $(LDFLAGS)

EtaPdf_cc.so : EtaPdf.cc EtaPdf.h
	root -b -l -q EtaPdf.cc+; ls # without the ls, make thinks the root command fails

#%.o : %.cc
#	$(CC) $(CCFLAGS) -c $< -o $@ $(LDFLAGS)

%.exe : %.cc $(OBJ)
	$(CC) $(CCFLAGS) -Wl,-rpath=$(PWD) $< -o $@ $(LDFLAGS)

clean :
	rm -f $(EXE) $(OBJ)
