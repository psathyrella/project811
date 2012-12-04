CC=g++ -O2  -g
CCFLAGS = `root-config --cflags` -I /afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms4/include
LDFLAGS = `root-config --libs` -L /afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms4/lib  -l RooFit -l RooFitCore
# note: RooClassFactory::makePdf("EtaPdf","theta","","1/(2*3.1415926535*sin(theta))");

EXE = test.exe
OBJ = MitStyleRemix.o EtaPdf_cc.so Simulator.o
#test.cc -o bin/test.exe
all : obj exe
exe : $(EXE)
obj : $(OBJ)

MitStyleRemix.o : MitStyleRemix.cc MitStyleRemix.h
	$(CC) $(CCFLAGS) -c MitStyleRemix.cc -o MitStyleRemix.o $(LDFLAGS)

Simulator.o : Simulator.cc Simulator.h
	$(CC) $(CCFLAGS) -c Simulator.cc -o Simulator.o $(LDFLAGS)

EtaPdf_cc.so : EtaPdf.cc EtaPdf.h
	root -b -l -q EtaPdf.cc+; ls # without the ls, make thinks the root command fails

#%.o : %.cc
#	$(CC) $(CCFLAGS) -c $< -o $@ $(LDFLAGS)

%.exe : %.cc MitStyleRemix.o EtaPdf_cc.so Simulator.o
	$(CC) $(CCFLAGS) -Wl,-rpath=$(PWD) $< -o $@ $(LDFLAGS) MitStyleRemix.o EtaPdf_cc.so Simulator.o

clean :
	rm -f test.exe MitStyleRemix.o EtaPdf_cc.so Simulator.o
