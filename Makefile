CXXFLAGS := -Wall -g -O2 $(shell root-config --cflags) 
#-I$(shell cd $$CMSSW_BASE; scram tool tag boost include)
LDFLAGS := -g $(shell root-config --libs) -lJetMETObjects -LJetMETObjects/lib -Wl,-rpath,$(shell readlink -f JetMETObjects/lib)

dummy := $(shell [ -d obj ] || mkdir obj)
src := $(wildcard *.cpp)
obj := $(patsubst %.cpp,obj/%.o,$(src))

all: zz

obj/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $+

zz: $(obj)
	$(CXX) $(LDFLAGS) -o $@ $+


clean:
	rm -rf obj zz
