CXX = g++
CXXFLAGS = -I$$GUROBI_HOME/include -O3 
LDLIBS = -lgmp -lboost_filesystem -lmpfr -lgurobi_c++ -lgurobi110 -lboost_program_options
LDFLAGS = -L$$GUROBI_HOME/lib

OBJS = main.o helpers.o network.o circle_intersection.o solver.o graph_helpers.o arrangement.o

main: $(OBJS)
	$(CXX) $(CXXFLAGS) -o main $(LDFLAGS) $(OBJS) $(LDLIBS)
	
main.o: main.cpp network.hpp

helpers.o: helpers.hpp helpers.cpp

network.o: network.hpp helpers.hpp circle_intersection.hpp solver.hpp graph_helpers.hpp arrangement.hpp network.cpp

circle_intersection.o: global.hpp circle_intersection.hpp circle_intersection.cpp

solver.o: global.hpp solver.hpp graph_helpers.hpp solver.cpp

graph_helpers.o: global.hpp graph_helpers.hpp graph_helpers.cpp

arrangement.o: arrangement.hpp arrangement.cpp

clean:
	rm -f main $(OBJS) generate generate_networks.o

# program to generate instances
generate: generate_networks.o graph_helpers.o
	g++ -Wall -o generate generate_networks.o graph_helpers.o -lgmp -lboost_filesystem -lmpfr -lboost_program_options
	
generate_networks.o: generate_networks.cpp global.hpp graph_helpers.hpp
