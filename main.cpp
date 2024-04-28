#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>

#include "network.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#define EPS 0.000001
#define DEFAULT_DISASTER_INCREMENT 0.1
#define DEFAULT_CSV_OUT "out.csv"
#define DEFAULT_INSTANCE_FILE "all-instances.txt"

void single_instance(std::string filename, bool log);

std::vector<std::string> process_instance_file(std::string filename);

int main(int argc, char** argv)
{
    // define and parse command line options, adapted from 
    // https://www.boost.org/doc/libs/1_84_0/libs/program_options/example/first.cpp
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("single-instance,s", po::value<std::string>(), "filename of a single instance to run")
        ("log,l", po::bool_switch()->default_value(false), "log progress to stdout. not advised for a batch")
        ("instance-file,i", po::value<std::string>()->default_value(DEFAULT_INSTANCE_FILE), "file containing names of batch instances to run")
        ("out-file,o", po::value<std::string>()->default_value(DEFAULT_CSV_OUT), "file to store csv results in for a batch")
        ("disaster-increment,d", po::value<double>()->default_value(DEFAULT_DISASTER_INCREMENT), "disaster incremement for batch testing")
        ("max-disaster,m", po::value<double>()->default_value(-1), "optional max disaster size")

    ;

    po::variables_map vm;        
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    bool log = vm["log"].as<bool>();

    // run a single instance
    if (vm.count("single-instance")) {
        single_instance(vm["single-instance"].as<std::string>(), log);
        return 0;
    } 

    // run a batch of instances
    std::vector<std::string> network_names = process_instance_file(vm["instance-file"].as<std::string>());
    std::ofstream outfile(vm["out-file"].as<std::string>());
    outfile << (CSV_HEADER + "total_time\n");

    bool instance_ran = true;
    double disaster_radius = EPS;
    double prev_disaster_radius = 0;

    bool max_disaster_set = vm["max-disaster"].as<double>() > 0;

    while (instance_ran) {
        if (max_disaster_set) {
            if (disaster_radius > vm["max-disaster"].as<double>()) {
                std::cout << "Hit max disaster " << vm["max-disaster"].as<double>() << std::endl;
                break;
            }
        }

        instance_ran = false;
        for (auto f : network_names) {
            Network network(f, log);
            double max_disaster_radius = network.get_max_disaster();
            // disaster too big
            if (prev_disaster_radius > max_disaster_radius) {
                continue;
            } 
            // disaster just too big, do eps less than max
            if (disaster_radius > max_disaster_radius) {
                network.set_disaster(max_disaster_radius - EPS);
            }
            // disaster is fine
            if (disaster_radius < max_disaster_radius) {
                network.set_disaster(disaster_radius);
            }

            auto start = std::chrono::high_resolution_clock::now();

            network.solve();

            auto stop = std::chrono::high_resolution_clock::now();
            

            instance_ran = true;
            outfile << network.report_statistics();
            outfile << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()) + ",\n";
            outfile << std::flush;

        }
        // increment disaster radius for next round
        if (instance_ran) std::cout << "completed disaster radius up to " << disaster_radius << std::endl;
        if (prev_disaster_radius == 0) {
            disaster_radius = vm["disaster-increment"].as<double>();
            prev_disaster_radius = EPS;
        } else {
            prev_disaster_radius = disaster_radius;
            disaster_radius += vm["disaster-increment"].as<double>();
        }
    }

    return 0;
}

// run a single instance
void single_instance(std::string filename, bool log) {
    Network test(filename, log);
    test.solve();
    test.construct_tikz_figures();

    std::cout << CSV_HEADER + "\n";
    std::cout << test.report_statistics() << std::endl << std::endl;
    
    std::cout << "n_nodes: " << test.n_nodes << std::endl;
    std::cout << "r_protect: " << test.r_protect << std::endl;

    std::cout << "orig edges size: " << test.orig_edges.size() << std::endl;



    std::cout << "edge overlaps size: " << test.edge_overlaps.size() << std::endl;
    std::cout << "feas edge size: " << test.feas_edges.size() << std::endl;

}

// process a file of input network names
std::vector<std::string> process_instance_file(std::string filename) {
    std::ifstream in(filename);
    std::string line;

    std::vector<std::string> network_names;
    if (in.fail()) {
        throw std::runtime_error("Instance file not found. Ensure the path is correct.");
    }
    while (in.peek() != EOF && in.peek() != '\n') {
        getline(in, line);
        network_names.push_back(line);
    }
    return network_names;
}