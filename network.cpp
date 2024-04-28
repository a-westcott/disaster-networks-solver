#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cassert>
#include <cmath>
#include <iterator>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/property_map.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/enum.h>

#include "helpers.hpp"
#include "network.hpp"
#include "solver.hpp"
#include "graph_helpers.hpp"
#include "arrangement.hpp"


typedef Traits_2::General_polygon_2                                                            Polygon_2;
typedef CGAL::Convex_hull_traits_adapter_2<Kernel, CGAL::Pointer_property_map<Point_2>::type > Convex_hull_traits_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                                                     Polygon_with_holes_2;
using namespace std;

// constructor
Network::Network(string filename_, bool log_) {
    filename = filename_;
    log = log_;
    if (log) std::cout << "Parsing input ..." << endl;
    parse_input(filename);
    if (log) std::cout << "Done! Computing segments ..." << endl;
    compute_segments();
    if (log) std::cout << "Done! Validating network ..." << endl;
    validate_network();
    if (log) std::cout << "Done! Determining feasibility ..." << endl;
    determine_max_disaster();
    if (log) std::cout << "Done!" << endl;

}

// solve instance
void Network::solve() {
    if (!disaster_feasible()) {
        throw runtime_error("Instance is infeasible. Disaster radius of " + to_string(r_disaster) + " given. Max radius is " + to_string(max_disaster_radius) + "\n");
    }

    if (log) std::cout << "Computing overlaps ..." << endl;

    auto start = chrono::high_resolution_clock::now();
    compute_overlaps();
    auto stop = chrono::high_resolution_clock::now();
    statistic_time_compute_overlaps = chrono::duration_cast<chrono::milliseconds>(stop - start);
    if (log) std::cout << "Done! Determining vul pieces ..." << endl;
    vul_pieces = determine_vul_pieces(node_points, feas_edges, r_protect);

    if (log) std::cout << "Done! Simultaneous destruction via arrangements ..." << endl;

    start = chrono::high_resolution_clock::now();
    
    auto tuple = determine_simultaneous_destruction(vul_pieces, r_disaster, log);
    vul_zone_intersection_nonempty = tuple.first;
    statistic_num_faces = tuple.second;


    stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    statistic_time_arrangement_sim_destruction = duration;
    if (log) std::cout << "Done! Simultaneous destruction determination duration: " << duration.count() << " milliseconds" << endl << "Calling solver." << endl << endl;
    
    call_solver(log);
}


// parse input network file
void Network::parse_input(string filename) {
    ifstream in(filename);
    string line;

    if (in.fail()) {
        throw runtime_error("File not found. Ensure the path is correct.");
    }

    // ignore "Nodes: ", read value
    getline(in, line, ' ');
    getline(in, line);
    n_nodes = stoi(line);

    // ignore "Protection radius: ", read value
    getline(in, line, ' ');
    getline(in, line, ' ');
    getline(in, line);
    r_protect = stod(line);

    // ignore "Disaster radius: ", read value
    getline(in, line, ' ');
    getline(in, line, ' ');
    getline(in, line);
    r_disaster = stod(line);

    // parse coordinates
    node_points = vector<Point_2>(n_nodes);
    prot_zones = vector<Circle_2>(n_nodes);
    for (int i=0; i<n_nodes; i++) {
        if (in.peek() == 'E') {
            throw runtime_error("The number of nodes specified does not match the number found\n");
        }
        int node;
        Point_2 p;
        // get node number
        getline(in, line, ':');
        node = stoi(line);
        
        // CGAL can read in the coordinates exactly if done in this manner
        in >> p;
        node_points[node] = p;
        // initialise the protection zones at the same time
        prot_zones[node] = Circle_2(p, r_protect*r_protect); // it wants squared radius for some reason

        // clear newline
        getline(in, line);

    }
    if (isdigit(in.peek())) {
        throw runtime_error("The number of nodes specified does not match the number found\n");
    }

    // collect edges
    set<Edge> edge_collection;

    // clear "Edges", possibly a newline (was previously happening here, I think incorrectly)
    getline(in, line);
    if (in.peek() != '(') getline(in, line);

    while (in.peek() != EOF) {
        int u,v;
        // first (
        getline(in, line, '(');
        // read first coord
        getline(in, line, ',');
        u = stoi(line);
        //read second coord
        getline(in, line, ')');
        v = stoi(line);
        // clear newline
        getline(in, line, '\n');


        if (u < v) {
            edge_collection.insert(pair(u,v));
        } else if (v < u) {
            edge_collection.insert(pair(v,u));
        } else {
            if (log) cout << "Self loop on " << u << " ignored" << endl;
        }

    }

    orig_edges = edge_collection;
}

// store edges as CGAL segments
void Network::compute_segments() {
    for (int u=0; u<n_nodes; u++) {
        for (int v = u+1; v < n_nodes; v++) {
            Edge e = pair(u,v);
            Point_2 u_c = node_points[u];
            Point_2 v_c = node_points[v];
            edge_segments[e] = Segment_2(u_c, v_c);
            
            // also collect all edges in complete graph here
            feas_edges.insert(e);
        }
    }
}

// ensure input valid
void Network::validate_network() {
    // ensure protection zones don't overlap
    for (int u=0; u<n_nodes; u++) {
        for (int v = u+1; v < n_nodes; v++) {
            assert(distance(node_points[u], node_points[v]) > 2*r_protect);
        }
    }
    // ensure original edges don't cross
    int n_edges = (int)orig_edges.size();
    Segment_2* segments = new Segment_2[n_edges];
    int i=0;
    for (auto e : orig_edges)  {
        segments[i++] = edge_segments[e];
    }

    assert(!do_curves_intersect(segments, segments+n_edges));

    delete[] segments;

    // ensure original graph is connected
    assert(determine_components(orig_edges, n_nodes).size() == 1);
}

// determine the maximum size of a disaster, such that a feasible solution exists
// subject to rounding errors
void Network::determine_max_disaster() {
    vector<Point_2> alt_out;
    CGAL::convex_hull_2(node_points.begin(), node_points.end(), std::back_inserter(alt_out));

    int n_ch_nodes = (int)alt_out.size();

    double min_angle = 1000;

    for (int i=0; i<n_ch_nodes; i++) {
        Point_2 p1,p2,p3;
        p1 = alt_out[i];
        p2 = alt_out[(i+1)%n_ch_nodes];
        p3 = alt_out[(i+2)%n_ch_nodes];
        double angle = atan2(CGAL::to_double(p3.y() - p2.y()), CGAL::to_double(p3.x() - p2.x()))
                     - atan2(CGAL::to_double(p1.y() - p2.y()), CGAL::to_double(p1.x() - p2.x()));
        
        // put angle in range (-pi, pi]
        if (angle > M_PI) angle -= 2*M_PI;
        else if (angle <= - M_PI) angle += 2*M_PI;
        // then take its absolute value, because we don't care which way the angle faces (or whatever the correct geometric interpretation is)
        angle = abs(angle);

        if (angle < min_angle) min_angle = angle;
    }

    max_disaster_radius = r_protect*sin(min_angle/2);
}

// determine whether the instance is feasible for the given disaster radius
bool Network::disaster_feasible() {
    return r_disaster < max_disaster_radius;
}

double Network::get_max_disaster() {
    return max_disaster_radius;
}

void Network::set_disaster(double disaster_radius) {
    r_disaster = disaster_radius;
}

// determine which edges in the complete graph overlap
void Network::compute_overlaps() {
    // first determine which augmentation edges are infeasible
    set<Edge> infeasible;
    for (auto aug : feas_edges) {
        if (orig_edges.count(aug)) continue;
        for (auto orig : orig_edges) {
            auto result = intersection(edge_segments[orig], edge_segments[aug]);
            if (result) {
                if (const Point_2* p = boost::get<Point_2 >(&*result)) {
                    // ignore intersection if its at a shared endpoint only
                    if ((*p == node_points[orig.first] || *p == node_points[orig.second]) && (*p == node_points[aug.first] || *p == node_points[aug.second])) continue;
                }
                infeasible.insert(aug);
                break;
            }

        }
    }
    // remove the infeasible edges
    filter_infeasible(infeasible);

    // determine pairwise intersections of feasible augmentation edges
    for (auto e1 = feas_edges.begin(); e1 != feas_edges.end(); e1++) {
        for (auto e2 = e1; e2 != feas_edges.end(); e2++) {
            if (e1 == e2) continue; // inelegant but I can't see how to add one to the pointer without changing e1

            auto result = intersection(edge_segments[*e1], edge_segments[*e2]);
            if (result) {
                if (const Point_2* p = boost::get<Point_2 >(&*result)) {
                    // ignore intersection if its at a shared endpoint only
                    if ((*p == node_points[(*e1).first] || *p == node_points[(*e1).second]) && (*p == node_points[(*e2).first] || *p == node_points[(*e2).second])) continue;
                    edge_overlaps.insert(pair(*e1, *e2));
                }
            }
        }
    }
}


void Network::filter_infeasible(set<Edge> infeasible) {
    // remove any possible new edge which overlaps an original edge
    set<Edge> difference;
    set_difference(feas_edges.begin(), feas_edges.end(), infeasible.begin(), infeasible.end(), inserter(difference, difference.end()));
    feas_edges = difference;
}

// call gurobi solver
void Network::call_solver(bool log) {
    // first prepare the necessary pieces
    vector<Edge> feas_edge_vector(feas_edges.begin(), feas_edges.end());
    vector<double> edge_lengths(feas_edge_vector.size());
    // compute lengths
    for (long unsigned int i = 0; i < feas_edge_vector.size(); i++) {
        int u = feas_edge_vector[i].first;
        int v = feas_edge_vector[i].second;
        edge_lengths[i] = distance(node_points[u],node_points[v]);
    }
    // construct overlap vector
    set<set<Edge>> intersection_edge_sets;
    for (auto const& [edge_set, flag] : vul_zone_intersection_nonempty) {
        if (flag) {
            intersection_edge_sets.insert(edge_set);
        }
    }

    statistic_sim_destruct_before_size = (int)intersection_edge_sets.size();
    
    original_bridges = determine_bridges(n_nodes, orig_edges);

    SolverData solver_data;
    solver_data = gurobi_solve(log, n_nodes, orig_edges, feas_edge_vector, edge_lengths, edge_overlaps, intersection_edge_sets, original_bridges);    

    statistic_sim_destruct_after_size = solver_data.statistic_sim_destruct_after_size;
    statistic_sim_destruct_callback_size = solver_data.statistic_sim_destruct_callback_size;
    statistic_time_solver_preprocess = solver_data.statistic_time_solver_preprocess;
    statistic_time_gurobi_total = solver_data.statistic_time_gurobi_total;
    statistic_time_gurobi_solve = solver_data.statistic_time_gurobi_solve;
    statistic_max_components = solver_data.statistic_max_components;


    final_edges = solver_data.final_edges;
}

// euclidean length of a set of edges
double Network::length(std::set<Edge> edges) {
    double len = 0;
    for (auto e : edges) {
        len += distance(node_points[e.first], node_points[e.second]);
    }
    return len;
}


// format edges in a way that's easy to parse when reading the results csv
string format_edges_for_csv(set<Edge> edges) {
    string out = "";
    for (Edge e : edges) {
        out += to_string(e.first) + "-" + to_string(e.second) + ";";
    }
    // remove last semicolon
    out.pop_back();
    return out;
}

// header for csv file
string CSV_HEADER = "filename,n_nodes,n_orig_edges,n_feas_edges,r_protect,r_disaster,num_arr_faces,sim_destruct_before_size,sim_destruct_after_size,sim_destruct_callback_size,max_components,compute_overlaps_time,compute_arrangement_time,solver_preprocess_time,gurobi_time_total,gurobi_time_solve,n_aug_edges,orig_length,total_length,aug_edge_length,final_edges,";

// give a csv string of statistics 
string Network::report_statistics() {
    string line = "";
    line += filename + ",";
    line += to_string(n_nodes) + ",";
    line += to_string(orig_edges.size()) + ",";
    line += to_string(feas_edges.size()) + ",";
    line += to_string(r_protect) + ",";
    line += to_string(r_disaster) + ",";


    line += to_string(statistic_num_faces) + ",";
    line += to_string(statistic_sim_destruct_before_size) + ",";
    line += to_string(statistic_sim_destruct_after_size) + ",";
    line += to_string(statistic_sim_destruct_callback_size) + ",";
    line += to_string(statistic_max_components) + ",";


    line += to_string(statistic_time_compute_overlaps.count()) + ",";
    line += to_string(statistic_time_arrangement_sim_destruction.count()) + ",";
    line += to_string(statistic_time_solver_preprocess.count()) + ",";
    line += to_string(statistic_time_gurobi_total.count()) + ",";
    line += to_string(statistic_time_gurobi_solve) + ",";


    line += to_string(final_edges.size() - orig_edges.size()) + ",";
    line += to_string(length(orig_edges)) + ",";
    line += to_string(length(final_edges)) + ",";
    line += to_string(length(final_edges) - length(orig_edges)) + ",";

    line += format_edges_for_csv(final_edges) + ",";

    return line;
}

// construct a tikz figure of the network and its augmentation
void Network::construct_tikz_figures(string outpath) {
    if (outpath == "") {
        outpath = "out/";
    }


    double scale = min(1.0, 18/CGAL::to_double((*(node_points.end() - 1)).x()));

    ofstream out_radii(outpath + "radii.tex");
    ofstream out_original(outpath + "original_network.tex");
    ofstream out_augmented(outpath + "augmented_network.tex");

    out_radii << "Protection zone radius: " + to_string(r_protect) + "\n";
    out_radii << "\nDisaster radius: " + to_string(r_disaster);
    out_radii << ". Sample scale disaster: \n";
    out_radii << "\\begin{tikzpicture}[scale=" + to_string(scale) + "]\n    \\draw[color=red!30, very thick](0,0) circle (" + to_string(r_disaster) + ");\n\\end{tikzpicture}\n";


    string preamble = "\\begin{tikzpicture}\n    [scale = " + to_string(scale) + ", terminal/.style={circle,draw=black,fill=black,fill opacity = 1,inner sep=0pt,minimum size=4pt}]\n";

    out_original << preamble;
    out_augmented << preamble;

    for (long unsigned int i=0; i<node_points.size(); i++) {
        Point_2 node = node_points[i];
        string node_line = "    \\node[terminal,label={above left: \\small " + to_string(i) + "}] (" + to_string(i) + ") at (" + to_string(CGAL::to_double(node.x())) + "," + to_string(CGAL::to_double(node.y())) + "){};\n";
        string prot_line = "    \\draw[color=black!30, very thick](" + to_string(CGAL::to_double(node.x())) + "," + to_string(CGAL::to_double(node.y())) + ") circle (" + to_string(r_protect) + ");\n";
        out_original << node_line << prot_line;
        out_augmented << node_line << prot_line;
    }

    for (auto e : orig_edges) {
        string edge_line = "    \\draw[-] (" + to_string(e.first) + ") to (" + to_string(e.second) + ");\n";
        out_original << edge_line;
        out_augmented << edge_line;
    }

    for (auto e : final_edges) {
        if (orig_edges.count(e)) continue;
        string edge_line = "    \\draw[-,blue] (" + to_string(e.first) + ") to (" + to_string(e.second) + ");\n";
        out_augmented << edge_line;
    }

    string last_bit = "\\end{tikzpicture}\n";
    out_original << last_bit;
    out_augmented << last_bit;

}
