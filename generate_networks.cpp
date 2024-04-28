#define DEFAULT_PROT_RADIUS 1
#define DEFAULT_DIS_RADIUS 0.5
#define DEFAULT_INITIAL_GRID 5
#define DEFAULT_RETRIES_LIMIT 20
#define SEED 12345

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef CGAL::Delaunay_triangulation_2<Kernel>              Triangulation;
typedef Triangulation::Edge_iterator                        Edge_iterator;
typedef CGAL::Point_2<Kernel>                               Point;
typedef Kernel::Segment_2                                   Segment_2;


// just rolls off the tongue
typedef CGAL::Triangulation_2<Kernel, CGAL::Triangulation_data_structure_2<CGAL::Triangulation_vertex_base_2<Kernel>, CGAL::Triangulation_face_base_2<Kernel>>>::Finite_faces_iterator FaceIterator;


#include <cstdlib>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <set>
#include <random>

#include "global.hpp"
#include "graph_helpers.hpp"

using namespace std;

void generate_instance(int instance_number, 
                       int n_nodes, 
                       double square_size, 
                       string outfile_head,
                       double protection_radius,
                       double disaster_radius,
                       bool random_tree,
                       double edge_prop,
                       int index_offset);

Point random_coord(double grid_size);

set<Edge> identify_edges(Triangulation t, FaceIterator f, vector<Point> points);

Edge choose_from_set(set<Edge> edges);

pair<int, int> choose_two(int size);

vector<Point> generate_points(int n_nodes, double square_size, double protection_radius);

set<Edge> generate_edges_delaunay(Triangulation triangulation, vector<Point> points, int n_nodes, double edge_prop);

set<Edge> generate_edges_random_tree(vector<Point> points, int n_nodes, int n_tri_edges, double edge_prop);



int main(int ac, char* av[]) {
    // define and parse command line options, adapted from 
    // https://www.boost.org/doc/libs/1_84_0/libs/program_options/example/first.cpp
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("n-instances,i", po::value<int>(), "number of instances to generate")
        ("n-nodes,n", po::value<int>(), "number of nodes in each instance")
        ("square-size,s", po::value<double>()->default_value(DEFAULT_INITIAL_GRID), "initial grid dimension")
        ("protection-radius,p", po::value<double>()->default_value(DEFAULT_PROT_RADIUS), "protection zone radius")
        ("outfile-head,o", po::value<string>()->default_value("instances/network"), "head of output file path")
        ("disaster-radius,d", po::value<double>()->default_value(DEFAULT_DIS_RADIUS), "disaster radius")
        ("random-tree-based,r", po::bool_switch()->default_value(false), "generate edges starting with a random tree, as opposed to a Deluauny triangulation")
        ("prop-of-triangulation,t", po::value<double>(), "number of edges to have beyond a tree, as a proportion between a tree and a Deluauny triangulation.")
        ("index-offset", po::value<int>()->default_value(0), "filename number offset")

    ;

    po::variables_map vm;        
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << "\n";
        return 0;
    }

    if (!vm.count("n-instances")) {
        cout << "n-instances was not set.\n";
        return 0;
    } 

    if (!vm.count("n-nodes")) {
        cout << "n-nodes was not set.\n";
        return 0;
    } 


    // set random seed and get distribution
    srand(SEED*vm["n-nodes"].as<int>() + vm["random-tree-based"].as<bool>() + (int)(vm["prop-of-triangulation"].as<double>()*100));



    for (int i=0; i<vm["n-instances"].as<int>(); i++) {
        generate_instance(i, 
                          vm["n-nodes"].as<int>(), 
                          vm["square-size"].as<double>(), 
                          vm["outfile-head"].as<string>(), 
                          vm["protection-radius"].as<double>(), 
                          vm["disaster-radius"].as<double>(), 
                          vm["random-tree-based"].as<bool>(),
                          vm["prop-of-triangulation"].as<double>(),
                          vm["index-offset"].as<int>());
    }

    return 0;
}

// generate the network of a problem instance
void generate_instance(int instance_number, int n_nodes, double square_size, string outfile_head, double protection_radius, double disaster_radius, bool random_tree, double edge_prop, int index_offset) {
    vector<Point> accepted_points = generate_points(n_nodes, square_size, protection_radius);
    
    
    // triangulate, and identify which edges are in the triangulation
    Triangulation triangulation(accepted_points.begin(), accepted_points.end());
    set<Edge> edges;

    if (random_tree) {
        int n_tri_edges = n_nodes + (triangulation.number_of_faces() + 1) - 2; // method gives number of finite faces
        edges = generate_edges_random_tree(accepted_points, n_nodes, n_tri_edges, edge_prop);
    } else {
        edges = generate_edges_delaunay(triangulation, accepted_points, n_nodes, edge_prop);
    }


    // write out to file
    string filename = outfile_head + (random_tree ? "-tree" : "-deluauny") + "-n-" + to_string(n_nodes) + "-i-" + to_string(instance_number + index_offset) + ".txt";
    ofstream outfile(filename);
    outfile << "Nodes: " + to_string(n_nodes) + "\n";
    outfile << "Protection radius: " + to_string(protection_radius) + "\n";
    outfile << "Disaster radius: " + to_string(disaster_radius) + "\n";
    
    for (int i=0; i<n_nodes; i++) {
        outfile << to_string(i) + ": " + to_string(CGAL::to_double(accepted_points[i].x())) + " " + to_string(CGAL::to_double(accepted_points[i].y())) + "\n";
    }
    
    outfile << "Edges\n";
    for (Edge e : edges) {
        outfile << "(" + to_string(e.first) + "," + to_string(e.second) + ")\n";
    }


    // store filename
    ofstream network_names;
    network_names.open("all-instances.txt", ios_base::app);
    network_names << filename << "\n";

}

// give a random Point in a grid of `grid_size`
Point random_coord(double grid_size) {
    double x = (static_cast<double>(std::rand()) / RAND_MAX)*grid_size;
    double y = (static_cast<double>(std::rand()) / RAND_MAX)*grid_size;

    return Point(x, y);
}

// given a face of the triangulation, extract edges from it encoded by the indices 
// of its vertices in the original vector
set<Edge> identify_edges(Triangulation t, FaceIterator f, vector<Point> points) {
    set<Edge> edges;
    for (int i=0; i<=2; i++) {
        auto e = t.segment(f, i);

        Point u = e.point(0);
        Point v = e.point(1);
        
        int u_i = find(points.begin(), points.end(), u) - points.begin();
        int v_i = find(points.begin(), points.end(), v) - points.begin();

        Edge edge = u_i < v_i ? pair(u_i, v_i) : pair(v_i, u_i);
        edges.insert(edge);
    }
    return edges;
}

Edge choose_from_set(set<Edge> collection) {
    int n = rand() % collection.size();
    auto it = collection.begin();
    // 'advance' the iterator n times
    advance(it, n);
    return *it;
}

pair<int, int> choose_two(int size) {
    int i1 = (rand() % (size) );
    int i2 = (rand() % (size-1));

    if (i2 >= i1) i2++;
    return pair(i1, i2);
}

// generate coordinates of nodes within a grid
vector<Point> generate_points(int n_nodes, double square_size, double protection_radius) {
    vector<Point> accepted_points;
    int retries = 0;

    while ((int)accepted_points.size() < n_nodes) {
        // if we can't fit a point in, make the grid bigger
        if (retries >= DEFAULT_RETRIES_LIMIT) {
            accepted_points.clear(); 
            square_size+=1;
            retries = 0;
            continue;

        }

        Point candidate = random_coord(square_size);
        // check if our new point is too close to an existing one
        bool overlap = false;
        for (Point accepted : accepted_points) {
            double c_x = CGAL::to_double(candidate.x());
            double a_x = CGAL::to_double(accepted.x());            
            double c_y = CGAL::to_double(candidate.y());
            double a_y = CGAL::to_double(accepted.y());

            double distance = sqrt((c_x - a_x)*(c_x - a_x) + (c_y - a_y)*(c_y - a_y));
            if (distance <= 2*protection_radius) {
                retries++;
                overlap = true;
                break;
            }
        }
        if (overlap) continue;

        accepted_points.push_back(candidate);
    }

    // sort points by x coord so we can make sense of diagrams if we want
    sort(accepted_points.begin(), accepted_points.end(), [](Point u, Point v) {return u.x() < v.x(); });

    return accepted_points;
}

// determine edges for a triangulation based instance
set<Edge> generate_edges_delaunay(Triangulation triangulation, vector<Point> points, int n_nodes, double edge_prop) {
    set<Edge> edges;
    for (FaceIterator f : triangulation.finite_face_handles()) {
        set<Edge> new_edges = identify_edges(triangulation, f, points);
        edges.insert(new_edges.begin(), new_edges.end());
    }
    
    // determine which edges to cull
    int tree_edges = n_nodes - 1;
    int triangulated_edges = edges.size();
    int diff = triangulated_edges - tree_edges;
    int to_remove = diff - edge_prop*diff;

    for (int i=0; i<to_remove; i++) {
        set<Edge> bridges = determine_bridges(n_nodes, edges);
        set<Edge> removal_candidates;
        set_difference(edges.begin(), edges.end(), bridges.begin(), bridges.end(), inserter(removal_candidates, removal_candidates.end()));
        edges.erase(choose_from_set(removal_candidates));
    }
    return edges;
}

// determine edges for a tree based instance
set<Edge> generate_edges_random_tree(vector<Point> points, int n_nodes, int n_tri_edges, double edge_prop) {
    set<Edge> instance_edges, infeasible, cyclic;
    Segment_2* segments = new Segment_2[n_tri_edges];
    int n_components = n_nodes;

    while ((int)instance_edges.size() < n_nodes - 1) {
        int i1, i2;
        // pick our two points
        auto indices = choose_two(n_nodes);
        i1 = indices.first;
        i2 = indices.second;

        Edge e = i1 < i2 ? Edge(i1, i2) : Edge(i2, i1);
        if (infeasible.count(e) || cyclic.count(e)) continue;

        Point p1 = points[i1];
        Point p2 = points[i2];

        // check whether our new edge intersects the current ones
        Segment_2 seg(p1, p2);
        segments[((int)instance_edges.size())] = seg;
        if (do_curves_intersect(segments, segments + ((int)instance_edges.size()) + 1)) {
            infeasible.insert(e);
            continue;
        } 
        // check whether we intersect a node
        for (auto p : points) {
            if (p == p1 || p == p2) continue;
            if (seg.has_on(p)) {
                infeasible.insert(e);
                continue;
            }
        }
        // check whether we're joining two previously disconnected components
        instance_edges.insert(e);
        if (determine_n_components(instance_edges, n_nodes) == n_components) {
            instance_edges.erase(e);
            cyclic.insert(e);
            continue;
        }
        // if we're here we're happy with the edge to add, and can leave it in `instance_edges`
        n_components = determine_n_components(instance_edges, n_nodes);
    }

    // determine how many more edges to add
    int tree_edges = n_nodes - 1;
    int diff = n_tri_edges - tree_edges;
    int to_add = edge_prop*diff;

    while ((int)instance_edges.size() < n_nodes - 1 + to_add) {
        int i1, i2;
        // pick our two points
        auto indices = choose_two(n_nodes);
        i1 = indices.first;
        i2 = indices.second;

        Edge e = i1 < i2 ? Edge(i1, i2) : Edge(i2, i1);
        if (infeasible.count(e)) continue;

        Point p1 = points[i1];
        Point p2 = points[i2];

        // check whether our new edge intersects the current ones
        Segment_2 seg(p1, p2);
        segments[((int)instance_edges.size())] = seg;
        if (do_curves_intersect(segments, segments + ((int)instance_edges.size()) + 1)) {
            infeasible.insert(e);
            continue;
        } 
        // check whether we intersect a node
        for (auto p : points) {
            if (p == p1 || p == p2) continue;
            if (seg.has_on(p)) {
                infeasible.insert(e);
                continue;
            }
        }
        // otherwise we're happy to keep it
        instance_edges.insert(e);
    }


    return instance_edges;
}

