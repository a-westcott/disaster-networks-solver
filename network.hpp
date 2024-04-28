#ifndef NETWORK_HPP
#define NETWORK_HPP

#include <vector>
#include <utility>
#include <set>
#include <string>
#include <map>
#include <chrono>
#include <unordered_set>

#include "circle_intersection.hpp"

#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Segment_2.h>
// #include <CGAL/Arr_dcel_base.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef CGAL::Point_2<Kernel>                               Point_2;
typedef CGAL::Circle_2<Kernel>                              Circle_2;
typedef CGAL::Gps_circle_segment_traits_2<Kernel>           Traits_2;
typedef CGAL::General_polygon_set_2<Traits_2>               Polygon_set_2;
typedef Kernel::Segment_2                                   Segment_2;



class Network {
private:
    std::string filename;
    double r_disaster;
    std::map<Edge, Segment_2> edge_segments;
    std::vector<Circle_2> prot_zones;
    std::set<Edge> final_edges;
    std::set<Edge> original_bridges;
    std::map<Edge, std::set<Edge>> proximity_map;
    double max_disaster_radius;
    bool log;

    int statistic_num_faces;
    int statistic_sim_destruct_before_size;
    int statistic_sim_destruct_after_size;
    int statistic_sim_destruct_callback_size;
    std::chrono::milliseconds statistic_time_compute_overlaps;
    std::chrono::milliseconds statistic_time_arrangement_sim_destruction;
    std::chrono::milliseconds statistic_time_solver_preprocess;
    std::chrono::milliseconds statistic_time_gurobi_total;
    double statistic_time_gurobi_solve;
    int statistic_max_components;




    void parse_input(std::string filename);
    void validate_network();
    void compute_segments();
    void determine_max_disaster();
    void compute_overlaps();
    void filter_infeasible(std::set<Edge> infeasible);
    void determine_complete_vul_zones();
    void determine_vul_zone_overlaps();
    void call_solver(bool log);
    bool disaster_feasible();
    double length(std::set<Edge> edges);


public:
    Network(std::string filename_, bool log_);
    // TODO: probably should make these private
    int n_nodes;
    double r_protect;
    std::vector<Point_2> node_points;
    std::set<Edge> orig_edges;
    std::set<Edge> feas_edges;
    std::unordered_set<std::pair<Edge,Edge>> edge_overlaps; 
    std::map<Edge, std::set<std::pair<Point_2, Point_2>>> vul_pieces;
    std::map<Edge, Polygon_set_2> vul_zones;
    std::map<std::set<Edge>, Polygon_set_2> vul_zone_intersections;
    std::map<std::set<Edge>, bool> vul_zone_intersection_nonempty;


    double get_max_disaster();
    void set_disaster(double disaster_radius);
    void solve();
    std::string report_statistics();
    void construct_tikz_figures(std::string outpath="");
};

extern std::string CSV_HEADER;

#endif