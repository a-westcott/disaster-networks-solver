#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include <set>
#include <chrono>
#include <unordered_set>

#include "global.hpp"

struct SolverData {
    int statistic_sim_destruct_after_size;
    int statistic_sim_destruct_callback_size;
    std::chrono::milliseconds statistic_time_solver_preprocess;
    std::chrono::milliseconds statistic_time_gurobi_total;
    std::set<Edge> final_edges;
    double objective;
    int statistic_max_components;
    double statistic_time_gurobi_solve;
};

SolverData gurobi_solve(
    bool log,
    int n_nodes, 
    std::set<Edge> orig_edges,
    std::vector<Edge> feas_edges,
    std::vector<double> edge_lengths,
    std::unordered_set<std::pair<Edge,Edge>> edge_overlaps,
    std::set<std::set<Edge>> intersection_edge_sets,
    std::set<Edge> original_bridges
);  

#endif