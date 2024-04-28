#ifndef GRAPH_HELPERS_HPP
#define GRAPH_HELPERS_HPP

#include <set>

#include "global.hpp"

std::set<std::set<int>> determine_components(std::set<Edge>, int n_nodes);
std::set<Edge> determine_bridges(int n_nodes, std::set<Edge> edges);
int determine_n_components(std::set<Edge>, int n_nodes);


#endif