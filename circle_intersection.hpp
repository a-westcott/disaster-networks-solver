#ifndef CIRCLE_INTERSECTION_HPP
#define CIRCLE_INTERSECTION_HPP

#include <vector>
#include <utility>
#include <set>
#include <map>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include "global.hpp"

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;


std::map<Edge, std::set<std::pair<Point_2, Point_2>>> determine_vul_pieces(std::vector<Point_2> node_points, std::set<Edge> feas_edges, double r_protect);

#endif