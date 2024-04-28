#ifndef ARRANGEMENT_HPP
#define ARRANGEMENT_HPP

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gps_circle_segment_traits_2.h>

#include <vector>
#include <map>

#include "global.hpp"



typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef CGAL::Point_2<Kernel>                               Point_2;

std::pair<std::map<std::set<Edge>, bool>, int> determine_simultaneous_destruction(std::map<Edge, std::set<std::pair<Point_2, Point_2>>> vul_pieces, double r_disaster, bool log);



#endif