#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <set>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include "global.hpp"

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::Circle_2                                        Circle_2;
typedef CGAL::Gps_circle_segment_traits_2<Kernel>               Traits_2;
typedef Traits_2::General_polygon_2                             Polygon_2;

Polygon_2 construct_polygon (const Circle_2& circle);
Polygon_2 construct_polygon (const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4);
Polygon_2 construct_vulnerable_box(Point_2 p1, Point_2 p2, double r_protect);
double distance(Point_2 u, Point_2 v);
void print_edge_set(std::set<Edge> s);

#endif