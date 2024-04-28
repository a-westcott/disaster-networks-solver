#include "helpers.hpp"

#include <list>
#include <cassert>

// #include <CGAL/Lazy_exact_nt.h>
#include <CGAL/General_polygon_set_2.h>


typedef Traits_2::Curve_2                                       Curve_2;
typedef Traits_2::X_monotone_curve_2                            X_monotone_curve_2;


// Construct a polygon from a circle.
// from https://doc.cgal.org/4.10.2/Boolean_set_operations_2/Boolean_set_operations_2_2circle_segment_8cpp-example.html
Polygon_2 construct_polygon (const Circle_2& circle)
{
  // Subdivide the circle into two x-monotone arcs.
  Traits_2 traits;
  Curve_2 curve (circle);
  std::list<CGAL::Object>  objects;
  traits.make_x_monotone_2_object() (curve, std::back_inserter(objects));
  CGAL_assertion (objects.size() == 2);
  // Construct the polygon.
  Polygon_2 pgn;
  X_monotone_curve_2 arc;
  std::list<CGAL::Object>::iterator iter;
  for (iter = objects.begin(); iter != objects.end(); ++iter) {
    CGAL::assign (arc, *iter);
    pgn.push_back (arc);
  }
  return pgn;
}


// Construct a polygon from a rectangle.
// modified from https://doc.cgal.org/4.10.2/Boolean_set_operations_2/Boolean_set_operations_2_2circle_segment_8cpp-example.html
Polygon_2 construct_polygon (const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4)
{
  Polygon_2 pgn;
  X_monotone_curve_2 s1(p1, p2);    pgn.push_back(s1);
  X_monotone_curve_2 s2(p2, p3);    pgn.push_back(s2);
  X_monotone_curve_2 s3(p3, p4);    pgn.push_back(s3);
  X_monotone_curve_2 s4(p4, p1);    pgn.push_back(s4);

  return pgn;
}

// construct a rectangular region around two points representing a vulnerable segment of an edge
Polygon_2 construct_vulnerable_box(Point_2 p1, Point_2 p2, double r_disaster) {
    double x1, x2, y1, y2;
    x1 = CGAL::to_double(p1.x());
    y1 = CGAL::to_double(p1.y());
    x2 = CGAL::to_double(p2.x());
    y2 = CGAL::to_double(p2.y());

    // find the offset of the line in the perpendicular direction
    double dy, dx, length;
    dx = (y2 - y1);
    dy = (x1 - x2);
    length = std::sqrt(dx*dx + dy*dy);
    dx=dx*r_disaster/length;
    dy=dy*r_disaster/length;

    // get the points of the rectangle
    Point_2 a,b,c,d;
    a = Point_2(x1 + dx, y1 + dy);
    b = Point_2(x2 + dx, y2 + dy);
    c = Point_2(x2 - dx, y2 - dy);
    d = Point_2(x1 - dx, y1 - dy);

    return construct_polygon(a, b, c, d);

}

// euclidean distance between points
double distance(Point_2 u, Point_2 v) {
    double x1 = CGAL::to_double(u.x());
    double y1 = CGAL::to_double(u.y());
    double x2 = CGAL::to_double(v.x());
    double y2 = CGAL::to_double(v.y());

    return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));

}

void print_edge_set(std::set<Edge> s) {
    for (auto e : s) {
      std::cout << e.first << "-" << e.second << " ";
    }
    std::cout << std::endl;
}