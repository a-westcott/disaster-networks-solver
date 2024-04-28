#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/basic.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>


#include <iostream>

#include "arrangement.hpp"

typedef Kernel::FT                                          Number_type;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>           Traits;
typedef Traits::Curve_2                                     Curve;
typedef Traits::Rational_point_2                            Rational_point;
typedef Traits::Rational_segment_2                          Segment;
typedef Traits::Rational_circle_2                           Circle;

using namespace std;

// Define a functor for creating a label from sets of edges
struct Overlay_label {
  set<Edge> operator()(set<Edge> s1, set<Edge> s2) const
  {
    set<Edge> s;
    for (auto x : s1) {
        s.insert(x);
    }
    for (auto x : s2) {
        s.insert(x);
    }
    return s;
  }
};

typedef CGAL::Arr_face_extended_dcel<Traits, std::set<Edge>>                                               Dcel_res;
typedef CGAL::Arrangement_2<Traits, Dcel_res>                                                              Arrangement;
typedef CGAL::Arr_face_overlay_traits<Arrangement, Arrangement, Arrangement, Overlay_label>                Overlay_traits;


Arrangement construct_arrangement(Edge e, set<pair<Point_2, Point_2>> vul_segments, double r_disaster);

Arrangement recursive_determine(vector<Edge> edges, map<Edge, set<pair<Point_2, Point_2>>> *vul_pieces, double r_disaster, bool log);

// determine the sets of edges which can be simultaneously destroyed by some disaster
pair<map<set<Edge>, bool>,int> determine_simultaneous_destruction(map<Edge, set<pair<Point_2, Point_2>>> vul_pieces, double r_disaster, bool log) {
    map<set<Edge>, bool> vul_zone_intersection_nonempty;

    if (log) cout << vul_pieces.size() << " edges to process" << endl;

    Arrangement arr3;
    vector<Edge> edges;
    for (auto const& [edge, vul_segments] : vul_pieces) {
        edges.push_back(edge);
    }

    arr3 = recursive_determine(edges, &vul_pieces, r_disaster, log);


    // all the arrangements are overlaid, so extract the facial information 
    int n_faces = 0;
    for (auto fit = arr3.faces_begin(); fit != arr3.faces_end(); ++fit) {
        if (!fit->is_unbounded()) {
            vul_zone_intersection_nonempty[fit->data()] = true;
            n_faces++;
        }
    }

    return pair(vul_zone_intersection_nonempty, n_faces);
}


// construct arrangement of vul zone, labelling bounded faces with the edge
Arrangement construct_arrangement(Edge e, set<pair<Point_2, Point_2>> vul_segments, double r_disaster) {
    Arrangement arr;
    Point_2 p1, p2;

    for (auto pair : vul_segments) {
        p1 = pair.first;
        p2 = pair.second;

        // circle centers
        Rational_point p1_rat(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()));
        Rational_point p2_rat(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));

        // determine vulnerable box
        double x1, x2, y1, y2;
        x1 = CGAL::to_double(p1.x());
        y1 = CGAL::to_double(p1.y());
        x2 = CGAL::to_double(p2.x());
        y2 = CGAL::to_double(p2.y());

        double dy, dx, length;
        dx = (y2 - y1);
        dy = (x1 - x2);
        length = std::sqrt(dx*dx + dy*dy);
        dx=dx*r_disaster/length;
        dy=dy*r_disaster/length;

        Rational_point a,b,c,d;
        a = Rational_point(x1 + dx, y1 + dy);
        b = Rational_point(x2 + dx, y2 + dy);
        c = Rational_point(x2 - dx, y2 - dy);
        d = Rational_point(x1 - dx, y1 - dy);

        // insert circles and box
        insert(arr, Curve(Circle(p1_rat, Number_type(r_disaster*r_disaster))));
        insert(arr, Curve(Circle(p2_rat, Number_type(r_disaster*r_disaster))));

        Segment s1(a, b);
        Segment s2(b, c);
        Segment s3(c, d);
        Segment s4(d, a);
        insert(arr, s1);
        insert(arr, s2);
        insert(arr, s3);
        insert(arr, s4);
    }
    
    // label the bounded faces
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (fit->is_unbounded()) continue;
        set<Edge> s;
        s.insert(e);
        fit->set_data(s);
    }

    return arr;

}

// recursively determine the overlay of the arrangements of vulnerable zones
Arrangement recursive_determine(vector<Edge> edges, map<Edge, set<pair<Point_2, Point_2>>> *vul_pieces, double r_disaster, bool log) {
    if (edges.size() == 1) {
        return construct_arrangement(edges[0], (*vul_pieces)[edges[0]], r_disaster);
    } else if (edges.size() == 2) {
        Arrangement arr1 = construct_arrangement(edges[0], (*vul_pieces)[edges[0]], r_disaster);
        Arrangement arr2 = construct_arrangement(edges[1], (*vul_pieces)[edges[1]], r_disaster);
        Overlay_traits overlay_traits;
        Arrangement arr3;
        CGAL::overlay(arr1, arr2, arr3, overlay_traits);
        return arr3;
    }
    int midpoint = ((int)edges.size())/2;
    vector<Edge> first(edges.begin(), edges.begin() + midpoint);
    vector<Edge> second(edges.begin() + midpoint, edges.end());
    Arrangement arr1 = recursive_determine(first, vul_pieces, r_disaster, log);
    Arrangement arr2 = recursive_determine(second, vul_pieces, r_disaster, log);
    Overlay_traits overlay_traits;
    Arrangement arr3;
    CGAL::overlay(arr1, arr2, arr3, overlay_traits);

    if (log && (int)edges.size() > 20) cout << "completed arrangement of size " << edges.size() << endl;

    return arr3;

}