#include "circle_intersection.hpp"

#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

typedef CGAL::Exact_circular_kernel_2                        Kernel_c;
typedef CGAL::Point_2<Kernel_c>                              Point_2_c;
typedef CGAL::Circle_2<Kernel_c>                             Circle_2_c;
typedef CGAL::Circular_arc_2<Kernel_c>                       Circular_arc_2_c;
typedef CGAL::Line_arc_2<Kernel_c>                           Line_arc_2_c;
typedef CGAL::Gps_circle_segment_traits_2<Kernel_c>          Traits_2_c;

using namespace std;

// crazy types to get intersection of line with circle
typedef CGAL::CK2_Intersection_traits<Kernel_c, Circle_2_c, Line_arc_2_c>::type CircleIntersectionResult;
typedef pair<CGAL::Circular_arc_point_2<Kernel_c>, unsigned int> CircleIntersectPair;


// determine the sections of edges which are vulnerable
map<Edge, set<pair<Point_2, Point_2>>> determine_vul_pieces(vector<Point_2> node_points, set<Edge> feas_edges, double r_protect) {
    
    // initialise protection zones and edge arcs in circular kernel 
    int n_nodes = node_points.size();
    vector<Point_2_c> node_points_c = vector<Point_2_c>(n_nodes);

    vector<Circle_2_c> prot_zones = vector<Circle_2_c>(n_nodes);
    int i=0;
    for (Point_2 p : node_points) {
        double x,y; 
        x = CGAL::to_double(p.x());
        y = CGAL::to_double(p.y());
        node_points_c[i] = Point_2_c(x, y);
        prot_zones[i] = Circle_2_c(node_points_c[i], r_protect*r_protect);
        i++;
    }

    std::map<Edge, Line_arc_2_c> edge_arcs;
    for (auto e : feas_edges) {
            int u = e.first;
            int v = e.second;
            Point_2_c u_c = node_points_c[u];
            Point_2_c v_c = node_points_c[v];
            edge_arcs[e] = Line_arc_2_c(u_c, v_c);
    }


    vector<CircleIntersectionResult> pts;
    map<Edge, set<pair<Point_2_c, Point_2_c>>> vul_pieces_c;
    for (auto e : feas_edges) {
        int u = e.first;
        int v = e.second;

        // trim edge to section between endpoint protection zones.
        // by construction sure to have one intersection with each 
        pts.clear();
        intersection(prot_zones[u], edge_arcs[e], back_inserter(pts));
        intersection(prot_zones[v], edge_arcs[e], back_inserter(pts));

        double u_x, u_y, v_x, v_y;
        CircleIntersectPair sq_coord;
        CGAL::assign (sq_coord, pts[0]);
        u_x = CGAL::to_double(sq_coord.first.x());
        u_y = CGAL::to_double(sq_coord.first.y());
        CGAL::assign (sq_coord, pts[1]);
        v_x = CGAL::to_double(sq_coord.first.x());
        v_y = CGAL::to_double(sq_coord.first.y());

        vul_pieces_c[e].insert(pair(Point_2_c(u_x, u_y), Point_2_c(v_x, v_y)));

        // check all other protection zones
        for (int prot_zone_i=0; prot_zone_i<n_nodes; prot_zone_i++) {
            if (prot_zone_i == u || prot_zone_i == v) {continue;} // already done endpoints
            
            // check all pieces our edge is currently in.
            // a protection zone can intersect at most one piece
            for (auto piece : vul_pieces_c[e]) {
                Line_arc_2_c segment = Line_arc_2_c(piece.first, piece.second);
                pts.clear();
                intersection(prot_zones[prot_zone_i], segment, back_inserter(pts));

                if (pts.size() < 2) {
                    pts.clear();
                    continue;
                } else {
                    double x1, y1, x2, y2;
                    CircleIntersectPair sq_coord;
                    CGAL::assign (sq_coord, pts[0]);
                    x1 = CGAL::to_double(sq_coord.first.x());
                    y1 = CGAL::to_double(sq_coord.first.y());
                    CGAL::assign (sq_coord, pts[1]);
                    x2 = CGAL::to_double(sq_coord.first.x());
                    y2 = CGAL::to_double(sq_coord.first.y());

                    Point_2_c end1 = Point_2_c(x1, y1);
                    Point_2_c end2 = Point_2_c(x2, y2);

                    pair<Point_2_c, Point_2_c> new_piece1, new_piece2;
                    if (CGAL::has_smaller_distance_to_point(piece.first, end1, end2)) {
                        new_piece1 = pair(piece.first, end1);
                        new_piece2 = pair(piece.second, end2);
                    } else {
                        new_piece1 = pair(piece.first, end2);
                        new_piece2 = pair(piece.second, end1);
                    }

                    vul_pieces_c[e].erase(piece);
                    vul_pieces_c[e].insert(new_piece1);
                    vul_pieces_c[e].insert(new_piece2);

                    // this is the only piece that can be hit by this prot zone,
                    // so we can move on to the next prot zone
                    break; 
                }
            }
        }
    }
    // convert vul_pieces_c to vul_pieces and return
    map<Edge, set<pair<Point_2, Point_2>>> vul_pieces;
    for (auto const& [edge, pieces] : vul_pieces_c) {
        for (auto piece : pieces) {
            Point_2 u, v;
            double x1,y1,x2,y2;
            x1 = CGAL::to_double(piece.first.x());
            y1 = CGAL::to_double(piece.first.y());
            x2 = CGAL::to_double(piece.second.x());
            y2 = CGAL::to_double(piece.second.y());
            vul_pieces[edge].insert(pair(Point_2(x1, y1), Point_2(x2, y2)));
        }
    }

    return vul_pieces;

}