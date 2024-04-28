#include <algorithm>
#include <array>
#include <functional>
#include <iomanip>
#include <iostream>
#include <cassert>

#include <gurobi_c++.h>

#include "solver.hpp"
#include "graph_helpers.hpp"

using namespace std;

#define CALLBACK 1

// adapted from https://www.gurobi.com/documentation/current/examples/cb_cpp_cpp.html 
class mycallback: public GRBCallback
{
public:
    int numvars;
    GRBVar* vars;
    int n_nodes;
    vector<Edge> feas_edges;
    vector<vector<bool>> simultaneous_destruction_membership;
    vector<vector<set<int>>> callback_simultaneous_destruction_components;
    mycallback(int xnumvars, GRBVar* xvars, int n_nodes_, vector<Edge> feas_edges_, vector<vector<bool>> simultaneous_destruction_membership_, vector<vector<set<int>>> callback_simultaneous_destruction_components_) {
        numvars = xnumvars;
        vars = xvars;
        n_nodes = n_nodes_;
        feas_edges = feas_edges_;
        simultaneous_destruction_membership = simultaneous_destruction_membership_;
        callback_simultaneous_destruction_components = callback_simultaneous_destruction_components_;
    }
protected:
    // callback function for solver, determining if an omitted constraint should be added
    void callback () {

    try {
        if (where == GRB_CB_MIPSOL) {
            double* x = getSolution(vars, numvars);

            // collect edges currently in use by solution
            set<Edge> active_edges;
            for (int i=0; i<numvars; i++) {
                if (x[i] > 0.5) {
                    active_edges.insert(feas_edges[i]);
                }
            }

            assert(callback_simultaneous_destruction_components.size() == simultaneous_destruction_membership.size());

            for (long unsigned int i=0; i<simultaneous_destruction_membership.size(); i++) {
                vector<bool> row = simultaneous_destruction_membership[i];
                for (auto X : callback_simultaneous_destruction_components[i]) {
                    set<Edge> delta_X_survive;
                    int j = 0;
                    bool cont = false;
                    for (auto e : feas_edges) {
                        if ((X.count(e.first) && !X.count(e.second)) || (!X.count(e.first) && X.count(e.second))) {
                            if (!row[j]) {
                                if (active_edges.count(e)) {
                                    cont = true;
                                    break;
                                }
                                delta_X_survive.insert(e);
                            }
                        }
                        j++;
                    }
                    if (cont) continue;

                    // if delta_X_active empty, we need the model to use at least one of the edges in delta_X_survive
                    GRBLinExpr lhs = 0;        
                    for (auto e : delta_X_survive) {
                        auto pos = find(feas_edges.begin(), feas_edges.end(), e);
                        int index = pos - feas_edges.begin(); // certain to be in there as it was checked above
                        lhs += vars[index];
                    }

                    addLazy(lhs, GRB_GREATER_EQUAL, 1);
                    delete[] x;
                    return;

                }
            }
                        


        } 
    } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during callback" << endl;
    }
    }
};

// determine relevant sets of edges which are able to be simultaneously destroyed
vector<vector<bool>> determine_simultaneous_destruction_membership(set<set<Edge>> intersection_edge_sets, set<Edge> orig_edges, set<Edge> original_bridges, vector<Edge> feas_edges) {
    
    vector<vector<bool>> simultaneous_destruction_membership;
    for (auto edge_set : intersection_edge_sets) {
        // disregard the set if it contains no edges of the original
        set<Edge> original_overlap;
        set_intersection(edge_set.begin(), edge_set.end(), orig_edges.begin(), orig_edges.end(), inserter(original_overlap, original_overlap.end()));
        if (original_overlap.size() == 0) continue;
        // also disregard if only original edge is a non bridge
        if ((original_overlap.size() == 1) && (original_bridges.count(*original_overlap.begin()) == 0)) continue;

        // disregard subsets
        bool subset = false;
        for (auto alt_edge_set : intersection_edge_sets) {
            if (edge_set.size() >= alt_edge_set.size()) continue;
            set<Edge> overlap;
            set_intersection(edge_set.begin(), edge_set.end(), alt_edge_set.begin(), alt_edge_set.end(), inserter(overlap, overlap.end()));
            if (overlap == edge_set) {
                subset = true;
                break;
            }
        }
        if (subset) continue;
        
        vector<bool> row = vector<bool>(feas_edges.size());
        for (long unsigned int j=0; j<feas_edges.size(); j++) {
            row[j] = edge_set.count(feas_edges[j]) ? 1 : 0;
        }
        simultaneous_destruction_membership.push_back(row);

    }
    return simultaneous_destruction_membership;

}

// create the model and solve it, recording statistics 
SolverData gurobi_solve(
        bool log,
        int n_nodes, 
        std::set<Edge> orig_edges,
        vector<Edge> feas_edges,
        vector<double> edge_lengths,
        unordered_set<std::pair<Edge,Edge>> edge_overlaps,
        set<set<Edge>> intersection_edge_sets,
        set<Edge> original_bridges
) {
    SolverData final_data;
    final_data.statistic_max_components = 1;

    auto start = chrono::high_resolution_clock::now();
    vector<vector<bool>> simultaneous_destruction_membership = determine_simultaneous_destruction_membership(intersection_edge_sets, orig_edges, original_bridges, feas_edges);
    auto stop = chrono::high_resolution_clock::now();
    final_data.statistic_time_solver_preprocess = chrono::duration_cast<chrono::milliseconds>(stop - start);
    final_data.statistic_sim_destruct_after_size = (int)simultaneous_destruction_membership.size();

    // membership info logging
    if (log) {
        cout << "simultaneous_destruction_membership size: " << simultaneous_destruction_membership.size() << endl;
        cout << "simultaneous destruction contents: \n";
        for (auto row : simultaneous_destruction_membership) {
            for (long unsigned int i=0; i<row.size(); i++) {
                if (row[i]) {
                    cout << feas_edges[i].first << "-" << feas_edges[i].second << " ";
                }
            }
            cout << endl;
        }
        cout << endl;
    }

    // only really important for the callback but oh well. Maybe it helps with the lazy constraints anyway 
    sort(simultaneous_destruction_membership.begin(), simultaneous_destruction_membership.end(), [](const vector<bool> &a, const vector<bool> &b)
    { 
        return count(a.begin(), a.end(), true) > count(b.begin(), b.end(), true); 
    });

    try {
        auto start = chrono::high_resolution_clock::now();
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "solver.log");
        if (!log) env.set(GRB_IntParam_OutputFlag, 0);
        env.start();
        GRBModel model = GRBModel(env);

        // model.set("Presolve", "2");
        model.set("IntegralityFocus", "1");
        model.set("Heuristics", "0");


        int n_edges = feas_edges.size();
        GRBVar* x = model.addVars(n_edges, GRB_BINARY);



        // set objective
        GRBLinExpr obj = 0;
        for (int i=0; i<n_edges; i++) {
            obj += x[i]*edge_lengths[i];
        }
        model.setObjective(obj, GRB_MINIMIZE);


        // constraint 1: original edges stay in
        for (int i=0; i<n_edges; i++) {
            if (orig_edges.count(feas_edges[i])) {
                model.addConstr(x[i], GRB_EQUAL, 1);
            }
        }

        // constraint 2: new edges can't overlap
        for (int i=0; i<n_edges; i++) {
            for (int j=i+1; j<n_edges; j++) {
                Edge e1 = feas_edges[i];
                Edge e2 = feas_edges[j];
                if (!(orig_edges.count(e1) + orig_edges.count(e2))) {
                    if (edge_overlaps.count(pair(e1, e2))) {
                        model.addConstr(x[i] + x[j], GRB_LESS_EQUAL, 1);
                    }
                }
            }
        }

        // constraint 3: connectivity

        vector<vector<bool>> callback_simultaneous_destruction_membership;
        vector<vector<set<int>>> callback_simultaneous_destruction_components;
        
        for (auto row : simultaneous_destruction_membership) {
            set<Edge> row_set, original_overlap, original_less_overlap;
            for (int i=0; i<n_edges; i++) {
                if (row[i]) {
                    row_set.insert(feas_edges[i]);
                }
            }
            set_intersection(orig_edges.begin(), orig_edges.end(), row_set.begin(), row_set.end(), inserter(original_overlap, original_overlap.end()));
            set_difference(orig_edges.begin(), orig_edges.end(), original_overlap.begin(), original_overlap.end(), inserter(original_less_overlap, original_less_overlap.end()));

            set<set<int>> components = determine_components(original_less_overlap, n_nodes);

            if (components.size() == 1) continue;
            if (components.size() == 2) {
                // check if we have a new biggest component to record
                if ((int)components.size() > final_data.statistic_max_components) {
                    final_data.statistic_max_components = components.size();
                }
                set<int> X, X_comp;
                X = *components.begin();
                X_comp = *components.rbegin();

                assert(X.size() + X_comp.size() == (long unsigned int)n_nodes);

                GRBLinExpr lhs = 0;
                for (auto u : X) {
                    for (auto v : X_comp) {
                        Edge e = u < v ? pair(u,v) : pair(v,u);
                        auto pos = find(feas_edges.begin(), feas_edges.end(), e);
                        // our edge may not exist if it's infeasible
                        if (pos != feas_edges.end()) {
                            int index = pos - feas_edges.begin();
                            lhs += (1 - row[index])*x[index];
                        }
                    }
                }
                model.addConstr(lhs, GRB_GREATER_EQUAL, 1);

            } 
            // more than two components
            else {
                // check if we have a new biggest component to record
                if ((int)components.size() > final_data.statistic_max_components) {
                    final_data.statistic_max_components = components.size();
                }

                // take note of the work we have to do in the callback
                callback_simultaneous_destruction_membership.push_back(row);
                vector<set<int>> possible_components;
                vector<set<int>> all_components(components.begin(), components.end());
                // note the possible relevant components (including unions) for the row
                for (long unsigned int k=1; k<=components.size(); k++) { 
                    vector<bool> v(components.size());
                    fill(v.end() - k, v.end(), true);
                    do {
                        set<int> component_union;
                        // why is it seemingly so hard to take unions in c++
                        for (long unsigned int i=0; i<v.size(); i++) {
                            if (!v[i]) continue;
                            for (auto node : all_components[i]) {
                                component_union.insert(node);
                            }
                        }
                        if (component_union.size() > (long unsigned int)(n_nodes/2)) continue;
                        possible_components.push_back(component_union);
                        
                    // move to the next choice of k components
                    } while (std::next_permutation(v.begin(), v.end()));
                }
                callback_simultaneous_destruction_components.push_back(possible_components);

                // still add constraints for each of the "bits"
                for (auto X : components) {
                    set<int> X_comp;
                    for (int i=0; i<n_edges; i++) {
                        if (!X.count(i)) X_comp.insert(i);
                    }

                    GRBLinExpr lhs = 0;
                    for (auto u : X) {
                        for (auto v : X_comp) {
                            Edge e = u < v ? pair(u,v) : pair(v,u);
                            auto pos = find(feas_edges.begin(), feas_edges.end(), e);
                            // our edge may not exist if it's infeasible
                            if (pos != feas_edges.end()) {
                                int index = pos - feas_edges.begin();
                                lhs += (1 - row[index])*x[index];
                            }
                        }
                    }
                    model.addConstr(lhs, GRB_GREATER_EQUAL, 1);
                }
            }
        }

        final_data.statistic_sim_destruct_callback_size = (int)callback_simultaneous_destruction_membership.size();
        if (log) {
            cout << "simultaneous_destruction_membership size for callback: " << callback_simultaneous_destruction_membership.size() << endl;
            cout << "simultaneous destruction contents: \n";
            for (auto row : callback_simultaneous_destruction_membership) {
                for (long unsigned int i=0; i<row.size(); i++) {
                    if (row[i]) {
                        cout << feas_edges[i].first << "-" << feas_edges[i].second << " ";
                    }
                }
                cout << endl;
            }
            cout << endl;
        }

        mycallback cb = mycallback(n_edges, x, n_nodes, feas_edges, callback_simultaneous_destruction_membership, callback_simultaneous_destruction_components);


        model.set("LazyConstraints", "1");
        model.setCallback(&cb);

        model.optimize();

        auto stop = chrono::high_resolution_clock::now();
        final_data.statistic_time_gurobi_total = chrono::duration_cast<chrono::milliseconds>(stop - start);

        // process solution once obtained
        final_data.objective = model.get(GRB_DoubleAttr_ObjVal);
        final_data.statistic_time_gurobi_solve = 1000*model.get(GRB_DoubleAttr_Runtime); // converting s to ms

        set<Edge> final_edges;
        for (int i=0; i<n_edges; i++) {
            if (x[i].get(GRB_DoubleAttr_X) > 0.5) {
                final_edges.insert(feas_edges[i]);
                if (log) cout << feas_edges[i].first << "-" << feas_edges[i].second << " ";
            }
        }
        if (log) cout << endl << endl;


        final_data.final_edges = final_edges;

        return final_data;

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    SolverData s;
    return s; 
}