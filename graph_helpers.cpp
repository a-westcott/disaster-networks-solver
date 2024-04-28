#include <vector>
#include <deque>
#include <algorithm>
#include <iostream>
#include <functional>

#include "graph_helpers.hpp"

using namespace std;

// determine the set of components of a graph
set<set<int>> determine_components(set<Edge> edges, int n_nodes) {
    set<set<int>> components;
    vector<bool> visited(n_nodes, false);

    while (find(visited.begin(), visited.end(), false) != visited.end()) {
        deque<int> queue;
        int start_node = find(visited.begin(), visited.end(), false) - visited.begin();
        queue.push_front(start_node);
        set<int> component;
        component.insert(start_node);
        while (queue.size()) {
            int node = queue.front();
            queue.pop_front();
            if (visited[node]) continue;
            visited[node] = true;

            for (Edge e : edges) {
                if (e.first == node && !visited[e.second]) {
                    queue.push_back(e.second);
                    component.insert(e.second);
                }
                if (e.second == node && !visited[e.first]) {
                    queue.push_back(e.first);
                    component.insert(e.first);
                }
            }
        }
        components.insert(component);   
    }

    return components;

}

// determine which edges of a graph are bridges
// adapted from https://cp-algorithms.com/graph/bridge-searching.html
set<Edge> determine_bridges(int n_nodes, set<Edge> edges) {
    vector<bool> visited;
    vector<int> tin, low;
    int timer;
    set<Edge> bridges;

    function<void(int,int)> dfs = [&] (int v, int p = -1) {
        visited[v] = true;
        tin[v] = low[v] = timer++;

        vector<int> adj;
        for (Edge e : edges) {
            if (e.first == v) {
                adj.push_back(e.second);
            }
            if (e.second == v) {
                adj.push_back(e.first);
            }
        }

        for (int to : adj) {
            if (to == p) continue;
            if (visited[to]) {
                low[v] = min(low[v], tin[to]);
            } else {
                dfs(to, v);
                low[v] = min(low[v], low[to]);
                if (low[to] > tin[v]) {
                    Edge bridge = (v < to) ? pair(v, to) : pair(to, v);
                    bridges.insert(bridge);
                }
            }
        }
    };

    timer = 0;
    visited.assign(n_nodes, false);
    tin.assign(n_nodes, -1);
    low.assign(n_nodes, -1);
    for (int i=0; i<n_nodes; ++i) {
        if (!visited[i])
            dfs(i, -1);
    }
    return bridges;
}

// determine the number of components in a graph
int determine_n_components(set<Edge> edges, int n_nodes) {
    return (int)determine_components(edges, n_nodes).size();
}