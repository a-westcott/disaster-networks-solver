#ifndef GLOBAL_HPP
#define GLOBAL_HPP

#include <utility>
#include <set>

typedef std::pair<int,int> Edge;

namespace std {
    template <> struct hash<std::pair<Edge,Edge>>
    {
        size_t operator()(const std::pair<Edge,Edge> & x) const
        {
            return (
                (((51 + std::hash<int>()(x.first.first)) * 51 
                     + std::hash<int>()(x.first.second)) * 51 
                     + std::hash<int>()(x.second.first)) * 51 
                     + std::hash<int>()(x.second.second)
            );    
        }
    };
}

#endif