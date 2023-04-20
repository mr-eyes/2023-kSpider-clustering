#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>

using namespace boost;
using namespace std;

typedef adjacency_list<vecS, vecS, directedS> Graph;

void process_edge_line(const std::string &line, Graph &graph) {
    std::istringstream iss(line);
    int node1, node2;
    double weight;
    char separator;

    iss >> node1 >> separator >> node2 >> separator >> weight;

    while (node1 >= num_vertices(graph) || node2 >= num_vertices(graph)) {
        add_vertex(graph);
    }

    add_edge(node1, node2, graph);
}

int main(int argc, char *argv[]) {
    std::string edges_file_path = argv[1];
    std::ifstream input_file(edges_file_path);
    std::string line;

    Graph graph;
    // skip first line
    std::getline(input_file, line);
    cout << "parsing edges ..." << endl;
    while (std::getline(input_file, line)) {
        process_edge_line(line, graph);
    }
    input_file.close();
    cout << "parsing edges done" << endl;

    // Compute strongly connected components
    cout << "computing strongly connected components ..." << endl;
    std::vector<int> component(num_vertices(graph));
    int num_components = strong_components(graph, make_iterator_property_map(component.begin(), get(vertex_index, graph)));

    // Create a hash set of strongly connected nodes
    cout << "creating hash set of strongly connected nodes ..." << endl;
    std::unordered_set<int> strong_nodes;
    for (std::size_t i = 0; i < component.size(); ++i) {
        if (strong_nodes.find(component[i]) != strong_nodes.end()) {
            strong_nodes.insert(static_cast<int>(i));
        } else {
            strong_nodes.erase(static_cast<int>(i));
        }
    }

    // Write edges belonging to strongly connected components
    cout << "writing edges belonging to strongly connected components ..." << endl;
    std::ofstream output_file("strong_edges.txt");
    input_file.open(edges_file_path);
    while (std::getline(input_file, line)) {
        std::istringstream iss(line);
        int node1, node2;
        double weight;
        char separator;

        iss >> node1 >> separator >> node2 >> separator >> weight;

        if (strong_nodes.find(node1) != strong_nodes.end() && strong_nodes.find(node2) != strong_nodes.end()) {
            output_file << node1 << "," << node2 << "," << weight << std::endl;
        }
    }
    input_file.close();
    output_file.close();

    return 0;
}
