#include "../sub/ConsoleInterface/ConsoleInterface.h"
#include "../sub/matrix/blas.h"
#include "winograd.h"
#include "winograd_parallel.h"


using namespace s21;


class Interface : public ConsoleInterface {
 public:
  Interface();
  Graph &GetGraph() { return graph_; }

 private:
  Graph graph_;
};

Interface::Interface() {
    auto home_menu = AddMenu(
      {"Input Matrices",
       "Random Matrices",
       "Print Matrices",
       "graph traversal in depth",
       "searching for the shortest path between any two vertices",
       "searching for the shortest paths between all pairs of vertices in the "
       "graph",
       "searching for the minimal spanning tree in the graph",
       "solving the salesman problem",
       "Comparison of methods for solving the traveling salesman problem"});

}