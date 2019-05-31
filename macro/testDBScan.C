#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "ITStracking/Graph.h"
#include "ITStracking/DBScan.h"
#include <random>
#include <vector>

typedef std::pair<int, unsigned char> State;

void dumpClusters(o2::its::Graph<o2::its::Centroid>* graph, size_t range)
{
  for (size_t i{ 0 }; i < range; ++i) {
    auto cluster = graph->getCluster(i);
    std::cout << "cluster starting from " << i << ": ";
    for (auto id : cluster) {
      std::cout << id << " ";
    }
    std::cout << std::endl;
  }
}

void dumpEdges(const std::vector<std::vector<o2::its::Edge>>& edgesVector)
{
  for (size_t i{ 0 }; i < edgesVector.size(); ++i) {
    for (auto& edge : edgesVector[i])
      std::cout << edge.first << " -> " << edge.second << std::endl;
  }
}

void dumpStates(const std::vector<State>& states)
{
  for (auto& state : states) {
    std::cout << "Centroid: " << state.first << std::endl;
    switch (static_cast<int>(state.second)) {
      case 0:
        std::cout << "\tnoise\n";
        break;
      case 1:
        std::cout << "\tborder\n";
        break;
      case 2:
        std::cout << "\tcore\n";
        break;
      default:
        std::cout << "\terror!\n";
    }
  }
}

void checkEdgeConsistency(const std::vector<std::vector<std::pair<int, int>>>& edges_s, const std::vector<std::vector<std::pair<int, int>>>& edges_p, bool debug = false)
{
  std::cout<<"\tChecking Edges consistency...";
  for (size_t i{ 0 }; i < edges_p.size(); ++i) {
    if (edges_s[i].size() != 0 && edges_p[i].size() != 0) {
      if (edges_s[i].size() == edges_p[i].size()) {
        for (size_t j{ 0 }; j < edges_p[i].size(); ++j) {
          if ((edges_s[i][j].first != edges_p[i][j].first) ||
              (edges_s[i][j].second != edges_p[i][j].second)) {
            std::cout << "\tserial\t\tparallel\n"
                      << edges_s[i][j].first << " -> " << edges_s[i][j].second << "\t\t" << edges_p[i][j].first << " -> " << edges_p[i][j].second << std::endl;
          }
        }
      } else {
        std::cout << "Not even the same size!" << std::endl;
      }
    } else {
      if (debug)
        std::cout << "Edge[" << i << "] size: serial is " << edges_s[i].size() << "; parallel is " << edges_p[i].size() << std::endl;
    }
  }
  std::cout<<" done\n";
}

void checkStateConsistency(const std::vector<State>& states_s, const std::vector<State>& states_p)
{
  std::cout<<"\tChecking States consistency...";
  for (size_t i{ 0 }; i < states_p.size(); ++i) {
    if ((states_s[i].first != states_p[i].first) ||
        (states_s[i].second != states_p[i].second)) {
      std::cout << "\tserial\t\tparallel\n"
                << states_s[i].first << " -> " << states_s[i].second << "\t\t" << states_p[i].first << " -> " << states_p[i].second << std::endl;
    }
  }
  std::cout<<" done\n";
}

void checkClusterConsistency(o2::its::Graph<o2::its::Centroid>* graph_serial, o2::its::Graph<o2::its::Centroid>* graph_parallel, const size_t range)
{
  std::cout<<"\tChecking Clusters consistency...";
  for (size_t i{ 0 }; i < range; ++i) {
    auto cluster_s = graph_serial->getCluster(i);
    auto cluster_p = graph_parallel->getCluster(i);
    if (cluster_s.size() != cluster_p.size()) {
      std::cout << "Not even the same lenght for clusters: " << i << std::endl;
      std::cout << "\tSerial cluster starting from " << i << ": ";
      for (auto id : cluster_s) {
        std::cout << id << " ";
      }
      std::cout << std::endl;
      std::cout << "\tParallel cluster starting from " << i << ": ";
      for (auto id : cluster_p) {
        std::cout << id << " ";
      }
    } else {
      std::sort(cluster_s.begin(), cluster_s.end());
      std::sort(cluster_p.begin(), cluster_p.end());
      for (size_t c{ 0 }; c < cluster_p.size(); ++c) {
        if (cluster_p[c] != cluster_s[c]) {
          std::cout << "Mismatch found in cluster: " << i << std::endl;
          std::cout << "\tindex " << c << ": " << cluster_p[c] << " is not " << cluster_s[c] << std::endl;
        }
      }
    }
  }
  std::cout<<" done\n";
}

void testDBScan()
{
  std::cout << " --- Test Graph Class ---" << std::endl;

  std::default_random_engine generator;
  std::uniform_int_distribution<int> id_dist(1, 10000);
  std::uniform_real_distribution<float> coord_dist(-10.f, 10.f);
  std::vector<o2::its::Centroid> centroids;
  for (size_t i{ 0 }; i < 1000; ++i) {
    int indices[2] = { id_dist(generator), id_dist(generator) };
    float coordinates[3] = { coord_dist(generator), coord_dist(generator), coord_dist(generator) };
    centroids.emplace_back(indices, coordinates);
  }

  const float radius{ 3.5f };
  o2::its::Graph<o2::its::Centroid> centroids_graph_serial(1);
  o2::its::Graph<o2::its::Centroid> centroids_graph_parallel(4);

  centroids_graph_serial.init(centroids);
  centroids_graph_parallel.init(centroids);

  centroids_graph_serial.computeEdges([radius](const o2::its::Centroid& c1, const o2::its::Centroid& c2) { return o2::its::Centroid::ComputeDistance(c1, c2) < radius; });
  centroids_graph_parallel.computeEdges([radius](const o2::its::Centroid& c1, const o2::its::Centroid& c2) { return o2::its::Centroid::ComputeDistance(c1, c2) < radius; });

  // dumpEdges(centroids_graph_serial.getEdges());
  // dumpEdges(centroids_graph_parallel.getEdges());

  checkEdgeConsistency(centroids_graph_serial.getEdges(), centroids_graph_parallel.getEdges());

  // dumpClusters(&centroids_graph_serial, centroids.size());
  // dumpClusters(&centroids_graph_parallel, centroids.size());

  checkClusterConsistency(&centroids_graph_serial, &centroids_graph_parallel, centroids.size());

  std::cout << " --- Test DBScan Class ---" << std::endl;
  o2::its::DBScan<o2::its::Centroid> dbscan_serial{ 1 };
  o2::its::DBScan<o2::its::Centroid> dbscan_parallel{ 4 };

  dbscan_serial.init(centroids, [radius](const o2::its::Centroid& c1, const o2::its::Centroid& c2) { return o2::its::Centroid::ComputeDistance(c1, c2) < radius; });
  dbscan_parallel.init(centroids, [radius](const o2::its::Centroid& c1, const o2::its::Centroid& c2) { return o2::its::Centroid::ComputeDistance(c1, c2) < radius; });

  const int nContribs{ 3 };

  // dbscan_serial.classifyVertices([nContribs](std::vector<o2::its::Edge>& edges) { return edges.size() == 0 ? 0 : edges.size() > nContribs ? 2 : 1; });
  // dbscan_parallel.classifyVertices([nContribs](std::vector<o2::its::Edge>& edges) { std::cout<<"> "<<edges.size()<<std::endl; return edges.size() == 0 ? 0 : edges.size() > nContribs ? 2 : 1; });

  dbscan_serial.classifyVertices([nContribs](std::vector<o2::its::Edge>& edges) { return edges.size() == 0 ? 0 : edges.size() > nContribs ? 2 : 1; },
                                 [](State& s1, State& s2) { return static_cast<int>(s1.second) > static_cast<int>(s2.second); });
  dbscan_parallel.classifyVertices([nContribs](std::vector<o2::its::Edge>& edges) { return edges.size() == 0 ? 0 : edges.size() > nContribs ? 2 : 1; },
                                   [](State& s1, State& s2) { return static_cast<int>(s1.second) > static_cast<int>(s2.second); });

  // dumpStates(dbscan_serial.getStates());
  // std::cout<<" --- --- --- \n";
  // dumpStates(dbscan_parallel.getStates());

  checkStateConsistency(dbscan_serial.getStates(), dbscan_parallel.getStates());
}
#endif
