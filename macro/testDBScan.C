#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "ITStracking/Graph.h"
#include "ITStracking/DBScan.h"
#include <random>
#include <vector>

typedef std::pair<int, unsigned char> State;

void dumpDBSCANClusters(o2::its::DBScan<o2::its::Centroid>& dbscan)
{
  int counter{ 0 };
  auto clusters = dbscan.computeClusters();

  for (auto& cluster : clusters) {
    std::cout << "cluster " << counter << " size is " << cluster.size() << " dumping: " << std::endl;
    std::copy(cluster.begin(), cluster.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
    ++counter;
  }
}

void dumpRealClusters(o2::its::Graph<o2::its::Centroid>* graph, size_t range)
{
  for (size_t i{ 0 }; i < range; ++i) {
    auto cluster = graph->getClusterIndices(i);
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
      std::cout << i << " -> " << edge << std::endl;
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

void checkEdgeConsistency(const std::vector<std::vector<int>>& edges_s, const std::vector<std::vector<int>>& edges_p, bool debug = false)
{
  std::cout << "\tChecking Edges consistency...";
  for (size_t i{ 0 }; i < edges_p.size(); ++i) {
    if (edges_s[i].size() != 0 && edges_p[i].size() != 0) {
      if (edges_s[i].size() == edges_p[i].size()) {
        for (size_t j{ 0 }; j < edges_p[i].size(); ++j) {
          if (edges_s[i][j] != edges_p[i][j]) {
            std::cout << "\tserial\t\tparallel\n"
                      << i << " -> " << edges_s[i][j] << "\t\t" << i << " -> " << edges_p[i][j] << std::endl;
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
  std::cout << " done\n";
}

void checkStateConsistency(const std::vector<State>& states_s, const std::vector<State>& states_p)
{
  std::cout << "\tChecking States consistency...";
  for (size_t i{ 0 }; i < states_p.size(); ++i) {
    if ((states_s[i].first != states_p[i].first) ||
        (states_s[i].second != states_p[i].second)) {
      std::cout << "\tserial\t\tparallel\n"
                << states_s[i].first << " -> " << states_s[i].second << "\t\t" << states_p[i].first << " -> " << states_p[i].second << std::endl;
    }
  }
  std::cout << " done\n";
}

void checkRawClusterConsistency(o2::its::Graph<o2::its::Centroid>* graph_serial, o2::its::Graph<o2::its::Centroid>* graph_parallel, const size_t range)
{
  std::cout << "\tChecking raw clusters consistency...";
  for (size_t i{ 0 }; i < range; ++i) {
    auto cluster_s = graph_serial->getClusterIndices(i);
    auto cluster_p = graph_parallel->getClusterIndices(i);
    if (cluster_s.size() != cluster_p.size()) {
      std::cout << "Clusters size is different for cluster starting on " << i << std::endl;
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
  std::cout << " done\n";
}

void checkClusterConsistency(o2::its::DBScan<o2::its::Centroid>* graph_serial, o2::its::DBScan<o2::its::Centroid>* graph_parallel)
{
  std::cout << "\tChecking DBScan clusters consistency...";
  auto cluster_s = graph_serial->computeClusters();
  auto cluster_p = graph_parallel->computeClusters();
  if (cluster_s.size() != cluster_p.size()) {
    std::cout << "\n\t /!\\Clusters size is different!" << std::endl;
    return;
  }
  std::sort(cluster_s.begin(), cluster_s.end(), [](const std::vector<int> c_1, const std::vector<int> c_2) { return c_1.size() < c_2.size(); });
  std::sort(cluster_p.begin(), cluster_p.end(), [](const std::vector<int> c_1, const std::vector<int> c_2) { return c_1.size() < c_2.size(); });
  for (size_t i{ 0 }; i < cluster_s.size(); ++i) {
    std::sort(cluster_s[i].begin(), cluster_s[i].end(), std::greater<int>());
    std::sort(cluster_p[i].begin(), cluster_p[i].end(), std::greater<int>());
    for (size_t j{ 0 }; j < cluster_s[i].size(); ++j) {
      if (cluster_s[i][j] == !cluster_p[i][j]) {
        std::cout << "Mismatch found in cluster: " << i << std::endl;
        std::cout << "\tindex " << j << ": " << cluster_p[i][j] << " is not " << cluster_s[i][j] << std::endl;
      }
    }
  }
  std::cout << " done\n";
}

void testDBScan(const size_t nVert = 10, const float radius = 2.f, const int nContribs = 2)
{

  std::cout << " --- Test Graph Class ---" << std::endl;

  std::default_random_engine generator;
  std::uniform_int_distribution<int> id_dist(1, nVert);
  std::uniform_real_distribution<float> coord_dist(-10.f, 10.f);
  std::vector<o2::its::Centroid> centroids;
  for (size_t i{ 0 }; i < nVert; ++i) {
    int indices[2] = { id_dist(generator), id_dist(generator) };
    float coordinates[3] = { coord_dist(generator), coord_dist(generator), coord_dist(generator) };
    centroids.emplace_back(indices, coordinates);
  }

  o2::its::Graph<o2::its::Centroid> centroids_graph_serial(1);
  o2::its::Graph<o2::its::Centroid> centroids_graph_parallel(4);

  centroids_graph_serial.init(centroids);
  centroids_graph_parallel.init(centroids);

  centroids_graph_serial.computeEdges([radius](const o2::its::Centroid& c1, const o2::its::Centroid& c2) { return o2::its::Centroid::ComputeDistance(c1, c2) < radius; });
  centroids_graph_parallel.computeEdges([radius](const o2::its::Centroid& c1, const o2::its::Centroid& c2) { return o2::its::Centroid::ComputeDistance(c1, c2) < radius; });

  // dumpEdges(centroids_graph_serial.getEdges());
  // dumpEdges(centroids_graph_parallel.getEdges());

  checkEdgeConsistency(centroids_graph_serial.getEdges(), centroids_graph_parallel.getEdges());

  // dumpRealClusters(&centroids_graph_serial, centroids.size());
  // dumpRealClusters(&centroids_graph_parallel, centroids.size());

  checkRawClusterConsistency(&centroids_graph_serial, &centroids_graph_parallel, centroids.size());

  std::cout << " --- Test DBScan Class ---" << std::endl;
  o2::its::DBScan<o2::its::Centroid> dbscan_serial{ 1 };
  o2::its::DBScan<o2::its::Centroid> dbscan_parallel{ 4 };

  dbscan_serial.init(centroids, [radius](const o2::its::Centroid& c1, const o2::its::Centroid& c2) { return o2::its::Centroid::ComputeDistance(c1, c2) < radius; });
  dbscan_parallel.init(centroids, [radius](const o2::its::Centroid& c1, const o2::its::Centroid& c2) { return o2::its::Centroid::ComputeDistance(c1, c2) < radius; });

  checkStateConsistency(dbscan_serial.getStates(), dbscan_parallel.getStates());

  dbscan_serial.classifyVertices(nContribs);
  dbscan_parallel.classifyVertices(nContribs);

  // dumpDBSCANClusters(dbscan_serial);
  // dumpDBSCANClusters(dbscan_parallel);
  checkClusterConsistency(&dbscan_serial, &dbscan_parallel);
}
#endif
