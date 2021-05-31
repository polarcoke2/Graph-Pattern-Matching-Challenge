/**
 * @file backtrack.cc
 *
 */

#include "backtrack.h"

Backtrack::Backtrack() {}
Backtrack::~Backtrack() {}

namespace {
/* Time complexity: O(|V(query)|) */
Vertex find_root(const Graph &query, const CandidateSet &cs) {
  // root is a vertex in query which has the minimum |C(u)|/deg(u)
  Vertex root = 0; // Index starts with 0
  double best_score = cs.GetCandidateSize(root) / query.GetDegree(root);

  for (Vertex i=1; i<query.GetNumVertices(); i++) {
    double new_score = cs.GetCandidateSize(i) / query.GetDegree(i);
    if (new_score < best_score) {
      root = i;
      best_score = new_score;
    }
  }
  return root;
}

/*
 * In all cases, u is a query graph vertex
 */
std::vector<Vertex> get_neighbors(Vertex u, const Graph &query) {
  std::vector<Vertex> neighbors;
  for (int i=query.GetNeighborStartOffset(u); i<query.GetNeighborEndOffset(u); i++) {
    Vertex neighbor = query.GetNeighbor(i);
    neighbors.push_back(neighbor);
  }
  return neighbors;
}

/*
 * Parents of a vertex is its neighbors which are in partial embedding
 * In all cases, u is a query graph vertex
 */
std::vector<Vertex> get_parents(Vertex u, const Graph &query,
                                const std::unordered_map<Vertex, Vertex> &embedding) {
  std::vector<Vertex> parents;
  for (Vertex neighbor : get_neighbors(u, query)) {
    auto itr = embedding.find(neighbor);
    if (itr != embedding.end()) {
      parents.push_back(itr->first);
    }
  }
  return parents;
}

/* 
 * Get the initial extendable candidates of a data graph vertex
 * Some vertices in the candidate set may not be valid
 */
std::vector<Vertex> get_extendable_candidates(Vertex u, const Graph &query, 
                                          const Graph &data, const CandidateSet &cs,
                                          const std::unordered_map<Vertex, Vertex> &embedding) {
  std::vector<Vertex> candidates;
  std::vector<Vertex> parents = get_parents(u, query, embedding);
  // push back in reverse order
  for (int i=cs.GetCandidateSize(u)-1; i>=0; i--) {
    Vertex candidate = cs.GetCandidate(u, i);
    // for every parent of vertex u, if M[u_parent] and candidate are not connected, fail
    bool is_connected = true;
    for (Vertex parent : parents) {
      // parent is guaranteed to be included in embedding
      if (!data.IsNeighbor(embedding.find(parent)->second, candidate)) {
        is_connected = false;
        break;
      }
    }

    if (is_connected) {
      candidates.push_back(candidate);
    }
  }
  return candidates;
}

/*
 * This function determines matching order
 * Any vertex in query_next can be next extendable vertext
 * It is assumed that query_next is not empty
 */
Vertex get_extendable_vertex(const std::unordered_set<Vertex> &query_next, const CandidateSet &cs) {
  Vertex next = *query_next.begin(); // the first element in the set
  int candidate_size = cs.GetCandidateSize(next);
  for (auto itr = ++query_next.begin(); itr != query_next.end(); ++itr) {
    Vertex v = *itr;
    // choose an element with the smallest candidate size
    // TODO: the standard should be the size of "extendable candidates", not candidate size in CS
    if (cs.GetCandidateSize(v) < candidate_size) {
      next = v;
      candidate_size = cs.GetCandidateSize(v);
    }
  }

  return next;
}

void print_embedding(const std::unordered_map<Vertex, Vertex> &embedding) {
  std::vector<std::pair<Vertex, Vertex>> M; // an embedding is a mapping, thus M
  std::copy(embedding.begin(), embedding.end(), std::back_inserter(M));
  std::sort(M.begin(), M.end()); // the default comparator would do

  std::cout << "{";
  std::cout << "(" << M[0].first << ", " << M[0].second << ")";
  for (int i=1; i < M.size(); i++)
    std::cout << ", (" << M[i].first << ", " << M[i].second << ")";
  std::cout << "}" << std::endl;
}
}  // namespace

void Backtrack::PrintAllMatches(const Graph &data, const Graph &query,
                                const CandidateSet &cs) {
  std::cout << "t " << query.GetNumVertices() << "\n";

  // implement your code here.

  /*
   * query_next is candidates for extendable vertex
   * an unvisited vertex which is adjancent to some vertex in partial embedding
   * should be included in query_next
   */
  std::unordered_set<Vertex> query_next; 
  // Question: do we need both query_visited and data_visited when we know whether  
  // a vertex is visited by checking if the vertex is in partial_embedding?
  std::vector<bool> query_visited(query.GetNumVertices(), false);
  std::unordered_set<Vertex> data_visited; // choose set over vector to reduce space
  /* 
   * Let p be a pair_to_visit object.
   * If p.second is -1 (NULL) then p.first is an extendable vertex
   * If p.second is not -1 then (p.first, p.second) is an extendable candidate
   */ 
  std::stack<std::pair<Vertex, Vertex>> pair_to_visit;
  std::unordered_map<Vertex, Vertex> partial_embedding;
  
  Vertex root = find_root(query, cs);
  pair_to_visit.push(std::pair<Vertex, Vertex>(root, -1)); // -1 means NULL

  while (!pair_to_visit.empty()) {
    // get current vertex
    std::pair<Vertex, Vertex> current = pair_to_visit.top(); pair_to_visit.pop();

    // current.first is an extendable vertex
    if (current.second == -1) {
      // backtrack
      if (query_visited[current.first]) {
        query_visited[current.first] = false; // mark unvisited
        query_next.erase(current.first);
        // remove adjacent vertices to current.first
        for (Vertex neighbor : get_neighbors(current.first, query)) {
          bool is_connected = false;
          for (Vertex v : get_neighbors(neighbor, query)) {
            // If query_visited[v] is true, then partial_embedding[v] must exist
            if (query_visited[v]) {
              is_connected = true;
              break;
            }
          }
          if (!is_connected) {
            query_next.erase(neighbor);
          }
        }

        /* duplicate code */
        // (current.first, old_data_vertex) is in partial_embedding
        if (partial_embedding.find(current.first) != partial_embedding.end()) {
          Vertex old_data_vertex = partial_embedding[current.first];
          data_visited.erase(old_data_vertex);
        }
        /* end duplicate code */
      } 
      // grow
      else {
        pair_to_visit.push(current); // re-push current
        // get extendable candidates of current.first
        std::vector<Vertex> candidates = get_extendable_candidates(current.first, data, query, cs, partial_embedding);
        for (Vertex v : candidates) {
          if (data_visited.find(v) == data_visited.end()) {
            pair_to_visit.push(std::pair<Vertex, Vertex>(current.first, v));
          }
        }
      }
    }
    // (current.first, current.second) is an extendable candidate
    else {
      /* duplicate code */
      // (current.first, old_data_vertex) is in partial_embedding
      if (partial_embedding.find(current.first) != partial_embedding.end()) {
        Vertex old_data_vertex = partial_embedding[current.first];
        data_visited.erase(old_data_vertex);
      }
      /* end duplicate code */
      // Here is where a new pair is added to the embedding
      partial_embedding[current.first] = current.second; // update partial_embedding
      if (partial_embedding.size() == query.GetNumVertices()) {
        print_embedding(partial_embedding);
        continue;
      }
      data_visited.insert(current.second);

      // Find next extendable candidate
      query_next.insert(current.first);
      Vertex next = get_extendable_vertex(query_next, cs);
      query_visited[next] = true;
      query_next.erase(current.first);
      pair_to_visit.push(std::pair<Vertex, Vertex>(next, -1));
    }
  }
}
