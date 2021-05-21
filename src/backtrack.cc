/**
 * @file backtrack.cc
 *
 */

#include "backtrack.h"

Backtrack::Backtrack() {}
Backtrack::~Backtrack() {}

namespace {
// Time complexity: O(|V(query)|)
Vertex find_root(const Graph &query, const CandidateSet &cs) {
  // root is a vetex in query which has the minimum |C(u)|/deg(u)
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

std::vector<Vertex> extendable_candidates(Vertex u, const CandidateSet &cs) {
  std::vector<Vertex> candidates;
  // push back in reverse order
  for (int i=cs.GetCandidateSize(u)-1; i>=0; i--) {
    candidates.push_back(cs.GetCandidate(u, i));
  }
  return candidates;
}

// This function determines matching order
// It is assumed that query_univisited is not empty
Vertex extendable_vertex(const std::unordered_set<Vertex> &query_unvisited, const CandidateSet &cs) {
  Vertex next = *query_unvisited.begin(); // the first element in the set
  int candidate_size = cs.GetCandidateSize(next);
  for (auto itr = ++query_unvisited.begin(); itr != query_unvisited.end(); itr++) {
    Vertex v = *itr;
    // choose an element with the smallest candidate size
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

  std::unordered_set<Vertex> query_unvisited;
  for (Vertex i=0; i<query.GetNumVertices(); i++)
    query_unvisited.insert(i);
  // Question: do we need query_visited and data_visited when we know whether  
  // a vertex is visited by checking if the vertex is in parital_embedding?
  std::vector<bool> query_visited(query.GetNumVertices(), false);
  std::vector<bool> data_visited(data.GetNumVertices(), false);
  /* 
   * Let p be a pair_to_visit object.
   * If p.second is -1 (NULL) then p.first is an extendable vertex
   * If p.second is not -1 then (p.first, p.second) is an extendable candidate
   */ 
  std::stack<std::pair<Vertex, Vertex>> pair_to_visit;
  std::unordered_map<Vertex, Vertex> parital_embedding;
  
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
        query_unvisited.insert(current.first);
        /* duplicate code */
        // (current.first, old_data_vertex) is in partial_embedding
        if (parital_embedding.find(current.first) != parital_embedding.end()) {
          Vertex old_data_vertex = parital_embedding[current.first];
          data_visited[old_data_vertex] = false;
        }
        /* end duplicate code */
      } 
      // grow
      else {
        pair_to_visit.push(current); // re-push current
        query_visited[current.first] = true; // mark visited
        query_unvisited.erase(current.first);
        // get extendable candidates of current.first
        std::vector<Vertex> candidates = extendable_candidates(current.first, cs);
        for (Vertex v : candidates) {
          if (!data_visited[v]) {
            pair_to_visit.push(std::pair<Vertex, Vertex>(current.first, v));
          }
        }
      }
    }
    // (current.first, current.second) is an extendable candidate
    else {
      /* duplicate code */
      // (current.first, old_data_vertex) is in partial_embedding
      if (parital_embedding.find(current.first) != parital_embedding.end()) {
        Vertex old_data_vertex = parital_embedding[current.first];
        data_visited[old_data_vertex] = false;
      }
      /* end duplicate code */
      // Here is where a new pair is added to the embedding
      parital_embedding[current.first] = current.second; // update parital_embedding
      if (parital_embedding.size() == query.GetNumVertices()) {
        // Or if (query_unvisited.size() == 0)
        print_embedding(parital_embedding);
        continue;
      }
      data_visited[current.second] = true;

      // Find next extendable candidate
      Vertex next = extendable_vertex(query_unvisited, cs);
      query_visited[next] = true;
      query_unvisited.erase(current.first);
      pair_to_visit.push(std::pair<Vertex, Vertex>(next, -1));
    }
  }
}
