/**
 * @file backtrack.cc
 *
 */

#include "backtrack.h"
#include <fstream>

Backtrack::Backtrack() {}
Backtrack::~Backtrack() {}

namespace {
/* Print debugging messages if DEBUG is true */
bool DEBUG = false;
/* Verify if an embedding is correct if VERIFY is true  */
bool VERIFY = true;
bool ALL_CORRECT = true;

/* Print in a debug mode if PRINT_DEBUG is true */
bool PRINT_DEBUG = true;
/* Print in a submission format if PRINT_SUBMISSION is true */
bool PRINT_SUBMISSION = true;
/* Print in a file if PRINT_FILE is true */
bool PRINT_FILE = true;
std::ofstream out;
/* Print the number of embeddings if PRINT_CNT is true  */
bool PRINT_CNT = true;

/* Time complexity: O(|V(query)|) */
Vertex find_root(const Graph &query, const CandidateSet &cs) {
  // root is a vertex in query which has the minimum |C(u)|/deg(u)
  Vertex root = 0; // Index starts with 0
  double best_score = cs.GetCandidateSize(root) / query.GetDegree(root);

  for (int i=1; i<(int)query.GetNumVertices(); i++) {
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
  for (int i=query.GetNeighborStartOffset(u); i<(int)query.GetNeighborEndOffset(u); i++) {
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
int get_num_extendable_candidates(Vertex u, const Graph &query,
                                  const Graph &data, const CandidateSet &cs,
                                  const std::unordered_map<Vertex, Vertex> &embedding) {
  int count = 0;
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
      count++;
    }
  }
  return count;
}

std::vector<Vertex> get_extendable_candidates(Vertex u, const Graph &query,
                                              const Graph &data, const CandidateSet &cs,
                                              const std::unordered_map<Vertex, Vertex> &embedding) {

    DEBUG && std::cout << "[DEBUG] Get extendable candidates for: " << u << "\n";
    std::vector<Vertex> candidates;
    std::vector<Vertex> parents = get_parents(u, query, embedding);
    // push back in reverse order
    DEBUG && std::cout << "[DEBUG] candidate:";
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
            DEBUG && std::cout << " " << candidate;
            candidates.push_back(candidate);
        }
    }
    DEBUG && std::cout << "\n";
    return candidates;
}

/*
 * This function determines matching order
 * Any vertex in query_next can be next extendable vertext
 * It is assumed that query_next is not empty
 */
Vertex get_extendable_vertex(const std::unordered_set<Vertex> &query_next,
                             const Graph &query, const Graph &data, const CandidateSet &cs,
                             const std::unordered_map<Vertex, Vertex> &embedding) {
  Vertex next = *query_next.begin(); // the first element in the set
  DEBUG && std::cout << "[DEBUG] query_next: " << next;
  int candidate_size = get_num_extendable_candidates(next, query, data, cs, embedding);

  for (auto itr = ++query_next.begin(); itr != query_next.end(); ++itr) {
    Vertex v = *itr;
    DEBUG && std::cout << ", " << v;
    // choose an element with the smallest extendable candidate size
    int temp_candidate_size = get_num_extendable_candidates(v, query, data, cs, embedding);
    if (temp_candidate_size < candidate_size) {
      next = v;
      candidate_size = temp_candidate_size;
    }
  }
  DEBUG && std::cout << "\n";
  return next;
}

bool verify_embedding(const std::vector<std::pair<Vertex, Vertex>> M) {
  std::unordered_set<Vertex> data_visited;

  for (auto p : M) {
    // Check if M is injective (u != v -> M(u) != M(v))
    // It is assumed that there is no duplicate p.first 
    if (data_visited.find(p.second) == data_visited.end()) {
      data_visited.insert(p.second);
    } else {
      return false; // Wrong
    }
  }
  return true;
}

void print_embedding(const std::unordered_map<Vertex, Vertex> &embedding) {
  std::vector<std::pair<Vertex, Vertex>> M; // an embedding is a mapping, thus M
  std::copy(embedding.begin(), embedding.end(), std::back_inserter(M));
  std::sort(M.begin(), M.end()); // the default comparator would do

  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  // PRINT_FILE should be false in a submission version
  if (PRINT_FILE) {
    std::cout.rdbuf(out.rdbuf()); // redirect cout to out
  }

  // Print in a debug mode
  if (PRINT_DEBUG) {
    std::cout << "{";
    std::cout << "(" << M[0].first << ", " << M[0].second << ")";
    for (int i=1; i < (int)M.size(); i++)
      std::cout << ", (" << M[i].first << ", " << M[i].second << ")";
    std::cout << "}" << "\n";
  }  

  // Print in a submission format. It should be true in a submission version
  if (PRINT_SUBMISSION) {
    std::cout << "a";
    for (int i=0; i < (int)M.size(); i++)
      std::cout << " " << M[i].second;
    std::cout << "\n";
  }

  if (PRINT_FILE) {
    std::cout.rdbuf(coutbuf); // reset cout 
  }

  if (VERIFY) {
    bool result = verify_embedding(M);
    if (!result) ALL_CORRECT = false;
    // std::string result_string = verify_embedding(M) ? "correct" : "wrong";
    // std::cout << "Embedding is " << result_string << "\n";
  }
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

  if (PRINT_FILE) {
    out.open("out/out.txt");
  }
  
  Vertex root = find_root(query, cs);
  pair_to_visit.push(std::pair<Vertex, Vertex>(root, -1)); // -1 means NULL

  int cnt = 0; // the number of embedding
  int limit = 100000; // 10^5

  while (!pair_to_visit.empty()) {
    // get current vertex
    std::pair<Vertex, Vertex> current = pair_to_visit.top(); pair_to_visit.pop();

    // current.first is an extendable vertex
    if (current.second == -1) {
      // backtrack
      if (query_visited[current.first]) {
        DEBUG && std::cout << "[DEBUG] BACKTRACK vertex: " << current.first << "\n";
        query_visited[current.first] = false; // mark unvisited
        query_next.insert(current.first); // insert current.first back to query_next

        // remove adjacent vertices to current.first
        for (Vertex neighbor : get_neighbors(current.first, query)) {
          // neighbor is either in partial_embedding or query_next
          if (partial_embedding.find(neighbor) != partial_embedding.end()) continue;
          
          bool is_connected = false;
          for (Vertex v : get_neighbors(neighbor, query)) {
            // query_visited[v] can be true when partial_embedding[v] does not exist
            if (partial_embedding.find(v) != partial_embedding.end()) {
              is_connected = true;
              break;
            }
          }
          if (!is_connected) {
            query_next.erase(neighbor);
          }
        }

        // (current.first, old_data_vertex) is in partial_embedding
        if (partial_embedding.find(current.first) != partial_embedding.end()) {
          Vertex old_data_vertex = partial_embedding[current.first];
          data_visited.erase(old_data_vertex);
          partial_embedding.erase(current.first);
        }
      } 
      // grow
      else {
        DEBUG && std::cout << "[DEBUG] GROW vertex: " << current.first << "\n";
        pair_to_visit.push(current); // re-push current
        query_visited[current.first] = true; // mark visited

        // current.first is no longer in query_next
        query_next.erase(current.first);

        // Add unvisited query vertices to query_next
        // which are adjacent to current.first and not in query_next
        DEBUG && std::cout << "[DEBUG] new query_next for " << current.first << ": ";
        for (Vertex v : get_neighbors(current.first, query)) {
          if (!query_visited[v] && query_next.find(v) == query_next.end()) {
            DEBUG && std::cout << v << ", ";
            query_next.insert(v);
          }
        }
        DEBUG && std::cout << "\n";

        // get extendable candidates of current.first
        std::vector<Vertex> candidates = get_extendable_candidates(current.first, query, data, cs, partial_embedding);
        DEBUG && std::cout << "[DEBUG] new extendable candidates for " << current.first << ": ";
        for (Vertex v : candidates) {
          if (data_visited.find(v) == data_visited.end()) {
            DEBUG && std::cout << "(" << current.first << ", " << v << "), ";
            pair_to_visit.push(std::pair<Vertex, Vertex>(current.first, v));
          }
        }
        DEBUG && std::cout << "\n";
      }
    }
    // (current.first, current.second) is an extendable candidate
    else {
      // Here is where a new pair is added to the embedding
      if (partial_embedding.find(current.first) != partial_embedding.end()) {
        Vertex old_data_vertex = partial_embedding[current.first];
        data_visited.erase(old_data_vertex);
      }
      partial_embedding[current.first] = current.second; // update partial_embedding
      // std::cout << "[DEBUG] partial embedding size: " << partial_embedding.size() << "\n";
      if (DEBUG) print_embedding(partial_embedding); // for debugging

      if (partial_embedding.size() == query.GetNumVertices()) {
        print_embedding(partial_embedding);
        if (++cnt == limit) break;
        continue;
      }
      data_visited.insert(current.second);

      // Find next extendable vertex
      Vertex next = get_extendable_vertex(query_next, query, data, cs, partial_embedding);
      // std::cout << "[DEBUG] extendable vertex: " << next << "\n";
      pair_to_visit.push(std::pair<Vertex, Vertex>(next, -1));
    }
  }

  if (PRINT_FILE) {
    out.close();
  }

  if (VERIFY && ALL_CORRECT) {
    std::cout << "[DEBUG] All embeddings are correct" << "\n";
  }
  if (PRINT_CNT) {
    std::cout << "[DEBUG] # embeddings: " << cnt << "\n";
  }
}
