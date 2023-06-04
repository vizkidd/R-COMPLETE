#include <Rcpp.h>
using namespace Rcpp;

// C++ function to process a single chunk
List processChunk(List chunk, std::string score_col) {
   int n = chunk.size();
  List result(n);

  Environment wis_env = Environment::namespace_env("wisard");
  Environment gr_env = Environment::namespace_env("GenomicRanges");
  Function get_wis = wis_env["get_wis"];
  Function GRanges = gr_env["GRanges"];
  
  for (int i = 0; i < n; i++) {
    List x = chunk[i];

    // Call the R functions with the appropriate arguments
    List tmp_result = get_wis(GRanges(x), Rcpp::Named("max_score") = score_col, Rcpp::Named("overlap") = 0);
    //print(x);
    //print(tmp_result);
    //Rcout << "here" << std::endl;
    //int max_score = tmp_result["max_score"];
    //List alignments = tmp_result["alignments"];

    // Check conditions and store the result in the C++ list
    //if (tmp_result$max_score > 0 && tmp_result$alignments.size() > 0) {
    result[i] = tmp_result;
    //}
}
return result;
}

// Rcpp function that splits the input and processes chunks in parallel
// [[Rcpp::export]]
List run_WISARD_on_BLAST(List blast_hits, int chunk_size,  std::string score_col, int num_threads) {
  int n = blast_hits.size();
  int num_chunks = ceil(n / chunk_size);
  List result(num_chunks);
  // Split the blast_hits into chunks
  List chunks(num_chunks);
  for (int i = 0; i < num_chunks; i++) {
    int start_index = i * chunk_size;
    int end_index = std::min((i + 1) * chunk_size, n);
    chunks[i] = blast_hits[Range(start_index, end_index - 1)];
  }
  
  // Process chunks in parallel
  #pragma omp parallel for num_threads(num_threads)
  for (int i = 0; i < num_chunks; i++) {
    List chunk = chunks[i];
    result[i] = processChunk(chunk, score_col);
  }

  return result;
}