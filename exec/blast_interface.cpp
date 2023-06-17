#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sys/mman.h>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <omp.h>
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <string_view>
#include <ncbi-tools++/algo/blast/api/local_blast.hpp>
#include <ncbi-tools++/algo/blast/api/blast_results.hpp>
#include <ncbi-tools++/objtools/blast/seqdb_reader/seqdb.hpp>

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// Define the FastaSequenceData class/struct
struct FastaSequenceData {
  std::string header;
  std::string seq;
};

// Function to process a single FASTA block
void processFastaBlock(const FastaSequenceData& data, std::ostringstream& outputStream) {
  // Perform processing on the FastaSequenceData object
  outputStream << "Header: " << data.header << std::endl;
  outputStream << "Sequence: " << data.seq << std::endl;
  outputStream << std::endl;
  outputStream.flush();
}

// Function to stream FASTA file using mmap()
void streamFastaFile(const std::string& filename, const int num_threads) {
  
  omp_set_num_threads(num_threads);
  std::ostringstream* outputStreams = new std::ostringstream[num_threads];
  // Open the FASTA file using fopen()
  FILE* file = fopen(filename.c_str(), "rb");
  if (!file) {
    std::cerr << "Error: Failed to open file: " << filename << std::endl;
    return;
  }
  
  // Get the file size
  fseek(file, 0, SEEK_END);
  long fileSize = ftell(file);
  rewind(file);
  
  // Memory map the file
  char* fileData = static_cast<char*>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fileno(file), 0));
  if (fileData == MAP_FAILED) {
    std::cerr << "Error: Failed to map file: " << filename << std::endl;
    fclose(file);
    return;
  }
  
  // Process the FASTA file
  //char* p = fileData;
  //bool loopControl=true;
  int rec_count=0;
  #pragma omp parallel for num_threads(num_threads)
  for (char* p = fileData; p < fileData + fileSize; p = p + 1) {
    if(p != nullptr){
      // Read the header line
      char* headerEnd=strchr(p, '\n');
      //std::cout << "here1" << std::endl;
      std::string header(p + 1, headerEnd - p - 1);
      p += header.size() + 1;
      //std::cout << "here2" << std::endl;
      // Read the sequence
      char* seqEnd=strchr(p, '>');
      if (!seqEnd)
      {
        //break;
        //Last record, so read till \n
        seqEnd=fileData + fileSize; //strchr(p, '\n');
        //std::cout << seqEnd << std::endl;
      }
      std::string seq(p, seqEnd - p);
      //std::cout << "here3" << std::endl << (long)(fileData + fileSize) << std::endl;
      FastaSequenceData data;
      data.header= header;
      data.seq = seq;
      if(!data.header.empty() && !data.seq.empty()){
        rec_count++;
        //#pragma omp parallel
        //{
        int threadID = omp_get_thread_num();
        processFastaBlock(data, outputStreams[threadID]);
        //}
      }
      p = seqEnd - 1; //strchr(p, '>') - 1;
      //std::cout << "here4" << std::endl;
    }
  }
  
#pragma omp barrier
  
  for (int i = 0; i < num_threads; i++) {
    std::cout << outputStreams[i].str();
  }
  
  // Cleanup
  delete[] outputStreams;
  
  
  // Unmap the file
  munmap(fileData, fileSize);
  //std::cout << rec_count << std::endl;
  // Close the file using fclose()
  fclose(file);
}

//int main() {
//  std::string fastaFile = "ungrouped.cds";

// Stream FASTA file using mmap()
//  streamFastaFile(fastaFile);

//  return 0;
//}

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
