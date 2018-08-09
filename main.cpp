#include "bwt_se.hpp"

int main(int argc, char **argv) {

  if(argc != 3) {
    fprintf(stderr, "Usage: %s <input sequence> <delta in bits>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  string input_file(argv[1]);
  uint32_t delta = atoi(argv[2]);
  
  bwt_se<uint8_t>(input_file, delta);
  
  return EXIT_SUCCESS;
}
