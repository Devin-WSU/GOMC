#include "CheckpointUtilities.h"

void CheckpointUtilities::readVector3DDouble(FILE * stream, std::vector< std::vector< std::vector <double> > > & data){
      // read size of data
  ulong size_x = read_uint64_binary(stream);
  ulong size_y = read_uint64_binary(stream);
  ulong size_z = read_uint64_binary(stream);

  // read array
  data.resize(size_x);
  for(int i = 0; i < (int) size_x; i++) {
    data[i].resize(size_y);
    for(int j = 0; j < (int) size_y; j++) {
      data[i][j].resize(size_z);
      for(int k = 0; k < (int) size_z; k++) {
        data[i][j][k] = read_double_binary(stream);
      }
    }
  }
}
void CheckpointUtilities::readVector3DUint(FILE * stream, std::vector< std::vector< std::vector <uint> > > & data){
      // read size of data
  ulong size_x = read_uint64_binary(stream);
  ulong size_y = read_uint64_binary(stream);
  ulong size_z = read_uint64_binary(stream);

  // read array
  data.resize(size_x);
  for(int i = 0; i < (int) size_x; i++) {
    data[i].resize(size_y);
    for(int j = 0; j < (int) size_y; j++) {
      data[i][j].resize(size_z);
      for(int k = 0; k < (int) size_z; k++) {
        data[i][j][k] = read_uint32_binary(stream);
      }
    }
  }
}
void CheckpointUtilities::readVector2DUint(FILE * stream, std::vector< std::vector< uint > > & data){  // read size of data
  ulong size_x = read_uint64_binary(stream);
  ulong size_y = read_uint64_binary(stream);

  // read array
  data.resize(size_x);
  for(int i = 0; i < (int) size_x; i++) {
    data[i].resize(size_y);
    for(int j = 0; j < (int) size_y; j++) {
      data[i][j] = read_uint32_binary(stream);
    }
  }}
void CheckpointUtilities::readVector1DDouble(FILE * stream, std::vector< double > & data){
    // read size of data
  ulong size_x = read_uint64_binary(stream);

  // read array
  data.resize(size_x);
  for(int i = 0; i < (int) size_x; i++) {
    data[i] = read_double_binary(stream);
  }
}
double CheckpointUtilities::read_double_binary(FILE * stream){
  if(stream == NULL) {
    printf("Error opening checkpoint output file\n");
    exit(EXIT_FAILURE);
  }
  dbl_input_union temp;
  int ret = fscanf(stream, "%c%c%c%c%c%c%c%c",
                   &temp.bin_value[0],
                   &temp.bin_value[1],
                   &temp.bin_value[2],
                   &temp.bin_value[3],
                   &temp.bin_value[4],
                   &temp.bin_value[5],
                   &temp.bin_value[6],
                   &temp.bin_value[7]);
  if(ret != 8) {
    std::cerr << "CheckpointSetup couldn't read required data from binary!\n";
    exit(EXIT_FAILURE);
  }
  return temp.dbl_value;
}
int8_t CheckpointUtilities::read_uint8_binary(FILE * stream){
  if(stream == NULL) {
    printf("Error opening checkpoint output file\n");
    exit(EXIT_FAILURE);
  }
  int8_input_union temp;
  int ret = fscanf(stream, "%c",
                   &temp.bin_value[0]);
  if(ret != 1) {
    // We could add this back if we REQUIRE the PT flag as output, but
    // this would break all previous checkpoint files generated by legacy code.
    //std::cerr << "CheckpointSetup couldn't read required data from binary!\n";
    //exit(EXIT_FAILURE);
    // If we ran out of binary, then return 0, which evaluates to false
    // when casted to boolean.  This way legacy checkpoints can be used
    // as input to MPI builds.
    return 0;
  }
  return temp.int_value;
}
uint32_t CheckpointUtilities::read_uint32_binary(FILE * stream){
  if(stream == NULL) {
    printf("Error opening checkpoint output file\n");
    exit(EXIT_FAILURE);
  }
  uint32_input_union temp;
  int ret = fscanf(stream, "%c%c%c%c",
                   &temp.bin_value[0],
                   &temp.bin_value[1],
                   &temp.bin_value[2],
                   &temp.bin_value[3]);
  if(ret != 4) {
    std::cerr << "CheckpointSetup couldn't read required data from binary!\n";
    exit(EXIT_FAILURE);
  }
  // Fix endianness, implementation in lib/Endian.h
  temp.uint_value = ftoh32(temp.uint_value);
  return temp.uint_value;
}
uint64_t CheckpointUtilities::read_uint64_binary(FILE * stream){
  if(stream == NULL) {
    printf("Error opening checkpoint output file\n");
    exit(EXIT_FAILURE);
  }
  uint64_input_union temp;
  int ret = fscanf(stream, "%c%c%c%c%c%c%c%c",
                   &temp.bin_value[0],
                   &temp.bin_value[1],
                   &temp.bin_value[2],
                   &temp.bin_value[3],
                   &temp.bin_value[4],
                   &temp.bin_value[5],
                   &temp.bin_value[6],
                   &temp.bin_value[7]);
  if(ret != 8) {
    std::cerr << "CheckpointSetup couldn't read required data from binary!\n";
    exit(EXIT_FAILURE);
  }
  // Fix endianness, implementation in lib/Endian.h
  temp.uint_value = ftoh64(temp.uint_value);
  return temp.uint_value;
}