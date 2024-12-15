#include <argparse.h>

void get_opts(int argc, char** argv, struct options_t *opts) {
  if (argc == 1) {
    std::cout << "Usage:" << std::endl;
    std::cout << "\t-i <file_path>" << std::endl;
    std::cout << "\t-o <file_path>" << std::endl;
    std::cout << "\t-s <steps>" << std::endl;
    std::cout << "\t-t <theta>" << std::endl;
    std::cout << "\t-d <timestep>" << std::endl;
    exit(0);
  }

  struct option l_opts[] = {
    {"in_file", required_argument, NULL, 'i'},
    {"out_file", required_argument, NULL, 'o'},
    {"steps", required_argument, NULL, 's'},
    {"theta", required_argument, NULL, 't'},
    {"timestep", no_argument, NULL, 'd'}
  };

  int ind, c;
  while ((c = getopt_long(argc, argv, "i:o:s:t:d:", l_opts, &ind)) != -1) {
    switch (c) {
    case 0:
      break;
    case 'i':
      opts->in_file = (char *)optarg;
      break;
    case 'o':
      opts->out_file = (char *)optarg;
      break;
    case 's':
      opts->steps = atoi((char *)optarg);
      break;
    case 't':
      opts->theta = std::stod((char *)optarg);
      break;
    case 'd':
      opts->timestep = std::stod((char *)optarg);
      break;
    case ':':
      std::cerr << argv[0] << ": option -" << (char)optopt << "requires an argument." << std::endl;
      exit(1);
    }
  }
}
