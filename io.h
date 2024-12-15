#ifndef _IO_H
#define _IO_H

#include "particle.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>

std::vector<Particle*> read_file(char* file_name, int* num_particles);

void write_file(char* file_name, std::vector<Particle*> particles, int num_particles);

#endif