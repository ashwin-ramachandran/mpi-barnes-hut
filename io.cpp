#include <io.h>

std::vector<Particle*> read_file(char* file_name, int* num_particles) {
  // Open file
  std::ifstream in_file;
  in_file.open(file_name);

  in_file >> *num_particles;

  std::vector<Particle*> particles;

  for (int i = 0; i < *num_particles; i++) {
    Particle* p = new Particle();
    in_file >> p->index >> p->x_pos >> p->y_pos >> p->mass >> p->x_vel >> p->y_vel;
    particles.push_back(p);
  }

  in_file.close();

  return particles;
}

void write_file(char* file_name, std::vector<Particle*> particles, int num_particles) {
  std::ofstream out_file;
  out_file.open(file_name, std::ofstream::trunc);

  out_file << num_particles << std::endl;
  
  // Set precision for writing doubles
  out_file << std::setprecision(std::numeric_limits<double>::max_digits10) << std::fixed;

  for (int i = 0; i < num_particles; i++) {
    out_file << particles[i]->index << " " 
             << particles[i]->x_pos << " " 
             << particles[i]->y_pos << " " 
             << particles[i]->mass << " "
             << particles[i]->x_vel << " " 
             << particles[i]->y_vel << std::endl;
  }

  out_file.flush();
  out_file.close();
}