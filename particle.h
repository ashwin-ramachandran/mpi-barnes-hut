#ifndef POINT_H
#define POINT_H

#include <string>
#include <sstream>

struct Particle {
  // identifier of the particle
  int index;
  
  // (x, y) coordinates of the particle's position
  double x_pos, y_pos;

  // mass of the particle
  double mass;

  // x and y components of particle's velocity
  double x_vel, y_vel;

  // Custom toString method
  std::string toString() const {
    std::ostringstream oss;
    oss << "\n" << "index: " << index 
        << "\n" << "x_pos: " << x_pos 
        << "\n" << "y_pos: " << y_pos 
        << "\n" << "mass: " << mass 
        << "\n" << "x_vel: " << x_vel 
        << "\n" << "y_vel: " << y_vel 
        << "\n";
    return oss.str();
  }
};

#endif