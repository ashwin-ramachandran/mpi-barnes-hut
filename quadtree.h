#ifndef QUADTREE_H
#define QUADTREE_H

#include "particle.h"
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>

#define RLIMIT 0.03
#define G 0.0001

struct QuadTreeNode {
  // Paticle value (null if the node is an internal node)
  Particle* particle;

  // Boundary of the region repented by the internal node (N/A if the node is a leaf)
  // Use this value to determine the s-value (width of this region)
  double x_min, x_max, y_min, y_max;

  // Total mass of all child particles (N/A if the node is a leaf)
  double total_mass;

  // (x, y) coordinates of the center of mass (N/A if the node is a leaf)
  double center_of_mass_x_numerator, center_of_mass_y_numerator;

  // Pointers to the 4 children of this node
  QuadTreeNode* NW;
  QuadTreeNode* NE; 
  QuadTreeNode* SW; 
  QuadTreeNode* SE;

  // Constructor to default values
  QuadTreeNode() :
    particle(NULL),
    x_min(0.0), x_max(0.0), y_min(0.0), y_max(0.0),
    total_mass(0.0),
    center_of_mass_x_numerator(0.0), center_of_mass_y_numerator(0.0),
    NW(NULL), NE(NULL), SW(NULL), SE(NULL) {}

  // Constructor with boundaries
  QuadTreeNode(double x_min, double x_max, double y_min, double y_max) :
    particle(NULL),
    x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max),
    total_mass(0.0),
    center_of_mass_x_numerator(0.0), center_of_mass_y_numerator(0.0),
    NW(NULL), NE(NULL), SW(NULL), SE(NULL) {}

  // Custom toString method
  std::string toString() const {
    std::ostringstream oss;
    oss << "\n" << "particle: " << (particle == NULL ? -1 : particle->index) 
        << "\n" << "x_min: " << x_min 
        << "\n" << "x_max: " << x_max 
        << "\n" << "y_min: " << y_min 
        << "\n" << "y_max: " << y_max 
        << "\n" << "total_mass: " << total_mass 
        << "\n" << "center_of_mass_x_numerator: " << center_of_mass_x_numerator
        << "\n" << "center_of_mass_y_numerator: " << center_of_mass_y_numerator
        << "\n" << "center_of_mass: (" << center_of_mass_x_numerator / total_mass << ", " << center_of_mass_y_numerator / total_mass << ")"
        << "\n" << "NW: " << NW 
        << "\n" << "NE: " << NE 
        << "\n" << "SW: " << SW 
        << "\n" << "SE: " << SE 
        << "\n";
    return oss.str();
  }
};


bool particle_in_bounds(Particle* p, double x_min, double x_max, double y_min, double y_max);

void delete_tree(QuadTreeNode* node);

void insert_particle(QuadTreeNode* node, Particle* p);

QuadTreeNode* create_new_child_quadrant(QuadTreeNode* quadrant, double x_min, double x_max, double y_min, double y_max);

double compute_distance(double x1, double y1, double x2, double y2);

void compute_net_force(QuadTreeNode* node, Particle* p, double theta, double& force_x, double& force_y);

std::string get_quad_tree_string(const QuadTreeNode* node, int level = 0);

#endif