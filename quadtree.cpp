#include "quadtree.h"

bool particle_in_bounds(Particle* p, double x_min, double x_max, double y_min, double y_max) {
  return p->x_pos >= x_min && 
         p->x_pos <= x_max &&
         p->y_pos >= y_min &&
         p->y_pos <= y_max;
}

QuadTreeNode* create_new_child_quadrant(QuadTreeNode* quadrant, double x_min, double x_max, double y_min, double y_max) {
  if (quadrant == NULL) {
    return new QuadTreeNode(x_min, x_max, y_min, y_max);
  }
  return quadrant;
}

void delete_tree(QuadTreeNode* node) {
  if (node == NULL) {
    return;
  }

  delete_tree(node->NW);
  delete_tree(node->NE);
  delete_tree(node->SW);
  delete_tree(node->SE);

  delete node;
}

void insert_particle(QuadTreeNode* node, Particle* p) {
  // Base case: If the node does not contain a particle or children (is an empty node/quadrant), just put the particle here
  if (node->particle == NULL && 
      node->NW == NULL && 
      node->NE == NULL && 
      node->SW == NULL && 
      node->SE == NULL) {

        node->particle = p;
        return;
  }
  // Recursive case: If the node is currently an internal node or a leaf node (already directly contains a particle), 
  // then we must sub-divide the node into 4 children and/or determine the appropriate quadrant, then recursively 
  // call insert to insert the existing particle and current particle (if applicable) into appropriate quadrants
  // Calculate the 4 childrent/quadrant's new boundaries
  double new_NW_x_min = node->x_min;
  double new_NW_x_max = (node->x_min + node->x_max) / 2;
  double new_NW_y_min = (node->y_min + node->y_max) / 2;
  double new_NW_y_max = node->y_max;

  double new_NE_x_min = (node->x_min + node->x_max) / 2;
  double new_NE_x_max = node->x_max;
  double new_NE_y_min = (node->y_min + node->y_max) / 2;
  double new_NE_y_max = node->y_max;

  double new_SW_x_min = node->x_min;
  double new_SW_x_max = (node->x_min + node->x_max) / 2;
  double new_SW_y_min = node->y_min;
  double new_SW_y_max = (node->y_min + node->y_max) / 2;

  double new_SE_x_min = (node->x_min + node->x_max) / 2;
  double new_SE_x_max = node->x_max;
  double new_SE_y_min = node->y_min;
  double new_SE_y_max = (node->y_min + node->y_max) / 2;

  // Figure out appropriate quadrant for new particle
  if (particle_in_bounds(p, new_NW_x_min, new_NW_x_max, new_NW_y_min, new_NW_y_max)) {
    node->NW = create_new_child_quadrant(node->NW, new_NW_x_min, new_NW_x_max, new_NW_y_min, new_NW_y_max);
    insert_particle(node->NW, p);
  }
  else if (particle_in_bounds(p, new_NE_x_min, new_NE_x_max, new_NE_y_min, new_NE_y_max)) {
    node->NE = create_new_child_quadrant(node->NE, new_NE_x_min, new_NE_x_max, new_NE_y_min, new_NE_y_max);
    insert_particle(node->NE, p);
  }
  else if (particle_in_bounds(p, new_SW_x_min, new_SW_x_max, new_SW_y_min, new_SW_y_max)) {
    node->SW = create_new_child_quadrant(node->SW, new_SW_x_min, new_SW_x_max, new_SW_y_min, new_SW_y_max);
    insert_particle(node->SW, p);
  }
  else if (particle_in_bounds(p, new_SE_x_min, new_SE_x_max, new_SE_y_min, new_SE_y_max)) {
    node->SE = create_new_child_quadrant(node->SE, new_SE_x_min, new_SE_x_max, new_SE_y_min, new_SE_y_max);
    insert_particle(node->SE, p);
  }
  node->total_mass += p->mass;
  node->center_of_mass_x_numerator += (p->x_pos * p->mass);
  node->center_of_mass_y_numerator += (p->y_pos * p->mass);

  if (node->particle != NULL) {
    // Figure out appropriate quadrant for current particle
    if (particle_in_bounds(node->particle, new_NW_x_min, new_NW_x_max, new_NW_y_min, new_NW_y_max)) {
      node->NW = create_new_child_quadrant(node->NW, new_NW_x_min, new_NW_x_max, new_NW_y_min, new_NW_y_max);
      insert_particle(node->NW, node->particle);
    }
    else if (particle_in_bounds(node->particle, new_NE_x_min, new_NE_x_max, new_NE_y_min, new_NE_y_max)) {
      node->NE = create_new_child_quadrant(node->NE, new_NE_x_min, new_NE_x_max, new_NE_y_min, new_NE_y_max);
      insert_particle(node->NE, node->particle);
    }
    else if (particle_in_bounds(node->particle, new_SW_x_min, new_SW_x_max, new_SW_y_min, new_SW_y_max)) {
      node->SW = create_new_child_quadrant(node->SW, new_SW_x_min, new_SW_x_max, new_SW_y_min, new_SW_y_max);
      insert_particle(node->SW, node->particle);
    }
    else if (particle_in_bounds(node->particle, new_SE_x_min, new_SE_x_max, new_SE_y_min, new_SE_y_max)) {
      node->SE = create_new_child_quadrant(node->SE, new_SE_x_min, new_SE_x_max, new_SE_y_min, new_SE_y_max);
      insert_particle(node->SE, node->particle);
    }
    node->total_mass += node->particle->mass;
    node->center_of_mass_x_numerator += (node->particle->x_pos * node->particle->mass);
    node->center_of_mass_y_numerator += (node->particle->y_pos * node->particle->mass);

    // NULL out the node's particle attribute since it is now an internal node and no longer a leaf node
    node->particle = NULL;
  }
}

double compute_distance(double x1, double y1, double x2, double y2) {
  double distance_x = x2 - x1;
  double distance_y = y2 - y1;
  double distance = std::sqrt((distance_x * distance_x) + (distance_y * distance_y));
  return std::max(RLIMIT, distance);
}

void compute_net_force(QuadTreeNode* node, Particle* p, double theta, double& force_x, double& force_y) {
  // Cover the cases that one of the quadrants/child nodes that is passed in does not exist
  if (node == NULL) {
    return;
  }

  // If the node is a leaf node, return the force exerted on the particle p by the particle in the leaf node
  if (node->particle != NULL) {
    double distance_x = node->particle->x_pos - p->x_pos;
    double distance_y = node->particle->y_pos - p->y_pos;
    double euclidean_distance = compute_distance(node->particle->x_pos, node->particle->y_pos, p->x_pos, p->y_pos);
    force_x += ((G * node->particle->mass * p->mass * distance_x) / (euclidean_distance * euclidean_distance * euclidean_distance));
    force_y += ((G * node->particle->mass * p->mass * distance_y) / (euclidean_distance * euclidean_distance * euclidean_distance));
    return;
  }

  // If the node is an internal node, calculate l and d and then check if the l/d ratio is less than theta
  // If l/d < theta, then the internal node/body is "far enough away" so just compute the force using the
  // body's center of mass
  double quadrant_length = node->x_max - node->x_min;  // Equavaliently could use y_max - y_min
  double center_of_mass_x = node->center_of_mass_x_numerator / node->total_mass;
  double center_of_mass_y = node->center_of_mass_y_numerator / node->total_mass;
  double euclidean_distance = compute_distance(center_of_mass_x, center_of_mass_y, p->x_pos, p->y_pos);
  if ((quadrant_length / euclidean_distance) < theta) {
    double distance_x = center_of_mass_x - p->x_pos;
    double distance_y = center_of_mass_y - p->y_pos;
    force_x += ((G * node->total_mass * p->mass * distance_x) / (euclidean_distance * euclidean_distance * euclidean_distance));
    force_y += ((G * node->total_mass * p->mass * distance_y) / (euclidean_distance * euclidean_distance * euclidean_distance));
    return;
  }

  // Otherwise, recursively keep traversing to the body's/internal node's children/quadrants
  compute_net_force(node->NW, p, theta, force_x, force_y);
  compute_net_force(node->NE, p, theta, force_x, force_y);
  compute_net_force(node->SW, p, theta, force_x, force_y);
  compute_net_force(node->SE, p, theta, force_x, force_y);
}

std::string get_quad_tree_string(const QuadTreeNode* node, int level) {
    if (node == nullptr) {
        return "";  // Return an empty string for null nodes
    }

    std::ostringstream oss;
    std::string indent(level * 4, ' ');  // 4 spaces per level

    // Set the precision for floating point values
    oss << std::fixed << std::setprecision(std::numeric_limits<double>::digits10); // Max precision for double

    if (node->particle != nullptr) {
        // Leaf node: contains a particle
        oss << indent << "Leaf Node - Particle Index: " << node->particle->index
            << ", Position: (" << node->particle->x_pos << ", "
            << node->particle->y_pos << "), Mass: " << node->particle->mass
            << "\n";
    } else {
        // Internal node: contains total mass and center of mass
        oss << indent << "Internal Node - Total Mass: " 
            << node->total_mass << " Center of Mass: (" 
            << node->center_of_mass_x_numerator / node->total_mass << ", " 
            << node->center_of_mass_y_numerator / node->total_mass << ")\n";

        // Build strings for child nodes
        oss << indent << "  NW: ";
        if (node->NW != nullptr) {
            oss << get_quad_tree_string(node->NW, level + 1);
        } else {
            oss << "null\n";
        }

        oss << indent << "  NE: ";
        if (node->NE != nullptr) {
            oss << get_quad_tree_string(node->NE, level + 1);
        } else {
            oss << "null\n";
        }

        oss << indent << "  SW: ";
        if (node->SW != nullptr) {
            oss << get_quad_tree_string(node->SW, level + 1);
        } else {
            oss << "null\n";
        }

        oss << indent << "  SE: ";
        if (node->SE != nullptr) {
            oss << get_quad_tree_string(node->SE, level + 1);
        } else {
            oss << "null\n";
        }
    }

    return oss.str();
}
