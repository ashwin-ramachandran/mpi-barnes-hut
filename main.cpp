#include <argparse.h>
#include "particle.h"
#include "io.h"
#include "quadtree.h"
#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <cmath>

int get_start_end_index(int num_part, int num_processes, int process_id, int& start, int&end) {
  int particles_per_process = num_part / num_processes;
  int extra_processes = num_part % num_processes; // This is the number of processes that will have 1 extra particle (particles_per_process + 1)
  int num_processes_with_particles_per_process = num_processes - extra_processes;  // This is the number of processes that will have a normal number of particles (particles_per_process)
  if (process_id < num_processes_with_particles_per_process) {
    start = process_id * particles_per_process;
    end = std::min(start + particles_per_process, num_part);
  }
  else {
    int offset = particles_per_process * num_processes_with_particles_per_process;
    particles_per_process++;
    start = offset + (process_id - num_processes_with_particles_per_process) * particles_per_process;
    end = std::min(start + particles_per_process, num_part);
  }
  return particles_per_process;
}

int main(int argc, char* argv[]) {
  // MPI set up
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Parse command line args
  struct options_t opts;
  get_opts(argc, argv, &opts);

  // Read and parse input file to get particles
  int num_particles;
  std::vector<Particle*> particles = read_file(opts.in_file, &num_particles);

  int start_index, end_index;
  int particles_per_process = get_start_end_index(num_particles, size, rank, start_index, end_index);

  // Use a barrier to ensure that all processes reach here before marking start time
  MPI_Barrier(MPI_COMM_WORLD);
  double start_time = MPI_Wtime();

  // Steps represents "frames"
  for (int step = 0; step < opts.steps; step++) {
    // All processes will build the whole tree based off of the current position of all points
    QuadTreeNode* root = new QuadTreeNode(0.0, 4.0, 0.0, 4.0);
    for (int i = 0; i < num_particles; i++) {
      if (particles[i]->mass > 0) {
        insert_particle(root, particles[i]);
      }
    }

    // Compute all the forces acting on each particle (each process only does this for the particles that it was assigned)
    std::vector<double> new_forces(particles_per_process * 2);
    for (int i = start_index; i < end_index; i++) {
      double force_x = 0.0;
      double force_y = 0.0;
      if (particles[i]->mass > 0) {
        compute_net_force(root, particles[i], opts.theta, force_x, force_y);
      }
      new_forces[(i - start_index) * 2] = force_x;
      new_forces[((i - start_index) * 2) + 1] = force_y;
    }

    // Free/delete the tree
    delete_tree(root);

    std::vector<double> all_new_forces(num_particles * 2);
    if (rank == 0) {
      // Store results from this process itself
      for (int i = start_index; i < end_index; i++) {
        all_new_forces[i * 2] = new_forces[i * 2];
        all_new_forces[(i * 2) + 1] = new_forces[(i * 2) + 1];
      }

      // Receive results from all other processes
      for (int i = 1; i < size; i++) {
        int start, end;
        get_start_end_index(num_particles, size, i, start, end);
        MPI_Recv(&(all_new_forces[start * 2]), (end - start) * 2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      // Broadcast all_new_forces to all other processes
      MPI_Bcast(all_new_forces.data(), num_particles * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else {
      // Send results to process rank 0
      MPI_Send(new_forces.data(), particles_per_process * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

      // Wait to receive updated forces from process 0 and store it in all_new_forces
      MPI_Bcast(all_new_forces.data(), num_particles * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Update all the particle's positions and velocities
    for (int i = 0; i < num_particles; i++) {
      if (particles[i]->mass > 0) {
        double acceleration_x = all_new_forces[i * 2] / particles[i]->mass;
        double acceleration_y = all_new_forces[(i * 2) + 1] / particles[i]->mass;

        particles[i]->x_pos += ((particles[i]->x_vel * opts.timestep) + (0.5 * acceleration_x * opts.timestep * opts.timestep));
        particles[i]->y_pos += ((particles[i]->y_vel * opts.timestep) + (0.5 * acceleration_y * opts.timestep * opts.timestep));

        particles[i]->x_vel += (acceleration_x * opts.timestep);
        particles[i]->y_vel += (acceleration_y * opts.timestep);

        if (!particle_in_bounds(particles[i], 0.0, 4.0, 0.0, 4.0)) {
          particles[i]->mass = -1;
        }
      }
    }
  }

  // Use a barrier to ensure that all processes reach here before marking end time
  MPI_Barrier(MPI_COMM_WORLD);

  // Write results to output file and record compute time
  if (rank == 0) {
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    std::cout << elapsed_time << std::endl;
    write_file(opts.out_file, particles, num_particles);
  }

  MPI_Finalize();

  return 0;
}