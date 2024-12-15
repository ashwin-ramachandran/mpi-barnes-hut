# Solving an N-body problem using the Barnes-Hut algorithm & MPI

The code in this repository solves an N-body problem using the Barnes-Hut algorithm. It deal with dynamic changing data as the bodies move during the simulation.

## Background

The particular N-body problem that this code solves is an astrophysical simulation that models the continuous interaction of each body with every other body. In this case the interactions are caused by the gravitational force that each body exerts on every other body, where each body could represent a star, a planet, or a galaxy. As the simulation proceeds, the forces will cause the bodies to move, so this is a dynamic computation.

### Barnes-Hut Algorithm

There are several algorithms for solving the N-Body problem. The simplest solution directly computes the n*n interactions that exist between every pair of bodies. The complexity of this approach is O(n^2), where n is the number of bodies in the system. Instead we can use a superior solution, the Barnes-Hut method, which takes O(n log n) time; this solution computes direct interactions for pairs of bodies that are sufficiently close; for the effects of far away bodies, the algorithm approximates a group of bodies by a single body whose center of mass is the same as the center of mass of the group of bodies.

[This article](http://arborjs.org/docs/barnes-hut) provides a great introduction and explanation of the Barnes-Hut algorithm.

The algorithm consists of two phases. The first phase constructs a tree that represents a spatial partitioning of the bodies in the system. Each interior node in the tree represents the center of mass of the bodies in the subtree rooted at that node. The second phase computes the force interactions for each body in the system by using the MAC (multipole acceptance criteria), which determines whether the bodies in a subtree are sufficiently far away that their net force on a given body can be approximated by the center of mass of the bodies in a subtree. The velocity and position of the body is then calculated via a Leapfrog-Verlet integrator. Since we are approximating the effect of group of bodies by their center of mass, this algorithm takes O(n log n) steps. These two phases of the algorithm will occur iteratively, for every single step/frame, for however many steps/frames are specified.

### MAC (Multipole Acceptance Criteria)

Whether a body is far away from a subtree depends on the distance d between the body and the center of mass of the subtree, the length l of the side of the subtree, and a threshold theta. The values for d and l are computed during run-time, and theta is a command line argument to the program.

### Force

The gravitational force between two bodies (or one body and a subtree) are given as:

```
F = G*M0*M1*d / d^3
```

where M0 and M1 are the mass of the two bodies (or one body and a subtree) and d is the distance between the two bodies (or one body and the center of mass of a subtree). However, when two bodies are too close, the force computed using the above formula can be infinitely large. To avoid infinity values, an rlimit value is set such that if the distance between two bodies are shorter than rlimit, then rlimit is used as the distance between the two bodies. Note that the force computed here is a vector, which means forces from different bodies (or subtrees) cannot be added directly. What is needed is to project the force on x and y axis and then the projected forces on each axis are addable. The projected force on each axis can be computed as:

```
Fx = G*M0*M1*dx / d^3
Fy = G*M0*M1*dy / d^3
```

where dx and dx is the distance between the two bodies (or one body and the center of mass of a subtree) on the x and y axis.

### Leapfrog-Verlet Integration

Using the total force on a body, the new position (Px', Py') and velocity (Vx', Vy') of the body given its current position (Px, Py) and velocity (Vx, Vy) can be computed as:

```
ax = Fx/M0
ay = Fy/M0
Px' = Px + Vx*dt + 0.5*ax*dt^2
Py' = Py + Vy*dt + 0.5*ay*dt^2
Vx' = Vx + ax*dt
Vy' = Vy + ay*dt
```

where dt is the timestep.

### Parameters

The following values are used for the above parameters in this implementation:

* `G`: 0.0001
* `rlimit`: 0.03
* `dt`: 0.005

### Flags

The program accepts the following four command line parameters, with types indicated in parentheses:

* `-i inputfilename` (char *): input file name
* `-o outputfilename` (char *): output file name
* `-s steps` (int): number of iterations
* `-t \theta` (double): threshold for MAC
* `-d dt`(double): timestep

Both the input and output files will have the following format:

* The number of bodies (n)
* A list of n tuples:
  * index
  * x position
  * y position
  * mass
  * x velocity
  * y velocity

Note that the index of each particle is considered as the particle's identifier.

### Performance

The program also outputs the elapsed time Iin seconds) for the simulation.

### Domain Boundary

When the position of a particle is out of the boundary, the particle is considered to be lost. For lost particles, the particle's mass is set to -1, and it is not included in the force calculation any more. The x and y coordinates of each particle are within the range (0, 4). Therefore, the domain is the square (0,0) -- (4,4).

## Message Passing Interface (MPI)

[MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) is standard designed for C/C++ to enable parallel computing. The code in this repo further optimizes beyond the advantages of the Barnes-Hut algorithm, by using MPI to have multiple processes partion the input work and "divide and conquer" for the most computationally expensive part of the algorithm, which is force computation on each particle. So, for example, if there are 100,000 input particles and 4 processes, then each process will be given 25,000 of the particles, and all 4 processes will parallelly compute the forces on its assigned particles, rather than having a single process sequentially compute this for 100,000 particles.

## Compling & Running

This project requires MPICH. You can download source from: [https://www.mpich.org/](https://www.mpich.org/) and install it with the following commands:

```
 ./configure; make; sudo make install
```

Compile:

`make`

Run:

`mpirun -np <num_processes> ./nbody -i <input_file_path> -o <output_file_path> -s <num_steps (num_frames)> -t <theta_threshold_value> -d <timestep_value>`
