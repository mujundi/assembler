# DNA Assembly using Ant Colony Optimization

This repository hosts 2 programs: 
The first (compiled from the source files in the "Reads Generator" folder) randomly generates reads from a FASTA file.
The second (from the "Assembler" folder) assembles the fragments into the original sequence.

# GENERAL DESCRIPTION

This program is the result of a very simple genome assembly challenge I found online at the following address:

https://www.reddit.com/r/dailyprogrammer/comments/46km7n/20160219_challenge_254_hard_dna_shotgun_sequencing/

My first solution was a brute force algorithm that, while sufficient for solving the challenge with only 38 fragments, scaled horribly. I set out to approach the problem from a new perspective that would challenge my programming skills further. After reading the Wikipedia entry on the Traveling Salesman Problem, I sought to adapt this assembly problem into a TSP, and use an optimization algorithm to solve it. The basic idea is that each of the reads represents a node, and the amount of overlap the read shares with other reads is the weight of each edge that connects that node to all the other reads. This forms a Hamiltonian Path Problem, which can be converted to a TSP by establishing an imaginary starting node that is equidistant to all nodes in the graph. To adapt the algorithm for this particular application, it is searching to maximize the distance traversed by the ants, not minimize. By maximizing the distance it is finding the path with most overlap, giving us the shortest length found for the assembled sequence.


For more info on Ant Colony Optimization, visit:

https://en.wikipedia.org/wiki/Ant_colony_optimization_algorithms



# COMPILATION INSTRUCTIONS

1.) Edit the code files to use appropriate file addresses for the FASTA file and the text file containing the reads.

2.) Compile the "reads_generator.cpp" file and run the executable.

3.) Move the reads file to the desired location.

4.) Compile all source files in the "Assembler" folder into one executable and run it.




# Notes: 

- I used only the first 5,000 base pairs of the E. coli genome to keep computation times relatively low. The original genome is included in the "Reads Generator" folder as "ecoli_K12_MG165(complete).fasta".

- The parameters of the assembler are still being tested to find optimal values. 

- The program will output a text file containing the assembled sequence, as well as the best solution found in each of the 3 trials in each iteration, allowing examination of how quickly the program converges to the optimal solution.

- OpenMP was used to parallelize the execution of the algorithm so that 3 trials may be run concurrently.

- An example of the code's output is included in the "Example" folder.
