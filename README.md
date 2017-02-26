# assembler

This repository hosts 2 programs: 
The first (compiled from the source files in the "Reads Generator" folder) randomly generates reads from a FASTA file.
The second (from the "Assembler" folder) assembles the fragments into the original sequence.


COMPILATION INSTRUCTIONS

1.) Edit the code files to use appropriate file addresses for the FASTA file and the text file containing the reads.
2.) Compile the "reads_generator.cpp" file and run the executable.
3.) Move the reads file to the desired location.
4.) Compile all source files in the "Assembler" folder into one executable and run it.




Notes: 

- I used only the first 5,000 base pairs of the E. coli genome to keep computation times relatively low. The original genome is included in the "Reads Generator" folder as "ecoli_K12_MG165(complete).fasta".

- The parameters of the assembler are still being tested to find optimal values. 

- The program will output a text file containing the assembled sequence, as well as the best solution found in each of the 3 trials in each iteration, allowing examination of how quickly the program converges to the optimal solution.

- OpenMP was used to parallelize the execution of the algorithm so that 3 trials may be run concurrently.

- An example of the code's output is included in the "Example" folder.
