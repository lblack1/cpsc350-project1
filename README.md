# cpsc350-project1
Lloyd Black
2295968
CPSC350, Section 2
Rene German

A program for basic file reading, statistic computation, and file writing sans data structures (excepting strings). Specifically, operates on nucleotides.

1. Source Files - DNAReader.h, DNAReader.cpp, main.cpp

2. Issues - If you say yes to entering another text file, there's no graceful way to exit the program without doing so.

3. Resources - cplusplus.com/reference, stackoverflow, Data Structures and Algorithms in C (by Goodrich, Tamassia, and Mount).

4. Description of Program - Reads in text files and computes a number of statistics based upon the number of A's, C's, T's, and G's (nucleotides) per line, such as mean line length, variance of line length, relative probability of each nucleotide, relative probability of each nucleotide bigram, et cetera. Writes the results to a file "lloydblack.out" in the same directory, then generates 1000 lines composed using the same mean, variance and relative nucleotide probability in the same file. Allows for multiple files to be read (one at a time).

5. After compiling - ./YOUROUTFILE.out [TEXT FILE TO READ]
