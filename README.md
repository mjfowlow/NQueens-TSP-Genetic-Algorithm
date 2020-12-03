# NQueens-TSP-Genetic-Algorithm
A Genetic Algorithm that solves the NQueens and Travelling Sales Person problems. 
Written using java.

The Problems:

The NQueens problem involves placing N queens on a NxN chess board such that none of the queens check eachother. That is none are in the same row, column, or diagonal. 

For the Travelling Salesperson Problem (TSP) you are given a list of cities and the distances between each pair of cities, the goal is to find the shortest possible route that visits each city only once and returns to the starting city. 

To see the TSP the algorithm solves:
See TSPGraph.PNG

The Solution:

For a 12x12 board there are over 100 quadrillion ways to arrange the queens. So testing every possible arrangment takes to much work and time.

Instead of testing every possible solution this algorithm uses an evolutionary approach. Where an intial population of candidate solutions is created and then through parent selection, recombination, mutation and survivor selection a whole new population is created which should be more "fit" then the last.

The goal is to find strong candidate solutions and create even stronger solutions from the candidate solutions selected. Leading to an optimal solution. It follows the theory of "Survival of the fittest".

To Run:

- DownLoad the files in src folder.

- Navigate to the folder.

- Compile the code with javac *.java

- Run:  java GeneticAlgorithm  prob ps cps mr sp mg

       prob = problem            [NQ = NQueens, TSP = Travelling sales person]
       ps = populationSize       [Recommended Respectively: 100, 5]
       cps = childPopulationSize [Recommended Respectively: 300, 10]
       mr = mutationRate         [Recommended: 0.8]
       sp = selectionPressure    [Recommended: 1.7] 
       mg = maxGenerations       [1000]

       Example: java GeneticAlgorithm NQ 100 300 0.8 1.7 1000

                                OR

       Example: java GeneticAlgorithm TSP 5 10 0.8 1.7 1000

To-Come:
Create a GUI representation.
Allow for different methods of selection, recombination, and mutation.
Give the choice to solve more problems.