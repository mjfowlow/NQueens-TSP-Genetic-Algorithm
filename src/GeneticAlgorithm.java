import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/*
*Mark Fowlow
*Start: November 22, 2020
*Last Updated: December, 2020
*To Run:
*       java GeneticAlgorithm  prob ps cps mr sp mg
*
*       prob = problem            [NQ = NQueens, TSP = Travelling sales person]
*       ps = populationSize       [Recommended Respectively: 100, 5]
*       cps = childPopulationSize [Recommended Respectively: 300, 10]
*       mr = mutationRate         [0.8]
*       sp = selectionPressure    [1.8] 
*       mg = maxGeneration        [1000, 100]
*/

public class GeneticAlgorithm {

    private static Integer[][] population;

    //Parameter Values (Given as input)
    private static int populationSize;
    private static int childPopulationSize;
    private static double mutationRate;
    private static double selectionPressure;
    private static String problem;

    //termination Criteria
    private static int maxGenerations;
    private static int maxFitness;
    
    //Intialization
    private static int generation;
    private static int maxFitnessSeen;
    private static Integer [] bestSolution;
    private static int diversity;
    private static int boardSize;
    private static int numTowns;

    public static void main(String[] args) {
        
        //Get parameters from user input
        problem = args[0];
        populationSize = Integer.parseInt(args[1]);
        childPopulationSize = Integer.parseInt(args[2]);
        mutationRate = Double.parseDouble(args[3]);
        selectionPressure = Double.parseDouble(args[4]);
        maxGenerations = Integer.parseInt(args[5]);

        //Set parameters
        generation = 0;
        maxFitnessSeen = 0;

        //Intialize population and maxfitness depending on problem
        if(problem.equals("NQ")) {
            maxFitness = QueensModel.getMaxFitness();
            //Get boardSize
            boardSize = QueensModel.boardSize;
            //Use boardSize to intialize population
            population = new Integer [populationSize][boardSize];
            for(int i = 0; i < populationSize; i++) {
                population[i] =  QueensModel.getRandomGenotype();
            }
        }
        else {
            //MaxFitness for TSP is currently just a given
            maxFitness = 38;
            //Get numTowns
            numTowns = TSPModel.towns;
            //Use numTowns to intialize population
            population = new Integer [populationSize][TSPModel.towns];
            for(int i = 0; i < populationSize; i++) {
                population[i] = TSPModel.getRandomGenotype();
            }
        }

        //Set bestSolution as first genotype in population and calculate divsersity of population
        bestSolution = population[0];
        diversity = fitnessDiversity(population);

        //Print intial stats
        printStats(0, bestSolution);

        //Run the simulation as long as termination criteria is not met 
        while(generation < maxGenerations && maxFitnessSeen < maxFitness && diversity > 1) {
            runGeneration();  
        }
        
        //Once terminated, print the best solution
        System.out.println("Evolution Ended");
        System.out.println("Best fitness seen: " + maxFitnessSeen); 
        System.out.print("Genotype: ");
        for (int i = 0; i < bestSolution.length; i++)
        {
            System.out.print(" " + bestSolution[i]);
        }
        System.out.println();
    }
    
    //This is where a full generation happens.
    //parent selection, child creation, mutation, and survivor selection
    //
    //
    //Arguements:
    //          None
    //Return:
    //          None
    private static void runGeneration(){
        
        //Intialize the children populaton depending on the problem
        Integer [][] children;
        if (problem.equals("NQ")) {
            children = new Integer [childPopulationSize][boardSize];
        }
        else {
            children = new Integer [childPopulationSize][numTowns];
        }

        //Produce the children
        for (int i = 0; i < childPopulationSize; i = i + 2)
            {
                //Select two parents
                Integer [][] parents = parentSelection(population, selectionPressure);
                //Then recombine the parents to produce to children
                Integer [][] tempKids = recombine(parents[0],parents[1]);
                //Mutate the children
                mutation(tempKids[0], mutationRate);
                mutation(tempKids[1], mutationRate);
                //Add the new children to the new population of children
                children[i] = tempKids[0];
                children[i+1] = tempKids[1];
            }
            
        // perform survivor selection
        population = survivorSelection(children, populationSize);
    
        //Get the highest fitness of this generation. 
        //First in population since survivorSelection picks the best fitnesses first
        Integer[] bestGeno = population[0];
        //Get fitness of the gene
        int fitness;
        if(problem.equals("NQ")) {
            fitness = QueensModel.fitnessFunction(bestGeno);
        }
        else{
            fitness = TSPModel.fitnessFuncton(bestGeno);
        }
        
        //Update maxFitness seen if the max this generation is greater
        if (fitness > maxFitnessSeen)
        {
            maxFitnessSeen = fitness;
            bestSolution = bestGeno;
        }
            
        // check for diversity and increment the generation
        diversity = fitnessDiversity(population);
        generation++;

        //print stats after each generation     
        printStats(generation, bestGeno);
    }

    
    //performs linear ranking parent selection 
    //
    //Arguements:
    //          Integer[][] population   - All the genotypes in current population
    //          double selectionPressure - selectionPressure is value used by probability function
    //Return:
    //          Integer[][] parents      - returns two parents selected from population
    public static Integer[][] parentSelection(Integer[][] population, double selectionPressure) {
        
        //Intiate parents
        Integer [][] parents;
        if(problem.equals("NQ")){
            parents = new Integer [2][boardSize];
        }
        else {   
            parents = new Integer [2][TSPModel.towns];
        }
        
        //Determine the fitness of each member of the population
        Integer [] fitness = getFitnesses(population);
        //Determine the ranking of each member by its fitness
        double [] rank = rankFitnesses(fitness);
        //Determine the cumulative probability of selecting each member, using the linear ranking formula
        double [] cumulative = new double [populationSize];
        cumulative[0] = linearRankingProb(rank[0], selectionPressure, populationSize);
        for (int i = 1; i < populationSize; i ++)
        {
            cumulative[i] = cumulative[i-1] + linearRankingProb(rank[i], selectionPressure, populationSize);
        }
        // Select two parents, based on the cumulative probabilities
        //Get random number between 0 and 1
        double randNum = Math.random();
        //Get first parent
        int first = 0;
        while(cumulative[first] < randNum) {
            first++;
        }
        //Get second parent and make sure it is not the same as first parent
        int second = first;
        while(second == first) {
            randNum = Math.random();
            second = 0;
            while (cumulative[second] < randNum) {
                second++;
            }       
        }
        //Select the two parents from population
        parents[0] = population[first];
        parents[1] = population[second];

        return parents;
    }
    
    //Returns the rank of fitnesses
    //Worst = 0; best = values.length-1
    //Ranking is shared if tied:
    //For input values [30, 70, 50, 50, 40, 80] the rankings are [0, 4, 2.5, 2.5, 1, 5]
    public static double[] rankFitnesses(Integer values[])
    {
        double rank[] = new double[values.length];
        
        for (int i = 0; i < values.length; i++) {
            //check if the value already has a ranking
            //This could happen when dealing with duplicates
            if(rank[i] != 0) {
                continue;
            }
            //Get currentValue to find ranking for
            int currentValue = values[i];
            //Intialize ranking to 0
            double ranking = 0;
            //Create a duplicates list
            ArrayList<Integer> duplicates = new ArrayList<Integer>();
            duplicates.add(i);

            for (int j =0; j < values.length; j++) {
                if(values[j] < currentValue) {
                    ranking++;
                }
                if(values[j] == currentValue) {
                    if (i != j) {
                        duplicates.add(j);
                    }    
                }
            }
            if (duplicates.size() > 1){
                double nextRanking = ranking + 1;
                for(int k = 1; k < duplicates.size(); k++ ) {
                    ranking += nextRanking;
                    nextRanking++;
                }
                ranking  = ranking/duplicates.size();
                for(int q = 0; q < duplicates.size(); q++ ) {
                    rank[duplicates.get(q)] = ranking; 
                }
            }else{
                rank[i]  = ranking;
            }
        }
        return rank;
    }

    //Calculates the linear ranking probability of a genotype
    public static double linearRankingProb(double rank, double s, int populationSize)
    {
        double probability = ((2-s)/populationSize)+(((2*rank)*(s-1))/(populationSize*(populationSize-1)));
        return probability;
    }

    // returns an array of fitnesses for a population
    private static Integer[] getFitnesses(Integer[][] population)
    {
        Integer[] fitness = new Integer[populationSize];
        
        for (int i = 0; i < populationSize; i++){
            //Calculate fitness depening on problem
            if(problem.equals("NQ")) {
                fitness[i] = QueensModel.fitnessFunction(population[i]);
            }
            else{    
                fitness[i] = TSPModel.fitnessFuncton(population[i]);
        
            }
        }
        return fitness;
    }
    
    // creates 2 child genotypes using the cut-and-crossfill method
    public static Integer [][] recombine(Integer[] parent0, Integer[] parent1)
    {
        //Intiate children depending on problem
        Integer [][] children; 
        if(problem.equals("NQ")) {
            children = new Integer[2][boardSize];
        }
        else{
            children = new Integer[2][TSPModel.towns];
        }

        //Calculate children
        Integer[] child0 = cutAndCrossFill(parent0, parent1);
        Integer[] child1 = cutAndCrossFill(parent1, parent0);
        
        //Set children equal to parents
        children[0] = child0;
        children[1] = child1;
        return children;
    }

    private static Integer[] cutAndCrossFill(Integer[] parent0, Integer[] parent1) {
        //Set child variables
        Integer[] child0 = parent0.clone();
        Integer[] child1 = parent1.clone();

        //Get length and crossover point
        int lenGeno = child0.length;
        int crossover = lenGeno/2;  
        for (int i = crossover; i < lenGeno; i++ ) {
            //Starting at crossover in parent1 to crossover-1 
            //check if duplicate and swap when find not duplicate
            for (int j = crossover; j < crossover + lenGeno; j++) {
                //get next possible value to swap
                int swapValue = child1[j%lenGeno];
                //Check if duplicate
                boolean duplicate = false; 
                for (int k = 0; k < i; k++) {
                    if (swapValue == child0[k]) {
                        duplicate = true;
                    }
                }
                //If duplicate increment positon in parent 1
                if (duplicate) {
                    continue;
                }
                //If not duplicate swap values 
                else {
                    child0[i] = swapValue;
                    break;
                }   
            }
        }
        return child0;
    }

    // inverts the order of a series of genes in the genotype
    public static Integer[] mutation(Integer[] genotype, double mR) {
        Random rand = new Random();
        //Create a random number between 1 (inclusive) and 11 (exclsuive)
        int randNum = rand.nextInt(10) + 1;         

        if (randNum <= mR*10) {
            
            //Get two random values. Range depends on the problem
            int randomInt1;
            int randomInt2;
            if(problem.equals("NQ")) {
                //Create to random numbers within the range on the genoptype size
                randomInt1 = rand.nextInt(genotype.length);
                randomInt2 = rand.nextInt(genotype.length);
            }
            else {
                //Create to random numbers within the range of 1 (inclusive) and 7 (exclusive)
                //This way the first value in genotype is kept as 1
                randomInt1 = rand.nextInt(genotype.length-1) + 1;
                randomInt2 = rand.nextInt(genotype.length-1) + 1;
            }
            //Find the lower and upper alleles chosen
            int firstAllele;
            int secondAllele;
            if (randomInt1 < randomInt2) {
                firstAllele = randomInt1;
                secondAllele = randomInt2;
            } else {
                firstAllele = randomInt2;
                secondAllele = randomInt1;
            }
            //Create 2 temp values
            int temp1;
            int temp2;
            while(firstAllele < secondAllele) {
                //Have the temp values hold the orginial values
                temp1 = genotype[firstAllele];
                temp2 = genotype[secondAllele];
                //Swap the values
                genotype[secondAllele] = temp1;
                genotype[firstAllele] = temp2;
                //Increase the lower allele and decrease the upper allele
                firstAllele++;
                secondAllele--;
            } 
        }
        return genotype;
    }

    // creates a new population through (λ, μ) survivor selection
    //
    // this method will return the n fittest children as the new population
    public static Integer[][] survivorSelection(Integer[][] children, int n) {
        //Intiate the new populatio acorroding to the problem
        Integer [][] newPop;
        if(problem.equals("NQ")){
            newPop = new Integer [n][boardSize];
        }
        else{
            newPop = new Integer [n][numTowns];
        }

        //Get an array of the fitness values
        Integer fitnesses[] = getFitnesses(children);

        //Find the maximum child in children n times
        for (int i = 0; i < n; i++) {
            //Create two values to keep track of the max fitness and where it was located
            int maximumLocation = -1;
            int maximumValue  = -1;  
            //Find max
            for (int j = 0; j < fitnesses.length; j++) {
                if (fitnesses[j] > maximumValue) {
                    maximumLocation = j;
                    maximumValue = fitnesses[j];
                }
            }
            //Add the child to the survivor list
            newPop[i] = children[maximumLocation];
            //Set the maximum value in fitness to negative so it is ot selected again
            fitnesses[maximumLocation] = -1;
        }
        return newPop;
    }

    // counts the number of unique fitness values seen in the population
    public static int fitnessDiversity(Integer[][] population){

        //Get an array of the fitness values
        Integer fitnesses[] = getFitnesses(population);
        //Sort array so all equivalent fitness values are together
        Arrays.sort(fitnesses);
        int uniqueFitnessValues = 0;
        for (int i = 0; i < fitnesses.length; i++) {
            //Iterate through fitness list and skip over duplicates
            while(i < fitnesses.length - 1 && fitnesses[i] == fitnesses[i+1]) {
                i++;
            }
            uniqueFitnessValues++;
        }
        return uniqueFitnessValues;
    }

    //Prints stats so progres can been seen
    public static void printStats(int generation, Integer[] bestFitness){
        //Calculate fitness of gene past
        int fitness;
        if(problem.equals("NQ")) {
            fitness = QueensModel.fitnessFunction(bestFitness);
        }
        else{
            fitness = TSPModel.fitnessFuncton(bestFitness);
        }
        
        //print stats
        System.out.println("Generation: " + generation);
        System.out.println("Fitness: Max = " + fitness);
        System.out.println("Diversity: " + fitnessDiversity(population));
        System.out.println();
    }

}