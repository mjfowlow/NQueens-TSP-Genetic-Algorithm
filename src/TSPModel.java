import java.util.Arrays;
import java.util.Collections;

/*
Mark Fowlow
November 30, 2020
The program will represent the travelling sales man problem 
The graph it will use in saved with the folder ("TSPGraph.png")
This program will be used by GeneticAlgorithm
*/

public class TSPModel {
    
    public static int towns = 7;
    private static int[][] roadLengths = {{0,5,7,7,11,6,6},
                                   {5,0,4,8,14,11,9},
                                   {7,4,0,5,6,9,10},
                                   {7,8,5,0,6,8,11},
                                   {11,14,6,6,0,4,6},
                                   {6,11,9,8,4,0,5},
                                   {6,9,10,11,6,5,0}                                    
                                  }; 

    //The genotype is going to be a permutation 
    //Values 1-towns
    //No repeats
    public static Integer [] getRandomGenotype() {
        
        Integer [] genotype = new Integer [towns];
        
        //Make first value 1 to represent intial town
        genotype[0] = 1;

        //Create an array with values between 1 and number of town
        for (int i = 1; i < towns; i++) {
            genotype[i] = i+1;
        }
        //Shuffle the list 
        Collections.shuffle(Arrays.asList(genotype));
        Collections.swap(Arrays.asList(genotype), 0, Arrays.asList(genotype).indexOf(1));
        
        return genotype;
    }

    //Subtract from worst path
    //So if worth past is selected will get fitness of 0
    //Higher Value = Greater Fitness
    public static int fitnessFuncton(Integer[] genotype) {
        
        int totalCost = 0;
        for(int i = 0; i < genotype.length-1; i++) {
            int startingPoint = genotype[i];
            int finishingPoint = genotype[i+1];
            int cost = calculateCost(startingPoint,finishingPoint);
            totalCost += cost;
        }

        //Calculate the cost to return back to intial town (town 1)
        int lastTownVisited = genotype[genotype.length-1];
        totalCost += calculateCost(1, lastTownVisited);
        
        
        return 73-totalCost;
    }

    private static int calculateCost(int startingPoint, int finishingPoint) {
        int cost = roadLengths[startingPoint-1][finishingPoint-1];
        return cost; 
    }

}
