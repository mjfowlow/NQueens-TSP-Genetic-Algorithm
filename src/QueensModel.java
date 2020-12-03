import java.util.Arrays;
import java.util.Collections;

/*
*MarK Fowlow
*Started: November 22, 2020
*Last Updated: December 2, 2020
*
*This class is a model of the NQueens problem
*Which will be used by GeneticAlgorithm class
*/

public class QueensModel {

    public static int boardSize = 12;

    //This method creates a random genotype for NQueens problem
    //The genotype is a permutation with range [1 to boardSize] no repetition 
    public static Integer [] getRandomGenotype() {
        
        Integer [] genotype = new Integer [boardSize];
        //Create an array with values between 0 and boardSize
        for (int i = 0; i < boardSize; i++) {
            genotype[i] = i+1;
        }
        //Shuffle the list 
        Collections.shuffle(Arrays.asList(genotype));
        
        return genotype;
    }    

    //Method to calculate MaxFitness
    //What this method calculates is how amy pairs of queens can be check 
    public static int getMaxFitness() {
        int maxFitness = (int) (0.5 * boardSize * (boardSize - 1));
        return maxFitness;
    }

    
    //This method takes the maximum pair of possible queens in check
    //and subtracts 1 for every pair queens in check
    //So if get value of 0 it means every pair of queens is in check
    //Higher Value - Higher fitness
    public static int fitnessFunction(Integer [] genotype) {
        
        //Calculate maximum pair of queens in check
        int fitness = (int) (0.5 * boardSize * (boardSize - 1));
        for (int col = 0; col < boardSize - 1; col++ ) {
            int row = genotype[col];
            int column = col;
            for(int y = row + 1; y < boardSize+1; y++) {
                column += 1;
                if (column == boardSize) {
                    break;
                }
                if(genotype[column] == y) {
                    fitness -= 1;
                }
            }
            column = col;
            for(int y = row - 1; y > 0; y--) {
                column += 1;
                if (column == boardSize) {
                    break;
                }
                if(genotype[column] == y){
                    fitness -= 1;
                }
            }
        }
        return fitness;
    }

}
