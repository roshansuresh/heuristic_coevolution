package seakers.datastructure;

/**
 * Data structure to store the heuristic weights at each NFE of the main coevolutionary algorithm
 */

public class HeuristicWeightsData {

    public int NFE;
    public double[] heuristicWeights;
    public double[] fitness;
    public int numberOfFeasibleSolutions;

    public HeuristicWeightsData(int numberOfHeuristics, int numberOfOptimizationObjectives) {
        this.NFE = 0;
        this.fitness = new double[numberOfOptimizationObjectives]; // Number of optimization objectives is used only if CMOEA is true
        this.numberOfFeasibleSolutions = 0;
        this.heuristicWeights = new double[numberOfHeuristics];
    }

    public int getNFE() {
        return this.NFE;
    }

    public void setNFE(int nfe) {
        this.NFE = nfe;
    }

    public double[] getHeuristicWeights() {
        return this.heuristicWeights;
    }

    public void setHeuristicWeights(double[] heuristicWeights) {
        this.heuristicWeights = heuristicWeights;
    }

    public double[] getFitness() {
        return this.fitness;
    }

    public void setFitness(double[] fitness) {
        this.fitness = fitness;
    }

    public int getNumberOfFeasibleSolutions() {return this.numberOfFeasibleSolutions; }

    public void setNumberOfFeasibleSolutions(int numberOfFeasibleDesigns) {this.numberOfFeasibleSolutions = numberOfFeasibleDesigns; }
}
