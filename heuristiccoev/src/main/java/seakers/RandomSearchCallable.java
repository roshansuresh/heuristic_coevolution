package seakers;

import org.moeaframework.algorithm.RandomSearch;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.problem.AbstractProblem;
import seakers.datastructure.HeuristicWeightsData;
import seakers.problem.CoevolutionaryProblem;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

public class RandomSearchCallable implements Callable<Algorithm> {

    public AbstractProblem coevolutionaryProblem;
    public List<Solution> initialWeightsSolutions;
    private int numberOfOptimizationObjectives;
    public final int weightsMaxEvaluations;
    public final boolean evolvePopulation;
    private Result result;
    private final int runNumber;

    public RandomSearchCallable(AbstractProblem coevolutionaryProblem, List<Solution> initialWeightsSolutions, int numberOfOptimizationObjectives, int weightsMaxEvaluations, boolean evolvePopulation, Result result, int runNumber) {
        this.coevolutionaryProblem = coevolutionaryProblem;
        this.initialWeightsSolutions = initialWeightsSolutions;
        this.numberOfOptimizationObjectives = numberOfOptimizationObjectives;
        this.weightsMaxEvaluations = weightsMaxEvaluations;
        this.evolvePopulation = evolvePopulation;
        this.result = result;
        this.runNumber = runNumber;
    }


    @Override
    public Algorithm call() throws Exception {
        // Run Random Search
        System.out.println("Starting Random Search");

        int numberOfWeightsObjectives = 1;
        if (((CoevolutionaryProblem) coevolutionaryProblem).getCMOEA()) {
            numberOfWeightsObjectives = numberOfOptimizationObjectives;
        }

        long startTime = System.currentTimeMillis();

        // Evaluate initial weights population
        int NFE = 0;
        for (Solution initialWeightsSolution : initialWeightsSolutions) {
            coevolutionaryProblem.evaluate(initialWeightsSolution);
            NFE++;
        }

        // Sample and evaluate random weights solutions until max NFE is reached
        while (NFE < weightsMaxEvaluations) {

            // Create new weights solution and update evolvePopulation variable accordingly (since newSolution initializes the evolvePopulation variable as False)
            Solution newWeightsSolution = coevolutionaryProblem.newSolution();
            BinaryVariable var = new BinaryVariable(1);
            EncodingUtils.setBoolean(var, evolvePopulation);
            newWeightsSolution.setVariable(newWeightsSolution.getNumberOfVariables()-1, var);

            coevolutionaryProblem.evaluate(newWeightsSolution);

            NFE++;
        }

        long endTime = System.currentTimeMillis();
        System.out.println("Random Search completed in " + ((endTime - startTime)/1000.0) + "s");

        ArrayList<HeuristicWeightsData> allWeights = ((CoevolutionaryProblem) coevolutionaryProblem).getAllWeights();

        result.saveHeuristicWeights(allWeights, numberOfWeightsObjectives, ((CoevolutionaryProblem) coevolutionaryProblem).getOptimizationProblem().getNumberOfConstraints(), "run " + runNumber + File.separator + "random_search_heuristic_weights.csv", true);

        return null;
    }
}
