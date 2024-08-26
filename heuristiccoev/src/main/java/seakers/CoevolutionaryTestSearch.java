package seakers;

import org.apache.commons.math3.genetics.GeneticAlgorithm;
import org.moeaframework.algorithm.AbstractAlgorithm;
import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.core.comparator.DominanceComparator;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;
import seakers.datastructure.HeuristicWeightsData;
import seakers.problem.CoevolutionaryProblem;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.*;

public class CoevolutionaryTestSearch implements Callable<Algorithm> {

    private Algorithm coevolutionaryAlgorithm;
    //private String saveDirectory;
    //private int numberOfHeuristics;
    private int numberOfOptimizationObjectives;
    private int maxCoevolutionEvaluations;
    private Result result;
    private DominanceComparator coevolutionComparator;
    private boolean cMOEA;
    private final int runNumber;
    private final boolean isEOSS;
    private final boolean isPartitioning;
    private final boolean evolvePopulation;
    private final boolean integerWeights;
    private final boolean weightOfWeights;
    private final boolean periodicZeroInjection;

    public CoevolutionaryTestSearch(Algorithm algorithm, int maxCoevolutionEvaluations, int numberOfOptimizationObjectives, Result result, DominanceComparator coevolutionComparator, boolean cMOEA, boolean isEOSS, boolean isPartitioning, boolean evolvePopulation, boolean integerWeights, boolean weightOfWeights, boolean periodicZeroInjection, int runNumber) {
        this.coevolutionaryAlgorithm = algorithm;
        //this.saveDirectory = saveDirectory;
        this.maxCoevolutionEvaluations = maxCoevolutionEvaluations;
        //this.numberOfHeuristics = numberOfHeuristics;
        this.numberOfOptimizationObjectives = numberOfOptimizationObjectives;
        this.result = result;
        this.coevolutionComparator = coevolutionComparator;
        this.cMOEA = cMOEA;
        this.isEOSS = isEOSS;
        this.isPartitioning = isPartitioning;
        this.runNumber = runNumber;
        this.evolvePopulation = evolvePopulation;
        this.integerWeights = integerWeights;
        this.weightOfWeights = weightOfWeights;
        this.periodicZeroInjection = periodicZeroInjection;
    }

    @Override
    public Algorithm call() throws Exception {
        // Run Outer MOEA (evolving population of heuristic weights)
        System.out.println("Starting Coevolutionary Algorithm");

        int numberOfWeightsObjectives = 1;
        if (((CoevolutionaryProblem) coevolutionaryAlgorithm.getProblem()).getCMOEA()) {
            numberOfWeightsObjectives = numberOfOptimizationObjectives;
        }

        long startTime = System.currentTimeMillis();
        coevolutionaryAlgorithm.step();

        System.out.println("NFE = 0");

        while (!coevolutionaryAlgorithm.isTerminated() && (coevolutionaryAlgorithm.getNumberOfEvaluations() < maxCoevolutionEvaluations)) {
            if (periodicZeroInjection) {
                Solution zeroSolution = ((CoevolutionaryProblem) coevolutionaryAlgorithm.getProblem()).newSolution();
                for (int i = 0; i < zeroSolution.getNumberOfVariables(); i++) {
                    if (zeroSolution.getVariable(i) instanceof RealVariable) {
                        RealVariable zeroVar = null;
                        if (weightOfWeights) {
                            if (i == zeroSolution.getNumberOfVariables()-1) { // last variable in weightOfWeights run is the weight of weights
                                zeroVar = new RealVariable(0, ((RealVariable) zeroSolution.getVariable(i)).getLowerBound(), ((RealVariable) zeroSolution.getVariable(i)).getUpperBound());
                            }
                        } else {
                            if (integerWeights) {
                                zeroVar = new RealVariable(0, ((RealVariable) zeroSolution.getVariable(i)).getLowerBound(), ((RealVariable) zeroSolution.getVariable(i)).getUpperBound());
                            } else {
                                zeroVar = new RealVariable(-3, ((RealVariable) zeroSolution.getVariable(i)).getLowerBound(), ((RealVariable) zeroSolution.getVariable(i)).getUpperBound());
                            }
                        }
                        zeroSolution.setVariable(i, zeroVar);
                    }
                }
                BinaryVariable evolvePopulationVariable = new BinaryVariable(1);
                if (evolvePopulation) {
                    EncodingUtils.setBoolean(evolvePopulationVariable, true);
                } else {
                    EncodingUtils.setBoolean(evolvePopulationVariable, false);
                }
                zeroSolution.setVariable(zeroSolution.getNumberOfVariables()-1, evolvePopulationVariable);
                coevolutionaryAlgorithm.getProblem().evaluate(zeroSolution); // NOTE: The zero solution will also be part of the allWeights arraylist, leading to more than the usual number of weights evaluated
                //coevolutionaryAlgorithm.getProblem().evaluate(zeroSolution);

                // Remove least fit solution and add zero solution
                Population currentPopulation = ((AbstractEvolutionaryAlgorithm) coevolutionaryAlgorithm).getPopulation();
                Solution lowestFitnessSolution = currentPopulation.get(0);
                for (int i = 1; i < currentPopulation.size(); i++) {
                    if (coevolutionComparator.compare(lowestFitnessSolution, currentPopulation.get(i)) == -1) { // high objective means low fitness (solution objectives must be minimized)
                        lowestFitnessSolution = currentPopulation.get(i);
                    }
                }
                ((AbstractEvolutionaryAlgorithm) coevolutionaryAlgorithm).getPopulation().remove(lowestFitnessSolution);
                ((AbstractEvolutionaryAlgorithm) coevolutionaryAlgorithm).getPopulation().add(zeroSolution);
            }
            coevolutionaryAlgorithm.step();
        }

        coevolutionaryAlgorithm.terminate();
        long endTime = System.currentTimeMillis();

        if (isEOSS) {
            if (isPartitioning) {
                ((PartitioningProblem) ((CoevolutionaryProblem) coevolutionaryAlgorithm.getProblem()).getOptimizationProblem()).getEvaluationManager().clear();
            } else {
                ((AssigningProblem) ((CoevolutionaryProblem) coevolutionaryAlgorithm.getProblem()).getOptimizationProblem()).getEvaluationManager().clear();
            }
        }

        System.out.println("Coevolutionary Algorithm completed in " + ((endTime - startTime)/1000.0) + "s");

        ArrayList<HeuristicWeightsData> allWeights = ((CoevolutionaryProblem) coevolutionaryAlgorithm.getProblem()).getAllWeights();

        Population finalPopulation = ((AbstractEvolutionaryAlgorithm) coevolutionaryAlgorithm).getPopulation();
        result.saveFinalPopulationOrArchive(finalPopulation, "run " + runNumber + File.separator + "coevolutionary_algorithm_heuristic_weights_finalpop.csv", numberOfWeightsObjectives, ((CoevolutionaryProblem) coevolutionaryAlgorithm.getProblem()).getOptimizationProblem().getNumberOfConstraints());

        if (cMOEA) {
            NondominatedPopulation finalArchive = ((AbstractEvolutionaryAlgorithm) coevolutionaryAlgorithm).getArchive();
            result.saveFinalPopulationOrArchive(finalArchive, "run " + runNumber + File.separator + "coevolutionary_algorithm_heuristic_weights_finalarchive.csv", numberOfWeightsObjectives, ((CoevolutionaryProblem) coevolutionaryAlgorithm.getProblem()).getOptimizationProblem().getNumberOfConstraints());
        }

        result.saveHeuristicWeights(allWeights, numberOfWeightsObjectives, ((CoevolutionaryProblem) coevolutionaryAlgorithm.getProblem()).getOptimizationProblem().getNumberOfConstraints(), "run " + runNumber + File.separator + "coevolutionary_algorithm_heuristic_weights.csv", true);

        return coevolutionaryAlgorithm;
    }

}
