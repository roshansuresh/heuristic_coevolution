package seakers;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.AggregateConstraintComparator;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.DominanceComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.indicator.Hypervolume;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.core.operator.TournamentSelection;
import org.moeaframework.core.operator.real.PM;
import org.moeaframework.core.operator.real.SBX;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.problem.CDTLZ.C1_DTLZ1;
import seakers.datastructure.HeuristicWeightsData;
import seakers.problem.CoevolutionaryProblem;
import seakers.problem.ModifiedC_DTLZProblem;

import java.util.*;
import java.util.stream.IntStream;

public class FullFactorialIntegerWeights {

    public static void main(String[] args) {
        // Save location
        String saveDir = "C:\\SEAK Lab\\SEAK Lab Github\\Coevolution-based Heuristic Incorporation\\results";

        boolean cMOEA = false; // if the coevolutionary optimization is single or multi objective
        boolean evolveInternalPopulation = false; // whether to update internal population for subsequent heuristic weights
        boolean integerWeights = true; // whether to use integer formulation for heuristic weights
        boolean weightOfWeights = false; // if an additional weight multiplicative parameter design decision for the outer GA is used
        boolean periodicZeroInjection = false; // if zero solution is added for evaluation to population at different NFE

        // Internal Test problem
        AbstractProblem problem = new C1_DTLZ1(2); // number of variables is (n_objs + k - 1) with recommended k = 5
        int numberOfHeuristics = 2; // 2 heuristics for C_DTLZ problems
        AbstractProblem heuristicProblem = new ModifiedC_DTLZProblem(problem, numberOfHeuristics, integerWeights);

        String[] variableNames = new String[problem.getNumberOfVariables()];
        String[] objectiveNames = new String[problem.getNumberOfObjectives()];
        String[] constraintNames = new String[problem.getNumberOfConstraints()];
        String[] heuristicNames = new String[]{"Heuristic 1", "Heuristic 2"};

        for (int i = 0; i < variableNames.length; i++) {
            variableNames[i] = "Variable " + i;
        }

        for (int i = 0; i < objectiveNames.length; i++) {
            objectiveNames[i] = "True Objective " + i;
        }

        for (int i = 0; i < constraintNames.length; i++) {
            constraintNames[i] = "Constraint " + i;
        }

        int internalPopulationSize = 25;
        int internalMaxEvaluations = 200;

        List<Solution> internalInitialSolutions = new ArrayList<>();
        for (int i = 0; i < internalPopulationSize; i++) {
            Solution initialSolution = problem.newSolution();
            for (int j = 0; j < initialSolution.getNumberOfVariables(); j++) {
                initialSolution.getVariable(j).randomize();
            }
            heuristicProblem.evaluate(initialSolution);
            internalInitialSolutions.add(initialSolution);
        }
        //Initialization initialization = new InjectedInitialization(problem, internalPopulationSize, internalInitialSolutions);
        Population initialInternalPopulation = new Population(internalInitialSolutions);
        DominanceComparator internalComparator = new ChainedComparator(new AggregateConstraintComparator(), new ParetoObjectiveComparator());

        // Internal Operators
        double crossoverProbability = 1.0;
        double mutationProbability = 1.0/problem.getNumberOfVariables();
        Variation crossover = new SBX(crossoverProbability, 0.5);
        Variation mutation = new PM(mutationProbability, 0.5);
        Selection internalSelection = new TournamentSelection(2, internalComparator);
        Variation internalVariation = new CompoundVariation(crossover, mutation);

        double[] epsilonBox = new double[heuristicProblem.getNumberOfObjectives()];
        Arrays.fill(epsilonBox, 0.0001);
        EpsilonBoxDominanceArchive internalArchive = new EpsilonBoxDominanceArchive(epsilonBox);

        int numberOfWeightsObjectives = 1;
        if (cMOEA) {
            numberOfWeightsObjectives = heuristicProblem.getNumberOfObjectives();
        }

        double[] internalObjectivesMinimum = new double[heuristicProblem.getNumberOfObjectives()];
        Arrays.fill(internalObjectivesMinimum, 0.0);
        double[] internalObjectivesMaximum = new double[heuristicProblem.getNumberOfObjectives()];
        Arrays.fill(internalObjectivesMaximum, 1.0);
        Hypervolume internalHypervolume = new Hypervolume(heuristicProblem, internalObjectivesMinimum, internalObjectivesMaximum);

        Result result;
        if (weightOfWeights) {
            result = new Result(saveDir, numberOfHeuristics+1);
        } else {
            result = new Result(saveDir, numberOfHeuristics);
        }

        // Initialize full factorial integer weights
        ArrayList<ArrayList<Integer>> weightsCombinations = new ArrayList<>();
        int[] integerWeightOptions = IntStream.range(0, 11).toArray();
        Iterator<int[]> iterator = CombinatoricsUtils.combinationsIterator(integerWeightOptions.length, numberOfHeuristics);
        while (iterator.hasNext()) {
            int[] combination = iterator.next();
            ArrayList<Integer> weightCombination = new ArrayList<>();
            for (int value : combination) {
                weightCombination.add(integerWeightOptions[value]);
            }
            weightsCombinations.add(weightCombination);
        }

        // Add repeated weights since apache commons combinations iterator does not add them
        for (int integerWeightOption : integerWeightOptions) {
            ArrayList<Integer> repeatedWeights = new ArrayList<>();
            for (int j = 0; j < numberOfHeuristics; j++) {
                repeatedWeights.add(integerWeightOption);
            }
            weightsCombinations.add(repeatedWeights);
        }

        AbstractProblem coevolutionaryProblem = new CoevolutionaryProblem(heuristicProblem, saveDir, weightsCombinations.size(), internalMaxEvaluations, initialInternalPopulation, internalArchive, internalComparator, internalSelection, internalVariation, internalHypervolume, cMOEA, evolveInternalPopulation, integerWeights, weightOfWeights, variableNames, objectiveNames, constraintNames, heuristicNames, heuristicNames, numberOfWeightsObjectives, result, false, false, 0);

        // Evaluate each weight combination and save the results
        ArrayList<HeuristicWeightsData> generatedSolutions = new ArrayList<>();
        Iterator<ArrayList<Integer>> weightsIterator = weightsCombinations.iterator();
        while (weightsIterator.hasNext()) {
            ArrayList<Integer> currentWeights = weightsIterator.next();
            Solution weightsSolution = coevolutionaryProblem.newSolution();

            // Set weights as variables to solution
            for (int i = 0; i < currentWeights.size(); i++) {
                RealVariable var = new RealVariable(currentWeights.get(i), 0, 10);
                weightsSolution.setVariable(i, var);
            }

            // Evaluate solution
            coevolutionaryProblem.evaluate(weightsSolution);

            // Save results of evaluation into HeuristicWeightsData instance
            HeuristicWeightsData currentData = new HeuristicWeightsData(numberOfHeuristics, numberOfWeightsObjectives);
            result.fillHeuristicWeightsData(currentData, weightsSolution, 0, numberOfWeightsObjectives, problem.getNumberOfConstraints());
            generatedSolutions.add(currentData);
        }

        // Save data into csv
        result.saveHeuristicWeights(generatedSolutions, numberOfWeightsObjectives, problem.getNumberOfConstraints(),"integer_weights_full_factorial_results.csv", false);
    }

}
