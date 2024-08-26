package seakers;

import org.moeaframework.core.*;
import org.moeaframework.core.comparator.*;
import org.moeaframework.core.indicator.Hypervolume;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.core.operator.TournamentSelection;
import org.moeaframework.core.operator.real.PM;
import org.moeaframework.core.operator.real.SBX;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.problem.CDTLZ.C1_DTLZ1;
import seakers.problem.CoevolutionaryProblem;
import seakers.problem.ModifiedC_DTLZProblem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ZeroWeightsInternalMOEA {

    public static void main (String[] args) {
        // Save location
        String saveDir = "C:\\SEAK Lab\\SEAK Lab Github\\Coevolution-based Heuristic Incorporation\\results";

        boolean cMOEA = false; // if the coevolutionary optimization is single or multi objective
        boolean evolveInternalPopulation = false; // whether to update internal population for subsequent heuristic weights
        boolean integerWeights = false; // whether to use integer formulation for heuristic weights
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

        int internalPopulationSize = 20;
        int internalMaxEvaluations = 100;

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

        double[] internalObjectivesMinimum = new double[heuristicProblem.getNumberOfObjectives()];
        Arrays.fill(internalObjectivesMinimum, 0.0);
        double[] internalObjectivesMaximum = new double[heuristicProblem.getNumberOfObjectives()];
        Arrays.fill(internalObjectivesMaximum, 1.0);
        Hypervolume internalHypervolume = new Hypervolume(heuristicProblem, internalObjectivesMinimum, internalObjectivesMaximum);

        // Initialize zero weights solution
        int numberOfObjectives = 1;
        if (cMOEA) {
            numberOfObjectives = 2;
        }
        Solution zeroWeightsSolutions = new Solution(numberOfHeuristics+1, numberOfObjectives, 0);
        for (int i = 0; i < numberOfHeuristics ; i++) {
            RealVariable var = new RealVariable(0.0, 0.0, 1.0);
            zeroWeightsSolutions.setVariable(i, var);
        }
        BinaryVariable binaryVar = new BinaryVariable(1);
        EncodingUtils.setBoolean(binaryVar, false);
        zeroWeightsSolutions.setVariable(numberOfHeuristics, binaryVar);

        Result result;
        if (weightOfWeights) {
            result = new Result(saveDir, numberOfHeuristics+1);
        } else {
            result = new Result(saveDir, numberOfHeuristics);
        }

        AbstractProblem internalProblem = new CoevolutionaryProblem(heuristicProblem, saveDir, 1, internalMaxEvaluations, initialInternalPopulation, internalArchive, internalComparator, internalSelection, internalVariation, internalHypervolume, cMOEA, evolveInternalPopulation, integerWeights, weightOfWeights, variableNames, objectiveNames, constraintNames, heuristicNames, heuristicNames, numberOfObjectives, result, false, false, 0);
        internalProblem.evaluate(zeroWeightsSolutions);

        System.out.println("Zero weights solution fitness = " + Arrays.toString(zeroWeightsSolutions.getObjectives()));
    }
}