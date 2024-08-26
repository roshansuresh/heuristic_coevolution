package seakers;

import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.algorithm.single.AggregateObjectiveComparator;
import org.moeaframework.algorithm.single.GeneticAlgorithm;
import org.moeaframework.algorithm.single.LinearObjectiveComparator;
import org.moeaframework.algorithm.single.MinMaxObjectiveComparator;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.*;
import org.moeaframework.core.indicator.Hypervolume;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.real.PM;
import org.moeaframework.core.operator.real.SBX;
import org.moeaframework.core.operator.real.UM;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.problem.CDTLZ.C1_DTLZ1;
import seakers.operators.SimulatedBinaryWeightsCrossover;
import seakers.operators.UniformIntegerWeightsMutation;
import seakers.operators.UniformWeightsCrossover;
import seakers.problem.CoevolutionaryProblem;
import seakers.problem.ModifiedC_DTLZProblem;
import seakers.problem.SimpleDTLZ1_2Problem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.*;

public class CoevolutionaryAlgorithmTest {

    public static ExecutorService pool;
    public static CompletionService<Algorithm> cs;

    public static void main(String[] args) {

        // Save location
        //String saveDir = System.getProperty("user.dir") + File.separator + "results"; // File not found error!
        String saveDir = "C:\\SEAK Lab\\SEAK Lab Github\\Coevolution-based Heuristic Incorporation\\results";

        boolean cMOEA = false; // if the coevolutionary optimization is single or multi objective
        boolean evolveInternalPopulation = true; // if initial internal population for subsequent heuristic weights should be updated or not
        boolean integerWeights = true; // if formulation for the heuristic weights should be integer or real
        boolean weightOfWeights = false; // if an additional weight multiplicative parameter design decision for the outer GA is used
        boolean periodicZeroInjection = false; // if zero solution is added for evaluation to population at different NFE

        // Internal Test problem
        //AbstractProblem problem = new C1_DTLZ1(2); // number of variables is (n_objs + k - 1) with recommended k = 5
        AbstractProblem problem = new SimpleDTLZ1_2Problem();
        int numberOfHeuristics = 2; // 2 heuristics for C_DTLZ problems
        AbstractProblem heuristicProblem = new ModifiedC_DTLZProblem(problem, numberOfHeuristics, integerWeights);
        int internalPopulationSize = 20;
        int internalMaxEvaluations = 200;
        //int numberOfInternalCores = 1; // ExecutorCompletionService no longer used in CoevolutionaryProblem class

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

        List<Solution> internalInitialSolutions = new ArrayList<>();
        for (int i = 0; i < internalPopulationSize; i++) {
            Solution initialSolution = problem.newSolution();
            for (int j = 0; j < initialSolution.getNumberOfVariables(); j++) {
                initialSolution.getVariable(j).randomize();
            }
            heuristicProblem.evaluate(initialSolution);
            internalInitialSolutions.add(initialSolution);
        }
        Population initialInternalPopulation = new Population(internalInitialSolutions);
        DominanceComparator internalComparator;
        if (heuristicProblem.getNumberOfConstraints() == 0) {
            internalComparator = new ParetoObjectiveComparator();
        } else {
            internalComparator = new ChainedComparator(new AggregateConstraintComparator(), new ParetoObjectiveComparator());
        }

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

        Result result;
        if (weightOfWeights) {
            result = new Result(saveDir, numberOfHeuristics+1);
        } else {
            result = new Result(saveDir, numberOfHeuristics);
        }

        int numberOfCoevolutionaryObjectives = 1;
        if (cMOEA) {
            numberOfCoevolutionaryObjectives = problem.getNumberOfObjectives();
        }

        // Initialize MOEA
        int coevolutionPopulationSize = 10;
        int coevolutionMaxEvaluations = 50;

        // Initialize Coevolutionary Problem class
        double[] internalObjectivesMinimum = new double[heuristicProblem.getNumberOfObjectives()];
        Arrays.fill(internalObjectivesMinimum, 0.0);
        double[] internalObjectivesMaximum = new double[heuristicProblem.getNumberOfObjectives()];
        Arrays.fill(internalObjectivesMaximum, 1.0);
        Hypervolume internalHypervolume = new Hypervolume(heuristicProblem, internalObjectivesMinimum, internalObjectivesMaximum);
        AbstractProblem coevolutionaryProblem = new CoevolutionaryProblem(heuristicProblem, saveDir, coevolutionPopulationSize, internalMaxEvaluations, initialInternalPopulation, internalArchive, internalComparator, internalSelection, internalVariation, internalHypervolume, cMOEA, evolveInternalPopulation, integerWeights, weightOfWeights, variableNames, objectiveNames, constraintNames, heuristicNames, heuristicNames, numberOfCoevolutionaryObjectives, result, false, false, 0);

        //Initialization coevolutionInitialization = new RandomInitialization(coevolutionaryProblem, coevolutionPopulationSize);
        //Population coevolutionPopulation = new Population(coevolutionInitialization.initialize());

        // NOTE: Injected initialization is used since RandomInitialization will randomize all variable in a new solution (including the last binary variable
        // signifying whether to update initial internal population or not, which we don't want. The initial population of weights must use the same initial internal population
        List<Solution> initialCoevolutionSolutions = new ArrayList<>();
        for (int i = 0; i < coevolutionPopulationSize; i++) {
            initialCoevolutionSolutions.add(coevolutionaryProblem.newSolution()); // This new solution by default generates a random set of weights with the internal population not to be updated

        }
        Initialization coevolutionInitialization = new InjectedInitialization(coevolutionaryProblem, coevolutionPopulationSize, initialCoevolutionSolutions);
        Population coevolutionPopulation = new Population();

        Variation coevolutionaryCrossover;
        Variation coevolutionaryMutation;
        if (integerWeights) {
            coevolutionaryCrossover = new UniformWeightsCrossover(1.0);
            coevolutionaryMutation = new UniformIntegerWeightsMutation(1.0/numberOfHeuristics);
        } else {
            coevolutionaryCrossover = new SimulatedBinaryWeightsCrossover(1.0, 0.5);
            coevolutionaryMutation = new PM(1.0/numberOfHeuristics, 0.5); // UM modifies only the RealVariable decisions
        }

        double epsilonDouble = 0.0001;
        double[] epsilonBoxDouble = new double[problem.getNumberOfObjectives()];
        Arrays.fill(epsilonBoxDouble, epsilonDouble);
        EpsilonBoxDominanceArchive coevolutionArchive = new EpsilonBoxDominanceArchive(epsilonBoxDouble);
        DominanceComparator coevolutionComparator;
        if (cMOEA) {
            coevolutionComparator = new ParetoObjectiveComparator();
        } else {
            coevolutionComparator = new LinearObjectiveComparator();
        }

        Variation coevolutionVariation = new CompoundVariation(coevolutionaryCrossover, coevolutionaryMutation);

        Selection coevolutionSelection = new TournamentSelection(2, coevolutionComparator);

        Algorithm coevolutionMOEA;
        if (cMOEA) {
            coevolutionMOEA = new EpsilonMOEA(coevolutionaryProblem, coevolutionPopulation, coevolutionArchive, coevolutionSelection, coevolutionVariation, coevolutionInitialization, coevolutionComparator);
        } else {
            coevolutionMOEA = new GeneticAlgorithm(coevolutionaryProblem, (AggregateObjectiveComparator) coevolutionComparator, coevolutionInitialization, coevolutionSelection, coevolutionVariation);
        }


        int numberOfCores = 1;
        pool = Executors.newFixedThreadPool(numberOfCores);
        cs = new ExecutorCompletionService<>(pool);

        cs.submit(new CoevolutionaryTestSearch(coevolutionMOEA, coevolutionMaxEvaluations, heuristicProblem.getNumberOfObjectives(), result, coevolutionComparator, cMOEA, false, false, evolveInternalPopulation, integerWeights, weightOfWeights, periodicZeroInjection, 0));

        try {
            cs.take().get();
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        pool.shutdown();

    }

}

