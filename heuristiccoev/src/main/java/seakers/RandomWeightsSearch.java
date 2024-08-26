package seakers;

import com.mathworks.engine.EngineException;
import com.mathworks.engine.MatlabEngine;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.AggregateConstraintComparator;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.DominanceComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.indicator.Hypervolume;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.core.operator.OnePointCrossover;
import org.moeaframework.core.operator.TournamentSelection;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;
import seakers.problem.CoevolutionaryProblem;
import seakers.problem.metamaterial.ModifiedArteryProblem;
import seakers.problem.metamaterial.ModifiedEqualNormalStiffnessProblem;
import seakers.trussaos.initialization.SynchronizedMersenneTwister;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.*;

public class RandomWeightsSearch {

    public static ExecutorService pool;
    public static CompletionService<Algorithm> ecs;
    /**
     * Matlab Engine for function evaluation
     */
    private static MatlabEngine engine;

    public static void main(String[] args) throws InterruptedException, ExecutionException {

        // Save location
        //String saveDir = System.getProperty("user.dir") + File.separator + "results"; // File not found error!
        String saveDir = "C:\\SEAK Lab\\SEAK Lab Github\\Coevolution-based Heuristic Incorporation\\results";

        engine = MatlabEngine.startMatlab();

        String matlabScriptsLocation = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS";
        engine.feval("addpath", matlabScriptsLocation); // Add location of MATLAB scripts used to compute objectives, constraints and heuristics to MATLAB's search path

        PRNG.setRandom(new SynchronizedMersenneTwister());

        int numCPU = 4;
        int numRuns = 30;
        pool = Executors.newFixedThreadPool(numCPU);
        ecs = new ExecutorCompletionService<>(pool);

        boolean cMOEA = false; // if the coevolutionary optimization is single or multi objective
        boolean useDE = false; // Use Differential Evolution for outer optimization (otherwise GA is used)
        boolean evolveInternalPopulation = true; // if initial internal population for subsequent heuristic weights should be updated or not
        boolean integerWeights = false; // if formulation for the heuristic weights should be integer or real
        boolean weightOfWeights = false; // if an additional weight multiplicative parameter design decision for the outer GA is used

        /**
         * modelChoice = 0 --> Fibre Stiffness Model
         *             = 1 --> Truss Stiffness Model
         *             = 2 --> Beam Model
         */
        int modelChoice = 1; // Fibre stiffness model cannot be used for the artery problem

        boolean arteryProblem = false; // Solve the artery optimization (otherwise the original truss problem is solved)
        double targetStiffnessRatio = 1;
        if (arteryProblem) {
            targetStiffnessRatio = 0.421;
        }

        // Heuristic Enforcement Methods
        /**
         * partialCollapsibilityConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint, AHS]
         * nodalPropertiesConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint, AHS]
         * orientationConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint, AHS]
         * intersectionConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint, AHS]
         *
         * heuristicsConstrained = [partialCollapsibilityConstrained, nodalPropertiesConstrained, orientationConstrained, intersectionConstrained]
         */
        boolean[] partialCollapsibilityConstrained = {false, false, false, false, false, false, false};
        boolean[] nodalPropertiesConstrained = {false, false, false, false, false, false, false};
        boolean[] orientationConstrained = {false, false, false, false, false, false, false};
        boolean[] intersectionConstrained = {false, false, false, false, false, false, false};

        boolean[][] heuristicsConstrained = new boolean[4][7];
        for (int i = 0; i < 7; i++) {
            heuristicsConstrained[0][i] = partialCollapsibilityConstrained[i];
            heuristicsConstrained[1][i] = nodalPropertiesConstrained[i];
            heuristicsConstrained[2][i] = orientationConstrained[i];
            heuristicsConstrained[3][i] = intersectionConstrained[i];
        }

        // Dimensions for printable solutions
        double printableRadius = 250e-6; // in m
        double printableSideLength = 10e-3; // in m
        double printableModulus = 1.8162e6; // in Pa
        double sideNodeNumber = 3.0D;
        int nucFactor = 3; // Not used if PBC model is used

        String[] constraintNames;
        if (arteryProblem) {
            constraintNames = new String[]{"FeasibilityViolation", "ConnectivityViolation"};
        } else {
            constraintNames = new String[]{"FeasibilityViolation", "ConnectivityViolation", "StiffnessRatioViolation"};
        }
        String[] allHeuristicNames = new String[]{"PartialCollapsibilityViolation","NodalPropertiesViolation","OrientationViolation","IntersectionViolation"};

        boolean[] heuristicIncorporated = new boolean[]{true, true, true, true}; // Decides which heuristics to enforce in the optimization framework, allows for enforcement of selective heuristic combinations
        ArrayList<String> heuristicNames = new ArrayList<>();
        for (int i = 0; i < heuristicIncorporated.length; i++) {
            if (heuristicIncorporated[i]) {
                heuristicNames.add(allHeuristicNames[i]);
            }
        }
        int numberOfHeuristics = heuristicNames.size();

        int totalNumberOfMembers;
        if (sideNodeNumber >= 5) {
            int sidenumSquared = (int) (sideNodeNumber*sideNodeNumber);
            totalNumberOfMembers =  sidenumSquared * (sidenumSquared - 1)/2;
        }
        else {
            totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sideNodeNumber*sideNodeNumber))/(CombinatoricsUtils.factorial((int) ((sideNodeNumber*sideNodeNumber) - 2)) * CombinatoricsUtils.factorial(2)));
        }
        int numberOfRepeatableMembers = (int) (2 * (CombinatoricsUtils.factorial((int) sideNodeNumber)/(CombinatoricsUtils.factorial((int) (sideNodeNumber - 2)) * CombinatoricsUtils.factorial(2))));
        int numVariables = totalNumberOfMembers - numberOfRepeatableMembers;

        int internalPopulationSize = 30;
        int internalMaxEvaluations = 120;

        AbstractProblem metamaterialProblem;
        if (arteryProblem) {
            metamaterialProblem = new ModifiedArteryProblem(integerWeights, numberOfHeuristics, modelChoice, numVariables, printableRadius, printableSideLength, printableModulus, sideNodeNumber, nucFactor, targetStiffnessRatio, engine, heuristicsConstrained, heuristicNames.toArray(new String[0]));
        } else {
            metamaterialProblem = new ModifiedEqualNormalStiffnessProblem(integerWeights, numberOfHeuristics, modelChoice, numVariables, printableRadius, printableSideLength, printableModulus, sideNodeNumber, nucFactor, targetStiffnessRatio, engine, heuristicsConstrained, heuristicNames.toArray(new String[0]));
        }

        String[] variableNames = new String[metamaterialProblem.getNumberOfVariables()];
        String[] objectiveNames = new String[metamaterialProblem.getNumberOfObjectives()]; // Attributes names for the unpenalized objectives as recorded in solutions

        for (int i = 0; i < variableNames.length; i++) {
            variableNames[i] = "Variable" + i;
        }

        for (int i = 0; i < objectiveNames.length; i++) {
            objectiveNames[i] = "TrueObjective" + (i+1);
        }

        int numberOfCoevolutionaryObjectives = 1;
        if (cMOEA) {
            numberOfCoevolutionaryObjectives = metamaterialProblem.getNumberOfObjectives();
        }

        DominanceComparator internalComparator = new ChainedComparator(new AggregateConstraintComparator(), new ParetoObjectiveComparator());

        // Internal Operators
        double crossoverProbability = 1.0;
        double mutationProbability = 1.0/numVariables;
        Variation crossover = new OnePointCrossover(crossoverProbability);
        Variation mutation = new BitFlip(mutationProbability);
        Selection internalSelection = new TournamentSelection(2, internalComparator);
        Variation internalVariation;

        double[] epsilonDouble = new double[]{0.001, 0.001};
        EpsilonBoxDominanceArchive internalArchive;

        // Set random search parameters
        int weightsPopulationSize = 1;
        int weightsMaxEvaluations = 1;

        // Initialize Coevolutionary Problem class
        double[] internalObjectivesMinimum = new double[]{-1, 0}; // Both metamaterial problems involve maximizing the first objective and minimizing the second
        double[] internalObjectivesMaximum = new double[]{0, 1};
        Hypervolume internalHypervolume = new Hypervolume(metamaterialProblem, internalObjectivesMinimum, internalObjectivesMaximum);

        Population initialInternalPopulation;
        AbstractProblem coevolutionaryProblem;

        for (int n = 0; n < numRuns; n++) {
            // Initialize random initial internal population and evaluate (to be used for initial fitness calculation)
            List<Solution> internalInitialSolutions = new ArrayList<>();
            for (int i = 0; i < internalPopulationSize; i++) {
                Solution initialSolution = metamaterialProblem.newSolution();
                for (int j = 0; j < initialSolution.getNumberOfVariables(); j++) {
                    initialSolution.getVariable(j).randomize();
                }
                metamaterialProblem.evaluate(initialSolution);
                for (double objective : initialSolution.getObjectives()) {
                    if (Double.isNaN(objective)) {
                        continue;
                    }
                }
                internalInitialSolutions.add(initialSolution);
            }
            initialInternalPopulation = new Population(internalInitialSolutions);

            Result result;
            if (weightOfWeights) {
                result = new Result(saveDir, numberOfHeuristics+1);
            } else {
                result = new Result(saveDir, numberOfHeuristics);
            }

            internalArchive = new EpsilonBoxDominanceArchive(epsilonDouble);

            internalVariation = new CompoundVariation(crossover, mutation);

            coevolutionaryProblem = new CoevolutionaryProblem(metamaterialProblem, saveDir, weightsPopulationSize, internalMaxEvaluations, initialInternalPopulation, internalArchive, internalComparator, internalSelection, internalVariation, internalHypervolume, cMOEA, evolveInternalPopulation, integerWeights, weightOfWeights, variableNames, objectiveNames, constraintNames, allHeuristicNames, heuristicNames.toArray(new String[0]), numberOfCoevolutionaryObjectives, result, false, false, n);

            List<Solution> initialWeightsSolutions = new ArrayList<>();

            // Initialize solution of all zero weights and add to the initial population
            Solution zeroSolution = coevolutionaryProblem.newSolution();
            double[] zeroSolutionWeights = new double[zeroSolution.getNumberOfVariables()-1];
            if (weightOfWeights) {
                RealVariable zeroWeightVariable = new RealVariable(0.0, 0.0, 1.0);
                zeroSolution.setVariable(0, zeroWeightVariable);
            } else {
                if (!integerWeights) {
                    Arrays.fill(zeroSolutionWeights, -3); // Real weights = 10^(decisions)
                }
                EncodingUtils.setReal(zeroSolution, 0, zeroSolution.getNumberOfVariables()-1, zeroSolutionWeights);
            }

            initialWeightsSolutions.add(zeroSolution);

            // Add other randomly generated solutions
            for (int i = 0; i < weightsPopulationSize-1; i++) {
                initialWeightsSolutions.add(coevolutionaryProblem.newSolution()); // This new solution by default generates a random set of weights with the internal population not to be updated
            }

            ecs.submit(new RandomSearchCallable(coevolutionaryProblem, initialWeightsSolutions, numberOfCoevolutionaryObjectives, weightsMaxEvaluations, evolveInternalPopulation, result, n));
        }

        for (int i = 0; i < numRuns; i++) {
            try {
                ecs.take().get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }

        pool.shutdown();
        engine.close();

    }

}
