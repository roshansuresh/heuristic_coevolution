package seakers;

import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.algorithm.single.AggregateObjectiveComparator;
import org.moeaframework.algorithm.single.DifferentialEvolution;
import org.moeaframework.algorithm.single.GeneticAlgorithm;
import org.moeaframework.algorithm.single.LinearObjectiveComparator;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.DominanceComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.indicator.Hypervolume;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.core.operator.real.DifferentialEvolutionSelection;
import org.moeaframework.core.operator.real.DifferentialEvolutionVariation;
import org.moeaframework.core.operator.real.PM;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;
import seakers.operators.DifferentialEvolutionWeightsVariation;
import seakers.operators.SimulatedBinaryWeightsCrossover;
import seakers.operators.UniformIntegerWeightsMutation;
import seakers.operators.UniformWeightsCrossover;
import seakers.problem.CoevolutionaryProblem;
import seakers.problem.eoss.ModifiedAssignmentProblem;
import seakers.problem.eoss.ModifiedPartitioningProblem;
import seakers.vassarexecheur.search.intialization.SynchronizedMersenneTwister;
import seakers.vassarexecheur.search.intialization.partitioning.RandomPartitioningReadInitialization;
import seakers.vassarexecheur.search.operators.partitioning.PartitioningCrossover;
import seakers.vassarexecheur.search.operators.partitioning.PartitioningMutation;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.evaluation.AbstractArchitectureEvaluator;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;

import java.util.*;
import java.util.concurrent.*;

public class EOSSCoevolutionaryMOEARun {

    public static ExecutorService pool;
    public static CompletionService<Algorithm> ecs;

    public static void main(String[] args) {
        // Save location
        //String saveDir = System.getProperty("user.dir") + File.separator + "results"; // File not found error!
        //String saveDir = "C:\\SEAK Lab\\SEAK Lab Github\\Coevolution-based Heuristic Incorporation\\results";
        String saveDir = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\Heuristic Coevolution\\results";

        boolean cMOEA = false; // if the coevolutionary optimization is single or multi objective
        boolean useDE = false; // Use Differential Evolution for outer optimization (otherwise GA is used)
        boolean evolveInternalPopulation = true; // if initial internal population for subsequent heuristic weights should be updated or not
        boolean integerWeights = false; // if formulation for the heuristic weights should be integer or real
        boolean weightOfWeights = false; // if an additional weight multiplicative parameter design decision for the outer GA is used
        boolean periodicZeroInjection = false; // if zero solution is added for evaluation to population at different NFE

        boolean assigningProblem = false; // True -> assigning problem, False -> partitioning problem

        // Heuristic Enforcement Methods
        /**
         * dutyCycleConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * instrumentOrbitRelationsConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * interferenceConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * packingEfficiencyConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * spacecraftMassConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * synergyConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * instrumentCountConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         *
         * if partitioning problem:
         * heuristicsConstrained = [dutyCycleConstrained, instrumentOrbitRelationsConstrained, interferenceConstrained, packingEfficiencyConstrained, spacecraftMassConstrained, synergyConstrained]
         * else:
         * heuristicsConstrained = [dutyCycleConstrained, instrumentOrbitRelationsConstrained, interferenceConstrained, packingEfficiencyConstrained, spacecraftMassConstrained, synergyConstrained, instrumentCountConstrained]
         */
        boolean[] dutyCycleConstrained = {false, false, false, false, false, false};
        boolean[] instrumentOrbitRelationsConstrained = {false, false, false, false, false, false};
        boolean[] interferenceConstrained = {false, false, false, false, false, false};
        boolean[] packingEfficiencyConstrained = {false, false, false, false, false, false};
        boolean[] spacecraftMassConstrained = {false, false, false, false, false, false};
        boolean[] synergyConstrained = {false, false, false, false, false, false};
        boolean[] instrumentCountConstrained = {false, false, false, false, false, false}; // only for assigning problem

        boolean[][] heuristicsConstrained;
        if (assigningProblem) {
            heuristicsConstrained = new boolean[7][6];
        } else {
            heuristicsConstrained = new boolean[6][6];
        }

        for (int i = 0; i < heuristicsConstrained[0].length; i++) {
            heuristicsConstrained[0][i] = dutyCycleConstrained[i];
            heuristicsConstrained[1][i] = instrumentOrbitRelationsConstrained[i];
            heuristicsConstrained[2][i] = interferenceConstrained[i];
            heuristicsConstrained[3][i] = packingEfficiencyConstrained[i];
            heuristicsConstrained[4][i] = spacecraftMassConstrained[i];
            heuristicsConstrained[5][i] = synergyConstrained[i];
            if (assigningProblem) {
                heuristicsConstrained[6][i] = instrumentCountConstrained[i];
            }
        }

        String[] constraintNames = new String[]{};
        String[] allHeuristicNames;
        if (assigningProblem) {
            allHeuristicNames = new String[]{"DCViolation","InstrOrbViolation","InterInstrViolation","PackEffViolation","SpMassViolation","SynergyViolation","InstrCountViolation"};
        } else {
            allHeuristicNames = new String[]{"DCViolation","InstrOrbViolation","InterInstrViolation","PackEffViolation","SpMassViolation","SynergyViolation"};
        }

        boolean[] heuristicIncorporated = new boolean[]{true, true, true, true, true, true, true}; // Decides which heuristics to enforce in the optimization framework, allows for enforcement of selective heuristic combinations
        ArrayList<String> heuristicNames = new ArrayList<>();
        for (int i = 0; i < allHeuristicNames.length; i++) {
            if (heuristicIncorporated[i]) {
                heuristicNames.add(allHeuristicNames[i]);
            }
        }
        int numberOfHeuristics = heuristicNames.size();

        int numCPU = 4;
        int numRuns = 30; // comment if line 126 is uncommented
        pool = Executors.newFixedThreadPool(numCPU);
        ecs = new ExecutorCompletionService<>(pool);

        // Next two lines for Partitioning problem (or any problem where runs need to be initialized with specific populations), populate with specific run numbers as needed
        //int[] runNumbers = {8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29}; // change numCPU accordingly
        //int numRuns = runNumbers.length; // If previous line is uncommented, uncomment this and line 219 and comment line 218

        //setup for epsilon MOEA
        double[] epsilonDouble = new double[]{0.0001, 0.0001};

        double crossoverProbability = 1.0;

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.7;
        double instrCountThreshold = 15; // only for assigning problem
        boolean considerFeasibility = true;

        int internalPopulationSize = 72;
        int internalMaxEvaluations = 216;

        //String resourcesPath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur"; // for lab system
        String resourcesPath = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur"; // for laptop

        double mutationProbability;
        BaseParams params;
        AbstractArchitectureEvaluator evaluator;
        HashMap<String, String[]> instrumentSynergyMap;
        HashMap<String, String[]> interferingInstrumentsMap;
        if (assigningProblem) {
            mutationProbability = 1. / 60.;
            params = new ClimateCentricAssigningParams(resourcesPath, "FUZZY-ATTRIBUTES","test", "normal");

            instrumentSynergyMap = getInstrumentSynergyNameMap(params);
            interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

            evaluator = new ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
        } else {
            mutationProbability = 1. / 24.; // Based on the 12 instruments for the ClimateCentric Problem
            params = new ClimateCentricPartitioningParams(resourcesPath, "FUZZY-ATTRIBUTES", "test", "normal");

            instrumentSynergyMap = getInstrumentSynergyNameMap(params);
            interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

            evaluator = new seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
        }

        Initialization internalInitialization;
        CompoundVariation internalVariation;

        ArchitectureEvaluationManager evaluationManager = null;

        // Problem class
        AbstractProblem satelliteProblem;

        EpsilonBoxDominanceArchive internalArchive;

        DominanceComparator internalComparator = new ParetoObjectiveComparator();
        TournamentSelection internalSelection = new TournamentSelection(2, internalComparator);

        PRNG.setRandom(new SynchronizedMersenneTwister());

        // Initialize MOEA
        int coevolutionPopulationSize = 5;
        int coevolutionMaxEvaluations = 50;

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

        DominanceComparator coevolutionComparator;
        if (cMOEA) {
            coevolutionComparator = new ParetoObjectiveComparator();
        } else {
            coevolutionComparator = new LinearObjectiveComparator();
        }

        Variation coevolutionVariation;
        Selection coevolutionSelection = new TournamentSelection(2, coevolutionComparator);

        int numberOfCoevolutionaryObjectives = 1;

        Population initialInternalPopulation;
        AbstractProblem coevolutionaryProblem;
        Initialization coevolutionInitialization;
        Algorithm coevolutionMOEA;

        for (int n = 0; n < numRuns; n++) {
            int currentRunNumber = n;
            //int currentRunNumber = runNumbers[n];

            evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
            evaluationManager.init(1);

            if (assigningProblem) {
                satelliteProblem = new ModifiedAssignmentProblem(integerWeights, numberOfHeuristics, new int[]{1}, params.getProblemName(), evaluationManager, (ArchitectureEvaluator) evaluator, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, instrCountThreshold, heuristicsConstrained, heuristicNames.toArray(new String[0]));
            } else {
                satelliteProblem = new ModifiedPartitioningProblem(integerWeights, numberOfHeuristics, params.getProblemName(), evaluationManager, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, heuristicsConstrained, heuristicNames.toArray(new String[0]));
            }
            if (cMOEA) {
                numberOfCoevolutionaryObjectives = satelliteProblem.getNumberOfObjectives();
            }

            // Initialize variable and objectives
            String[] variableNames;
            if (assigningProblem) {
                variableNames = getVariableNames(satelliteProblem.getNumberOfVariables() - 1); // Not counting the initial integer variable set in the AssigningArchitecture (which does not contribute to the architecture evaluation but present for legacy reasons)
            } else {
                variableNames = getVariableNames(satelliteProblem.getNumberOfVariables());
            }

            String[] objectiveNames = getObjectiveNames(satelliteProblem.getNumberOfObjectives());

            // Initialize Coevolutionary Problem class
            Hypervolume internalHypervolume = getHypervolume(satelliteProblem);

            double[] epsilonBoxDouble = new double[satelliteProblem.getNumberOfObjectives()];
            Arrays.fill(epsilonBoxDouble, 0.0001);
            EpsilonBoxDominanceArchive coevolutionArchive = new EpsilonBoxDominanceArchive(epsilonBoxDouble);

            // Initialize initial internal population and evaluate (to be used for initial fitness calculation)
            if (assigningProblem) {
                internalInitialization = new RandomInitialization(satelliteProblem, internalPopulationSize);
            } else {
                //initialization = new RandomPartitioningAndAssigning(popSize, (PartitioningProblem) satelliteProblem, params.getInstrumentList(), params.getOrbitList());
                //initialization = new RandomFeasiblePartitioning(popSize, (PartitioningProblem) satelliteProblem, params.getInstrumentList(), params.getOrbitList());

                internalInitialization = new RandomPartitioningReadInitialization(saveDir, currentRunNumber, internalPopulationSize, (PartitioningProblem) satelliteProblem, params.getInstrumentList(), params.getOrbitList());
            }
            Solution[] initialInternalSolutions = internalInitialization.initialize();

            initialInternalPopulation = new Population();
            for (int i = 0; i < internalPopulationSize; i++) {
                Solution initialSolution = initialInternalSolutions[i];
                satelliteProblem.evaluate(initialSolution);
                initialInternalPopulation.add(initialSolution);
            }

            Result result;
            if (weightOfWeights) {
                result = new Result(saveDir, numberOfHeuristics+1);
            } else {
                result = new Result(saveDir, numberOfHeuristics);
            }

            internalArchive = new EpsilonBoxDominanceArchive(epsilonDouble);

            if (assigningProblem) {
                internalVariation = new CompoundVariation(new OnePointCrossover(crossoverProbability), new BitFlip(mutationProbability));
            } else {
                internalVariation = new CompoundVariation(new PartitioningCrossover(crossoverProbability, params), new PartitioningMutation(mutationProbability, params));
            }

            coevolutionaryProblem = new CoevolutionaryProblem(satelliteProblem, saveDir, coevolutionPopulationSize, internalMaxEvaluations, initialInternalPopulation, internalArchive, internalComparator, internalSelection, internalVariation, internalHypervolume, cMOEA, evolveInternalPopulation, integerWeights, weightOfWeights, variableNames, objectiveNames, constraintNames, allHeuristicNames, heuristicNames.toArray(new String[0]), numberOfCoevolutionaryObjectives, result, !assigningProblem, true, currentRunNumber);

            // NOTE: Injected initialization is used since RandomInitialization will randomize all variable in a new solution (including the last binary variable
            // signifying whether to update initial internal population or not, which we don't want. The initial population of weights must use the same initial internal population
            List<Solution> initialCoevolutionSolutions = new ArrayList<>();

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

            initialCoevolutionSolutions.add(zeroSolution);

            // Add other randomly generated solutions
            for (int i = 0; i < coevolutionPopulationSize-1; i++) {
                initialCoevolutionSolutions.add(coevolutionaryProblem.newSolution()); // This new solution by default generates a random set of weights with the internal population not to be updated
            }
            coevolutionInitialization = new InjectedInitialization(coevolutionaryProblem, coevolutionPopulationSize, initialCoevolutionSolutions);

            coevolutionVariation = new CompoundVariation(coevolutionaryCrossover, coevolutionaryMutation);

            if (cMOEA) {
                coevolutionMOEA = new EpsilonMOEA(coevolutionaryProblem, coevolutionPopulation, coevolutionArchive, coevolutionSelection, coevolutionVariation, coevolutionInitialization, coevolutionComparator);
            } else {
                if (useDE) {
                    AggregateObjectiveComparator coevolutionDEComparator = new LinearObjectiveComparator();
                    DifferentialEvolutionVariation coevolutionDEVariation = new DifferentialEvolutionWeightsVariation(1.0, 0.5);
                    DifferentialEvolutionSelection coevolutionDESelection = new DifferentialEvolutionSelection();
                    coevolutionMOEA = new DifferentialEvolution(coevolutionaryProblem, coevolutionDEComparator, coevolutionInitialization, coevolutionDESelection, coevolutionDEVariation);
                } else {
                    coevolutionMOEA = new GeneticAlgorithm(coevolutionaryProblem, (AggregateObjectiveComparator) coevolutionComparator, coevolutionInitialization, coevolutionSelection, coevolutionVariation);
                }
            }
            ecs.submit(new CoevolutionaryTestSearch(coevolutionMOEA, coevolutionMaxEvaluations, satelliteProblem.getNumberOfObjectives(), result, coevolutionComparator, cMOEA, true, (!assigningProblem), evolveInternalPopulation, integerWeights, weightOfWeights, periodicZeroInjection, currentRunNumber));
        }

        for (int i = 0; i < numRuns; i++) {
            try {
                ecs.take().get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }
        //evaluationManager.clear(); // NullPointerException encountered
        pool.shutdown();
    }

    /**
     * Creates instrument synergy map used to compute the instrument synergy violation heuristic (only formulated for the
     * Climate Centric problem for now) (added by roshansuresh)
     * @param params
     * @return Instrument synergy hashmap
     */
    protected static HashMap<String, String[]> getInstrumentSynergyNameMap(BaseParams params) {
        HashMap<String, String[]> synergyNameMap = new HashMap<>();
        if (params.getProblemName().equalsIgnoreCase("ClimateCentric")) {
            synergyNameMap.put("ACE_ORCA", new String[]{"DESD_LID", "GACM_VIS", "ACE_POL", "HYSP_TIR", "ACE_LID"});
            synergyNameMap.put("DESD_LID", new String[]{"ACE_ORCA", "ACE_LID", "ACE_POL"});
            synergyNameMap.put("GACM_VIS", new String[]{"ACE_ORCA", "ACE_LID"});
            synergyNameMap.put("HYSP_TIR", new String[]{"ACE_ORCA", "POSTEPS_IRS"});
            synergyNameMap.put("ACE_POL", new String[]{"ACE_ORCA", "DESD_LID"});
            synergyNameMap.put("ACE_LID", new String[]{"ACE_ORCA", "CNES_KaRIN", "DESD_LID", "GACM_VIS"});
            synergyNameMap.put("POSTEPS_IRS", new String[]{"HYSP_TIR"});
            synergyNameMap.put("CNES_KaRIN", new String[]{"ACE_LID"});
        }
        else {
            System.out.println("Synergy Map for current problem not formulated");
        }
        return synergyNameMap;
    }

    /**
     * Creates instrument interference map used to compute the instrument interference violation heuristic (only formulated for the
     * Climate Centric problem for now)
     * @param params
     * @return Instrument interference hashmap
     */
    protected static HashMap<String, String[]> getInstrumentInterferenceNameMap(BaseParams params) {
        HashMap<String, String[]> interferenceNameMap = new HashMap<>();
        if (params.getProblemName().equalsIgnoreCase("ClimateCentric")) {
            interferenceNameMap.put("ACE_LID", new String[]{"ACE_CPR", "DESD_SAR", "CLAR_ERB", "GACM_SWIR"});
            interferenceNameMap.put("ACE_CPR", new String[]{"ACE_LID", "DESD_SAR", "CNES_KaRIN", "CLAR_ERB", "ACE_POL", "ACE_ORCA", "GACM_SWIR"});
            interferenceNameMap.put("DESD_SAR", new String[]{"ACE_LID", "ACE_CPR"});
            interferenceNameMap.put("CLAR_ERB", new String[]{"ACE_LID", "ACE_CPR"});
            interferenceNameMap.put("CNES_KaRIN", new String[]{"ACE_CPR"});
            interferenceNameMap.put("ACE_POL", new String[]{"ACE_CPR"});
            interferenceNameMap.put("ACE_ORCA", new String[]{"ACE_CPR"});
            interferenceNameMap.put("GACM_SWIR", new String[]{"ACE_LID", "ACE_CPR"});
        }
        else {
            System.out.println("Interference Map fpr current problem not formulated");
        }
        return interferenceNameMap;
    }

    static String[] getVariableNames(int numberOfVariables) {
        String[] variableNames = new String[numberOfVariables];

        for (int i = 0; i < variableNames.length; i++) {
            variableNames[i] = "Variable" + i;
        }
        return variableNames;
    }

    static String[] getObjectiveNames(int numberOfObjectives) {
        String[] objectiveNames = new String[numberOfObjectives]; // Attributes names for the unpenalized objectives as recorded in solutions

        for (int i = 0; i < objectiveNames.length; i++) {
            objectiveNames[i] = "TrueObjective" + (i + 1);
        }

        return objectiveNames;
    }

    private static Hypervolume getHypervolume(AbstractProblem satelliteProblem) {
        double[] internalObjectivesMinimum = new double[]{-1, 0}; // Both satellite problems involve maximizing the first objective and minimizing the secon
        double[] internalObjectivesMaximum = new double[]{0, 1};
        Hypervolume internalHypervolume = new Hypervolume(satelliteProblem, internalObjectivesMinimum, internalObjectivesMaximum);

        return internalHypervolume;
    }
}
