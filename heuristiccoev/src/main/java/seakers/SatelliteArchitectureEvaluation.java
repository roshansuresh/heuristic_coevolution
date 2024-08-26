package seakers;

import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.problem.AbstractProblem;
import seakers.architecture.util.IntegerVariable;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.evaluation.AbstractArchitectureEvaluator;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;
import seakers.vassarheur.problems.Assigning.ClimateCentricAssigningParams;
import seakers.vassarheur.problems.PartitioningAndAssigning.ClimateCentricPartitioningParams;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;

public class SatelliteArchitectureEvaluation {

    public static void main(String[] args) throws FileNotFoundException {

        boolean assigningProblem = true; // True -> assigning problem, False -> partitioning problem

        String saveDir = "C:\\SEAK Lab\\SEAK Lab Github\\Coevolution-based Heuristic Incorporation\\results";

        String resourcesPath = "C:\\SEAK Lab\\SEAK Lab Github\\VASSAR\\VASSAR_resources-heur"; // for lab system

        double dcThreshold = 0.5;
        double massThreshold = 3000.0; // [kg]
        double packEffThreshold = 0.7;
        double instrCountThreshold = 15; // only for assigning problem
        boolean considerFeasibility = true;

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

        int numberOfHeuristicConstraints = 0;
        int numberOfHeuristicObjectives = 0;
        for (int i = 0; i < heuristicsConstrained.length; i++) {
            if (heuristicsConstrained[i][5]) {
                numberOfHeuristicConstraints++;
            }
            if (heuristicsConstrained[i][4]) {
                numberOfHeuristicObjectives++;
            }
        }

        // Problem class
        AbstractProblem satelliteProblem;

        BaseParams params;
        AbstractArchitectureEvaluator evaluator;
        HashMap<String, String[]> instrumentSynergyMap;
        HashMap<String, String[]> interferingInstrumentsMap;

        if (assigningProblem) {
            params = new ClimateCentricAssigningParams(resourcesPath, "FUZZY-ATTRIBUTES","test", "normal");

            instrumentSynergyMap = getInstrumentSynergyNameMap(params);
            interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

            evaluator = new ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
        } else {
            params = new ClimateCentricPartitioningParams(resourcesPath, "FUZZY-ATTRIBUTES", "test", "normal");

            instrumentSynergyMap = getInstrumentSynergyNameMap(params);
            interferingInstrumentsMap = getInstrumentInterferenceNameMap(params);

            evaluator = new seakers.vassarheur.problems.PartitioningAndAssigning.ArchitectureEvaluator(considerFeasibility, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold);
        }

        ArchitectureEvaluationManager evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(1);

        if (assigningProblem) {
            satelliteProblem = new AssigningProblem(new int[]{1}, params.getProblemName(), evaluationManager, (ArchitectureEvaluator) evaluator, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, instrCountThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        } else {
            satelliteProblem = new PartitioningProblem(params.getProblemName(), evaluationManager, params, interferingInstrumentsMap, instrumentSynergyMap, dcThreshold, massThreshold, packEffThreshold, numberOfHeuristicObjectives, numberOfHeuristicConstraints, heuristicsConstrained);
        }

        // Create architecture to evaluate
        if (assigningProblem) {
            int[] architectureDecisions = {1,0,0,1,0,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,1,0,1,0,0,1,0,0,1,0,1,0,1,1,0,0,1,0,1,1,0,0,0,0,0};

            AssigningArchitecture architecture = (AssigningArchitecture) satelliteProblem.newSolution();
            for (int i = 1; i < architecture.getNumberOfVariables(); i++) {
                BinaryVariable var = new BinaryVariable(1);
                EncodingUtils.setBoolean(var, (architectureDecisions[i-1] == 1));
                architecture.setVariable(i, var);
            }
            satelliteProblem.evaluate(architecture);
        } else {
            int[] architectureDecisions = {0,0,0,1,1,2,3,3,4,4,5,6,1,1,3,2,0,4,2,-1,-1,-1,-1,-1};

            PartitioningArchitecture architecture = (PartitioningArchitecture) satelliteProblem.newSolution();
            for (int i = 0; i < architecture.getNumberOfVariables(); i++) {
                if (i < params.getNumInstr()) {
                    IntegerVariable var = new IntegerVariable(architectureDecisions[i], 0, params.getNumInstr()-1);
                    architecture.setVariable(i, var);
                }
                if (i >= params.getNumInstr()) {
                    IntegerVariable var = new IntegerVariable(architectureDecisions[i], -1, params.getNumOrbits()-1);
                    architecture.setVariable(i, var);
                }
            }
            satelliteProblem.evaluate(architecture);
        }

        evaluationManager.clear();
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
}
