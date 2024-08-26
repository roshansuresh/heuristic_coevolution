package seakers.problem.eoss;

import org.moeaframework.core.Solution;
import seakers.problem.AbstractInternalProblem;
import seakers.vassarexecheur.search.problems.assigning.AssigningArchitecture;
import seakers.vassarexecheur.search.problems.assigning.AssigningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;
import seakers.vassarheur.problems.Assigning.ArchitectureEvaluator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class ModifiedAssignmentProblem extends AssigningProblem implements AbstractInternalProblem {

    private double[] heuristicWeights; // Heuristic weights are the size of only the incorporated heuristics
    private final boolean integerWeights;
    private final String[] incorporatedHeuristicNames;

    public ModifiedAssignmentProblem(boolean integerWeights, int numberOfHeuristics, int[] alternativesForNumberOfSatellites, String eossProblem, ArchitectureEvaluationManager evaluationManager, ArchitectureEvaluator evaluator, BaseParams params, HashMap<String, String[]> interferenceMap, HashMap<String, String[]> synergyMap, double dcThreshold, double massThreshold, double packingEfficiencyThreshold, double instrumentCountThreshold, boolean[][] heuristicsEnforced, String[] incorporatedHeuristicNames) {
        super(alternativesForNumberOfSatellites, eossProblem, evaluationManager, evaluator, params, interferenceMap, synergyMap, dcThreshold, massThreshold, packingEfficiencyThreshold, instrumentCountThreshold, 0, 0, heuristicsEnforced);
        this.heuristicWeights = new double[numberOfHeuristics];
        this.integerWeights = integerWeights;
        this.incorporatedHeuristicNames = incorporatedHeuristicNames;
    }

    @Override
    public void evaluate(Solution solution) {
        // Evaluate design using the parent class
        super.evaluate((AssigningArchitecture) solution);
        double norm = 1.0;
        if (integerWeights) { // If real heuristic weights are used, then don't normalize
            norm = 10.0;
            double sum = Arrays.stream(heuristicWeights).sum();
            if (sum == 0.0) {
                norm = 1e-5;
            }
        }

        // Add heuristic penalization to the objectives (the parent class also computes the heuristic violations)
        //String[] heuristics = new String[]{"DCViolation","InstrOrbViolation","InterInstrViolation","PackEffViolation","SpMassViolation","SynergyViolation"};
        double totalHeuristicPenalty = 0.0;
        for (int i = 0; i < incorporatedHeuristicNames.length; i++) {
            double heuristicValue =  (double) solution.getAttribute(incorporatedHeuristicNames[i]);
            totalHeuristicPenalty += (heuristicWeights[i]/norm)*(heuristicValue/incorporatedHeuristicNames.length);
        }

        for (int i = 0; i < solution.getNumberOfObjectives(); i++) {
            double penalizedObjective = solution.getObjective(i) + totalHeuristicPenalty;
            solution.setObjective(i, penalizedObjective);
        }
    }

    @Override
    public void setHeuristicWeights(double[] heuristicWeights) {
        this.heuristicWeights = heuristicWeights;
    }

}
