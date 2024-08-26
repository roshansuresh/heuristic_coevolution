package seakers.problem.eoss;

import org.moeaframework.core.Solution;
import seakers.problem.AbstractInternalProblem;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningArchitecture;
import seakers.vassarexecheur.search.problems.partitioning.PartitioningProblem;
import seakers.vassarheur.BaseParams;
import seakers.vassarheur.evaluation.ArchitectureEvaluationManager;

import java.util.Arrays;
import java.util.HashMap;

public class ModifiedPartitioningProblem extends PartitioningProblem implements AbstractInternalProblem {

    private double[] heuristicWeights;
    private final boolean integerWeights;
    private final String[] incorporatedHeuristicNames;

    public ModifiedPartitioningProblem(boolean integerWeights, int numberOfHeuristics, String eossProblem, ArchitectureEvaluationManager evaluationManager, BaseParams params, HashMap<String, String[]> interferenceMap, HashMap<String, String[]> synergyMap, double dcThreshold, double massThreshold, double packingEfficiencyThreshold, boolean[][] heuristicsConstrained, String[] incorporatedHeuristicNames) {
        super(eossProblem, evaluationManager, params, interferenceMap, synergyMap, dcThreshold, massThreshold, packingEfficiencyThreshold, 0, 0, heuristicsConstrained);
        this.heuristicWeights = new double[numberOfHeuristics];
        this.integerWeights = integerWeights;
        this.incorporatedHeuristicNames = incorporatedHeuristicNames;
    }

    @Override
    public void evaluate(Solution solution) {
        // Evaluate design using the parent class
        super.evaluate((PartitioningArchitecture) solution);
        double norm = 1.0;
        if (integerWeights) { // If real heuristic weights are used, then don't normalize
            norm = 10.0;
            double sum = Arrays.stream(heuristicWeights).sum();
            if (sum == 0.0) {
                norm = 1e-5;
            }
        }

        // Add heuristic penalization to the objectives (the parent class also computes the heuristic violations,
        //String[] heuristics = new String[]{"DCViolation","InstrOrbViolation","InterInstrViolation","PackEffViolation","SpMassViolation","SynergyViolation"};
        double totalHeuristicPenalty = 0.0;
        for (int i = 0; i < incorporatedHeuristicNames.length; i++) {
            double heuristicValue = (double) solution.getAttribute(incorporatedHeuristicNames[i]);
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
