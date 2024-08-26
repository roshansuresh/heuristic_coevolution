package seakers.problem.metamaterial;

import com.mathworks.engine.MatlabEngine;
import org.moeaframework.core.Solution;
import seakers.problem.AbstractInternalProblem;
import seakers.trussaos.problems.ConstantRadiusArteryProblem;

import java.util.Arrays;

public class ModifiedArteryProblem extends ConstantRadiusArteryProblem implements AbstractInternalProblem {

    private double[] heuristicWeights;
    private final boolean integerWeights;
    private final String[] incorporatedHeuristicNames;

    public ModifiedArteryProblem(boolean integerWeights, int numberOfHeuristics, int modelSelection, int numberOfVariables, double memberRadius, double sideElementLength, double modulusYoungs, double sideNodeNumber, double nucFac, double targetCRatio, MatlabEngine eng, boolean[][] heuristicEnforcement, String[] incorporatedHeuristicNames) {
        super("", modelSelection, numberOfVariables, 0, 0, memberRadius, sideElementLength, modulusYoungs, sideNodeNumber, nucFac, targetCRatio, eng, heuristicEnforcement);
        this.integerWeights = integerWeights;
        this.heuristicWeights = new double[numberOfHeuristics];
        this.incorporatedHeuristicNames = incorporatedHeuristicNames;
    }

    @Override
    public void evaluate(Solution solution) {
        // Evaluate design using the parent class
        super.evaluate(solution);
        double norm = 1.0;
        if (integerWeights) { // If real heuristic weights are used, then don't normalize
            norm = 10.0;
            double sum = Arrays.stream(heuristicWeights).sum();
            if (sum == 0.0) {
                norm = 1e-5;
            }
        }

        // Add heuristic penalization to the objectives (the parent class also computes the heuristic violations,
        //String[] heuristics = new String[]{"PartialCollapsibilityViolation","NodalPropertiesViolation","OrientationViolation","IntersectionViolation"};
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

    public void setHeuristicWeights(double[] newWeights) {
        this.heuristicWeights = newWeights;
    }
}
