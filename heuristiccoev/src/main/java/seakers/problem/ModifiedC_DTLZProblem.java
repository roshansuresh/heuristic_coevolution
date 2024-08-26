package seakers.problem;

import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;

import java.util.Arrays;

/**
 * Any Cx_DTLZy problem modified with heuristic penalization
 */

public class ModifiedC_DTLZProblem extends AbstractProblem {

    private AbstractProblem cdltzProblem;
    private double[] heuristicPenaltyWeights;
    private boolean integerWeights;

    public ModifiedC_DTLZProblem(AbstractProblem cdltzProblem, int numberOfHeuristics, boolean integerWeights) {
        super(cdltzProblem.getNumberOfVariables(), cdltzProblem.getNumberOfObjectives(), cdltzProblem.getNumberOfConstraints());
        this.cdltzProblem = cdltzProblem;
        this.heuristicPenaltyWeights = new double[numberOfHeuristics]; // number of heuristics = 2
        this.integerWeights = integerWeights;
    }

    @Override
    public void evaluate(Solution solution) {
        // First evaluate solution using underlying C_DTLZ problem to obtain non-penalized objectives and constraints
        cdltzProblem.evaluate(solution);
        double norm = 1.0;
        if (integerWeights) { // If real heuristic weights are used, then don't normalize
            norm = 10.0;
            double sum = Arrays.stream(heuristicPenaltyWeights).sum();
            if (sum == 0.0) {
                norm = 1e-5;
            }
        }

        // Compute heuristic values (2 inequality heuristics are used: 1 quadratic and 1 linear)
        double quad = -0.5; // h_quad(x) = \sum(x_i^2) - 0.5 <= 0
        double lin = -0.5; // h_lin(x) = \sum(x_i) - 0.5 <= 0
        for (int i = 0; i < solution.getNumberOfVariables(); i++) {
            Variable var = solution.getVariable(i);
            quad += Math.pow(EncodingUtils.getReal(var), 2);
            lin += EncodingUtils.getReal(var);
        }
        quad = Math.max(quad, 0.0)/(solution.getNumberOfVariables() - 0.5);
        lin = Math.max(lin, 0.0)/(solution.getNumberOfVariables() - 0.5);

        // Add heuristic violation penalties to objectives and set modified objectives to solution
        for (int i = 0; i < solution.getNumberOfObjectives(); i++) {
            double penalizedObjective = 0.0;
            if (cdltzProblem instanceof SimpleDTLZ1_2Problem) {
                solution.setAttribute("TrueObjective " + (i+1), solution.getObjective(i)/450.0); // Normalized to enable Hypervolume computation for fitness of the heuristic weights
                penalizedObjective = solution.getObjective(i)/450.0 + (heuristicPenaltyWeights[0]/norm)*quad + (heuristicPenaltyWeights[1]/norm)*lin;
            } else {
                if (i == 0) {
                    solution.setAttribute("TrueObjective " + (i+1), solution.getObjective(i)); // Objectives are normalized for hypervolume computation
                }
                penalizedObjective = solution.getObjective(i) + (heuristicPenaltyWeights[0]/norm)*quad + (heuristicPenaltyWeights[1]/norm)*lin;
            }
            solution.setObjective(i, penalizedObjective);

        }
        solution.setAttribute("Heuristic 1", quad);
        solution.setAttribute("Heuristic 2", lin);
    }

    @Override
    public Solution newSolution() {
        return cdltzProblem.newSolution();
    }

    public void setHeuristicWeights(double[] newWeights) {
        this.heuristicPenaltyWeights = newWeights;
    }
}
