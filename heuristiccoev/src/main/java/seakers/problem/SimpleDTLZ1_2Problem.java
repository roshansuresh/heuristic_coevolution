package seakers.problem;

import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;

public class SimpleDTLZ1_2Problem extends AbstractProblem {

    public SimpleDTLZ1_2Problem() {
        super(6, 2); // Number of variables = n + k -1 (n = number of objectives, k = 5 recommended by Deb et al)
    }

    @Override
    public void evaluate(Solution solution) {
        double[] objectives = new double[2];
        double x1 = EncodingUtils.getReal(solution.getVariable(0));
        objectives[0] = x1;
        objectives[1] = (1 - x1);

        solution.setObjectives(objectives);
    }

    @Override
    public Solution newSolution() {
        Solution newSolution = new Solution(numberOfVariables, numberOfObjectives);
        for (int i = 0; i < numberOfVariables; i++) {
            RealVariable var = new RealVariable(PRNG.nextDouble(), 0.0, 1.0);
            newSolution.setVariable(i, var);
        }
        return newSolution;
    }
}
