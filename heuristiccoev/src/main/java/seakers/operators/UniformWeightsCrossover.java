package seakers.operators;

import jmetal.encodings.variable.Binary;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.Variation;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;

public class UniformWeightsCrossover implements Variation {

    private final double probability;

    public UniformWeightsCrossover(double probability) {
        this.probability = probability;
    }

    @Override
    public int getArity() {return 2;}

    @Override
    public Solution[] evolve(Solution[] parents) {
        Solution child1 = parents[0].copy();
        Solution child2 = parents[1].copy();
        if (PRNG.nextDouble() <= probability) {
            for (int i = 0; i < child1.getNumberOfVariables()-1; ++i) { // only the decisions corresponding to the heuristic weights
                if (PRNG.nextBoolean()) {
                    Variable temp = child1.getVariable(i);
                    child1.setVariable(i, child2.getVariable(i));
                    child2.setVariable(i, temp);
                }
            }
        }
        //BinaryVariable var = new BinaryVariable(1);
        //EncodingUtils.setBoolean(var, true);
        //BinaryVariable var2 = var.copy();
        //child1.setVariable(child1.getNumberOfVariables(), var);
        //child2.setVariable(child2.getNumberOfVariables(), var2);

        return new Solution[]{child1, child2};
    }

}
