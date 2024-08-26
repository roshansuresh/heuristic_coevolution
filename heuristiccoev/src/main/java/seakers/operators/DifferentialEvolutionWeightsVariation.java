package seakers.operators;

import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.operator.real.DifferentialEvolutionVariation;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;

/**
 * Differential Evolution Variation class for the Weights design decisions.
 * Same as the native DifferentialEvolutionVariation class in the moeaframework library, except the last binary variable
 * (associated with whether or not to update the initial internal population) is copied from any of the parents
 */

public class DifferentialEvolutionWeightsVariation extends DifferentialEvolutionVariation {
    private final double CR;
    private final double F;

    public DifferentialEvolutionWeightsVariation(double CR, double F) {
        super(CR, F);
        this.CR = CR;
        this.F = F;
    }

    @Override
    public Solution[] evolve(Solution[] parents) {
        Solution result = parents[0].copy();

        // Follow DE Variation procedure for the actual weight variables
        int jrand = PRNG.nextInt(result.getNumberOfVariables()-1);

        for(int j = 0; j < result.getNumberOfVariables()-1; ++j) {
            if (PRNG.nextDouble() <= this.CR || j == jrand) {
                RealVariable v0 = (RealVariable)result.getVariable(j);
                RealVariable v1 = (RealVariable)parents[1].getVariable(j);
                RealVariable v2 = (RealVariable)parents[2].getVariable(j);
                RealVariable v3 = (RealVariable)parents[3].getVariable(j);
                double y = v3.getValue() + this.F * (v1.getValue() - v2.getValue());
                if (y < v0.getLowerBound()) {
                    y = v0.getLowerBound();
                }

                if (y > v0.getUpperBound()) {
                    y = v0.getUpperBound();
                }

                v0.setValue(y);
            }
        }

        // If any of the parents have the final boolean variable as true, set to true otherwise set to false
        boolean decision = false;
        for (int i = 1; i < 4; i++) {
            if (EncodingUtils.getBoolean(parents[i].getVariable(result.getNumberOfVariables()-1))) {
                decision = true;
                break;
            }
        }
        EncodingUtils.setBoolean(result.getVariable(result.getNumberOfVariables()-1), decision);

        return new Solution[]{result};
    }
}
