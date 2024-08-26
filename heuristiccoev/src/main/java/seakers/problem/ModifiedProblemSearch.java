package seakers.problem;

import org.apache.commons.math3.util.Precision;
import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import seakers.Result;
import seakers.problem.eoss.ModifiedAssignmentProblem;
import seakers.problem.eoss.ModifiedPartitioningProblem;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.StringJoiner;
import java.util.concurrent.Callable;

public class ModifiedProblemSearch {

    private final Algorithm algorithm;
    private final String saveDirectory;
    private double[] heuristicWeights;
    private final int internalMaxEvaluations;
    private final int currentCoevolutionNFE;
    private final String[] variableNames;
    private final String[] objectiveNames;
    private final String[] constraintNames;
    private final String[] incorporatedHeuristicNames;
    private final String[] allHeuristicNames;
    private final int coevolutionaryRunNumber;
    private final boolean isPartitioning;
    private final boolean isEOSS;
    private final boolean weightOfWeights;

    public ModifiedProblemSearch(Algorithm algorithm, String saveDirectory, double[] heuristicWeights, int internalMaxEvaluations, int currentCoevolutionNFE, String[] variableNames, String[] objectiveNames, String[] constraintNames, String[] allHeuristicsNames, String[] incorporatedHeuristicNames, boolean isPartitioning, boolean isEOSS, boolean weightOfWeights, int coevolutionaryRunNumber) {
        this.algorithm = algorithm;
        this.saveDirectory = saveDirectory;
        this.heuristicWeights = heuristicWeights;
        this.internalMaxEvaluations = internalMaxEvaluations;
        this.currentCoevolutionNFE = currentCoevolutionNFE;
        this.variableNames = variableNames;
        this.objectiveNames = objectiveNames;
        this.constraintNames = constraintNames;
        this.incorporatedHeuristicNames = incorporatedHeuristicNames;
        this.allHeuristicNames = allHeuristicsNames;
        this.isPartitioning = isPartitioning;
        this.isEOSS = isEOSS;
        this.coevolutionaryRunNumber = coevolutionaryRunNumber;
        this.weightOfWeights = weightOfWeights;
    }

    public Population runMOEA() throws Exception {
        Result result;
        if (weightOfWeights) {
            result = new Result(saveDirectory, incorporatedHeuristicNames.length+1);
        } else {
            result = new Result(saveDirectory, incorporatedHeuristicNames.length);
        }

        // Run MOEA
        long startTime = System.currentTimeMillis();
        algorithm.step();

        ArrayList<Solution> allSolutions = new ArrayList<>();
        Population initialPopulation = ((AbstractEvolutionaryAlgorithm) algorithm).getPopulation();
        for (Solution currentSolution : initialPopulation) {
            currentSolution.setAttribute("NFE",0);
            allSolutions.add(currentSolution);
        }

        while (!algorithm.isTerminated() && (algorithm.getNumberOfEvaluations() < internalMaxEvaluations)) {
            try {
                algorithm.step();
            } catch (Exception e) {
                e.printStackTrace();
            }

            Population currentPopulation = ((AbstractEvolutionaryAlgorithm) algorithm).getPopulation();
            for (int i = 1; i < 3; i++) {
                Solution solution = currentPopulation.get(currentPopulation.size() - i);
                solution.setAttribute("NFE", algorithm.getNumberOfEvaluations());
                allSolutions.add(solution);
            }

            if ((algorithm.getProblem() instanceof ModifiedAssignmentProblem) || (algorithm.getProblem() instanceof ModifiedPartitioningProblem)) {
                if (algorithm.getProblem() instanceof ModifiedAssignmentProblem) {
                    if (algorithm.getNumberOfEvaluations() % 100 == 0) {
                        ((ModifiedAssignmentProblem) algorithm.getProblem()).getEvaluationManager().getResourcePool().poolClean();
                        System.out.println("Internal Rete clean initiated");
                    }
                } else {
                    if (algorithm.getNumberOfEvaluations() % 100 == 0) {
                        ((ModifiedPartitioningProblem) algorithm.getProblem()).getEvaluationManager().getResourcePool().poolClean();
                        System.out.println("Iternal Rete clean initiated");
                    }
                }
            }
        }

        algorithm.terminate();
        long endTime = System.currentTimeMillis();
        System.out.println("Internal Population Evolution Done in " + ((endTime - startTime)/1000.0) + " s");

        // Save solutions to CSV file
        Population finalPopulation = ((AbstractEvolutionaryAlgorithm) algorithm).getPopulation();
        NondominatedPopulation finalArchive =  ((AbstractEvolutionaryAlgorithm) algorithm).getArchive();

        StringJoiner weightsSJ = new StringJoiner("_");
        for (int i = 0; i < incorporatedHeuristicNames.length; i++) {
            weightsSJ.add("w" + i + "-" + Double.toString(Precision.round(heuristicWeights[i], 5)).replace(".",";"));
        }
        if (weightOfWeights) {
            weightsSJ.add("ww-" + Double.toString(Precision.round(heuristicWeights[incorporatedHeuristicNames.length], 5)).replace(".",";"));
        }
        String populationFilename = "run " + coevolutionaryRunNumber + File.separator + "Heuristic_Weights-" + weightsSJ.toString() + "_Coev_nfe-" + currentCoevolutionNFE + "_" + "finalpop" + ".csv";
        String archiveFilename = "run " + coevolutionaryRunNumber + File.separator + "Heuristic_Weights-" + weightsSJ.toString() + "_Coev_nfe-" + currentCoevolutionNFE + "_" + "finalarchive" + ".csv";
        String allSolutionsFilename = "run " + coevolutionaryRunNumber + File.separator + "Heuristic_Weights-" + weightsSJ.toString() + "_Coev_nfe-" + currentCoevolutionNFE + "_allSolutions.csv";

        result.saveAllInternalSolutions(allSolutionsFilename, allSolutions, variableNames, objectiveNames, constraintNames, allHeuristicNames, isPartitioning, isEOSS);
        result.saveInternalPopulationOrArchive(populationFilename, finalPopulation, variableNames, objectiveNames, constraintNames, allHeuristicNames, isPartitioning, isEOSS);
        result.saveInternalPopulationOrArchive(archiveFilename, finalArchive, variableNames, objectiveNames, constraintNames, allHeuristicNames, isPartitioning, isEOSS);

        return ((AbstractEvolutionaryAlgorithm) algorithm).getPopulation();
    }

    public Population getFinalPopulation() {
        // Call after internal population evolution is complete!
        return ((AbstractEvolutionaryAlgorithm) algorithm).getPopulation();
    }


}
