package OnLattice.utils;

import OnLattice.base.*;
import HAL.Gui.GifMaker;
import HAL.Gui.GridWindow;
import java.awt.GraphicsEnvironment;
import java.io.IOException;
import com.opencsv.CSVWriter;
import java.util.ArrayList;
import java.util.List;
import static OnLattice.base.CellType.celltype.TUMOR;
import static HAL.Util.*;
import static OnLattice.base.ExampleGrid.GAP_JUNCTION_MODULATION;
import static OnLattice.base.ModulationMapExporter.ExportTumorCellsModulation;
import static OnLattice.utils.SimulationManager.AnalysisType.CHEMO;

/**
 * The SimulationRunner class is responsible for executing the simulation loop
 * for a set of grids and simulation parameters. It advances the simulation timesteps,
 * updates cell behaviors, manages visualization (if available), and collects output
 * metrics for later analysis.
 */
public class SimulationRunner {

    /**
     * Runs the simulation for the provided grids and parameters.
     *
     * @param grids             Array of simulation grids.
     * @param effectAntiMet     Treatment effect parameter for anti-metastatic cells.
     * @param effectProMet      Treatment effect parameter for pro-metastatic cells.
     * @param conversionThreshold Threshold for cell-type conversion.
     * @param S2                Simulation parameter S2 (switchSensitivity) 
     * @param S4                Simulation parameter S4 (divisionSensitivity)
     * @param effectPerTumorCell Effect per tumor cell.
     * @param DRUG_DEATH        Drug-induced cell death rate.
     * @param DRUG_UPTAKE       Drug uptake rate.
     * @param timesteps         Total number of timesteps for the simulation.
     * @param sampleIndex       Index corresponding to the parameter set sample.
     * @param replicationIndex  Replication number for the simulation run.
     * @param writer            CSVWriter for output logging.
     * @param gifMaker          GifMaker for generating a visualization output.
     * @param win               GridWindow for live visualization (if not headless).
     * @throws IOException      When file operations fail.
     */
    public static void runSimulation(ExampleGrid[] grids,
                                     double effectAntiMet,
                                     double effectProMet,
                                     double conversionThreshold,
                                     double S2, //switchSensitivity
                                     double S4, //divisionSensitivity
                                     double effectPerTumorCell,
                                     double DRUG_DEATH,
                                     double DRUG_UPTAKE,
                                     int timesteps,
                                     int sampleIndex,
                                     int replicationIndex,
                                     CSVWriter writer,
                                     GifMaker gifMaker,
                                     GridWindow win) throws IOException {

        // Initialize simulation timing and data batch collection.
        int localTimeStep = 0;
        List<String[]> batchData = new ArrayList<>();
        boolean isHeadless = GraphicsEnvironment.isHeadless();
        int numGrids = grids.length;

        // Loop over simulation timesteps.
        for (int t = 0; t < timesteps; t++) {
            // Iterate through each grid.
            for (int i = 0; i < numGrids; i++) {
                ExampleGrid currentGrid = grids[i];

                // Cache grid dimensions locally.
                int gridX = currentGrid.xDim;
                int gridY = currentGrid.yDim;

                // Early timesteps: ensure a tumor cell exists at the grid center.
                if (t < 5 && currentGrid.GetPop() < 1 && currentGrid.IsEmpty(gridX / 2, gridY / 2)) {
                    Cell centralCell = currentGrid.NewAgentSQ(gridX / 2, gridY / 2);
                    centralCell.Init(TUMOR);
                }

                // Update cell behaviors by advancing the simulation for this timestep.
                currentGrid.StepCells(conversionThreshold, effectPerTumorCell, effectProMet, effectAntiMet,
                        DRUG_DEATH, localTimeStep, S2, S4);

                // Visualization step: only if not running in headless mode.
                if (!isHeadless) {
                    // Determine layout dimensions based on number of grids.
                    int rows, cols;
                    if (numGrids == 6) {
                        rows = 2;
                        cols = 3;
                    } else if (numGrids == 2) {
                        rows = 1;
                        cols = 2;
                    } else {
                        rows = 1;
                        cols = 1;
                    }
                    // Calculate the drawing offset based on the grid's index.
                    int offsetX = (i % cols) * gridX;
                    int offsetY = (i / cols) * gridY;
                    currentGrid.DrawModel(win, offsetX, offsetY);
                    gifMaker.AddFrame(win);
                }

                // Check whether chemotherapy was applied at the current timestep.
                String chemoApplied = currentGrid.chemoApplicationTimesteps.contains(localTimeStep) ? "1" : "0";

                // Optional (only relevant for chemo analysis) - At a specific timestep and grid, export the modulation map.
                if (t == 149 && i == 1) {
                    ExportTumorCellsModulation(grids[1], t, "tumorMod_" + GAP_JUNCTION_MODULATION + ".csv");
                    }

                // Collect output metrics at the final timestep; use placeholder values otherwise.
                if (t >= timesteps -1) { // change this if you want to collect morphological metrics at every timestep.
                    Outputs.isFrontCell(currentGrid);
                    double fractalDimension = Outputs.calculateFractalDimension(Outputs.frontCellX, Outputs.frontCellY, currentGrid);
                    double lacunarity = Outputs.calculateEdgeLacunarity(currentGrid, Outputs.frontCellX, Outputs.frontCellY);
                    double eccentricity = Outputs.calculateBulkEccentricity(currentGrid);
                    Outputs.frontCellX.clear();
                    Outputs.frontCellY.clear();

                    batchData.add(new String[]{
                            Integer.toString(sampleIndex),
                            Integer.toString(replicationIndex),
                            Integer.toString(t),
                            Integer.toString(i),
                            Double.toString(effectAntiMet),
                            Double.toString(effectProMet),
                            Double.toString(conversionThreshold),
                            Double.toString(S2),  // switchSensitivity
                            Double.toString(S4), // divisionSensitivity
                            Double.toString(effectPerTumorCell),
                            Integer.toString(currentGrid.GetPop()),
                            Double.toString(fractalDimension),
                            Double.toString(lacunarity),
                            Double.toString(eccentricity),
                            chemoApplied
                    });
                } else {
                    // For intermediate timesteps, record default "NA" placeholders for morphology metrics.
                    batchData.add(new String[]{
                            Integer.toString(sampleIndex),
                            Integer.toString(replicationIndex),
                            Integer.toString(t),
                            Integer.toString(i),
                            Double.toString(effectAntiMet),
                            Double.toString(effectProMet),
                            Double.toString(conversionThreshold),
                            Double.toString(S2),
                            Double.toString(S4),
                            Double.toString(effectPerTumorCell),
                            Integer.toString(currentGrid.GetPop()),
                            "NA", "NA", "NA",
                            chemoApplied
                    });
                }
                localTimeStep++;
            }
            // Write the batch data for this timestep and clear the list for the next.
            OnLattice.utils.OutputManager.writeBatch(writer, batchData);
            batchData.clear();
        }
    }
}
