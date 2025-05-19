package OnLattice.utils;

import OnLattice.analysis.DensityAnalysis;
import OnLattice.base.ExampleGrid;
import OnLattice.base.Cell;
import OnLattice.analysis.SpatialDistributions;
import OnLattice.utils.SimulationManager.AnalysisType;
import java.util.Properties;
import static OnLattice.base.CellType.celltype.TUMOR;

/**
 * The GridFactory class creates an array of ExampleGrid objects based on the specified
 * analysis mode and configuration properties.
 */
public class GridFactory {

    public static ExampleGrid[] createGrids(AnalysisType mode, Properties config) {
        // Retrieve the number of spatial distributions from the configuration.
        int spatialDistributions = Integer.parseInt(config.getProperty("spatial.distributions"));
        ExampleGrid[] sim;
        int defaultAst = 27000; // Default number of astrocytes for initialization.
        switch (mode) {
            case BASE:
            case BASE_PSA:
                sim = new ExampleGrid[spatialDistributions];
                for (int i = 0; i < spatialDistributions; i++) {
                    sim[i] = new ExampleGrid(SimulationManager.x, SimulationManager.y);
                    // Create grids and initialize astrocytes using a fixed spatial method.
                    SpatialDistributions.initRandomAstrocytes(sim[i], defaultAst);
                    // Set per-grid toggles read from config file.
                    boolean gridAstroSwitch = Boolean.parseBoolean(config.getProperty("grid." + i + ".enable.astrocyteSwitch"));
                    boolean gridChemotherapy = Boolean.parseBoolean(config.getProperty("grid." + i + ".enable.chemotherapy"));
                    sim[i].enableAstrocyteSwitch = gridAstroSwitch;
                    sim[i].enableChemotherapy = gridChemotherapy;
                    // Place the tumor cell at the center of the grid.
                    Cell centralCell = sim[i].NewAgentSQ(sim[i].xDim / 2, sim[i].yDim / 2);
                    centralCell.Init(TUMOR);
                }
                break;
            case DENSITY:
                sim = new ExampleGrid[spatialDistributions];
                // Create grids and initialize with a fixed number of neutral agents.
                for (int i = 0; i < spatialDistributions; i++) {
                    sim[i] = new ExampleGrid(SimulationManager.x, SimulationManager.y);
                    DensityAnalysis.initRandomNeutralAgents(sim[i], 45000);
                    // Read grid-specific astrocyte density and convert neutral cells.
                    double gridAstroDensity = Double.parseDouble(config.getProperty("grid." + i + ".astroGridDensity"));
                    int numAstrocytes = (int)(gridAstroDensity * 45000);
                    DensityAnalysis.ConvertNeutralCellsToAstrocytes(sim[i], numAstrocytes);
                    // Set per-grid toggles.
                    boolean gridAstroSwitch = Boolean.parseBoolean(config.getProperty("grid." + i + ".enable.astrocyteSwitch"));
                    boolean gridChemotherapy = Boolean.parseBoolean(config.getProperty("grid." + i + ".enable.chemotherapy"));
                    sim[i].enableAstrocyteSwitch = gridAstroSwitch;
                    sim[i].enableChemotherapy = gridChemotherapy;
                    // Place the tumor cell at the grid center.
                    Cell centralCell = sim[i].NewAgentSQ(sim[i].xDim / 2, sim[i].yDim / 2);
                    centralCell.Init(TUMOR);
                }
                break;
            case DENSITY_PSA:
                sim = new ExampleGrid[spatialDistributions];
                // Create grids and initialize with 45,000 neutral agents.
                for (int i = 0; i < spatialDistributions; i++) {
                    sim[i] = new ExampleGrid(SimulationManager.x, SimulationManager.y);
                    DensityAnalysis.initRandomNeutralAgents(sim[i], 45000);
                    // Read grid-specific astrocyte density and perform conversion.
                    double gridAstroDensity = Double.parseDouble(config.getProperty("grid." + i + ".astroGridDensity"));
                    int numAstrocytes = (int)(gridAstroDensity * 45000);
                    DensityAnalysis.ConvertNeutralCellsToAstrocytes(sim[i], numAstrocytes);
                    // Set grid toggles.
                    boolean gridAstroSwitch = Boolean.parseBoolean(config.getProperty("grid." + i + ".enable.astrocyteSwitch"));
                    boolean gridChemotherapy = Boolean.parseBoolean(config.getProperty("grid." + i + ".enable.chemotherapy"));
                    sim[i].enableAstrocyteSwitch = gridAstroSwitch;
                    sim[i].enableChemotherapy = gridChemotherapy;
                    // Place the tumor cell at the center.
                    Cell centralCell = sim[i].NewAgentSQ(sim[i].xDim / 2, sim[i].yDim / 2);
                    centralCell.Init(TUMOR);
                }
                break;
            case SPATIAL:
            case SPATIAL_PSA:
                sim = new ExampleGrid[spatialDistributions];
                // Create grids and initialize astrocytes using various spatial distribution methods.
                for (int i = 0; i < spatialDistributions; i++) {
                    sim[i] = new ExampleGrid(SimulationManager.x, SimulationManager.y);
                    // Read the distribution type for this grid from the configuration.
                    String distType = config.getProperty("distribution." + i, "random").toLowerCase();
                    switch (distType) {
                        case "uniform":
                            SpatialDistributions.initUniformAstrocytes(sim[i], defaultAst);
                            break;
                        case "random":
                            SpatialDistributions.initRandomAstrocytes(sim[i], defaultAst);
                            break;
                        case "clumped":
                            SpatialDistributions.initClumpedAstrocytes(sim[i], defaultAst, 6, 6, 23, 1235);
                            break;
                        case "inverseradial":
                            SpatialDistributions.initInverseRadialAstrocytes(sim[i], defaultAst, 1235, 0.3);
                            break;
                        case "radial":
                            SpatialDistributions.initRadialAstrocytes(sim[i], defaultAst, 1, 0.1);
                            break;
                        case "gradient":
                            SpatialDistributions.initGradientAstrocytes(sim[i], defaultAst, 2);
                            break;
                        default:
                            SpatialDistributions.initRandomAstrocytes(sim[i], defaultAst);
                            break;
                    }
                    // Set per-grid toggles.
                    boolean gridAstroSwitch = Boolean.parseBoolean(config.getProperty("grid." + i + ".enable.astrocyteSwitch"));
                    boolean gridChemotherapy = Boolean.parseBoolean(config.getProperty("grid." + i + ".enable.chemotherapy"));
                    sim[i].enableAstrocyteSwitch = gridAstroSwitch;
                    sim[i].enableChemotherapy = gridChemotherapy;
                    // Place the tumor cell at the grid center.
                    Cell centralCell = sim[i].NewAgentSQ(sim[i].xDim / 2, sim[i].yDim / 2);
                    centralCell.Init(TUMOR);
                }
                break;
            case CHEMO:
                sim = new ExampleGrid[spatialDistributions];
                // Create grids and initialize astrocytes using the random distribution.
                for (int i = 0; i < spatialDistributions; i++) {
                    sim[i] = new ExampleGrid(SimulationManager.x, SimulationManager.y);
                    SpatialDistributions.initRandomAstrocytes(sim[i], defaultAst);
                    // Set per-grid toggles.
                    boolean gridAstroSwitch = Boolean.parseBoolean(config.getProperty("grid." + i + ".enable.astrocyteSwitch"));
                    boolean gridChemotherapy = Boolean.parseBoolean(config.getProperty("grid." + i + ".enable.chemotherapy"));
                    sim[i].enableAstrocyteSwitch = gridAstroSwitch;
                    sim[i].enableChemotherapy = gridChemotherapy;
                    // Place the tumor cell at the grid center.
                    Cell centralCell = sim[i].NewAgentSQ(sim[i].xDim / 2, sim[i].yDim / 2);
                    centralCell.Init(TUMOR);
                }
                break;
            default:
                // Create a single grid with default settings.
                sim = new ExampleGrid[1];
                sim[0] = new ExampleGrid(SimulationManager.x, SimulationManager.y);
                SpatialDistributions.initRandomAstrocytes(sim[0], defaultAst);
                sim[0].enableAstrocyteSwitch = false;
                sim[0].enableChemotherapy = false;
                // Place the tumor cell at the center.
                Cell centralCell = sim[0].NewAgentSQ(sim[0].xDim / 2, sim[0].yDim / 2);
                centralCell.Init(TUMOR);
                break;
        }
            return sim;
    }
}
