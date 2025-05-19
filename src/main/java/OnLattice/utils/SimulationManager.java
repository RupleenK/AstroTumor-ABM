package OnLattice.utils;

import OnLattice.base.ExampleGrid;
import HAL.Gui.GifMaker;
import HAL.Gui.GridWindow;
import com.opencsv.CSVWriter;
import java.awt.GraphicsEnvironment;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Properties;

/**
 * The SimulationManager class loads configuration parameters, sets global values, and runs the simulation.
 */
public class SimulationManager {

    // Enum representing various analysis.
    public enum AnalysisType {
        BASE, BASE_PSA, DENSITY, DENSITY_PSA, SPATIAL, SPATIAL_PSA, CHEMO
    }

    // Global parameters (to be overridden by config)
    public static int x;
    public static int y;
    public static int visScale;
    public static int timesteps;
    public static double DRUG_UPTAKE;
    public static double DRUG_DEATH;
    public static int numSamples;
    public static int replications;
    public static int spatialDistributions;

    // Loads configuration file.
    private static Properties loadConfiguration(String filename) {
        Properties config = new Properties();
        try (InputStream input = SimulationManager.class.getClassLoader().getResourceAsStream(filename)) {
            if (input == null) {
                throw new IOException("Unable to find configuration file: " + filename);
            }
            config.load(input);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return config;
    }

    // Overrides defaults using the configuration file and returns the analysis mode.
    private static AnalysisType configureFromFile(Properties config) {
        String modeStr = config.getProperty("analysis.mode");
        AnalysisType mode = AnalysisType.valueOf(modeStr.toUpperCase());
        x = Integer.parseInt(config.getProperty("grid.width"));
        y = Integer.parseInt(config.getProperty("grid.height"));
        visScale = Integer.parseInt(config.getProperty("vis.scale"));
        timesteps = Integer.parseInt(config.getProperty("timesteps"));
        numSamples = Integer.parseInt(config.getProperty("num.samples"));
        replications = Integer.parseInt(config.getProperty("replications"));
        spatialDistributions = Integer.parseInt(config.getProperty("spatial.distributions"));
        return mode;
    }

    // Performs parameter sampling using the SobolSampler.
    public static double[][] performSobolSampling(int numSamples, ParamSpec[] specs) {
        return SobolSampler.sample(numSamples, specs);
    }

    // -------------------------------------------------------------------
    // Main method for launching the simulation.
    // -------------------------------------------------------------------
    
    public static void main(String[] args) {

        // Determine if the environment supports graphics.
        boolean isHeadless = GraphicsEnvironment.isHeadless();

        // Set default configuration values
        String configFile = "config_chemo_targeted.properties";
        int jobIndex = 0;

        if (args.length >= 2) {
            configFile = args[0] + ".properties";
            try {
                jobIndex = Integer.parseInt(args[1]) - 1; 
            } catch (Exception e) {
                System.err.println("Invalid job index provided, defaulting to 0.");
                jobIndex = 0;
            }
        } else if (args.length == 1) {
            try {
                jobIndex = Integer.parseInt(args[0]) - 1;
            } catch (Exception e) {
                System.err.println("Invalid job index provided, defaulting to 0.");
                jobIndex = 0;
            }
        }

        // Load and apply configuration.
        Properties config = loadConfiguration(configFile);
        AnalysisType mode = configureFromFile(config);

        // Create a new ConfigManager instance using the same config file.
        ConfigManager configMgr = new ConfigManager(configFile);
        ParamSpec[] specs = configMgr.getParamSpecs();

        // Update parameters from configuration file.
        numSamples = Integer.parseInt(config.getProperty("num.samples"));
        replications = Integer.parseInt(config.getProperty("replications"));
        spatialDistributions = Integer.parseInt(config.getProperty("spatial.distributions"));

        // Perform Sobol sampling for parameter sets.
        double[][] allSamples = performSobolSampling(numSamples, specs);

        // Validate job index.
        if (jobIndex < 0 || jobIndex >= allSamples.length) {
            System.out.println("Job index " + jobIndex + " is out of range. Exiting.");
            return;
        }
        double[] paramSet = allSamples[jobIndex];
        System.out.println("Running parameter set " + jobIndex + ": " + Arrays.toString(paramSet));

        // Set up visualization window if not headless.
        GridWindow win = null;
        GifMaker gifMaker = null;
        if (!isHeadless) {
            int rows, cols;
            if (spatialDistributions == 6) {
                rows = 2;
                cols = 3;
            } else if (spatialDistributions == 2) {
                rows = 1;
                cols = 2;
            } else {
                rows = 1;
                cols = 1;
            }
            int windowWidth = cols * x;
            int windowHeight = rows * y;
            win = new GridWindow("Model", windowWidth, windowHeight, visScale);
            gifMaker = new GifMaker("output.gif", 1, true);
        }

        // Open the CSV output file for writing simulation data.
        String filePath = "HAL_" + System.currentTimeMillis() + ".csv";
        try (CSVWriter writer = new CSVWriter(
                new FileWriter(filePath),
                CSVWriter.DEFAULT_SEPARATOR,
                CSVWriter.NO_QUOTE_CHARACTER,
                CSVWriter.DEFAULT_ESCAPE_CHARACTER,
                CSVWriter.DEFAULT_LINE_END)) {

            // Write CSV header information.
            OutputManager.writeHeader(writer);

            // Run simulation replications for the selected parameter set.
            for (int replicationIndex = 0; replicationIndex < replications; replicationIndex++) {
                System.out.printf("Running replication %d for parameter set %d%n", replicationIndex, jobIndex);
                // Create the simulation grids based on the chosen analysis mode.
                ExampleGrid[] grids = OnLattice.utils.GridFactory.createGrids(mode, config);
                // Run the simulation using the provided parameters and output tools.
                OnLattice.utils.SimulationRunner.runSimulation(
                        grids,
                        paramSet[0],  // effectAntiMet
                        paramSet[1],  // effectProMet
                        paramSet[2],  // conversionThreshold
                        paramSet[3],  // S2 (switchSensitivity)
                        paramSet[4],  // S4 (divisionSensitivity)
                        paramSet[5],  // effectPerTumorCell
                        DRUG_DEATH,
                        DRUG_UPTAKE,
                        timesteps,
                        jobIndex,         // sample index for output
                        replicationIndex, // replication index for output
                        writer,
                        gifMaker,
                        win
                );
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("End.");
    }
}
