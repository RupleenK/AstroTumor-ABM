package OnLattice.base;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Rand;
import static HAL.Util.*;
import java.util.ArrayList;

// The ExampleGrid class extends AgentGrid2D and represents a two-dimensional simulation grid.
// It manages cell placement, behavior (including tumor division and astrocyte switching),
// and the diffusion of chemotherapy agents.
public class ExampleGrid extends AgentGrid2D<OnLattice.base.Cell> {

    // Map storing modulation values for each grid site.
    public double[] modulationMap;

    // Static lists to record raw and normalized astrocyte influences.
    public static ArrayList<Double> astroLocalInfluences = new ArrayList<>();
    public static ArrayList<Double> astroNormalizedInfluences = new ArrayList<>();

    // Static lists to record raw and normalized tumor division influences.
    public static ArrayList<Double> tumorLocalInfluences = new ArrayList<>();
    public static ArrayList<Double> tumorNormalizedInfluences = new ArrayList<>();

    // Records timesteps when chemotherapy is applied.
    public ArrayList<Integer> chemoApplicationTimesteps = new ArrayList<>();

    // Random number generator for grid operations.
    public Rand rng = new Rand();

    // PDEGrid2D for handling chemotherapy diffusion.
    public PDEGrid2D chemotherapy;

    // Flags to control astrocyte switching and chemotherapy application.
    public boolean enableAstrocyteSwitch;
    public boolean enableChemotherapy;

    // Precomputed neighborhood indices for cell division.
    public int[] divHood = MooreHood(false);

    // Parameters for cell behavior
    public double effectProMet;
    public double effectAntiMet;
    public double S4;
    public double conversionThreshold;
    public double effectPerTumorCell;
    public double S2;

    // Constructs an ExampleGrid with the specified dimensions.
    public ExampleGrid(int xDim, int yDim) {
        super(xDim, yDim, Cell.class);
        chemotherapy = new PDEGrid2D(xDim, yDim);
        modulationMap = new double[xDim * yDim];
    }

    // Draws the current state of the grid onto the provided GridWindow.
    public void DrawModel(GridWindow win, int xOffset, int yOffset) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                int pixColor;
                Cell cell = GetAgent(x, y);
                if (cell != null) {
                    pixColor = cell.getColor();
                } else if (enableChemotherapy) {
                    double chemoVal = chemotherapy.Get(x, y);
                    pixColor = (chemoVal > 0) ? HeatMapParula(chemoVal) : BLACK;
                } else {
                    pixColor = BLACK;
                }
                if ((x + xOffset) < win.xDim && (y + yOffset) < win.yDim) {
                    win.SetPix(x + xOffset, y + yOffset, pixColor);
                }
            }
        }
    }

    // Returns the number of tumor cells present in the grid.
    public int GetPop() {
        int count = 0;
        for (int i = 0; i < length; i++) {
            OnLattice.base.Cell agent = GetAgent(i);
            if (agent != null && agent.isTumor()) {
                count++;
            }
        }
        return count;
    }

    // Chemo parameters
    public static double TIMESTEP = 17.0 / 24;
    public double SPACE_STEP = 25;
    public double DRUG_PERIOD = 41.5;  // Approximately 29.4 timesteps for a 21-day cycle.
    public double DRUG_DURATION = 1; // Duration of chemotherapy (1 day).
    public double DRUG_BOUNDARY_VAL = 0.095;
    public static double CHEMO_DECAY_RATE = 0.90;
    public static double DRUG_UPTAKE = -0.005 * TIMESTEP;
    public double DRUG_DIFF_RATE = 10 * 60 * 60 * 24 * (TIMESTEP / (SPACE_STEP * SPACE_STEP));
    public static double DRUG_DEATH = 1.0;
    public static double GAP_JUNCTION_MODULATION = 0.3;

    private static final int TUMOR_THRESHOLD = 3000;
    private int lastChemoCycle = -1;
    private boolean chemoStarted = false;
    private int chemoStartTime = -1;
    public int gridIndex = -1;

    // Manages chemotherapy administration over time.
    private void manageChemotherapy(int currentTimeStep) {
        if (!enableChemotherapy) return;

        // Start chemotherapy when tumor population exceeds threshold.
        if (!chemoStarted && GetPop() >= 3000) {
            chemoStarted = true;
            chemoStartTime = currentTimeStep;
              }

        if (chemoStarted) {
            double elapsedDays = (currentTimeStep - chemoStartTime) * TIMESTEP;
            int currentCycle = (int) (elapsedDays / DRUG_PERIOD);

            if (currentCycle > lastChemoCycle) {
                lastChemoCycle = currentCycle;
                 // Apply fresh drug: perform diffusion with boundary condition.
                chemotherapy.DiffusionADI(DRUG_DIFF_RATE, DRUG_BOUNDARY_VAL);
                chemotherapy.Update();
                if (!chemoApplicationTimesteps.contains(currentTimeStep)) {
                    chemoApplicationTimesteps.add(currentTimeStep);
                }
            } else {
                chemotherapy.DiffusionADI(DRUG_DIFF_RATE);
                chemotherapy.Update();
            }

            // Apply drug decay at each grid point.
            for (int x = 0; x < xDim; x++) {
                for (int y = 0; y < yDim; y++) {
                    double currentVal = chemotherapy.Get(x, y);
                    chemotherapy.Set(x, y, Math.max(currentVal * CHEMO_DECAY_RATE, 0));
                }
            }
            chemotherapy.Update();
        }
    }

    // Steps through each cell in the grid and updates its behavior.
    public void StepCells(double conversionThreshold,
                          double effectPerTumorCell,
                          double effectProMet,
                          double effectAntiMet,
                          double DRUG_DEATH,
                          int currentTimeStep,
                          double S2,
                          double S4) {

        int tumorCellCount = 0;

         // Iterate through each cell in the grid.
        for (Cell cell : this) {
            if (cell.isTumor()) {
                tumorCellCount++;
                cell.TumorCellDivision(effectProMet, effectAntiMet, S4);
                if (enableChemotherapy) {
                    cell.ChemoDeath();
                }
            } else if (cell.isAstrocyte()) {
                cell.AstrocyteSwitching(conversionThreshold, effectPerTumorCell, S2);
            }
            cell.updateColor();
        }

        if (enableChemotherapy) {
            manageChemotherapy(currentTimeStep);
        }
    }

    // Checks if the given coordinates are within grid bounds.
    public boolean InBounds(int x, int y) {
            return (x | (xDim - 1 - x) | y | (yDim - 1 - y)) >= 0;
    }

    // Checks if the grid cell at (x, y) is empty.
    public boolean IsEmpty(int x, int y) {
        return GetAgent(x, y) == null;
    }
}
