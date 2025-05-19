package OnLattice.base;

import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import java.util.ArrayList;
import static HAL.Util.*;
import static OnLattice.base.CellType.celltype.ASTROCYTE;
import static OnLattice.base.CellType.celltype.TUMOR;
import static OnLattice.base.ExampleGrid.*;

/**
 * The Cell class represents a cell agent on a 2D grid. It implements
 * behaviors such as tumor cell division, astrocyte switching, and chemotherapy-induced death.
 */
public class Cell extends AgentSQ2Dunstackable<ExampleGrid> {
    private int color;       // Visualization color
    private CellType type;   // The cell's type and state information

    // Initializes the cell with the specified primary cell type.
    public void Init(CellType.celltype cellType) {
        this.type = new CellType(cellType);
        updateColor();
    }

    // Returns the CellType of the cell.
    public CellType getCellTypeInfo() {
        return this.type;
    }

    // Checks if the cell is a tumor cell.
    public boolean isTumor() {
        return this.getCellTypeInfo().getCellTypeEnum() == CellType.celltype.TUMOR;
    }

    // Checks if the cell is an astrocyte.
    public boolean isAstrocyte() {
        return this.getCellTypeInfo().getCellTypeEnum() == CellType.celltype.ASTROCYTE;
    }

    // Checks if the cell is neutral.
    public boolean isNeutral() {
        return this.getCellTypeInfo().getCellTypeEnum() == CellType.celltype.NEUTRAL;
    }

    // Sets the cell's primary type.
    public void setCellType(CellType.celltype newType) {
        this.type.setType(newType);
    }

    // Returns the current visualization color of the cell.
    public int getColor() {
        return this.color;
    }

    // Implements tumor cell division. The method calculates the net influence from neighboring astrocytes
    // (using a weighted function based on distance) and determines if cell division should occur based on
    // a sigmoid function of the normalized influence. If division is successful, a new tumor cell is placed
    // in an adjacent empty location if there is one.
    public void TumorCellDivision(double effectProMet, double effectAntiMet, double S4) {
        int maxRadius = 3;
        double totalInfluence = 0.0;

        // Loop over a square neighborhood centered on this cell.
        for (int dx = -maxRadius; dx <= maxRadius; dx++) {
            for (int dy = -maxRadius; dy <= maxRadius; dy++) {
                int checkX = Xsq() + dx;
                int checkY = Ysq() + dy;
                if (G.InBounds(checkX, checkY)) {
                    Cell neighbor = G.GetAgent(checkX, checkY);
                    if ((neighbor != null) && (neighbor.isAstrocyte())) {
                        double distance = Math.sqrt(dx * dx + dy * dy);
                        double influenceScale = 1.0 / (1.0 + Math.exp(3 * (distance - 1.5)));
                        if (neighbor.getCellTypeInfo().getAstrocyteState() == CellType.AstrocyteState.PRO_METASTATIC) {
                            totalInfluence += influenceScale * effectProMet;
                        } else if (neighbor.getCellTypeInfo().getAstrocyteState() == CellType.AstrocyteState.ANTI_METASTATIC) {
                            totalInfluence -= influenceScale * effectAntiMet;
                        }
                    }
                }
            }
        }
        // Record the influence values for analysis.
        double normalizedInfluence = totalInfluence / 7.23;
        double rawVal = S4 * normalizedInfluence;
        double adjDivProb = 1.0 / (1.0 + Math.exp(-rawVal));

        // If a random draw is less than the adjusted division probability, perform cell division.
        if (G.rng.Double() < adjDivProb) {
            ArrayList<Integer> allowedSquares = new ArrayList<>();
            int options = G.MapEmptyHood(G.divHood, Xsq(), Ysq());
            for (int i = 0; i < options; i++) {
                Cell neighbor = G.GetAgent(G.divHood[i]);
                if (neighbor == null || !neighbor.isAstrocyte()) {
                    allowedSquares.add(G.divHood[i]);
                }
            }
            if (!allowedSquares.isEmpty()) {
                int chosenSq = allowedSquares.get(G.rng.Int(allowedSquares.size()));
                Cell newCell = G.NewAgentSQ(chosenSq);
                newCell.Init(TUMOR);
            }
        }
    }

    // Implements the astrocyte switching behavior. If the current cell is an astrocyte in the ANTI_METASTATIC state
    // and the local tumor influence (calculated from nearby tumor cells) exceeds a conversion threshold, the cell may
    // switch its state to PRO_METASTATIC.
    public void AstrocyteSwitching(double conversionThreshold, double effectPerTumorCell, double S2) {
        
        // Exit early if switching is disabled or if the cell is not an anti-metastatic astrocyte.
        if (!G.enableAstrocyteSwitch ||
                this.type.getCellTypeEnum() != ASTROCYTE ||
                this.type.getAstrocyteState() != CellType.AstrocyteState.ANTI_METASTATIC) {
            return;
        }

        double localInfluence = 0.0;
        int maxRadius = 3;
        // Sum the influence from nearby tumor cells.
        for (int dx = -maxRadius; dx <= maxRadius; dx++) {
            for (int dy = -maxRadius; dy <= maxRadius; dy++) {
                int checkX = Xsq() + dx;
                int checkY = Ysq() + dy;
                double distance = Math.sqrt(dx * dx + dy * dy);
                if (G.InBounds(checkX, checkY) && distance <= maxRadius) {
                    Cell neighbor = G.GetAgent(checkX, checkY);
                    if (neighbor != null && neighbor.isTumor()) {
                        localInfluence += effectPerTumorCell / (1.0 + Math.exp(3 * (distance - 1.5)));
                    }
                }
            }
        }

        // Record and normalize the influence for analysis.
        double normalizedInfluence = localInfluence / 7.23;
        double probabilityOfSwitching = 1.0 / (1.0 + Math.exp(-S2 * (normalizedInfluence - conversionThreshold)));
        // Switch state if probability criteria is met.
        if (localInfluence > 0.0 && G.rng.Double() < probabilityOfSwitching) {
            this.type.setAstrocyteState(CellType.AstrocyteState.PRO_METASTATIC);
            updateColor();
        }
    }

    // Simulates chemotherapy-induced cell death.
    // This method reduces the chemotherapy concentration at the cell's location, computes a death probability based on
    // that concentration, and disposes of the cell if a random number is less than the death probability.
    public void ChemoDeath() {
        // Retrieve current drug concentration at this cell's location.
        double current = G.chemotherapy.Get(Isq());
        double uptake = current * ExampleGrid.DRUG_UPTAKE;
        uptake = Math.min(uptake, current);  // Ensure uptake does not exceed available concentration.
        G.chemotherapy.Add(Isq(), -uptake);
        double localConcentration = current - uptake;

        // Adjust local concentration if this is a tumor cell and astrocyte switching is enabled.
        int proAstroCount = 0;
        if (this.isTumor() && G.enableAstrocyteSwitch) {
            int[] hood = MooreHood(false);  // Get positions around this cell.
            int nHood = MapHood(hood); // Find out how many neighbor positions there are.
            for (int i = 0; i < nHood; i++) {  // Loop over each neighbor position.
                Cell neighbor = G.GetAgent(hood[i]);   // Get the cell in this neighbor position.
                if (neighbor != null && neighbor.isAstrocyte() &&  // If there is a cell and it is a supportive cell (in the pro-tumor state),
                        neighbor.getCellTypeInfo().getAstrocyteState() == CellType.AstrocyteState.PRO_METASTATIC) {
                    proAstroCount++;
                }
            }

            // Adjust the remaining drug level based on how many prometastatic astrocytes are around.
            double modulationFactor = 1.0 - (ExampleGrid.GAP_JUNCTION_MODULATION * 1); //((double) proAstroCount / 8.0));
            localConcentration *= modulationFactor;
            G.modulationMap[Isq()] = modulationFactor;
        }

        // Compute cell death probability via a sigmoid function.
        double k = 5;
        double threshold = 0.5;
        double deathProb = 1.0 / (1.0 + Math.exp(-k * (localConcentration - threshold)));
        double rand = G.rng.Double();

        if (this.isTumor() && rand < deathProb && G.enableChemotherapy && localConcentration>0.001) {
            Dispose();
        }
    }

    // Updates the cell's visualization color based on its type and state.
    public void updateColor() {
        if (this.isTumor()) {
            this.color = CYAN;
        } else if (this.isAstrocyte()) {
            if (this.type.getAstrocyteState() == CellType.AstrocyteState.ANTI_METASTATIC) {
                this.color = MAGENTA;
            } else if (this.type.getAstrocyteState() == CellType.AstrocyteState.PRO_METASTATIC) {
                this.color = YELLOW;
            }
        } else {
            this.color = BLUE;
        }
    }

    // Counts the number of pro-metastatic astrocytes in the cell's neighborhood.
    public int countProAstrosAround() {
        int[] hood = MooreHood(false);
        int nHood = G.MapHood(hood, Xsq(), Ysq());
        int count = 0;
        
        for (int i = 0; i < nHood; i++) {
            Cell neighbor = G.GetAgent(hood[i]);
            if (neighbor != null && neighbor.isAstrocyte()
                    && neighbor.getCellTypeInfo().getAstrocyteState() == CellType.AstrocyteState.PRO_METASTATIC) {
                count++;
            }
        }
        return count;
    }
}
