package OnLattice.analysis;

import OnLattice.base.Cell;
import OnLattice.base.CellType;
import OnLattice.base.ExampleGrid;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.*;
import static OnLattice.base.CellType.celltype.ASTROCYTE;
import static OnLattice.base.CellType.celltype.NEUTRAL;

/**
 * The DensityAnalysis class provides methods for performing density-based analyses.
 * In particular, it offers methods to initialize neutral cells in a grid and then convert
 * a subset of those neutral cells into astrocytes. These methods are useful when you want
 * to simulate differential density of astrocytes by distributing cells uniformly at random
 * and then selectively converting them into astrocytes.
 */

public class DensityAnalysis {

    // Randomly initializes a specified number of neutral agents in the grid.
    public static void initRandomNeutralAgents(ExampleGrid grid, int neutralsum) {
        Random rand = new Random(12345);
        HashSet<String> placedCells = new HashSet<>(); 
        List<int[]> positions = new ArrayList<>();

        int centerX = grid.xDim / 2;
        int centerY = grid.yDim / 2;

        // Generate all possible positions in the grid, excluding the center.
        for (int x = 0; x < grid.xDim; x++) {
            for (int y = 0; y < grid.yDim; y++) {
                if (x == centerX && y == centerY) continue; // Skip center
                positions.add(new int[]{x, y});
            }
        }

        // Shuffle positions to randomize cell placement.
        Collections.shuffle(positions, rand);

        int placedNeutralCell = 0;

        // Place neutral cells until reaching the target or running out of positions.
        for (int i = 0; i < neutralsum && i < positions.size(); i++) {
            int[] pos = positions.get(i);
            int intX = pos[0], intY = pos[1];

            // Check if cell is empty before placing.
            if (grid.IsEmpty(intX, intY)) {
                Cell newCell = grid.NewAgentSQ(intX, intY);
                newCell.Init(NEUTRAL);
                placedCells.add(intX + "," + intY); // Track the placed cell
                placedNeutralCell++;
            }
        }
    }


    // Converts a subset of neutral cells in the grid into astrocytes.
    public static void ConvertNeutralCellsToAstrocytes(ExampleGrid grid, int numAstrocytes) {
        ArrayList<Cell> neutralCells = new ArrayList<>();

        // Iterate over all cells in the grid and collect neutral cells (excluding the center)
        for (Cell cell : grid) {
            if (cell.isNeutral() &&
                    !(cell.Xsq() == grid.xDim / 2 && cell.Ysq() == grid.yDim / 2)) {
                neutralCells.add(cell);
            }
        }

        Random rand = new Random(12345);
        int placedAstrocytes = 0;

        // Randomly convert neutral cells into astrocytes until the target is reached.
        while (placedAstrocytes < numAstrocytes && !neutralCells.isEmpty()) {
            int index = rand.nextInt(neutralCells.size());
            Cell cell = neutralCells.get(index);
            if (cell != null) {
                // Convert the cell's type to ASTROCYTE and update its astrocyte state and color.
                cell.setCellType(ASTROCYTE);
                cell.getCellTypeInfo().setAstrocyteState(CellType.AstrocyteState.ANTI_METASTATIC);
                cell.updateColor();
                neutralCells.remove(index);
                placedAstrocytes++;
            }
        }
    }
}
