package OnLattice.base;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class EdgeLacunarityTest {

    @Test
    public void testEdgeLacunarityForSparseTumorEdge() {
        ExampleGrid grid = new ExampleGrid(20, 20);
        // Place tumor cells spaced apart diagonally
        for (int i = 5; i <= 15; i += 2) {
            Cell cell = grid.NewAgentSQ(grid.I(i, i));
            cell.Init(CellType.celltype.TUMOR);
        }

        Outputs.isFrontCell(grid);
        double lacunarity = Outputs.calculateEdgeLacunarity(grid, Outputs.frontCellX, Outputs.frontCellY);

        assertTrue(lacunarity >= 0.0, "Lacunarity should be non-negative.");
        System.out.println("Edge lacunarity: " + lacunarity);
    }
}