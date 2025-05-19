package OnLattice.base;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class FractalDimensionTest {

    @Test
    public void testFractalDimensionForSimpleTumor() {
        ExampleGrid grid = new ExampleGrid(20, 20);
        // Create a 3x3 solid tumor block in the center
        for (int x = 8; x <= 10; x++) {
            for (int y = 8; y <= 10; y++) {
                Cell cell = grid.NewAgentSQ(grid.I(x, y));
                cell.Init(CellType.celltype.TUMOR);
            }
        }

        Outputs.isFrontCell(grid);
        double fd = Outputs.calculateFractalDimension(Outputs.frontCellX, Outputs.frontCellY, grid);

        assertTrue(fd > 0.0 && fd < 2.5, "Fractal dimension should be within valid range.");
        System.out.println("Fractal dimension: " + fd);
    }
}