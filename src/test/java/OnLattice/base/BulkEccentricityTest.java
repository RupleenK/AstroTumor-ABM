package OnLattice.base;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class BulkEccentricityTest {

    @Test
    public void testEccentricityCircularVsLinear() {
        ExampleGrid grid = new ExampleGrid(20, 20);
        // Create a circular cluster
        for (int dx = -2; dx <= 2; dx++) {
            for (int dy = -2; dy <= 2; dy++) {
                if (Math.sqrt(dx*dx + dy*dy) <= 2) {
                    int x = 10 + dx;
                    int y = 10 + dy;
                    Cell cell = grid.NewAgentSQ(grid.I(x, y));
                    cell.Init(CellType.celltype.TUMOR);
                }
            }
        }

        double ecc = Outputs.calculateBulkEccentricity(grid);
        assertTrue(ecc == 0, "Eccentricity should be [0].");
        System.out.println("Bulk eccentricity (circular tumor): " + ecc);
    }
}