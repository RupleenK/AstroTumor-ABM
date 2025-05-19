package OnLattice.base;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class TumorCellDivisionTest {

    @Test
    public void tumorDivisionProbabilityCurve() {
        double effectProMet = 1.0;
        double effectAntiMet = 1.0;
        double S4 = 8.0;

        for (int proCount = 0; proCount <= 8; proCount++) {
            ExampleGrid grid = new ExampleGrid(100, 100);
            Cell tumor = grid.NewAgentSQ(50, 50);
            tumor.Init(CellType.celltype.TUMOR);

            int placed = 0;
            for (int dx = -1; dx <= 1 && placed < proCount; dx++) {
                for (int dy = -1; dy <= 1 && placed < proCount; dy++) {
                    if (dx == 0 && dy == 0) continue;
                    Cell astro = grid.NewAgentSQ(50 + dx, 50 + dy);
                    astro.Init(CellType.celltype.ASTROCYTE);
                    astro.getCellTypeInfo().setAstrocyteState(CellType.AstrocyteState.PRO_METASTATIC);
                    placed++;
                }
            }

            double totalInfluence = 0.0;
            int maxRadius = 3;
            for (int dx = -maxRadius; dx <= maxRadius; dx++) {
                for (int dy = -maxRadius; dy <= maxRadius; dy++) {
                    int x = 50 + dx;
                    int y = 50 + dy;
                    if (grid.InBounds(x, y)) {
                        Cell neighbor = grid.GetAgent(x, y);
                        if (neighbor != null && neighbor.isAstrocyte()) {
                            double distance = Math.sqrt(dx * dx + dy * dy);
                            double influenceScale = 1.0 / (1.0 + Math.exp(3 * (distance - 1.5)));
                            if (neighbor.getCellTypeInfo().getAstrocyteState() == CellType.AstrocyteState.PRO_METASTATIC) {
                                totalInfluence += influenceScale * effectProMet;
                            }
                        }
                    }
                }
            }

            double normalized = totalInfluence / 7.23;
            double divisionProb = 1.0 / (1.0 + Math.exp(-S4 * normalized));

            System.out.printf("ProAstro: %d, Influence: %.4f, Normalized: %.4f, DivProb: %.4f\n",
                    proCount, totalInfluence, normalized, divisionProb);
        }
    }
}
