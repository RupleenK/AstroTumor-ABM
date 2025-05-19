package OnLattice.base;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class AstrocyteSwitchingTest {

    @Test
    public void astrocyteShouldNotSwitchWhenDisabled() {
        ExampleGrid grid = new ExampleGrid(100, 100);
        grid.enableAstrocyteSwitch = false;

        Cell astro = grid.NewAgentSQ(50, 50);
        astro.Init(CellType.celltype.ASTROCYTE);
        astro.getCellTypeInfo().setAstrocyteState(CellType.AstrocyteState.ANTI_METASTATIC);

        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                if (!(dx == 0 && dy == 0)) {
                    Cell tumor = grid.NewAgentSQ(50 + dx, 50 + dy);
                    tumor.Init(CellType.celltype.TUMOR);
                }
            }
        }

        astro.AstrocyteSwitching(0.1, 1.0, 1.0);

        assertEquals(CellType.AstrocyteState.ANTI_METASTATIC,
                astro.getCellTypeInfo().getAstrocyteState(),
                "Astrocyte should not switch when switching is disabled");
    }

    @Test
    public void astrocyteShouldSwitchWhenInfluenceIsHigh() {
        ExampleGrid grid = new ExampleGrid(100, 100);
        grid.enableAstrocyteSwitch = true;

        Cell astro = grid.NewAgentSQ(50, 50);
        astro.Init(CellType.celltype.ASTROCYTE);
        astro.getCellTypeInfo().setAstrocyteState(CellType.AstrocyteState.ANTI_METASTATIC);

        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                if (!(dx == 0 && dy == 0)) {
                    Cell tumor = grid.NewAgentSQ(50 + dx, 50 + dy);
                    tumor.Init(CellType.celltype.TUMOR);
                }
            }
        }

        for (int i = 0; i < 100; i++) {
            astro.AstrocyteSwitching(0.5, 1.0, 8.0); // Strong influence
            if (astro.getCellTypeInfo().getAstrocyteState() == CellType.AstrocyteState.PRO_METASTATIC) {
                break;
            }
        }

        assertEquals(CellType.AstrocyteState.PRO_METASTATIC,
                astro.getCellTypeInfo().getAstrocyteState(),
                "Astrocyte should switch when influence is high");
    }

    @Test
    public void astrocyteShouldNotSwitchWhenInfluenceIsTooLow() {
        ExampleGrid grid = new ExampleGrid(100, 100);
        grid.enableAstrocyteSwitch = true;

        Cell astro = grid.NewAgentSQ(50, 50);
        astro.Init(CellType.celltype.ASTROCYTE);
        astro.getCellTypeInfo().setAstrocyteState(CellType.AstrocyteState.ANTI_METASTATIC);

        // No tumor cells nearby
        astro.AstrocyteSwitching(0.5, 1.0, 8.0);

        assertEquals(CellType.AstrocyteState.ANTI_METASTATIC,
                astro.getCellTypeInfo().getAstrocyteState(),
                "Astrocyte should not switch when tumor influence is too low");
    }

    @Test
    public void astrocyteSwitchingInfluenceCurve() {
        double conversionThreshold = 0.5;
        double effectPerTumorCell = 1.0;
        double S2 = 8.0;

        for (int tumorCount = 0; tumorCount <= 8; tumorCount++) {
            ExampleGrid grid = new ExampleGrid(100, 100);
            grid.enableAstrocyteSwitch = true;

            Cell astro = grid.NewAgentSQ(50, 50);
            astro.Init(CellType.celltype.ASTROCYTE);
            astro.getCellTypeInfo().setAstrocyteState(CellType.AstrocyteState.ANTI_METASTATIC);

            int placed = 0;
            for (int dx = -1; dx <= 1 && placed < tumorCount; dx++) {
                for (int dy = -1; dy <= 1 && placed < tumorCount; dy++) {
                    if (dx == 0 && dy == 0) continue;
                    Cell tumor = grid.NewAgentSQ(50 + dx, 50 + dy);
                    tumor.Init(CellType.celltype.TUMOR);
                    placed++;
                }
            }

            double localInfluence = 0.0;
            int maxRadius = 3;
            for (int dx = -maxRadius; dx <= maxRadius; dx++) {
                for (int dy = -maxRadius; dy <= maxRadius; dy++) {
                    int x = 50 + dx;
                    int y = 50 + dy;
                    double distance = Math.sqrt(dx * dx + dy * dy);
                    if (grid.InBounds(x, y) && distance <= maxRadius) {
                        Cell neighbor = grid.GetAgent(x, y);
                        if (neighbor != null && neighbor.isTumor()) {
                            localInfluence += effectPerTumorCell / (1.0 + Math.exp(3 * (distance - 1.5)));
                        }
                    }
                }
            }

            double normalized = localInfluence / 7.23;
            double probability = 1.0 / (1.0 + Math.exp(-S2 * (normalized - conversionThreshold)));

            System.out.printf("Tumors: %d, Influence: %.4f, Normalized: %.4f, ProbSwitch: %.4f\n",
                    tumorCount, localInfluence, normalized, probability);
        }
    }
}
