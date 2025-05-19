package OnLattice.base;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class ChemoDeathTest {

    @Test
    public void tumorDiesWhenChemoIsHigh() {
        ExampleGrid grid = new ExampleGrid(100, 100);
        grid.enableChemotherapy = true;

        Cell tumor = grid.NewAgentSQ(50, 50);
        tumor.Init(CellType.celltype.TUMOR);

        grid.chemotherapy.Set(tumor.Isq(), 1.0); // Use Add to apply immediately

        for (int i = 0; i < 50; i++) {
            tumor.ChemoDeath();
            grid.chemotherapy.Update();
            if (grid.GetAgent(50, 50) == null) {
                break;
            }
        }
        assertNull(grid.GetAgent(50, 50), "Tumor should die when chemo is high");
    }

    @Test
    public void tumorSurvivesWhenChemoIsAbsent() {
        ExampleGrid grid = new ExampleGrid(100, 100);
        grid.enableChemotherapy = true;

        Cell tumor = grid.NewAgentSQ(50, 50);
        tumor.Init(CellType.celltype.TUMOR);

        grid.chemotherapy.Add(tumor.Isq(), 0.0); // Very low concentration

        for (int i = 0; i < 50; i++) {
            tumor.ChemoDeath();
            grid.chemotherapy.Update();
        }
        assertNotNull(grid.GetAgent(50, 50), "Tumor should survive at no chemo");
    }

    @Test
    public void chemoEffectModulatedByProAstrocytes() {
        ExampleGrid grid = new ExampleGrid(100, 100);
        grid.enableChemotherapy = true;
        grid.enableAstrocyteSwitch = true;

        Cell tumor = grid.NewAgentSQ(50, 50);
        tumor.Init(CellType.celltype.TUMOR);
        grid.chemotherapy.Add(tumor.Isq(), 1.0);

        // Add pro-metastatic astrocytes around
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                if (dx == 0 && dy == 0) continue;
                Cell astro = grid.NewAgentSQ(50 + dx, 50 + dy);
                astro.Init(CellType.celltype.ASTROCYTE);
                astro.getCellTypeInfo().setAstrocyteState(CellType.AstrocyteState.PRO_METASTATIC);
            }
        }
        boolean survived = false;
        for (int i = 0; i < 50; i++) {
            tumor.ChemoDeath();
            if (grid.GetAgent(50, 50) != null) {
                survived = true;
            } else {
                survived = false;
                break;
            }
        }
        assertTrue(survived, "Pro-metastatic astrocytes should reduce chemo effectiveness and allow tumor survival");
    }
}