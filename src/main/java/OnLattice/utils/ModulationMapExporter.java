package OnLattice.base;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;

/**
 * The ModulationMapExporter class exports a CSV file representing the spatial distribution of
 * modulation factors on the grid.
 */
public class ModulationMapExporter {
    // Exports tumor cell modulation factors with their coordinates in CSV format.
    public static void ExportTumorCellsModulation(ExampleGrid grid, int timeStep, String filename) {
        try (PrintWriter writer = new PrintWriter(new FileWriter(filename))) {
            // Write a header row for the CSV output.
            writer.println("x,y,modFactor");

            // Loop through every position in the grid using a linear index.
            for (int i = 0; i < grid.length; i++) {
                Cell c = grid.GetAgent(i);
                if (c != null && c.isTumor()) {
                    // Convert the 1D index 'i' into 2D coordinates (x, y).
                    int x = grid.ItoX(i);
                    int y = grid.ItoY(i);

                    // Compute the modulation factor.
                    int proAstroCount = c.countProAstrosAround();
                    double modFactor = 1.0 - (ExampleGrid.GAP_JUNCTION_MODULATION * ((double) proAstroCount / 8.0));

                    // Write a CSV row: x-coordinate, y-coordinate, and modulation factor.
                    writer.printf("%d,%d,%.3f\n", x, y, modFactor);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
