package OnLattice.analysis;

import OnLattice.base.Cell;
import OnLattice.base.CellType;
import OnLattice.base.ExampleGrid;
import HAL.Rand;
import java.util.ArrayList;
import java.util.*;

/**
 * Provides various methods for initializing spatial distributions of astrocytes in a grid.
 */
public class SpatialDistributions {

    // Determines whether the given (x, y) coordinate is the center of the grid.
    private static boolean isCenterCell(ExampleGrid grid, int x, int y) {
        return (x == grid.xDim / 2 && y == grid.yDim / 2);
    }

    // Initializes the grid with a uniform (approximately square) distribution of astrocytes.
    public static void initUniformAstrocytes(ExampleGrid grid, int numAstrocytes) {
        int numRows = (int) Math.sqrt(numAstrocytes);
        int numCols = (int) Math.ceil((double) numAstrocytes / numRows);
        double rowSpacing = (double) grid.xDim / numRows;
        double colSpacing = (double) grid.yDim / numCols;

        int placed = 0;
        for (int i = 0; i < numRows && placed < numAstrocytes; i++) {
            for (int j = 0; j < numCols && placed < numAstrocytes; j++) {
                int posX = (int) (i * rowSpacing + rowSpacing / 2);
                int posY = (int) (j * colSpacing + colSpacing / 2);
                posX = Math.min(posX, grid.xDim - 1);
                posY = Math.min(posY, grid.yDim - 1);
                // Skip the center cell.
                if (isCenterCell(grid, posX, posY)) continue;
                if (grid.IsEmpty(posX, posY)) {
                    Cell newCell = grid.NewAgentSQ(posX, posY);
                    newCell.Init(CellType.celltype.ASTROCYTE);
                    placed++;
                }
            }
        }
        System.out.println("✅ Uniformly placed " + placed + " astrocytes out of target " + numAstrocytes);
    }

    // Initializes the grid with a gradient distribution of astrocytes.
    public static void initGradientAstrocytes(ExampleGrid grid, int numAstrocytes, double gradientExponent) {
        ArrayList<Integer> siteIndices = new ArrayList<>();
        ArrayList<Double> siteWeights = new ArrayList<>();

        // Loop over grid positions, skipping the center and occupied cells.
        for (int x = 0; x < grid.xDim; x++) {
            for (int y = 0; y < grid.yDim; y++) {
                if (isCenterCell(grid, x, y)) continue;
                if (!grid.IsEmpty(x, y)) continue;
                double weight = Math.pow((double) (grid.yDim - y) / grid.yDim, gradientExponent);
                if (weight > 1e-12) {
                    int iSq = grid.I(x, y);
                    siteIndices.add(iSq);
                    siteWeights.add(weight);
                }
            }
        }
        if (numAstrocytes > siteIndices.size()) {
            System.out.println("Warning: Requested " + numAstrocytes + " astrocytes, but only " + siteIndices.size() + " available sites.");
            numAstrocytes = siteIndices.size();
        }
        ArrayList<Integer> chosenSites = WeightedSampleWithoutReplacement(siteIndices, siteWeights, numAstrocytes, grid.rng);
        for (int iSq : chosenSites) {
            int xx = iSq % grid.xDim;
            int yy = iSq / grid.xDim;
            Cell newCell = grid.NewAgentSQ(xx, yy);
            newCell.Init(CellType.celltype.ASTROCYTE);
        }
       System.out.println("✅ Placed " + numAstrocytes + " astrocytes in a gradient (exponent=" + gradientExponent + ").");
    }

    // Initializes the grid with a random distribution of astrocytes.
    public static void initRandomAstrocytes(ExampleGrid grid, int numAstrocytes) {
        Random rand = new Random(54321);
        HashSet<String> placedCells = new HashSet<>();
        int placedAstrocytes = 0;
        int centerX = grid.xDim / 2;
        int centerY = grid.yDim / 2;

        while (placedAstrocytes < numAstrocytes) {
            int intX = rand.nextInt(grid.xDim);
            int intY = rand.nextInt(grid.yDim);
            if (intX == centerX && intY == centerY) continue;
            String cellKey = intX + "," + intY;
            if (grid.IsEmpty(intX, intY) && !placedCells.contains(cellKey)) {
                Cell newCell = grid.NewAgentSQ(intX, intY);
                newCell.Init(CellType.celltype.ASTROCYTE);
                placedCells.add(cellKey);
                placedAstrocytes++;
            }
        }
       System.out.println("✅ Randomly placed " + placedAstrocytes + " astrocytes out of target " + numAstrocytes);
    }

    // Initializes the grid with an inverse radial distribution of astrocytes.
    public static void initInverseRadialAstrocytes(ExampleGrid grid, int numAstrocytes, long seed, double exponent) {
        ArrayList<Integer> siteIndices = new ArrayList<>();
        ArrayList<Double> siteWeights = new ArrayList<>();
        int centerX = grid.xDim / 2;
        int centerY = grid.yDim / 2;
        double maxDist = Math.hypot(centerX, centerY);

        // Iterate over every grid cell.
        for (int x = 0; x < grid.xDim; x++) {
            for (int y = 0; y < grid.yDim; y++) {
                if (isCenterCell(grid, x, y)) continue; // Always skip center cell
                if (!grid.IsEmpty(x, y)) continue;
                double d = Math.hypot(x - centerX, y - centerY);
                double weight = Math.pow(d / maxDist, exponent);
                siteIndices.add(grid.I(x, y));
                siteWeights.add(weight);
            }
        }

        if (siteIndices.size() < numAstrocytes) {
            System.out.println("Warning: Not enough available grid cells to place " + numAstrocytes + " astrocytes.");
            numAstrocytes = siteIndices.size();
        }

        // Use the helper method for weighted sampling without replacement.
        ArrayList<Integer> chosenSites = WeightedSampleWithoutReplacement(siteIndices, siteWeights, numAstrocytes, grid.rng);
        for (int iSq : chosenSites) {
            int xx = iSq % grid.xDim;
            int yy = iSq / grid.xDim;
            Cell newCell = grid.NewAgentSQ(xx, yy);
            newCell.Init(CellType.celltype.ASTROCYTE);
        }
        System.out.println("✅ Inverse radial: Placed " + numAstrocytes + " astrocytes (seed: " + seed + ", exponent: " + exponent + ").");
    }

    // Initializes the grid with clumped astrocyte groups.
    public static void initClumpedAstrocytes(ExampleGrid grid, int numAstrocytes, int numClumpsRow, int numClumpsCol, int clumpRadius, long seed) {
        Random rand = new Random(seed);
        int totalClumps = numClumpsRow * numClumpsCol;
        int astroPerClump = numAstrocytes / totalClumps;
        int remainder = numAstrocytes % totalClumps;

        int totalPlaced = 0;
        int blockWidth = grid.xDim / numClumpsRow;
        int blockHeight = grid.yDim / numClumpsCol;

        for (int i = 0; i < numClumpsRow; i++) {
            for (int j = 0; j < numClumpsCol; j++) {
                int clumpIndex = i * numClumpsCol + j;
                int astroCountHere = astroPerClump + (clumpIndex < remainder ? 1 : 0);
                int blockStartX = i * blockWidth;
                int blockStartY = j * blockHeight;
                int centerX = blockStartX + blockWidth / 2;
                int centerY = blockStartY + blockHeight / 2;
                // If the block center is the grid center, adjust to skip it.
                if (isCenterCell(grid, centerX, centerY)) {
                    centerX++;
                }

                int placedInClump = 0;
                int attempts = 0;
                int maxAttempts = astroCountHere * 10;
                while (placedInClump < astroCountHere && attempts < maxAttempts) {
                    attempts++;
                    double angle = 2 * Math.PI * rand.nextDouble();
                    double rFrac = Math.sqrt(rand.nextDouble());
                    double rActual = clumpRadius * rFrac;
                    int dx = (int) Math.round(rActual * Math.cos(angle));
                    int dy = (int) Math.round(rActual * Math.sin(angle));
                    int posX = centerX + dx;
                    int posY = centerY + dy;
                    if (!grid.InBounds(posX, posY)) continue;
                    if (isCenterCell(grid, posX, posY)) continue;
                    if (!grid.IsEmpty(posX, posY)) continue;
                    Cell newCell = grid.NewAgentSQ(posX, posY);
                    newCell.Init(CellType.celltype.ASTROCYTE);
                    placedInClump++;
                    totalPlaced++;
                }
            }
        }
        System.out.println("✅ Clumped: Placed " + totalPlaced + " astrocytes (target: " + numAstrocytes + ").");
    }

    // Initializes the grid with a radial gradient distribution of astrocytes.
    public static void initRadialAstrocytes(ExampleGrid grid, int numAstrocytes, double exponent, double centerWeight) {
        ArrayList<Integer> siteIndices = new ArrayList<>();
        ArrayList<Double> siteWeights = new ArrayList<>();
        int centerX = grid.xDim / 2;
        int centerY = grid.yDim / 2;
        double maxRadius = Math.sqrt(Math.pow(grid.xDim / 2.0, 2) + Math.pow(grid.yDim / 2.0, 2));

        for (int x = 0; x < grid.xDim; x++) {
            for (int y = 0; y < grid.yDim; y++) {
                if (isCenterCell(grid, x, y)) continue;
                if (!grid.IsEmpty(x, y)) continue;
                double dx = x - centerX;
                double dy = y - centerY;
                double dist = Math.sqrt(dx * dx + dy * dy);
                double ratio = dist / maxRadius;
                double weight = centerWeight * Math.pow(Math.max(0, 1.0 - ratio), exponent);
                if (weight > 1e-12) {
                    int iSq = grid.I(x, y);
                    siteIndices.add(iSq);
                    siteWeights.add(weight);
                }
            }
        }
        if (numAstrocytes > siteIndices.size()) {
            System.out.println("Warning: Requested more astrocytes than available sites! Clamping to " + siteIndices.size());
            numAstrocytes = siteIndices.size();
        }
        ArrayList<Integer> chosenSites = WeightedSampleWithoutReplacement(siteIndices, siteWeights, numAstrocytes, grid.rng);
        for (int iSq : chosenSites) {
            int xx = iSq % grid.xDim;
            int yy = iSq / grid.xDim;
            Cell newCell = grid.NewAgentSQ(xx, yy);
            newCell.Init(CellType.celltype.ASTROCYTE);
        }
        System.out.println("✅ Radial: Placed " + numAstrocytes + " astrocytes (exponent=" + exponent + ").");
    }

    // Performs weighted sampling without replacement.Each candidate site is assigned a weight based on distance from the center.
    // The method then randomly selects a specified number of these cells according to their weights,
    // ensuring that each chosen cell is unique (no replacement). This allows modeling non-uniform distributions—where
    // some areas are more likely to be chosen than others—while ensuring that you don’t place multiple astrocytes in the same grid site.
    public static ArrayList<Integer> WeightedSampleWithoutReplacement(ArrayList<Integer> siteIndices, ArrayList<Double> siteWeights, int k, Rand rng) {
        ArrayList<Integer> chosen = new ArrayList<>(k);
        double totalWeight = 0.0;
        for (double w : siteWeights) {
            totalWeight += w;
        }
        // Perform k picks
        for (int pick = 0; pick < k; pick++) {
            double r = rng.Double() * totalWeight;
            double cumulative = 0.0;
            int chosenIdx = -1;
            for (int i = 0; i < siteWeights.size(); i++) {
                cumulative += siteWeights.get(i);
                if (cumulative >= r) {
                    chosenIdx = i;
                    break;
                }
            }
            chosen.add(siteIndices.get(chosenIdx));
            totalWeight -= siteWeights.get(chosenIdx);
            siteWeights.set(chosenIdx, 0.0);
        }
        return chosen;
    }
}
