package OnLattice.base;

import java.util.ArrayList;
import java.util.List;
import static HAL.Util.*;

/**
 * The class provides static methods for computing spatial morphology metrics of the tumor.
 */
public class Outputs {

    // Lists to store the x and y coordinates of the front cells.
    public static ArrayList<Integer> frontCellX = new ArrayList<>();
    public static ArrayList<Integer> frontCellY = new ArrayList<>();

    /**
     * Determines whether the cell at (x, y) in the grid is a front cell.
     * A front cell is defined as a tumor cell that has at least one adjacent empty square.
     */

    public static void isFrontCell(ExampleGrid grid) {
        // Clear previous front cell coordinates.
        frontCellX.clear();
        frontCellY.clear();
        // Precompute the Moore neighborhood offsets.
        int[] hood = MooreHood(false);

        // Iterate over each cell in the grid.
        for (Cell cell : grid) {
            if (cell.isTumor()) {
                int x = cell.Xsq();  // Get the cell's x-coordinate.
                int y = cell.Ysq();  // Get the cell's y-coordinate.
                // Count empty neighboring cells using the precomputed Moore neighborhood.
                int emptyNeighbors = grid.MapEmptyHood(hood, x, y);
                if (emptyNeighbors > 0) {
                    frontCellX.add(x);
                    frontCellY.add(y);
                }
            }
        }
    }

    // Checks whether any front cell falls within a specified box.
    private static boolean containsBoundaryCell(int boxX, int boxY, int size, List<Integer> frontCellX, List<Integer> frontCellY) {
        for (int i = 0; i < frontCellX.size(); i++) {
            int cellX = frontCellX.get(i);
            int cellY = frontCellY.get(i);
            if (cellX >= boxX && cellX < boxX + size && cellY >= boxY && cellY < boxY + size) {
                return true;
            }
        }
        return false;
    }

    // Calculates the slope for the regression between log(box sizes) and log(box counts).
    private static double calculateSlope(ArrayList<Double> logBoxSizes, ArrayList<Double> logBoxCounts) {
        if (logBoxSizes.size() < 2) {
            throw new IllegalArgumentException("Not enough data for regression.");
        }
        double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
        int n = logBoxSizes.size();
        for (int i = 0; i < n; i++) {
            double x = logBoxSizes.get(i);
            double y = logBoxCounts.get(i);
            sumX += x;
            sumY += y;
            sumXY += x * y;
            sumXX += x * x;
        }
        // Return the negative slope computed from the regression.
        return -(sumXY - n * (sumX / n) * (sumY / n)) / (sumXX - n * Math.pow(sumX / n, 2));
    }

    // Calculates the fractal dimension of the tumor boundary based on the distribution of front cells.
    // The method uses a box-counting approach: it subdivides the grid into boxes of various sizes,
    // counts the number of boxes that contain at least one front cell, and then performs a linear regression
    // on the logarithm of the box sizes versus the logarithm of the counts.
    public static double calculateFractalDimension(List<Integer> frontCellX, List<Integer> frontCellY, ExampleGrid grid) {
        // Define the set of box sizes to use for the box-counting method
        int[] boxSizes = {2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256};

        ArrayList<Double> logBoxSizes = new ArrayList<>();
        ArrayList<Double> logBoxCounts = new ArrayList<>();

        // Loop over each box size and compute the number of boxes that contain at least one front cell.
        for (int size : boxSizes) {
            int count = 0;
            int numBoxesX = (int) Math.ceil((double) grid.xDim / size);
            int numBoxesY = (int) Math.ceil((double) grid.yDim / size);
            for (int x = 0; x < numBoxesX; x++) {
                for (int y = 0; y < numBoxesY; y++) {
                    if (containsBoundaryCell(x * size, y * size, size, frontCellX, frontCellY)) {
                        count++;
                    }
                }
            }
            logBoxSizes.add(Math.log(size));
            logBoxCounts.add(Math.log(count));
        }
        return calculateSlope(logBoxSizes, logBoxCounts);
    }

    // Calculates the edge lacunarity of the tumor boundary. This method first creates a binary representation of the front cells,
    // then subdivides the region into boxes of various sizes, and computes the lacunarity
    // based on the variance-to-mean ratio of the counts within those boxes.
    public static double calculateEdgeLacunarity(ExampleGrid grid, List<Integer> frontCellX, List<Integer> frontCellY) {

        // Create a binary grid indicating where front cells are located.
        boolean[][] isEdge = new boolean[grid.xDim][grid.yDim];
        for (int i = 0; i < frontCellX.size(); i++) {
            int fx = frontCellX.get(i);
            int fy = frontCellY.get(i);
            isEdge[fx][fy] = true;
        }

        // Determine the smallest box that contains all edge cells.
        int minX = grid.xDim, maxX = 0, minY = grid.yDim, maxY = 0;
        for (int x = 0; x < grid.xDim; x++) {
            for (int y = 0; y < grid.yDim; y++) {
                if (isEdge[x][y]) {
                    if (x < minX) minX = x;
                    if (x > maxX) maxX = x;
                    if (y < minY) minY = y;
                    if (y > maxY) maxY = y;
                }
            }
        }
        if (minX > maxX || minY > maxY) {
            return 0;
        }

        // Create a subgrid that bounds the tumor edge.
        int subWidth = maxX - minX + 1;
        int subHeight = maxY - minY + 1;
        int[][] subgrid = new int[subWidth][subHeight];
        for (int i = 0; i < subWidth; i++) {
            for (int j = 0; j < subHeight; j++) {
                int gridX = minX + i;
                int gridY = minY + j;
                subgrid[i][j] = isEdge[gridX][gridY] ? 1 : 0;
            }
        }

        // Divide the tumor edge into boxes of varying sizes and compute lacunarity.
        int[] boxSizes = {2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256};
        double lacunaritySum = 0;
        int countSizes = 0;
        for (int boxSize : boxSizes) {
            if (boxSize > subWidth || boxSize > subHeight) {
                continue;
            }
            ArrayList<Integer> boxCounts = new ArrayList<>();
            for (int i = 0; i < subWidth; i += boxSize) {
                for (int j = 0; j < subHeight; j += boxSize) {
                    int count = 0;
                    for (int dx = 0; dx < boxSize && (i + dx) < subWidth; dx++) {
                        for (int dy = 0; dy < boxSize && (j + dy) < subHeight; dy++) {
                            count += subgrid[i + dx][j + dy];
                        }
                    }
                    boxCounts.add(count);
                }
            }
            // Compute lacunarity for the current box size.
            double mean = boxCounts.stream().mapToInt(Integer::intValue).average().orElse(0.0);
            double variance = boxCounts.stream().mapToDouble(c -> Math.pow(c - mean, 2)).average().orElse(0.0);
            if (mean > 0) {
                lacunaritySum += variance / (mean * mean);
                countSizes++;
            }
        }
        return countSizes > 0 ? lacunaritySum / countSizes : 0;
    }

    // Calculates the bulk eccentricity of the tumor.
    // This method computes the covariance matrix of the tumor cell coordinates,
    // determines its eigenvalues, and calculates the eccentricity as a function of
    // the ratio of the minor to major axis lengths.
    public static double calculateBulkEccentricity(ExampleGrid grid) {
        ArrayList<Integer> tumorX = new ArrayList<>();
        ArrayList<Integer> tumorY = new ArrayList<>();

        // Populate the lists with tumor cell coordinates
        for (Cell cell : grid) {
            if (cell.isTumor()) {
                tumorX.add(cell.Xsq());
                tumorY.add(cell.Ysq());
            }
        }

        int nTumor = tumorX.size();
        if (nTumor < 2) {
            return Double.NaN;
        }

        double sumX = 0, sumY = 0;
        for (int i = 0; i < nTumor; i++) {
            sumX += tumorX.get(i);
            sumY += tumorY.get(i);
        }
        double centerX = sumX / nTumor;
        double centerY = sumY / nTumor;

        double covXX = 0, covXY = 0, covYY = 0;
        for (int i = 0; i < nTumor; i++) {
            double x = tumorX.get(i) - centerX;
            double y = tumorY.get(i) - centerY;
            covXX += x * x;
            covXY += x * y;
            covYY += y * y;
        }
        covXX /= nTumor;
        covXY /= nTumor;
        covYY /= nTumor;

        double T = covXX + covYY;
        double D = covXX * covYY - covXY * covXY;
        double eigenvalue1 = (T + Math.sqrt(T * T - 4 * D)) / 2.0;
        double eigenvalue2 = (T - Math.sqrt(T * T - 4 * D)) / 2.0;
        double a = (eigenvalue1 > 0) ? Math.sqrt(eigenvalue1) : 0;
        double b = (eigenvalue2 > 0) ? Math.sqrt(eigenvalue2) : 0;
        if (a == 0 || b == 0) {
            return Double.NaN;
        }
        return Math.sqrt(1 - (b * b) / (a * a));
    }
}
