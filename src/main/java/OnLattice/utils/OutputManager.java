package OnLattice.utils;

import com.opencsv.CSVWriter;
import java.io.IOException;
import java.util.List;

/**
 * The OutputManager class handles writing simulation results to CSV files.
 */
public class OutputManager {

    // Writes the CSV header.
    public static void writeHeader(CSVWriter writer) throws IOException {
        writer.writeNext(new String[]{
                "ParamSetIndex", "Replication", "TimeStep", "GridIndex",
                "effectAntiMet", "effectProMet", "conversionThreshold", "S2", "S4",
                "effectPerTumorCell", "TumorCellCount", "FractalDimension", "Lacunarity", "Eccentricity", "chemoApplied"
        });
    }

    // Writes a batch of data rows to the CSV.
    public static void writeBatch(CSVWriter writer, List<String[]> batchData) throws IOException {
        writer.writeAll(batchData);
    }
}
