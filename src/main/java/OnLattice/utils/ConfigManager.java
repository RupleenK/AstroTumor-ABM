package OnLattice.utils;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * The ConfigManager class loads and provides access to configuration parameters
 * from a properties file.
 */
public class ConfigManager {
    private Properties config;

    public ConfigManager(String filename) {
        config = new Properties();
        try (InputStream input = getClass().getClassLoader().getResourceAsStream(filename)) {
            if (input == null) {
                throw new IOException("Resource '" + filename + "' not found on classpath.");
            }
            config.load(input);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    // Returns the double value associated with the specified key.
    public double getDouble(String key) {
        return Double.parseDouble(config.getProperty(key));
    }

    // Returns a ParamSpec for a given parameter - discrete (list) vs continuous (min,max).
    public ParamSpec getParamSpec(String paramName) {
        String discreteKey = paramName + ".discrete";
        if (config.containsKey(discreteKey)) {
            String valuesStr = config.getProperty(discreteKey);
            String[] tokens = valuesStr.split(",");
            double[] values = new double[tokens.length];
            for (int i = 0; i < tokens.length; i++) {
                values[i] = Double.parseDouble(tokens[i].trim());
            }
            return new ParamSpec(values);
        } else {
            double min = getDouble(paramName + ".min");
            double max = getDouble(paramName + ".max");
            return new ParamSpec(min, max);
        }
    }

    // Returns an array of ParamSpec objects corresponding to your simulation parameters.
    public ParamSpec[] getParamSpecs() {
        ParamSpec[] specs = new ParamSpec[6];
        specs[0] = getParamSpec("effectAntiMet");
        specs[1] = getParamSpec("effectProMet");
        specs[2] = getParamSpec("conversionThreshold");
        specs[3] = getParamSpec("S2");
        specs[4] = getParamSpec("S4");
        specs[5] = getParamSpec("effectPerTumorCell");
        return specs;
    }
}
