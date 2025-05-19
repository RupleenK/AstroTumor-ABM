package OnLattice.utils;

/**
 * The ParamSpec class encapsulates a simulation parameter specification.
 * It supports both continuous (min, max) and discrete (enumerated values) parameters.
 */
public class ParamSpec {

    public boolean isDiscrete;      // true if this parameter is discrete
    public double[] discreteValues; // only used if isDiscrete is true
    public double min;              // minimum value (used if continuous)
    public double max;              // maximum value (used if continuous)

    // Constructor for continuous parameters.
    public ParamSpec(double min, double max) {
        this.isDiscrete = false;
        this.min = min;
        this.max = max;
    }

    // Constructor for discrete parameters.
    public ParamSpec(double[] discreteValues) {
        this.isDiscrete = true;
        this.discreteValues = discreteValues;
    }
}
