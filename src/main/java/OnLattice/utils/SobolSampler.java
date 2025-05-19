package OnLattice.utils;

import org.apache.commons.math3.random.SobolSequenceGenerator;

/**
 * The SobolSampler class generates parameter samples using a Sobol sequence.
 * Given a number of samples and an array of parameter specifications (ParamSpec),
 * the sample method returns a 2D array of generated samples.
 */
public class SobolSampler {

    // Generates a 2D array of parameter samples.
    public static double[][] sample(int numSamples, ParamSpec[] specs) {
        int dims = specs.length;                         // number of dimensions equals the number of parameters
        double[][] samples = new double[numSamples][dims]; // array to store the generated samples
        SobolSequenceGenerator sobol = new SobolSequenceGenerator(dims); // Sobol sequence generator

        // Loop over each sample to generate a sample point.
        for (int i = 0; i < numSamples; i++) {
            double[] point = sobol.nextVector();         // generate next point from the Sobol sequence
            // Loop over each parameter dimension.
            for (int d = 0; d < dims; d++) {
                ParamSpec spec = specs[d];               // get parameter specification for this dimension
                double u = point[d];                   
                if (!spec.isDiscrete) {
                    // Continuous parameter: map value linearly from [0,1] to [min, max]
                    samples[i][d] = spec.min + u * (spec.max - spec.min);
                } else {
                    // Discrete parameter: choose the corresponding discrete value.
                    int n = spec.discreteValues.length;
                    int index = (int) Math.floor(u * n);   // convert continuous value to a discrete index
                    if (index >= n) index = n - 1;           // ensure index is within array bounds
                    samples[i][d] = spec.discreteValues[index];
                }
            }
        }
        return samples;
    }
}
