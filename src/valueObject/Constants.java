package valueObject;

public class Constants {


	
	public static final String NUM_KEY = "NUM";
		
	public static final String OUTPUT_FOLDER_KEY = "OUTPUT_FOLDER";
	public static final String INPUT_FILE_KEY = "INPUT_FILE";
	
	public static final String CONVERT_FACTOR_KEY = "CONVERT_FACTOR";
	
	public static final String VERBOSE_KEY = "VERBOSE";
	
	public static final String CHR_UPPER_BOUND_ID_FILE_KEY = "CHR_UPPER_BOUND_ID_FILE";
	
	public static final String LEARNING_RATE_KEY = "LEARNING_RATE";
	
	public static final String MAX_ITERATION_KEY = "MAX_ITERATION";
	
	public static final String CHR_LENGTH_KEY = "CHROMOSOME_LENGTH";
	
	public static final String THRESHOLD_KEY = "CONTACT_THRESHOLD";
	
	
	//maximum number of threads should be used 
	public static final int MAX_NUM_THREAD = 120;
	
	//the starting learning rate for the line search
	public static double INITIAL_LEARNING_RATE = 0.001;		
	
	//maximum number of iterations
	public static int MAX_ITER = 200000;
	
	//this constant is used to check if the norm of the gradient is near zero
	public static final double NEAR_ZERO = 10e-6;

	//if the distance is larger than LARGE_DISTANCE_FOR_TANH, it will be scale down to this value
	public static final double SCALE_DISTANCE = 15.0;
	
	//average distance will be scaled down to this value
	public static final double AVG_DIST = 10.0;//test Sep 01, 2016
	
	public static final double WIDE_CURVE = 25;

}
