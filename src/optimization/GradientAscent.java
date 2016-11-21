package optimization;

import static valueObject.Constants.NEAR_ZERO;
import utility.Helper;
import valueObject.Constants;

/**
 * This class performs gradient descent using line search to find optimal step size
 * @author Tuan
 *
 */
public class GradientAscent {
	
	private OptimizedObject optimizedObject;
	
	//variables 
	private double[] variables;
	
	//this array stores updated variables, when trying different step size, used by line search
	private double[] triedVariables;
	
	//the derivatives
	private double[] derivatives;
	
	//the objective function to be optimized
	private double objectiveFn;
	
	//default step size
	private double step_size = 0.5;
	
	//gradient 2-norm
	private double gradientNorm;
	
	private boolean verbose;
	
	private double initialLearingRate = Constants.INITIAL_LEARNING_RATE;
	

	
	/**
	 * this should be used when calling gradient ascent multiple times, 
	 * set variables in very run but initialize derivative and tried_derivative once
	 * @param obj is optimized
	 * @param verbose
	 */
	public GradientAscent(OptimizedObject obj, boolean ver){
		this.optimizedObject = obj;
		this.verbose = ver;
	}
	/**
	 * this should be used when one calling is required
	 * @param obj
	 * @param x : variables
	 * @param verbose = true or false
	 */
	public GradientAscent(OptimizedObject obj, double[] x,boolean ver){
		this.optimizedObject = obj;
		this.variables = x;
		this.verbose = ver;
		
		//initialize gradient
		initialize();
	}

	/**
	 * initialization is performed if this same GradientAscent instance is called multiple times
	 */
	private void initialize(){
		if (variables != null){
			derivatives = new double[variables.length];			
			triedVariables = new double[variables.length];
		}
		
		gradientNorm = 0;
		objectiveFn = 0;
		step_size = 0.5;//default value		
	}
	
	/**
	 * The main steps of gradient descent is performed here
	 * @param: maximum number of iteration
	 * @throws InterruptedException 
	 */
	public void performGradientAscent(int max_iteration) throws Exception{
		int count = 0;
		initialize();
		
		//initial objective function value
		double oldObj;
		objectiveFn = optimizedObject.calGradientAndObjective(variables, derivatives);
		
		//perform line search and go up for the first step
//		if (isConvergence()){
//			return;
//		}
		
		step_size = lineSearch(variables, derivatives, objectiveFn);
		updateVariables(variables, derivatives, step_size);
		
		
		while(! isConvergence() && count < max_iteration){
			
			count++;			
			oldObj = objectiveFn;
			
			objectiveFn = optimizedObject.calGradientAndObjective(variables, derivatives);
			
			//if the objective doesn't increase, it is like taking a step backward, accept it and 
			//recalculate the step size
			if (objectiveFn < oldObj){		
				try{
					step_size = lineSearch(variables, derivatives, objectiveFn);
				}catch(Exception ex){
					break;
				}
				//continue;
			}
			
			if (verbose){
				System.out.printf("Iteration %1$d, objective function:%2$.5f, length of gradient %3$.7f\n",count,objectiveFn,gradientNorm);
			}
			
			//update variables with gradient
			updateVariables(variables, derivatives, step_size);					
			//step_size = initialLearingRate / (Math.sqrt(count));
		}
		
	}
	
	/**
	 * check if the size of derivatives/gradient is close to zero
	 */
	private boolean isConvergence(){
		double sum = 0;
		for(double d:derivatives){
			sum += Math.pow(d, 2);
		}
		sum = Math.sqrt(sum);
		
		gradientNorm = sum;
		
		if (sum < NEAR_ZERO * Math.abs(objectiveFn)){
			return true;
		}

		return false;
	}

	/**
	 * using Wolfe conditions to find inexact step size
	 * @param currentFn
	 * @param direction
	 * @return the approximately optimal step size for gradient ascent
	 * @throws Exception 
	 */
	private double lineSearch(double[] variables, double[] direction, double currentFn) throws Exception{
		
		double alpha = initialLearingRate * 2;

		double f1,f2;
		
		do{
			alpha = 0.5 * alpha;
			
			updateVariables(triedVariables,variables,direction,alpha);
			
			f1 = optimizedObject.calObjective(triedVariables);
			f2 = currentFn + 0.5 * alpha * gradientNorm * gradientNorm;
			
		}while ( f1 < f2 && alpha >= NEAR_ZERO);
		
		if (alpha <= NEAR_ZERO){
			System.err.println("Error in the objective function, try to reduce the initial learning rate !!! \n "
					+ "If it is already small, stop the program and call up Tuan! :( ");
			throw new Exception("Error in the objective function, try to reduce the initial learning rate !!! \n "
					+ "If it is already small, stop the program and call up Tuan! :( ");
		}
		
		if (verbose){
			System.out.println("Learning rate is:" + alpha);		
		}
		
		return alpha;
	
	}
	
	/**
	 * update variables directly into the source
	 * @param variables
	 * @param change
	 * @param stepSize
	 */
	
	private void updateVariables(double[] variables,double[] change,double stepSize){
		if (variables == null || change == null || variables.length != change.length){
			return;
		}
		for(int i = 0; i < variables.length; i++){
			variables[i] += stepSize * change[i];
		}
	}
	
	/**
	 * making update to the target variables from previous step variables ( this is mainly used in line search for optimal step size)
	 * @param target
	 * @param source
	 * @param change
	 * @param stepSize
	 */
	private void updateVariables(double[] target, double[] source, double[] change, double stepSize){
		if (target == null || source == null || change == null ||
				target.length != source.length || target.length != change.length){
			return;
		}
		
		for(int i = 0; i < target.length; i ++){
			target[i] = source[i] + stepSize * change[i];
		}
			
	}

	public double[] getVariables() {
		return variables;
	}

	public void setVariables(double[] variables) {
		this.variables = variables;
		if (this.variables != null && this.derivatives == null){
			initialize();
		}
	}
	public double getInitialLearingRate() {
		return initialLearingRate;
	}
	public void setInitialLearingRate(double initialLearingRate) {
		this.initialLearingRate = initialLearingRate;
	}
	
}
