package noisy_mds;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import optimization.GradientAscent;
import optimization.OptimizedObject;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import shrec3D.Evaluate;
import utility.Helper;
import valueObject.Constants;
import valueObject.Constraint;
import evaluation.CalRMSD;

public class StructureGeneratorLorentz implements OptimizedObject{

	private Helper helper = Helper.getHelperInstance();
	
	//number of structures will be generated
	private int NUM;
	
	//list of constraints, each contains position i,j, IF, dist
	private List<Constraint> lstCons;
	
	//list to map 0..n to pos1...posn
	private List<Integer> lstPos;
	
	//factor to convert IF to physical distance
	private double convertFactor = -1.0;
	
	//input file
	private String INPUT_FILE = null;
	
	//output folder
	private String OUTPUT_FOLDER = "";
	
	//structure
	private double[] str;
	
	//number of points
	private int n;
	
	//learning rate for optimization
	private double LEARNING_RATE = Constants.INITIAL_LEARNING_RATE;
	
	//print intermediate result during optimization or not
	private boolean VERBOSE = true;
	
	//maximum distance that will be scaled down to	
	//private double maxScale = Constants.SCALE_DISTANCE; 
	
	//list contains indices of sub constraints for parallel threads
	private ArrayList<Integer> lstSubDataSetId = new ArrayList<Integer>();

	//file prefix to name output file structure
	private String FILE_PREFIX;
	
	private int MAX_ITER = Constants.MAX_ITER;
	
	private String parameterFile;
	
	private double avgIF;
	
	private double totalIF;
	
	private double maxIF;	
	
	private double AVG_DIST = Constants.AVG_DIST;
	private double WIDE_CURVE = Constants.WIDE_CURVE;
			
	//interval to ignore when calculating Spearman correlation, a[i,i + interval] = 0
	private int interval = 5;
	
	private double contactThres;
	private int[] chrLens = null;
	private HashMap<Integer,Integer> idToChr = new HashMap<Integer,Integer>(); //to map index to chromosome
	
	public StructureGeneratorLorentz(String parameterFile){
		this.parameterFile = parameterFile;
	}
	
	/**
	 * read input contacts
	 * @throws Exception
	 */
	private void readInput() throws Exception{
		
		
		lstPos = new ArrayList<Integer>();
		
		lstCons = helper.readContactList(INPUT_FILE, lstPos,contactThres);
		//lstCons = helper.readContactMatrixAsList(INPUT_FILE);
		//n = helper.determineNbrOfPoints(INPUT_FILE);
		
		n = lstPos.size();		
		System.out.println("Number of points: " + n);
		//map position to absolute id
		Map<Integer,Integer> mapPosToID = new HashMap<Integer,Integer>();
		for(int i = 0; i < lstPos.size(); i++){
			mapPosToID.put(lstPos.get(i), i);
		}
		//
		
		//output the mapping of coordinate to id in the output structure
		PrintWriter pw = new PrintWriter(OUTPUT_FOLDER + "/" + FILE_PREFIX + "_coordinate_mapping.txt");
		for(int i = 0; i < lstPos.size(); i++){
			pw.println(lstPos.get(i) + "\t" + i);			
		}		
		pw.close();
		
		//calculate chromosome number for each index
		if (chrLens != null){
			
			for(int i = 1; i <chrLens.length; i++){
				chrLens[i] = chrLens[i - 1] + chrLens[i];
			}
			for(int j = 0; j < chrLens[0]; j++){
				idToChr.put(j, 0);
			}
			for(int i = 1; i < chrLens.length; i++){
				for(int j = chrLens[i-1]; j < chrLens[i]; j++){
					idToChr.put(j, i);
				}
			}
			
		}else{
			//if idToChr is null, make the whole as one chromosome
			for(int i = 0; i < n; i++){				
				idToChr.put(i, 0);				
			}
		}				

		
		//correct lstCons to remove gap, if there is no gap, this doesn't change anything
		avgIF = 0.0;
		for(Constraint con: lstCons){
//			if (Math.abs(con.getPos1() - con.getPos2()) == 1 && con.getIF() < 20.0){
//				System.out.println(con.getPos1() + "\t" + con.getPos2());
//			}
			con.setPos1(mapPosToID.get(con.getPos1()));
			con.setPos2(mapPosToID.get(con.getPos2()));
			avgIF += con.getIF();
		}
		avgIF /= lstCons.size();
		
		maxIF = 0.0;
		//scale average distance to AVG_DIST
		double avgDist = 0.0;
		double avgAdjIF = 0.0;
		int avgAdjCount = 0;
		for(Constraint con : lstCons){						
			
			con.setIF(con.getIF()/avgIF); //normalize IF by avgIF
			avgDist += (1.0 / Math.pow(con.getIF(),convertFactor));
			
			totalIF += con.getIF();
			if (con.getIF() > maxIF){
				maxIF = con.getIF();
			}
			
			if (Math.abs(con.getPos1() - con.getPos2()) == 1 && idToChr.get(con.getPos1()) == idToChr.get(con.getPos2())) {
				avgAdjCount++;
				avgAdjIF += con.getIF();
			}
		}
		avgDist /= lstCons.size();
		avgAdjIF /= avgAdjCount;
		
		
		//August 18 2016
		maxIF = Math.min(avgAdjIF, maxIF);
		
		Collections.sort(lstCons);
		addAdjacentContacts(avgAdjIF);
		//addNonContact();
		
		System.out.println("Number of constraints: " + lstCons.size());
		double max = 0;
		for(Constraint con : lstCons){			
			con.setDist(AVG_DIST / (Math.pow(con.getIF(),convertFactor) * avgDist ));
			if (AVG_DIST / (Math.pow(con.getIF(),convertFactor) * avgDist ) > max){
				max = AVG_DIST / (Math.pow(con.getIF(),convertFactor) * avgDist) ;
			}
		}
		
		System.out.println("Max distance is: " + max);
		
		
//		BufferedReader br = new BufferedReader(new FileReader("/Users/Tuan/workspace/DataAndEvaluation/output/SyntheticYeast/GMDS/poisson_chainDres25_dist20_4.txt"));
//		String line;
//		String[] st;
//		
//		double[][] dist = new double[n][n];
//		
//		while((line = br.readLine()) != null){
//			st = line.split("\\s+");
//			int x = Integer.parseInt(st[0]);
//			int y = Integer.parseInt(st[1]);
//			double z = Double.parseDouble(st[2]);
//			
//			dist[x-1][y-1] = z;
//		}
//		br.close();
//		
////		List<Constraint> tmp = helper.readContactList("/Users/Tuan/workspace/DataAndEvaluation/output/SyntheticYeast/GMDS/poisson_chainDres25_dist01_1.txt", 
////				new ArrayList<Integer>(),contactThres);
//		
//		
////		Collections.sort(lstCons);
//		
//		double[] lst1 = new double[lstCons.size()];
//		double[] lst2 = new double[lstCons.size()];
//		int t = 0;
//		for(int i = 0; i < lstCons.size(); i++){
//			Constraint con1 = lstCons.get(i);
//			
//			
//				//System.out.println();
//				
//			lst1[i] = con1.getDist();
//			lst2[i] = dist[con1.getPos1()][con1.getPos2()];
//			
//		}
//		
//		double rmsd = CalRMSD.rmse(lst2, lst1);
//		
//		System.out.println(Evaluate.calSpearmanCorrelation(lst1, lst2) + ", " + rmsd);
//		System.out.println();
//		
//		//to test
////		PrintWriter pw = new PrintWriter("interactions1.txt");
////		for(Constraint con : lstCons){			
////			con.setDist(AVG_DIST / (Math.pow(con.getIF(),convertFactor) * avgDist ));
////			pw.printf("%.2f,%.2f,%d,%d\n",con.getIF(),con.getDist(),idToChr.get(con.getPos1()), idToChr.get(con.getPos2()));
////		}		
////		pw.close();

		
	}
	
	
	//add adjacent contacts if not exist
	
	private void addAdjacentContacts(double IF){
		int id;
		ArrayList<Constraint> ctList = new ArrayList<Constraint>();
		Constraint ct;
		for(int i = 0; i < n - 1; i++){
			if (idToChr.get(i) == idToChr.get(i+ 1) ){	
				
				ct = new Constraint(i,i + 1,IF);				
				id = Collections.binarySearch(lstCons, ct);
				if (id < 0){
					ctList.add(ct);
				}else{					
					if (lstCons.get(id).getIF() < IF){
						lstCons.remove(id);
						ctList.add(ct);
					}
					
					
				}
			}
			
		}
		
		lstCons.addAll(ctList);
	}

	private void readParameters(String paraFile)throws Exception{
		File file = new File(paraFile);
		FileReader fr = null;
		BufferedReader br = null;
		String ln;
		String[] st;

		Pattern splitRegex = Pattern.compile("[=\\s#]+");
		
		try{
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			
			while((ln = br.readLine()) != null){
				if (ln.startsWith("#")){
					continue;
				}
				
				st = splitRegex.split(ln);
				if (st[0].equalsIgnoreCase(Constants.NUM_KEY)){
					NUM = Integer.parseInt(st[1]);
					
				}else if (st[0].equalsIgnoreCase(Constants.CONVERT_FACTOR_KEY)) {
					convertFactor = Double.parseDouble(st[1]);
					
				}else if (st[0].equalsIgnoreCase(Constants.OUTPUT_FOLDER_KEY)){
					OUTPUT_FOLDER = "";
					for (int i = 1; i < st.length; i++){
						OUTPUT_FOLDER += st[i];
					}
					
				}else if (st[0].equalsIgnoreCase(Constants.INPUT_FILE_KEY)){
					INPUT_FILE = "";
					for (int i = 1; i < st.length; i++){
						INPUT_FILE += st[i];
					}
					String[] tmp = INPUT_FILE.split("[\\/ \\. \\\\]");
					if (INPUT_FILE.contains(".")){
						FILE_PREFIX = tmp[tmp.length - 2];
					}else{
						FILE_PREFIX = tmp[tmp.length - 1];
					}
				
				}else if (st[0].equalsIgnoreCase(Constants.VERBOSE_KEY)){
					VERBOSE = Boolean.parseBoolean(st[1]);				
				
				}else if (st[0].equalsIgnoreCase(Constants.LEARNING_RATE_KEY)){
					LEARNING_RATE = Double.parseDouble(st[1]);
					
				}else if (st[0].equalsIgnoreCase(Constants.MAX_ITERATION_KEY)){
					MAX_ITER = Integer.parseInt(st[1]);
				
				}else if (st[0].equalsIgnoreCase(Constants.CHR_LENGTH_KEY)) {
					String[] lens = st[1].trim().split(",");
					chrLens = new int[lens.length];
					for(int i = 0; i < lens.length; i++){
						chrLens[i] = Integer.parseInt(lens[i]);
					}
				}else if (st[0].equalsIgnoreCase(Constants.THRESHOLD_KEY)){
					contactThres = Double.parseDouble(st[1]);
				}
			}
			

		}catch(Exception e){
			
			e.printStackTrace();
			throw e;			
			
		}finally{
			if (br != null){
				br.close();
			}
			if (fr != null){
				fr.close();
			}
		}

	}
	
	/**
	 * Initialize variables before running optimization
	 */
	private void initialize(){
				
		str = new double[n * 3];
		
		//get the number of processor available
		int numOfcores = Runtime.getRuntime().availableProcessors();
		
		if (numOfcores == 0){
			numOfcores = 2;// default number when this parameter cannot be detected
		}
		//each core can take care of 2 threads		
		//limit number of threads to avoid excessive communication cost
		numOfcores = Math.min(numOfcores * 2 , Constants.MAX_NUM_THREAD);
		
		String inEclipseStr = System.getProperty("runInEclipse");
		if ("true".equalsIgnoreCase(inEclipseStr)){
			numOfcores = 1;
		}
		
		//August,18, test to see if adjacent distance is better
		numOfcores = 1;
		
		System.out.println("Number of processors:" + numOfcores);
		//divide the set of points into equal subsets, each will be processed by one processor (thread)

		helper.divideDataSet(lstCons.size(), numOfcores, lstSubDataSetId);
		
		
		File outputFolder = new File(OUTPUT_FOLDER);
		if (!outputFolder.exists()) {
			outputFolder.mkdir();
		}
		
//		double[][] interMap = computeInterMap(lstCons, idToChr,chrLens.length);
//		//interMap = MakeInputLabData.standardNorm(interMap);
//		try{
//			PrintWriter pw = new PrintWriter("interMap.txt");
//			for(int i = 0; i < interMap.length; i++){
//				for(int j = 0; j < interMap.length - 1; j++){
//					pw.print(interMap[i][j]+",");
//				}
//				pw.print(interMap[i][interMap.length - 1]);
//				if (i < interMap.length - 1){
//					pw.println();
//				}
//			}
//			pw.close();
//		}catch(Exception ex){
//			
//		}
	}
	
//	private double[][] computeInterMap(List<Constraint> lst, Map<Integer,Integer> map, int n){
//		
//		double[][] a = new double[n][n];
//		int totalInter = 0, countInter = 0, totalIntra = 0, countIntra = 0;
//		double thres = 1.0;
//		
//		//SummaryStatistics interStat = new SummaryStatistics();
//		//SummaryStatistics intraStat = new SummaryStatistics();
//		for(Constraint ct : lst){
//			
//			if ( map.get(ct.getPos1()) == map.get(ct.getPos2()) ){
//				//intraStat.addValue(ct.getIF());
//				
//				totalIntra++;
//				if (ct.getIF() > thres){
//					countIntra++;
//				}
//			}else{
//				//interStat.addValue(ct.getIF());
//				
//				totalInter++;
//				if (ct.getIF() < thres){
//					countInter++;
//				}
//			}
//			
//			//a[map.get(ct.getPos1())][map.get(ct.getPos2())] += ct.getIF();
//			//a[map.get(ct.getPos2())][map.get(ct.getPos1())] += ct.getIF();
//		}
//		
//		System.out.printf("Inter less than 0.65: %.3f, \n Intra more than 0.65: %.3f", countInter*100.0/totalInter, countIntra*100.0/totalIntra);
//		/*
//		System.out.printf("Inter, mean: %.2f, min: %.2f, max: %.2f, sd: %.2f", interStat.getMean(), interStat.getMin(), interStat.getMax(), interStat.getStandardDeviation());
//		System.out.println();
//		System.out.printf("Intra, mean: %.2f, min: %.2f, max: %.2f, sd: %.2f", intraStat.getMean(), intraStat.getMin(), intraStat.getMax(), intraStat.getStandardDeviation());
//		*/
//		
//		return a;
//	}
	
	
	/**
	 * Initialize the genome structure, adjacent points are initialized to be closer together than the others
	 */
	private void initializeStructure() throws Exception{
		
		double chrX=0,chrY=0,chrZ=0,size = 0.1;
		
		for(int i = 0; i < n; i++){
			
//			//reset starting point for every chromosome
//			if (i == 0 || idToChr.get(i) != idToChr.get(i - 1)){			
				chrX = Math.random() * size;
				chrY = Math.random() * size;
				chrZ = Math.random() * size;
				
//			}else {
//				//extend in X,Y,Z coordinate
//				chrX += size * (Math.random() - 0.5);
//				chrY += size * (Math.random() - 0.5);
//				chrZ += size * (Math.random() - 0.5);
//			}

			str[i * 3] = chrX;
			str[i * 3 + 1] = chrY;
			str[i * 3 + 2] = chrZ;		

		}


	}

	public void generateStructure() throws Exception{
		
		readParameters(parameterFile);
		
		if (convertFactor == -1){
			convertFactor = 0.1;
			double cor, minCor = 1.0;
			double bestConvertFactor = -1;
			for(; convertFactor < 3.0; convertFactor += 0.1){
				cor = run(convertFactor + "");
				if (minCor > cor){
					minCor = cor;
					bestConvertFactor = convertFactor;
				}
			}
			
			PrintWriter pw = new PrintWriter(OUTPUT_FOLDER + "/" + "best_alpha_log.txt");
			pw.printf("\n\nBest convert factor: %.2f, pick models generated using this convert factor as your final models \n", bestConvertFactor);
			System.out.printf("\n\nBest convert factor: %.2f, pick models generated using this convert factor as your final models \n", bestConvertFactor);
			pw.close();
			
			//generate models after searching for best factor
			convertFactor = bestConvertFactor;
			run();
			
			
		}else{
			run();
		}
	}
	
	public double run(String... cFactor) throws Exception{
		String fileName;
		
		readInput();
		initialize();
		
		String logFileName = "";		
		PrintWriter logPW = null;
		double rmsd,cor,corDist;
		double avgRMSD = 0,avgCor = 0, avgCorDist = 0;
		boolean isOutput = false;
		
		int run_nbr = NUM;
		//if search for best alpha, just run 3 times for each alpha candidate
		if (cFactor != null && cFactor.length > 0){
			run_nbr = 3;
		}
		
		for(int i = 0; i < run_nbr; i++) {		
			initializeStructure();
			
			GradientAscent gradientAscent = new GradientAscent(this, str, VERBOSE);
			if (LEARNING_RATE != 0){
				gradientAscent.setInitialLearingRate(LEARNING_RATE);
			}
			
			gradientAscent.performGradientAscent(MAX_ITER);
			
			String currentTimeMillis = System.currentTimeMillis() + "";
			
			if (cFactor != null && cFactor.length > 0){
				fileName = FILE_PREFIX + "_" + currentTimeMillis + "_" + cFactor[0] + ".pdb" ;
				isOutput = false;
			}else{
				fileName = FILE_PREFIX + "_" + currentTimeMillis + ".pdb" ;
				
				isOutput = true;
			}
			
			//print out log			
			if (isOutput){
				
				helper.writeStructure(OUTPUT_FOLDER + "/" + fileName,str, idToChr, "Tuan test");
				
				logFileName =  FILE_PREFIX + "_log_" + currentTimeMillis + ".txt";
				logPW = new PrintWriter(OUTPUT_FOLDER + "/" + logFileName);
				logPW.println("Input file: " + INPUT_FILE);
				logPW.println("Convert factor: " + convertFactor);
				logPW.println("Learning rate: " + LEARNING_RATE);
				if (chrLens != null){
					logPW.print("Chromosome lengths:");
					for(int k = 0; k < chrLens.length; k++){
						logPW.print(chrLens[k] + " ");
					}
					logPW.println();
				}
	
				logPW.flush();
			
			}
			
			rmsd = CalRMSD.rmse(str, lstCons);
			interval = 0;
			cor = CalRMSD.correlationIFvsDist(str, lstCons, interval);
			corDist = CalRMSD.correlationWishDistvsDist(str, lstCons, interval);
			
			avgRMSD += rmsd;
			avgCor += cor;
			avgCorDist += corDist;
			
			if (isOutput){
				logPW.println("RMSE: " + rmsd);
				logPW.println("Spearman correlation IFs vs. Reconstructed Dist: " + cor);
				logPW.println("Spearman correlation WishDist vs. Reconstructed Dist: " + corDist);
				logPW.flush();
				logPW.close();
				
			}
			
//			System.out.println("RMSE: " + rmsd);
//			System.out.println("Spearman correlation IFs vs. Reconstructed Dist: " + cor);
//			System.out.println("Spearman correlation WishDist vs. Reconstructed Dist: " + corDist);
		}
				
		avgRMSD /= run_nbr;
		avgCor /= run_nbr;
		avgCorDist /= run_nbr;
		
		PrintWriter pw = new PrintWriter(OUTPUT_FOLDER + "/" + FILE_PREFIX + "_log.txt");
		pw.println("Input file: " + INPUT_FILE);
		pw.println("Convert factor: " + convertFactor);
		pw.println("Learning rate: " + LEARNING_RATE);
		if (chrLens != null){
			pw.print("Chromosome lengths:");
			for(int k = 0; k < chrLens.length; k++){
				pw.print(chrLens[k] + " ");
			}
			pw.println();
		}

		pw.println("AVG RMSE: " + avgRMSD);
		pw.println("AVG Spearman correlation IFs vs. Reconstructed Dist: " + avgCor);
		pw.println("AVG Spearman correlation WishDist vs. Reconstructed Dist: " + avgCorDist);
		pw.close();
		
		System.out.println("AVG RMSE: " + avgRMSD);
		System.out.println("AVG Spearman correlation IFs vs. Reconstructed Dist: " + avgCor);
		System.out.println("AVG Spearman correlation WishDist vs. Reconstructed Dist: " + avgCorDist);

		return avgCor;
	}
	
	
	
	/**
	 * Calculate objective function and gradient
	 */
	@Override
	public double calGradientAndObjective(double[] x, double[] der)
			throws InterruptedException {
		
		double cost = 0.0;
		
		GradientCaculator[] gradCalculator = new GradientCaculator[lstSubDataSetId.size()];
		
		gradCalculator[0] = new GradientCaculator(x, 0 , lstSubDataSetId.get(0), der != null);
		for(int i = 1; i < gradCalculator.length; i++){
			gradCalculator[i] = new GradientCaculator(x, lstSubDataSetId.get(i - 1) + 1, lstSubDataSetId.get(i), der != null);
		}
		
		//start threads
		for(int i = 0; i < gradCalculator.length; i++){
			gradCalculator[i].start();
		}
		
		//wait for all threads to finish
		try {
			for(int i=0; i< gradCalculator.length; i++){
				gradCalculator[i].join();
			}			
		} catch (InterruptedException e) {			
			e.printStackTrace();
			throw e;
		}

		//aggregate the cost
		for(int i = 0; i < gradCalculator.length; i++){
			cost += gradCalculator[i].getCost();
		}
	
		if (der != null){
			//aggregate the gradient
			for(int i = 0; i < der.length; i++){
				der[i] = 0;
				for(int k = 0; k < gradCalculator.length; k++){
					der[i] += gradCalculator[k].getChange()[i];
				} 				
			}			
		}		
		
		return cost;
		
	}

	/**
	 * calculate objective function only
	 */
	@Override
	public double calObjective(double[] x) throws InterruptedException {
		
		return calGradientAndObjective(x,null);
	}
	
	
	/**
	 * This class is used to calculate the gradient for a subset of datapoints
	 * in a single threaded program, all data points are i = 1 .. n and j = i+1 ... n
	 * to calculate the gradient in parallel, one thread will be in charged of calculating for i = begin ... end
	 * 
	 * for any modification of the objective function, this function will need to be modified accordingly
	 * @author Tuan
	 *
	 */
	class GradientCaculator extends Thread{
		//the first index to calculate the gradient
		private int beg;
		//the last index to calculate the gradient
		private int end;
		//gradient to be returned
		double[] change;
		//objective function
		double cost=0;
		//indicate if calculation for gradient is needed
		boolean isGradientNeeded;
		
		double[] structure;
		
		GradientCaculator(double[] str, int b, int e, boolean isGradient){
			this.beg = b;
			this.end = e;		
			this.structure = str;
			this.isGradientNeeded = isGradient;
			
			if (isGradientNeeded){
				this.change = new double[n * 3];
			}
		}
		
		public void run(){
			double dist,x,ez,tmp,z,ifr,tanh;
			int i,j;
			Constraint con;
			for(int k = beg; k <= end; k ++){
				
				con = lstCons.get(k);
				
				i = con.getPos1();
				j = con.getPos2();
				dist = con.getDist();
				
				
				ifr = con.getIF();
				
				if (ifr <= 0) continue;
					
//				ifr = ifr * ifr;
				
				//if (Math.abs(lstPos.get(i) - lstPos.get(j)) == 1 && idToChr.get(i) == idToChr.get(j)){
				if (Math.abs(i - j) == 1 && idToChr.get(i) == idToChr.get(j)){
					ifr = 1.0 * maxIF;
				}				
//				else {
//					ifr = 1.0;
//				}
				
				x = Math.sqrt(helper.calEuclidianDist(structure[i * 3], structure[i * 3 + 1], 
						structure[i * 3 + 2], structure[j * 3], structure[j * 3 + 1], structure[j * 3 + 2]));
				
				//WIDE_CURVE = dist;
				int pow = 2;
				z = Math.pow(x - dist,pow);		
				tmp = - ifr * WIDE_CURVE * pow * Math.pow(x - dist,pow - 1) / ((z + WIDE_CURVE)*(z + WIDE_CURVE) * x); 


				cost += ifr * WIDE_CURVE / (z + WIDE_CURVE);
				
				if ( isGradientNeeded){					 
					
					change[i * 3] += tmp * (structure[i * 3] - structure[j * 3]);
					change[i * 3 + 1] += tmp * (structure[i * 3 + 1] - structure[j * 3 + 1]);
					change[i * 3 + 2] += tmp * (structure[i * 3 + 2] - structure[j * 3 + 2]);
					
					change[j * 3] += tmp * (structure[j * 3] - structure[i * 3]);
					change[j * 3 + 1] += tmp * (structure[j * 3 + 1] - structure[i * 3 + 1]);
					change[j * 3 + 2] += tmp * (structure[j * 3 + 2] - structure[i * 3 + 2]);
					
				}

				
			}			
		}

		public double[] getChange() {
			return change;
		}
		
		public double getCost() {
			return cost;
		}
		
	}

	
	public static void main(String[] args) throws Exception{
		
		String paraFile = "parameters_gm12878_chr1_1mb.txt";
		if (args != null && args.length >= 1){
			paraFile = args[0];
		}
		//check if the parameter file exists in the current directory
		File file = new File(paraFile);
		if (!file.exists()){
			throw new Exception("The parameter file " + paraFile + "cannot be found in the current directory");
		}
		
		StructureGeneratorLorentz generator = new StructureGeneratorLorentz(paraFile);
		generator.generateStructure();
		

	}

	

}
