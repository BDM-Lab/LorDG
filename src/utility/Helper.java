package utility;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;
import java.util.regex.Pattern;

import valueObject.Constraint;

public class Helper {	
	private DecimalFormat df2 = new DecimalFormat("0.00");
	
	private static Helper helper = new Helper();
	
	private Helper(){		
	}
	
	public static Helper getHelperInstance(){
		if (helper == null){
			synchronized(Helper.class){
				if (helper == null){
					helper = new Helper();
				}
			}
		}
		
		return helper;
	}
	


	/**
	 * Read contact matrix file
	 * @param fileName
	 * @return
	 * @throws Exception
	 */
	public List<Constraint> readContactMatrixAsList(String fileName) throws Exception{
		File file = new File(fileName);
		FileReader fr = null;
		BufferedReader br = null;
		ArrayList<Constraint> lst = new ArrayList<Constraint>();
		
		try{
			Pattern splitRegex = Pattern.compile("\\s+");
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			
			String ln;
			String[] st;			
			int count = 0;		
			double d;
			while((ln = br.readLine()) != null){
				
				st = splitRegex.split(ln.trim());
								
				for(int i = count + 1; i < st.length; i++){
					d = Double.parseDouble(st[i]);
					if (!Double.isInfinite(d) && !Double.isNaN(d) && d > 0.0){
						lst.add(new Constraint(count,i,d));
					}
				}				
				count++;
			}
		}catch(Exception ex){
			ex.printStackTrace();
			throw ex;
		}finally{
			if (br != null){
				br.close();
			}
			if (fr != null){
				fr.close();
			}
		}		
		return lst;
	}

	/**
	 * Read contact list file, each line is a contact of the form: pos1 pos1 IF
	 * @param fileName
	 * @return
	 * @throws Exception
	 */
	public List<Constraint> readContactList(String fileName, List<Integer> lstPos,double...thres) throws Exception{
		File file = new File(fileName);
		FileReader fr = null;
		BufferedReader br = null;
		ArrayList<Constraint> lst = new ArrayList<Constraint>();
		Set<Integer> setPos = new HashSet<Integer>();
		double thr = thres.length == 0 ? 0.0 : thres[0];
		
		try{
			Pattern splitRegex = Pattern.compile("[:\\s]+");
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			
			String ln;
			String[] st;			

			StringBuilder sb = new StringBuilder();
			int count = 0;
			int x,y,nbr = -1; // number of elements in one line
			double f;
			int progress = 0;
			while((ln = br.readLine()) != null){
				
				if (ln.startsWith("#") || ln.trim().length() == 0){
					continue;
				}
				if (nbr == -1){
					nbr = splitRegex.split(ln).length;
				}
				
				//read every a thoudsand lines and split at once
				sb.append(ln).append(" ");
				count++;
				
				if (count == 200000){
					count = 0;
					st = splitRegex.split(sb.toString());
					sb = new StringBuilder();
					
					if (st.length % nbr != 0){
						throw new Exception("There is a line that doesn't contain exactly 3 numbers");
					}
					//each line contains 3 numbers
					for(int i = 0; i < st.length / nbr; i++){
						x = Integer.parseInt(st[i * nbr + 0]);
						
						//position2
						y = Integer.parseInt(st[i * nbr + 1]);
						
						//interaction frequency
						f = Double.parseDouble(st[i * nbr + 2]);
						
						if (x != y && !Double.isNaN(f) && f > thr){
							lst.add(new Constraint(x, y, f)) ;
							setPos.add(x);
							setPos.add(y);
						}
						 
					}
					progress++;
					System.out.println(progress * 200000 + " input lines have been read !");
				}		
			}
			
			//if sb is not empty
			if (sb.toString().trim().length() > 0){
				st = splitRegex.split(sb.toString());
				//sb = new StringBuilder();
				
				if (st.length % nbr != 0){
					throw new Exception("There is a line that doesn't contain exactly 3 numbers");
				}
				//each line contains 3 numbers
				for(int i = 0; i < st.length / nbr; i++){
					x = Integer.parseInt(st[i * nbr + 0]);
					
					//position2
					y = Integer.parseInt(st[i * nbr + 1]);
					
					//interaction frequency
					f = Double.parseDouble(st[i * nbr + 2]);
					
					//keeping absolute positions, so that later they can be recovered from indices
					
					if (x != y && !Double.isNaN(f) && f > thr){
						lst.add(new Constraint(x, y, f)) ;
						setPos.add(x);
						setPos.add(y);
					}
				}
				
			}
			
			
			lstPos.addAll(setPos);
			
				
			Collections.sort(lstPos);
			
		}catch(Exception ex){
			ex.printStackTrace();
			throw ex;
		}finally{
			if (br != null){
				br.close();
			}
			if (fr != null){
				fr.close();
			}
		}		
		return lst;
	}

	/**
	 * Read contact matrix file
	 * @param fileName
	 * @return
	 * @throws Exception
	 */
	public int determineNbrOfPoints(String fileName) throws Exception{
		File file = new File(fileName);
		Scanner scan = null;
		
		try{
			Pattern splitRegex = Pattern.compile("\\s+");
			scan = new Scanner(file);
			
			String ln;
			String[] st;			
			
			ln = scan.nextLine();			
			st = splitRegex.split(ln.trim());
			return st.length;
			
		}catch(Exception ex){
			ex.printStackTrace();
			throw ex;
		}finally{
			if (scan != null){
				scan.close();
			}
		}		
		
	}

	/**
	 * Make distance matrix from contact matrix
	 * @param contMT
	 * @return
	 */
	public double[][] makeDistanceMatrix(double[][] contMT, double factor){
		double[][] distMT = new double[contMT.length][contMT.length];
		for (int i = 0; i < distMT.length; i++){
			for (int j = 0; j < distMT.length; j++){
				if (i != j){
					distMT[i][j] = 1.0 / Math.pow(contMT[i][j],factor);
				}
				
			}
		}		
		return distMT;
	}
	

	
	/**
	 * 
	 * @param total: number of constraints
	 * @param k: number of processors
	 * @return an arraylist contains ending index for each processor
	 */
	public void divideDataSet(int total, int k, ArrayList<Integer> lstSubDataSetId){
		
		int size = total / k;		
		for(int i = 1; i < k; i++){
			lstSubDataSetId.add(i * size);
		}
		if (total % k != 0 || k == 1){
			lstSubDataSetId.add(total - 1);
		}
	}
	
	
	/**
	 * 
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 * @return square euclidean distance of (x1,y1) and (x2,y2) 
	 */
	public double calEuclidianDist(double x1, double y1, double z1, double x2, double y2, double z2){
		return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
	}
	
	/**
	 * log function to take care of the case x <= 0	
	 * @param x
	 * @return
	 */
	public double log(double x){
		if (x <= 0){
			return Integer.MIN_VALUE;
		}else{
			return Math.log(x);
		}
	}
	
	/**
	 * 
	 * @param 
	 * @return hyperbolic function tanh(x)
	 */
	public double tanh(double x){		
		
		//this is an approximation for the tanh function
		if (x > 15){
			return 1.0;
		}else if (x < -15){
			return -1;
		}
		double ex = Math.pow(Math.E,x);
		double emx = 1/ex;
		double y = (ex - emx)/(ex + emx) ;		
		
		return y;

	}
	
	
	/**
	 * Zoom the structure (str) by the factor
	 * @param str
	 * @param fator
	 */
	public void zoomStructure(double[] str, double factor){
		if (str == null){
			return;
		}
		for(int i = 0; i < str.length; i++){
			str[i] *= factor;
		}
	}
	
	/**
	 * Read a PDB file and return a coordinate matrix, general one
	 * @param fileName
	 * @return
	 */
	public double[][] loadPDBStructure(String fileName) throws Exception{
		File file = new File(fileName);
		FileReader fr = null;
		BufferedReader br = null;
		double[][] coor;
		ArrayList<Double> lstX = new ArrayList<Double>();
		ArrayList<Double> lstY = new ArrayList<Double>();
		ArrayList<Double> lstZ = new ArrayList<Double>();
		try{
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			

			String ln;
			String[] st;
			Pattern splitRegex = Pattern.compile("[\\s]+");
			while((ln = br.readLine()) != null){
				if (! ln.startsWith("ATOM")){
					continue;
				}
				st = splitRegex.split(ln);				
				if (st.length > 7){
					lstX.add(Double.parseDouble(st[5]));
					lstY.add(Double.parseDouble(st[6]));
					lstZ.add(Double.parseDouble(st[7]));
				}
			}		
			
			
		}catch(Exception ex){
			ex.printStackTrace();
			throw ex;			
			
		}finally{
			if (br != null){
				br.close();
			}
			if (fr != null){
				fr.close();
			}
		}
		
		coor = new double[lstX.size()][3];
		for(int i = 0; i < lstX.size(); i++){
			coor[i][0] = lstX.get(i);
			coor[i][1] = lstY.get(i);
			coor[i][2] = lstZ.get(i);
		}
		
		return coor;
	}

	/**
	 * Convert a coordinate matrix into a matrix ready for output to write out PDB file
	 * @param a
	 * @return
	 */
	public double[] makeMatrixForPDBOutput(double[][] a){
		double[] str = new double[a.length * 3];
		for(int i = 0; i < a.length; i++){
			str[i * 3 + 0] = a[i][0];
			str[i * 3 + 1] = a[i][1];
			str[i * 3 + 2] = a[i][2];
		}
		
		return str;
	}

	public double[][] convertArrayVariableToStr(double[] a){
		int n = a.length / 3;
		double[][] str = new double[n][3];
		for(int i = 0; i < n; i++){
			str[i][0] = a[i * 3];
			str[i][1] = a[i * 3 + 1];
			str[i][2] = a[i * 3 + 2];
		}
		
		return str;
	}
	
	/**
	 * To extract one chromosome from the whole genome structure
	 * @param a: coordinates of all points in the genome
	 * @param start: starting index
	 * @param end: ending index (inclusive)
	 * @return coordinates of points of the chromosome
	 */
	public double[][] extractChr(double[][] a, int start, int end){
		
		int n = end - start + 1;
		double[][] b = new double[n][n];
		for(int i = 0; i < n; i++){
			b[i][0] = a[i + start][0];
			b[i][1] = a[i + start][1];
			b[i][2] = a[i + start][2];
		}
		
		return b;		
	}

	/**
	 * Please refer to the pdb format file for detail of each column
	 * @param pathFilename: output file name
	 * @param structure: every 3 consecutive points is one point in 3D
	 * @param idToChr: to identify if 2 fragments belong to the same chromosome
	 * @param header for the pdb file
	 */
	public void writeStructure(String pathFilename, double[] structure, HashMap<Integer,Integer> idToChr, String header, boolean... isTranslate) throws IOException{

		//number of fragments
		int n = structure.length / 3;
		
		if (idToChr == null){
			idToChr = new HashMap<Integer,Integer>();
			//if idToChr is null, make the whole as one chromosome
			for(int i = 0; i < n; i++){				
				idToChr.put(i, 0);				
			}
		}
		/////////
		
		if (isTranslate != null && isTranslate.length > 0 && isTranslate[0]) {
			//  make the minimum x-, y-, and z-coordinates of the structure equal to one (1)
			double translationX = Double.MAX_VALUE;
			double translationY = Double.MAX_VALUE;
			double translationZ = Double.MAX_VALUE;
			
			//  find the minimum x,y,z coordinate values
			for(int i = 0; i < n; i++){
				if(structure[i * 3] < translationX){
					translationX = structure[i];
				}
				if(structure[i * 3 + 1] < translationY){
					translationY = structure[i + 1];
				}
				if(structure[i * 3 + 2] < translationZ){
					translationZ = structure[i + 2];
				}
			}
			
			//  subtract one (1.0) to each of the translations to leave the minimum point at coordinate one (1.0) after the translation
			translationX -= 1;
			translationY -= 1;
			translationZ -= 1;
			
			//  translate all of the points in the structure
			for(int i = 0; i < n; i++){
				structure[i * 3] -= translationX;
				structure[i * 3 + 1] -= translationY;
				structure[i * 3 + 2] -= translationZ;
			}
			
		}
		

		//  write the structure to file
		try{						
			PrintWriter pw = new PrintWriter(pathFilename);
			
			
			//  write the header block
			pw.println(header.toUpperCase());

			int atomSerial = 1;
			int resName = 1;
			
			String line;
			for(int i = 0; i < n; i++){
				line = "";
				line += getAtomString("ATOM");
				line += getSerialString(atomSerial + "") + " ";//12th space
				atomSerial++; // increase atom id
				
				line += getNameString("CA");
				line += getAltLocString(" ");
				//line += getResNameString(resName + "") + " ";//21st space
				line += getResNameString("MET" + "") + " ";//21st space
				
				
				line += getChainIDString((char)(idToChr.get(i) + 'A' ) + "");
				line += getResSeqString(resName + "");
				//line += getResSeqString("MET" + "");
				
				//if (atomSerial % 10 == 0){ // 10 atom in the same residue name
					resName++;
				//}
				
				line += " "; //27th space (iCode)
				line += "   "; //28,29,30th spaces
				
				line += getXString(structure[i * 3]);
				line += getYString(structure[i * 3 + 1]);
				line += getZString(structure[i * 3 + 2]);
				
				line += getOccupancyString(0.2);
				line += getTempFactorString(10.0);
				
				pw.println(line);				
				pw.flush();				
			}
			
			for(int i=1; i<atomSerial-1; i++){
				if(idToChr.get(i-1) == idToChr.get(i)){
					line = "";
				}else{
					line = "#";
				}
				
				line += getConnectString("CONECT");
				line += getSerial1String(i+"");
				line += getSerial2String((i+1) + "");
				pw.println(line);
				
				pw.flush();					
			}
			
			pw.println("END");			
			pw.flush();					
			pw.close();
		}catch (IOException e) {
			System.err.println("Error: could not output file: " + pathFilename);
			e.printStackTrace();
			throw e;
		}
	}
	

	
	//to format the string to follow the format of pdb file format
	private String getAtomString(String st){
		//1-6
		int length = 6;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in atom name, length exceeds " +  length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = st + " ";
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getSerialString(String st){
		//7-11
		int length = 5;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in serial, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " " + st;
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getNameString(String st){
		//13-16
		int length = 4;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in name, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " " +st;
		}
		
		return st;
	}

	//to format the string to follow the format of pdb file format
	private String getAltLocString(String st){
		//17
		int length = 1;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in alt loc, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = st + " ";
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getResNameString(String st){
		//18-20
		int length = 3;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in res name, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " " + st;
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getChainIDString(String st){
		//22
		int length = 1;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in chain id, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = st + " ";
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getResSeqString(String st){
		//23-26
		int length = 4;
		int currentLength = st.length();
		if (currentLength > length){
			//System.err.println("Error in res seq, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = st + " ";
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getXString(double x){
		//31-38
		int length = 8;
		String st = "";
		if (x > 10000 || x < -1000){
			st = String.format("%8.2f",x);
		}else{
			st = String.format("%8.3f",x);
		}
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in X, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " "+st;
		}
		
		return st;
	}

	//to format the string to follow the format of pdb file format
	private String getYString(double x){
		//39-46
		int length = 8;
		String st = "";
		if (x > 10000 || x < -1000){
			st = String.format("%8.2f",x);
		}else{
			st = String.format("%8.3f",x);
		}
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in Y, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " "+st;
		}
		
		return st;
	}
	//to format the string to follow the format of pdb file format
	private String getZString(double x){
		//47-54
		int length = 8;
		String st = "";
		if (x > 10000 || x < -1000){
			st = String.format("%8.2f",x);
		}else{
			st = String.format("%8.3f",x);
		}
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in Z, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " "+st;
		}
		
		return st;
	}
	//to format the string to follow the format of pdb file format
	private String getOccupancyString(double x){
		//55-60
		int length = 6;
		String st = df2.format(x);
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in occupancy, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " "+st;
		}
		
		return st;
	}
	//to format the string to follow the format of pdb file format
	private String getTempFactorString(double x){
		//61-66
		int length = 6;
		String st = df2.format(x);
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in tempFactor, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " "+st;
		}		
		return st;
	}
	//to format the string to follow the format of pdb file format
	private String getConnectString(String st){
		//1-6
		int length = 6;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in connect, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = st + " ";
		}
		
		return st;
	}
	//to format the string to follow the format of pdb file format
	private String getSerial1String(String st){
		//7-11
		int length = 5;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in serial 1, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " " + st;
		}		
		return st;
	}
	//to format the string to follow the format of pdb file format	
	private String getSerial2String(String st){
		//12-16
		int length = 5;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in serial 2, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " " + st;
		}		
		return st;
	}


	/**
	* Read contact matrix file
	* @param fileName
	* @return
	* @throws Exception
	*/
	public double[][] readContactMatrix(String fileName) throws Exception{
		File file = new File(fileName);
		Scanner scan = null;
		double[][] a = null;
		try{
			Pattern splitRegex = Pattern.compile("\\s+");
			scan = new Scanner(file);
			String ln = scan.nextLine();
			String[] st;
			st = splitRegex.split(ln.trim());
			int n = st.length;
			a = new double[n][n];
			for(int i = 0; i < n; i++){
				a[0][i] = Double.parseDouble(st[i]);
			}
			int count = 1;
			
			while(scan.hasNextLine()){
				ln = scan.nextLine();			
				st = splitRegex.split(ln.trim());
				
				for(int i = 0; i < n; i++){
					a[count][i] = Double.parseDouble(st[i]);
				}
				
				count++;
			}
		}catch(Exception ex){
			ex.printStackTrace();
			throw ex;
		}finally{
			if (scan != null){
				scan.close();
			}
		}
		
		return a;
	}

}



