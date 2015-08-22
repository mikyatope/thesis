import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Collections;
import java.util.Map;

/*
 * Calculate Site Frequency Spectrum for a set of aligned sequences
 */
class SFS {

	//general fields
	private int    fixNum;
	private String classSite;
	
	//per region fields
	private int    countSites;
	private int    countDiscarded;
	private int    countSegregating;
	private int[]  sfs;
	private int[]  ufold_sfs;
	private int[]  countCases;
	private int    Dpolarized;
	private int    Dout1;
	private int    Dout2;

	//per column fields
	private List<Character> column;
	private String          outgroups;


	/**
	 * Constructors
	 */
	SFS(int getFixNum){
		fixNum    = getFixNum;
		classSite = "total";
	}
	SFS(int getFixNum, String getClassSite){
		fixNum    = getFixNum;
		classSite = getClassSite;
	}

	
	/**
	 * Initialise region
	 */
	void setRegion(){
		sfs              = new int[(fixNum/2)+1];   //add +1 for alleles segregating at 50%/50%
		ufold_sfs        = new int[fixNum];
		countSites       = 0;
		countDiscarded   = 0;
		countSegregating = 0;
		Dpolarized       = 0;
		Dout1            = 0;
		Dout2            = 0;
		countCases       = new int[5];
	}
	
	
	/**
	 * check if valid available sites, get only fixnum number of sites and initialise column for SFS count
	 * @param lRow
	 * @param getOutgroups
	 * @return
	 */
	boolean isValid(List<Character> lRow, String getOutgroups){
		
		column = lRow;
		outgroups = getOutgroups;
				
		//Check if enough number of valid sites
		if (column.size() >= fixNum) {
			countSites++;
			//get only fixed number of sites
			column = column.subList(0, fixNum);
			return true;
		} else {
			countDiscarded++;
			return false;	
		}
	}
	
	
    /**
     * Main SFS method
     * @param lRow
     * @throws IOException
     */
	void doSFS() throws IOException {
		
		//Count nucleotides in fixedNum subset
		Map<String, Integer> countNuc = new HashMap<String, Integer>();
		countNuc.put("A", Collections.frequency(column, (char)'A'));
		countNuc.put("C", Collections.frequency(column, (char)'C'));
		countNuc.put("T", Collections.frequency(column, (char)'T'));
		countNuc.put("G", Collections.frequency(column, (char)'G'));
		
		//Order Hash by values using custom class MapUtil
		countNuc = MapUtil.sortByValue(countNuc);

		//SFS frequencies
		getFrequencies(countNuc);
		getUnfoldedFreqAndDiv(countNuc, outgroups.substring(2,3), outgroups.substring(4,5));  //yakuba & simulans
	}

 
    /**
     * Check if position is polymorphic. Checks if any key (A, T, G, C) in a Map has the total value of the fixNum.
     * @param map
     * @param fixNum
     * @return
     */
    boolean isPolymorphic (Map<String, Integer> map, int fixNum ){
    	if (map.containsValue(fixNum)) {
    		return false;
    	} else {
    		return true;
    	}
    }
    
    /**
     * Get second value from a Map with ordered (descending) counts of(A,T,G,C)
     * @param map
     * @return
     */
    Integer getSecondValue (Map<String, Integer> map){
		int i = 0;
		Integer result = null;
		for (Map.Entry<String, Integer> entry: map.entrySet()) {
			if (i++ == 1) {
				result = entry.getValue();
				break;
			}
		}
		return result;
    }
    
    /**
     * Get second value from a Map with ordered (descending) counts of(A,T,G,C)
     * @param map
     * @return
     */
    String getSecondKey (Map<String, Integer> map){
		int i = 0;
		String result = null;
		for (Map.Entry<String, Integer> entry: map.entrySet()) {
			if (i++ == 1) {
				result = entry.getKey();
				break;
			}
		}
		return result;
    }
  
	/**
	 * Get frequencies Folded
	 * @param map
	 * @return
	 */
    void getFrequencies (Map<String, Integer> map){

		if (isPolymorphic(map, fixNum)) {
			sfs[getSecondValue(map)]++;
			countSegregating++;
		} else {
			sfs[0]++; //No polymorphism
		}
    }
    
    /**
     * Get frequencies Unfolded : Assume that Outgroup2 is CLOSER to the sample
     * Cases:
     * Case1 -> Equal outgroups. One allele equal to outgroups. COUNT the allele that is different from outgroups
     * Case2 -> Different outgroups, 1 allele equal to outgroup2. COUNT the allele that is different from outgroup2
     * Case3 -> Equal outgroups. Alleles NOT equal to outgroups. COUNT minor allele.
     * Case4 -> Different outgroups. Alleles NOT equal to outgroups. COUNT minor allele.
     * Case5 -> Outgroup(s) with gaps/ambiguous bases.
     * 
     * @param map
     * @param outgroup1
     * @param outgroup2
     * @return
     */ 
    void getUnfoldedFreqAndDiv (Map<String, Integer> map, String outgroup1, String outgroup2 ){

    	String major = map.keySet().iterator().next();   //get first key
    	
    	if (isPolymorphic(map, fixNum)) {
    		
    		//Get major and minor alleles and frequencies
			String minor = getSecondKey(map);
			Integer major_count = map.values().iterator().next();  //get first value
			Integer minor_count = getSecondValue(map);
			
			//Check derivate alleles				
			if ( outgroup1.equals("-") || outgroup2.equals("-") || outgroup1.equals("N") || outgroup2.equals("N") ) {
				
				countCases[4]++;                                                             //Polymorphic but outgroup(s) with gaps/ambiguous
				
			} else {	
				
				if ( outgroup1.equals(outgroup2) ) {                                         //First: Equal outgroups	
					if (minor.equals(outgroup1)) ufold_sfs[major_count]++;                   //Major is derivate (Case1)
					if (major.equals(outgroup1)) ufold_sfs[minor_count]++;                   //Minor is derivate (Case1)
					if (minor.equals(outgroup1) || major.equals(outgroup1)) countCases[0]++;
			
					if (!major.equals(outgroup1) && !minor.equals(outgroup1))  {
						ufold_sfs[minor_count-1]++;                                          //Minor is derivate by default (Case3)
						countCases[2]++;
					}
			
				} else {                                                                     //Second: Divergent outgroups			
					if (minor.equals(outgroup2)) ufold_sfs[major_count]++;                   //Major is derivate as minor & outg2 are the same (Case2)
					if (major.equals(outgroup2)) ufold_sfs[minor_count]++;                   //Minor is derivate as minor & outg2 are the same (Case2)
					if (minor.equals(outgroup1) || major.equals(outgroup1)) countCases[1]++;
					
					if (!major.equals(outgroup1) && !minor.equals(outgroup1)) {
						ufold_sfs[minor_count-1]++;                                          //Minor is derivate by default (Case4)
						countCases[3]++;
					}
				}
			}

    	} else {                                      
    		//No polymorphism 
    		ufold_sfs[0]++; 
    		
    		//Calculate Divergences
		    	if(outgroup2.equals(outgroup1) && !outgroup2.equals("N") && !outgroup2.equals("-")){
		    		if(!major.equals(outgroup2)) Dpolarized++;   //sample is different from both outgroups and outgroups are the same
		    	} else if (major.equals(outgroup2) && !outgroup2.equals("N") && !outgroup2.equals("-")) {
		    		Dout1++;                                    //sample equal to outgroup2
		    	} else if (major.equals(outgroup1) && !outgroup1.equals("N") && !outgroup1.equals("-")){
		    		Dout2++;                                    //sample equal to outgroup1
		    	}

    	}   	
    }
    
	
	/**
	 * Return results
	 */
    
	int[] foldedSFS(){
		return sfs;
	}
	int[] unfoldedSFS(){
		return ufold_sfs;
	}
	int countSites(){
		return countSites;
	}
	int countDiscarded(){
		return countDiscarded;
	}
	int countSegregating(){
		return countSegregating;
	}
	String classSite(){
		return classSite;
	}
	int[] unfoldCases() {
		return countCases;
	}
	int Dpolarized() {
		return Dpolarized;
	}
	int Dout1() {
		return Dout1;
	}
	int Dout2() {
		return Dout2;
	}
		
    /**
     * Return Sum of the SFS categories
     * @param sfs
     * @return
     */
    int totalSites(int[] sfs) {
        int sum= 0; 
        for (int i=0; i < sfs.length; i++) 
            sum += sfs[i];
        return sum;
    }
    
    /**
     * Return sum of the SFS categories (segregating sites only)
     * @param sfs
     * @return
     */
    int totalSegregating(int[] sfs) {
        int sum= 0; 
        for (int i=1; i < sfs.length; i++)   //discard first element (total of non-segregating)
            sum += sfs[i];
        return sum;
    }
 	
}

