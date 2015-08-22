import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/*
 * Calculate selection regimes using SFS data for different class sites
 */
class CalculateMKT {
	
	//Intermediate variables
	private float[] sfs_4f ;         //4fold SFS as totally neutral
	private float[] sfs ;            //other class site SFS  
	private float   Pi_weakdel ;     //fraction of weakly deleterious sites in site class
	private float   Pi_neutral ;     //fraction of neutral sites in site class
	private float   Pi ;             //segregating sites in class and region
	private float   P4 ;             //segregating 4fold class sites in region
	private float   mi ;             //total class sites in region
	private float   m4 ;             //total 4folf class sites in region
	private float   Di ;             //Absolute divergence (D) for class sites (Polarized)
	private float   D4 ;             //Absolute divergence (D) 4f class sites (Polarized)
	private float   Ki ;             //Divergence K (D/m) for class sites (Polarized)
	private float   K4 ;             //Divergence K (D/m) 4f class sites (Polarized)
	
	//Results
	private float d ;      //Highly deletereous fraction
	private float b ;      //Weak deletereous fraction
	private float f ;      //Neutral fraction
	private float gamma ;  //New Neutral fraction
	private float alpha ;  //Alpha
	private float DoS ;    //DoS
	private float omega ;  //Omega
	
	//nature method
	private float[] sfs_4f_nat ;         //4fold SFS as totally neutral
	private float[] sfs_nat ;            //other class site SFS  
	private float   Pi_weakdel_nat ;     //fraction of weakly deleterious sites in site class
	private float   Pi_neutral_nat ;     //fraction of neutral sites in site class
	
	private float b_nat ;      //Weak deletereous fraction
	private float f_nat ;      //Neutral fraction
	private float gamma_nat ;  //New Neutral fraction
	private float alpha_nat ;  //Alpha
	private float DoS_nat ;    //DoS
	private float omega_nat ;  //Omega
	
		
	CalculateMKT(int[] getSFS, int[] getSFS_4f, int getPi, int getP4, int getDi, 
			     int getD4, int get_mi, int get_m4, int get_group, int threshold){
		
		//Initialize pre-SFSs Lists	
		List<Float> pre_sfs = new ArrayList<Float>();
		List<Float> pre_sfs_4f = new ArrayList<Float>();
		List<Float> pre_sfs_nat = new ArrayList<Float>();
		List<Float> pre_sfs_4f_nat = new ArrayList<Float>();


		//Sum SFSs in groups of elements. Exceptions:
    	//Leave singletons alone (also fixed, first element in array)
    	//Sum last group regardless of group size
        for (int i=0, j=0, sum=0, sum_4f=0; i < getSFS.length; i++) {
      	
        	if (i == 0 || i == 1) {                      // [0] fixed, [1] singletons
         		pre_sfs.add((float) getSFS[i]);   
        		pre_sfs_4f.add((float) getSFS_4f[i]); 
        	} else {
        		sum += getSFS[i];
        		sum_4f += getSFS_4f[i];
        		
        		j++;
        		
        		if (j == get_group || (j > 0 && i == (getSFS.length-1)) ) {
        			pre_sfs.add((float) sum);
        			pre_sfs_4f.add((float) sum_4f);
        			sum = 0;
        			sum_4f = 0;
        			j=0;
        		}
        	}	
        }
        
        //Nature method
		//Sum SFSs in 2 groups using 'threshold'
        //Discard fixed (i== 0)
        for (int i=0, sum=0, sum_4f=0; i < getSFS.length; i++) {
      	
        	if (i == 0) {                                // discard fixed
         		pre_sfs_nat.add((float) getSFS[0]);   
        		pre_sfs_4f_nat.add((float) getSFS_4f[0]); 
        		
        	} else if (i <= threshold) {                 // first bin
        		sum += getSFS[i];
        		sum_4f += getSFS_4f[i];
        		
        		if (i==threshold) {
        			pre_sfs_nat.add((float) sum);
        			pre_sfs_4f_nat.add((float) sum_4f);
        			sum = 0;
        			sum_4f = 0;
        		}
        		
        	} else {                                     // second bin    
        		sum += getSFS[i];
        		sum_4f += getSFS_4f[i];
        		
        		if (i == (getSFS.length-1)) {
        			pre_sfs_nat.add((float) sum);
        			pre_sfs_4f_nat.add((float) sum_4f);
        		}
    		}
        		
        }
        
        //Convert pre-SFS lists to float[]
        sfs = toFloatArray(pre_sfs); 
        sfs_4f = toFloatArray(pre_sfs_4f); 
        sfs_nat = toFloatArray(pre_sfs_nat); 
        sfs_4f_nat = toFloatArray(pre_sfs_4f_nat);   

		
		//Initialize Pi/P4 mi/m4 and divergences
		Pi = (float) getPi;
		P4 = (float) getP4;
		mi = (float) get_mi;
		m4 = (float) get_m4;
		Di = (float) getDi;
		D4 = (float) getD4;
		Ki = Di/mi;
		K4 = D4/m4;
		
		
		//STEP1
		//divide SFSs by total segregating sites
        for (int i=0; i < sfs.length; i++) {
        	if (Pi > 0) {
        		sfs[i] = sfs[i]/Pi;	
        		sfs_4f[i] = sfs_4f[i]/P4;
        	} else {
        		sfs[i] = 0;	
        		sfs_4f[i] = 0;
        	}
        }
        ////Nature method (the same with different array size)
        for (int i=0; i < sfs_nat.length; i++) {
        	if (Pi > 0) {
        		sfs_nat[i] = sfs_nat[i]/Pi;	
        		sfs_4f_nat[i] = sfs_4f_nat[i]/P4;
        	} else {
        		sfs_nat[i] = 0;	
        		sfs_4f_nat[i] = 0;
        	}
        }
        
        
        //STEP2
        //SFS less SFS_4f  || at this step get rid off first [0] column (num. fixed alleles)
        float[] workSFS = new float[(sfs.length-1)];
        for (int i=1; i < sfs.length; i++) {
        	float j = Math.abs(sfs[i] - sfs_4f[i]) ;
        	workSFS[i-1] = j;
        }
        /////Nature method
        float[] workSFS_nat = new float[(sfs_nat.length-1)];
        for (int i=1; i < sfs_nat.length; i++) {
        	float j = Math.abs(sfs_nat[i] - sfs_4f_nat[i]) ;
        	workSFS_nat[i-1] = j;
        }
        
        
        //STEP3
		//Sum everything
        for (int i=0; i < workSFS.length; i++) 
        	Pi_weakdel += workSFS[i];
       
        Pi_neutral = 1-Pi_weakdel;
        
        /////Nature method
        for (int i=0; i < workSFS_nat.length; i++) 
        	Pi_weakdel_nat += workSFS_nat[i];
       
        Pi_neutral_nat = 1-Pi_weakdel_nat;
        
        
                
        //RESULTS
        //Highly deletereous
        d = 1-( (Pi/mi) / (P4/m4) );
        
        //Weak deletereous
        b     = (m4/mi)*(Pi_weakdel*Pi/P4);
        b_nat = (m4/mi)*(Pi_weakdel_nat*Pi/P4);
        
        //Neutral fraction
        f     = (m4/mi)*(Pi_neutral*Pi/P4);
        f_nat = (m4/mi)*(Pi_neutral_nat*Pi/P4);
        
        //New Neutral fraction
        gamma     = (m4/mi)*((Pi_neutral/Pi) - (Di/D4));
        gamma_nat = (m4/mi)*((Pi_neutral_nat/Pi) - (Di/D4));
        
        
        //Alpha
        alpha     = 1-( (D4*Pi_neutral*Pi) / (Di*P4) );
        alpha_nat = 1-( (D4*Pi_neutral_nat*Pi) / (Di*P4) );
        
        //DoS
        DoS     = ( Ki/(Ki+K4) ) - ( ( Pi_neutral*(Pi/mi) ) / ( ( Pi_neutral*(Pi/mi) )+(P4/m4) ) );
        DoS_nat = ( Ki/(Ki+K4) ) - ( ( Pi_neutral_nat*(Pi/mi) ) / ( ( Pi_neutral_nat*(Pi/mi) )+(P4/m4) ) );
        
        //Omega
        omega     = (Ki/K4)*alpha;
        omega_nat = (Ki/K4)*alpha_nat;
        
	}
	
	//intermediate values
	float Di(){
		return Di;
	}
	float D4(){
		return D4;
	}
	float Ki(){
		return Ki;
	}
	float K4(){
		return K4;
	}
	
	//results
	float d(){
		return d;
	}
	float b(){
		return b;
	}	
	float f(){
		return f;
	}
	float gamma(){
		return gamma;
	}
	float alpha(){
		return alpha;
	}
	float DoS(){
		return DoS;
	}
	float omega(){
		return omega;
	}
	float pi_weakdel(){
		return Pi_weakdel;
	}
	
	//Nature returns
	float b_nat(){
		return b_nat;
	}	
	float f_nat(){
		return f_nat;
	}
	float gamma_nat(){
		return gamma_nat;
	}
	float alpha_nat(){
		return alpha_nat;
	}
	float DoS_nat(){
		return DoS_nat;
	}
	float omega_nat(){
		return omega_nat;
	}
	float pi_weakdel_nat(){
		return Pi_weakdel_nat;
	}
	
	
	/**
	* Converts an array of Integer objects to an array of integer primitives
	* 
	* @param integerList the integer list
	* 
	* @return an array of integer primitives
	*/
	static float[] toFloatArray(List<Float> floatList) {
		float[] floatArray = new float[floatList.size()];
		
		for (int i = 0; i < floatList.size(); i++) 
			floatArray[i] = floatList.get(i);
		
		return floatArray;
	}
	
	
}
