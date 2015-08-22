import java.io.*;
import java.util.Arrays;
import java.util.List;


public class IntegrativeMKT {

	/**
	 * MAIN
	 * Get slides from a verticalFasta file (each line is a position)
	 * Work is done by the class SFS
	 * @param args
	 * @throws IOException
	 */
	public static void main(String args[]) throws IOException  {
				
		String inputFile  = args[0];
		String outputFile = args[1];
		int swSize      = Integer.parseInt(args[2]);
		int fixNum      = Integer.parseInt(args[3]);
		int sfsGroup    = Integer.parseInt(args[4]);
		int threshold   = Integer.parseInt(args[5]);
		String Chr = args[6];
		
		BufferedWriter outputBuffer = new BufferedWriter(new FileWriter(outputFile));
		
		//Initialize sliding window controller
	    SlidingWindows window = new SlidingWindows(inputFile, swSize);
	    
	    //Initialize vFasta position reader
	    GetRegionVFasta position = new GetRegionVFasta(inputFile, fixNum);
	    	     	
    	//initialize SFS's
    	SFS sfs_whole = new SFS(fixNum);
    	SFS sfs_4f    = new SFS(fixNum, "4");
    	SFS sfs_0f    = new SFS(fixNum, "0");
    	SFS sfs_U     = new SFS(fixNum, "U");
    	SFS sfs_i     = new SFS(fixNum, "i");
    	SFS sfs_I     = new SFS(fixNum, "I");
    	SFS sfs_N     = new SFS(fixNum, "N");
    	SFS sfs_2f    = new SFS(fixNum, "2");
    	
    	SFS[] SFSs = {sfs_whole, sfs_4f, sfs_0f, sfs_U, sfs_i, sfs_I, sfs_N, sfs_2f};
    	
    	//control variables for complex 2fold codon cases
    	int complexCount ;
    	
    	/*
    	boolean firstLastIsPolymorphic = false;
    	List<String> complex_outgroups = new ArrayList<String>();
    	List<List<Character>> complex_polymorphic = new ArrayList<List<Character>>();
    	*/
    	
    	
	    //Window loop
	    while(window.hasNext()){
	    	window.next();
	    	
	    	//Initialize count for complex 2fold positions
	    	complexCount = 0;

	    	
	    	//set regions for SFS
	    	for (SFS sfs:SFSs)
	    		sfs.setRegion();
	    		
	    	//positions Loop
	    	position.setRegion(window.start(), window.stop());
	    	while(position.hasNext()){	    		
	    		position.next();
	    		
	   	    	String outgroups = position.getOutgroups();              //get Outgroups and recoded
    			List<Character> polymorphic = position.getRandColumn();  //get aligned polymorphic position
				String recoded = outgroups.substring(0,1);               //separate recoded from outgroups
    			

				
				//******** NEEDS REFACTORING **********//
				//Correct for 2fold variations
				
				//Simple cases
				if(recoded.matches("(3|2)")) {                           
					
					//Ile - Met codon     
					if (recoded.matches("3") ) {
						if (polymorphic.contains((char)'G') ) {
							recoded = "0";
						} else {
							recoded = "4";
						}
														
					//other simple 2fold
				    } else {
				    	
						if (polymorphic.contains((char)'T') || polymorphic.contains((char)'C') )  {
							
							if (polymorphic.contains((char)'C') && polymorphic.contains((char)'T') ) {
								recoded = "4";
							}
							if ( polymorphic.contains((char)'A') || polymorphic.contains((char)'G') )  {
								recoded = "0";    //possible redundancy with next if
							}
						}

						if (polymorphic.contains((char)'A') || polymorphic.contains((char)'G') )  {
							
							if (polymorphic.contains((char)'G') || polymorphic.contains((char)'A') ) {
								recoded = "4";
							}
							if ( polymorphic.contains((char)'T') || polymorphic.contains((char)'C') )  {
								recoded = "0";
							}
						}
					
					}
														
				} 
				
				//Complex cases
				if (recoded.matches("(E|e|F|f|D|d|P|p)") || 
						(complexCount == 1 && recoded.matches("(E|e|F|f|D|d|P|p|0)")) )
				{   
					
					complexCount++;
					
					/* 
					complex_outgroups.add(outgroups);
					complex_polymorphic.add(polymorphic);
					
					if (complexCount == 1 && checkPolymorphic(polymorphic)) {
						firstLastIsPolymorphic = true;
					} else {
						
					}
					
					*/
				}
				
				//*************//
				
				
    			//do each SFS count for each class site
    			for (SFS sfs:SFSs) {

    				//Non 2fold cases
    				if ( recoded.equals(sfs.classSite()) || sfs.classSite().equals("total") ) {
    					if (sfs.isValid(polymorphic, outgroups)){ 
    	    				sfs.doSFS();
    	    			} 
    				}

    			}

	    			
	    	}//end position loop
	    	
	    	//Write output  		    	
	    	//Prepare Header and final output string
	    	String header =  "#Start"   +
				    			"\t#Stop"  +
				    			"\t#Chr"  +
				    			"\t#class" +
				    			"\t#validSites"  +
				    			"\t#Discarded"   +
				    			"\t#Segregating" +
				    			"\t#Dpolarized"  +
				    			"\t#Dout1"  +
				    			"\t#Dout2"  +
				    			"\t#Di"  +
				    			"\t#D4"  +
				    			"\t#Ki"  +
				    			"\t#K4"  +
				    			"\t#d"      +
				    			"\t#b"      +
				    			"\t#f"      +
				    			"\t#gamma"      +
				    			"\t#alpha"  +
				    			"\t#DoS"    +
				    			"\t#omega"  +
				    			"\t#Pi_weakdel"  +
				    			"\t#b_nat"      +
				    			"\t#f_nat"      +
				    			"\t#gamma_nat"      +
				    			"\t#alpha_nat"  +
				    			"\t#DoS_nat"    +
				    			"\t#omega_nat"  +
				    			"\t#Pi_weakdel_nat"  +
				    			"\t#countComplex" +
				    			"\t#SFS"    +
				    			"\t#SFSunf" +
				    			"\t#SFSunfClasses" +
				    			"\n" ;
	    	
	    	String base = (window.start()+1) + "\t" + (window.stop()+1) + "\t" + Chr + "\t";
	    	String output = "";     	
	    	
    	    /////////////////////////
	    	//** GET MKT results **//
	    	/////////////////////////
	    	
	    	// 4fold MKT values
	    	int[] sfs4fold = sfs_4f.foldedSFS();
	    	int Seg4f      = sfs_4f.countSegregating();
	    	int D4f        = sfs_4f.Dout1() + sfs_4f.Dpolarized();   //total (unpolarized) divergence with outgroup1
	    	int m4f        = sfs_4f.countSites();
    	
	    	for (SFS sfs:SFSs) {

	    		Float Di;
	    		Float D4;
	    		Float Ki;
	    		Float K4;
	    		
	    		
	    		Float d     ;
	    		Float b     ;
	    		Float f     ;
	    		Float gamma ;
	    		Float alpha ;
	    		Float DoS   ;
	    		Float omega ;
	    		Float pi_weakdel;	 
	    		
	    		//Nature method
	    		Float b_nat     ;
	    		Float f_nat     ;
	    		Float gamma_nat ;
	    		Float alpha_nat ;
	    		Float DoS_nat   ;
	    		Float omega_nat ;
	    		Float pi_weakdel_nat;	
	    		
				if ( sfs.classSite().equals("4") || sfs.classSite().equals("total") ) {
	    			
		    		Di         = null;
		    		D4         = null;
		    		Ki         = null;
					K4         = null;
					
					d          = null;
		    		b          = null;
		    		f          = null;
		    		gamma      = null;
		    		alpha      = null;
		    		DoS        = null;
		    		omega      = null;
		    		pi_weakdel = null;
		    		
		    		b_nat          = null;
		    		f_nat           = null;
		    		gamma_nat       = null;
		    		alpha_nat       = null;
		    		DoS_nat         = null;
		    		omega_nat       = null;
		    		pi_weakdel_nat  = null;
		    		
	    		} else { 
	    			
			    	//selected class MKT values
			    	int[] sfsClass = sfs.foldedSFS();
			    	int SegClass   = sfs.countSegregating();
			    	int DClass     = sfs.Dout1() + sfs.Dpolarized();   //total (unpolarized) divergence with outgroup1
			    	int mClass     = sfs.countSites();
			    	
			    	//do MKT
			    	CalculateMKT intMKT = new CalculateMKT(sfsClass, sfs4fold, SegClass, Seg4f, 
												    	   DClass, D4f, mClass, m4f, sfsGroup, threshold);
			    	
			    	Di          = intMKT.Di();
			    	D4          = intMKT.D4();
			    	Ki          = intMKT.Ki();
			    	K4          = intMKT.K4();
			    	
    				d          = intMKT.d();
    				b          = intMKT.b();
    				f          = intMKT.f();
    				gamma      = intMKT.gamma();
    				alpha      = intMKT.alpha();
    				DoS        = intMKT.DoS();
    				omega      = intMKT.omega();
    				pi_weakdel = intMKT.pi_weakdel();
    				
    				//Nature method
    				b_nat          = intMKT.b_nat();
    				f_nat          = intMKT.f_nat();
    				gamma_nat      = intMKT.gamma_nat();
    				alpha_nat      = intMKT.alpha_nat();
    				DoS_nat        = intMKT.DoS_nat();
    				omega_nat      = intMKT.omega_nat();
    				pi_weakdel_nat = intMKT.pi_weakdel_nat();
	    		}

			    output += base + 
	    				sfs.classSite()         + "\t" +
	    				sfs.countSites()        + "\t" +
	    				sfs.countDiscarded()    + "\t" +
	    				sfs.countSegregating()  + "\t" +
	    				sfs.Dpolarized()        + "\t" +
	    				sfs.Dout1()             + "\t" +
	    				sfs.Dout2()             + "\t" +
	    				Di              + "\t" +
	    				D4              + "\t" +
	    				Ki              + "\t" +
	    				K4              + "\t" +
	    				d              + "\t" +
	    				b              + "\t" +
	    				f              + "\t" +
	    				gamma          + "\t" +
	    				alpha          + "\t" +
	    				DoS            + "\t" +
	    				omega          + "\t" +
	    				pi_weakdel     + "\t" +
	    				b_nat              + "\t" +
	    				f_nat              + "\t" +
	    				gamma_nat          + "\t" +
	    				alpha_nat          + "\t" +
	    				DoS_nat            + "\t" +
	    				omega_nat          + "\t" +
	    				pi_weakdel_nat     + "\t" ;
	    				
	    		if (sfs.classSite().equals("2")) {
	    			output += complexCount + "\t";
	    		} else {
	    			output += "NA" + "\t";
	    		}
	    				
	    		output +=	Arrays.toString( sfs.foldedSFS() )  + "\t" + 
	    				Arrays.toString( sfs.unfoldedSFS() )    + "\t" +
	    				Arrays.toString( sfs.unfoldCases() )    + "\n" ;  	
	    	}
	    	
	    	outputBuffer.write(header+output+"\n");
	    	//System.out.println(header+output+"\n");

	    	
	    }//end window loop
	    
	    outputBuffer.close();
	    
	}//end main method	
	
	
	
	/* fast check if polymorphic*/
	static boolean checkPolymorphic (List<Character> column){

		boolean result = false;
		
		//get some nucleotide. 
		//(returned columns are clean, don't have N's or other ambiguous bases)		
		Character nuc = column.get(0);
		
		if ( column.contains((char)'A') && !nuc.equals((char)'A') )
			result = true;
		
		if ( column.contains((char)'T') && !nuc.equals((char)'T') )
			result = true;
		
		if ( column.contains((char)'C') && !nuc.equals((char)'C') )
			result = true;
		
		if ( column.contains((char)'C') && !nuc.equals((char)'C') )
			result = true;
				
		return result;		
	}
	
	
}