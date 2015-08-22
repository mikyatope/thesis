import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Read vFasta line by line
 * @author mike
 *
 */

class GetRegionVFasta {

	private long   stop;
	private long   pos;               //position inside current region
	private long   row;               //Current reading row. Starts at 'getStart'
	private String currentRow;
	private String outgroups;
	
	private String      filename;
	private int         fixNum;       //Number of analysed nucleotides per position
	private long        rowLength;    //first row length
	private FileChannel inChan;
	private ByteBuffer  buf;
	
	
	private boolean positionStatus ;	
	
	
	/**
	 * Constructor
	 * @param getFilename
	 * @param getStart
	 * @param getStop
	 * @param getRowLength
	 * @param getFixNum
	 * @throws IOException
	 */
	GetRegionVFasta(String getFilename, int getFixNum) throws IOException{
		
		filename = getFilename;
		fixNum   = getFixNum;
		
		//get first row length
		//this program assumes fixed row length for a vFasta file!
		rowLength = byteRowLength(filename); 
		
		//check fixNum is less than row length
		if (rowLength < fixNum) fixNum = (int)rowLength;
		
		//Initialise randomAcces for the vFasta file
		RandomAccessFile inFile = new RandomAccessFile(filename, "r");
		inChan                  = inFile.getChannel();
	    buf                     = ByteBuffer.allocate(1);	
	}
	
	
	/**
	 * Initialise new Region
	 */
	void setRegion(long getStart, long getStop){
		positionStatus = true;
		row            = getStart;
		stop           = getStop;
	}
	

	/**
	 * Return true if next position available in region
	 * @return boolean
	 */
	boolean hasNext(){
		return positionStatus;
	}
	
	/**
	 * loop through multiple alignment positions as String in region
	 * @return
	 * @throws IOException 
	 */
	void next() throws IOException {
		//System.out.println(row + " " + stop);
		if (row <= stop) {
			readRow();
			if ( row == stop ) {              //Check if last loop
				positionStatus = false;
			} else {
				row++;
			}
		} else {
			positionStatus = false;
		}
		
	}
	
	/**
	 * Return current position inside region
	 * @return
	 */
	long currentPosition() {
		return pos;
	}
	
	/**
	 * Get all characters in a row (RandomAccess)
	 * Separate outgroups(+recoded) from polymorphic seqs.
	 * @throws IOException
	 */
	void readRow() throws IOException{
		currentRow = "";
		
		for (pos = 0; pos<=rowLength-1; pos++){        //Character by character RandomAccess buffer
			long absPos = ((rowLength+1)*row)+pos;     //Calculate position
			inChan.position(absPos);
			inChan.read(buf); 		                   //Read (buffer1)				
	    	buf.flip();                                //Get char into variable
	    	char c = (char) buf.get();    	
	    	currentRow += c;	                       //grow sequence
	    	buf.clear();                               //reset buffer and position
	    	inChan.position(0);
		}
		
		//Separate outgroups from polymorphic seqs.
		outgroups  = currentRow.substring(0, 5);
		currentRow = currentRow.substring(5);	
	}
	
	
	/**
	 * return current row / outgroups / site and discarded counts
	 * @return
	 */
	String getRawColumn(){
		return currentRow;
	}
	
	
	/**
	 * Shuffle column
	 * @return
	 */
	List<Character> getRandColumn(){
		
		currentRow = currentRow.replaceAll("N|M|R|W|S|Y|K|-", "");   //remove gaps and ambiguous bases
		
		List<Character> lRow = new ArrayList<Character>();           //create character list of remaining nucleotides 
		for (char c : currentRow.toCharArray()){                    
			lRow.add(c);
		}
		Collections.shuffle(lRow);                                   //build-in shuffle method

		return lRow;
	}
	
	String getOutgroups(){
		return outgroups;
	}
	


	/**
	 * Count number of bytes in first line of a file
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	int byteRowLength(String filename) throws IOException {
		LineNumberReader in = new LineNumberReader(new FileReader(filename));
		try {
			String s;
			s = in.readLine();
			int len = s.getBytes().length;
			return len;
		} finally {
			in.close();
		}
	}
	
	
}
