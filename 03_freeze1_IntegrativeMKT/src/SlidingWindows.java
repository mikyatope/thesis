import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Class to calculate sliding-windows coordinates
 * @author mike
 */

class SlidingWindows {

	private String filename;             
	private int swSize;                //Sliding Window size

	private long totalLength;          //total row number in file (total positions for a vFasta)

	private long start;
	private long stop;
	private boolean windowStatus ;
	
	
	/**
	 * Constructor
	 * @param getFilename
	 * @param getFixNum
	 * @param getSwSize
	 * @throws IOException 
	 */
	SlidingWindows(String getFilename, int getSwSize) throws IOException {
		
		//initialize variables
		filename = getFilename;
		swSize = getSwSize;

		windowStatus = true;
	    
		//get total row number (total positions for a vFasta)
		//returns 1 less row if last does not have \n!!!
		totalLength = countLines(filename);                             
	}

	/**
	 * Return true if next Sliding Window Available
	 * @return
	 */
	boolean hasNext(){
		return windowStatus;
	}
	
	/**
	 * Sliding Windows coordinates loop
	 */
	void next(){
		if (stop < totalLength) {
			if (stop == 0) {                                //Check if first loop, else do normal               
				start = 0;                                  //first start
				stop = stop + swSize -1;                    //first stop
			} else {
				start = start + swSize;	                    //normal start  
				if ( (stop + swSize) -1 > totalLength ) {   //if last loop is bigger than swsize
					stop = totalLength;
					windowStatus = false;
				} else {
					stop = stop + swSize;                   //normal stop
				}
			}
		} else {
			windowStatus = false;
		}
	}
	
	/**
	 * Return current Sliding Window Start
	 * @return
	 */
	long start(){
		return start;
	}
	
	/**
	 * Return current Sliding Window Stop
	 * @return
	 */
	long stop(){
		return stop;
	}
	

	/**
	 * Count lines in a file
	 * from: http://stackoverflow.com/questions/453018/number-of-lines-in-a-file-in-java
	 * ALERT! does not count last line if it doesn't end with \n 
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	static int countLines(String filename) throws IOException {
		InputStream is = new BufferedInputStream(new FileInputStream(filename));
		try {
		    byte[] c = new byte[1024];
		    int count = 0;
		    int readChars = 0;
		    while ((readChars = is.read(c)) != -1) {
		        for (int i = 0; i < readChars; ++i) {
		            if (c[i] == '\n')
		                ++count;
		        }
		    }
		    return count;
		} finally {
		    is.close();
		}
	}

}