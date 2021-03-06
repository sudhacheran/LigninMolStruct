package com.mol.util;


import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.mol.construct.childNode;
import org.apache.log4j.Logger;

/**
 * Generation the Adjacency and Connectivity matrices 
 * - Writing the matrices to CSV
 * @author Sudha
 *
 */

public class AdjacencyAndConnectivityMatrix {
	  static Logger logger = Logger.getLogger(AdjacencyAndConnectivityMatrix.class);
	  
	  
	  private boolean adjMatrix[][];
	  private char connMatrix[][];
	  private int numVertices;
	  public String[] nodeLabels;

	  // Initialize the matrix
	  public AdjacencyAndConnectivityMatrix(int numVertices) {
		
	    this.numVertices = numVertices;
	    adjMatrix = new boolean[numVertices][numVertices];
	    connMatrix = new char[numVertices][numVertices];
	    nodeLabels = new String[numVertices];
	    
	  }

	  // Add edges
	  public void addEdge(int i, int j, char[] bndpos) {
	    adjMatrix[i][j] = true;
	    //adjMatrix[j][i] = true;
	    connMatrix[i][j] = bndpos[0];
	    connMatrix[j][i] = bndpos[1];	    
	    
	  }

	  // Remove edges
	  public void removeEdge(int i, int j) {
	    adjMatrix[i][j] = false;
	    adjMatrix[j][i] = false;
	    connMatrix[i][j] = Character.MIN_VALUE;
	    connMatrix[j][i] = Character.MIN_VALUE;	
	  }
	  
	  
	  /*
	   * Methods generates the Adjacency matrix and connectivity matrix based on the Lignin graph generated
	   */
	  public void getAdjandConnMatrix(List<childNode> cNodeList)
		{
		    
		  	//logger.debug("Generating Adjacency and connectivity matrices"); 
			for (childNode cn : cNodeList)
			{
				int node1 = Integer.parseInt(cn.getpType().substring(0,cn.getpType().length()-1));
				int node2 = Integer.parseInt(cn.getcType().substring(0,cn.getcType().length()-1));	
				
				nodeLabels[node1-1]=cn.getpType();
				nodeLabels[node2-1]=cn.getcType();				
				//System.out.println(nodeLabels[node1-1]+"--"+nodeLabels[node2-1]);
				
				String bndType = cn.getBondType();
				char b1 = bndType.charAt(0);
				char b2 = bndType.charAt(bndType.length()-1);
				
				String [] onward = new String[] {Constants.BO4, Constants.B5, Constants.BB, Constants._4O5, Constants._55};
				String [] bckward = new String[] {Constants._4OB, Constants._5B, Constants._5O4};
				boolean fwd = Arrays.stream(onward).anyMatch(bndType::equals);
				boolean bck = Arrays.stream(bckward).anyMatch(bndType::equals);
				//System.out.println(cn.getpType()+"-"+cn.getcType()+"bndType="+bndType+"fwd="+fwd+"bck="+bck);
				
				if (fwd)
					addEdge(node1-1, node2-1,new char[]{b1,b2});		  
				if (bck)
					addEdge(node2-1, node1-1,new char[]{b2,b1});			
				
			}
			
		}
	  
	  
	  // Print the matrix
	  public String toString() {
	    StringBuilder s = new StringBuilder();
	    s.append("Adjacency matrix\n");
	    s.append("   ");
	    for (int i = 0; i < numVertices; i++) {
	    	s.append((nodeLabels[i]!=""?nodeLabels[i]: i+1)+ " ");
	    	//System.out.println(nodeLabels[i]);
	    }
	    s.append("\n");
	    
	    for (int i = 0; i < numVertices; i++) {
	      s.append((nodeLabels[i]!=""?nodeLabels[i]: i+1) + "  ");
	      for (boolean j : adjMatrix[i]) {
	        s.append((j ? 1 : 0) + "  ");
	      }
	      s.append("\n");
	    }	    
	    	
	    s.append("\n");
	    s.append("Connectivity matrix\n");
	    s.append("   ");
	    for (int i = 0; i < numVertices; i++) {
	    	s.append((nodeLabels[i]!=""?nodeLabels[i]: i+1) + " ");
	    }
	    s.append("\n");
	    
	    for (int i = 0; i < numVertices; i++) {
	      s.append((nodeLabels[i]!=""?nodeLabels[i]: i+1) + "  ");
	      for (int j = 0; j < numVertices; j++) {
	        s.append(connMatrix[i][j]=='\0'? '0'+"  " : connMatrix[i][j] + "  ");
	      }
	      s.append("\n");
	    }	    
	    
	    return s.toString();
	  }
	  
	  
	  public String getCSVformat(String ligninid) {
		    StringBuilder s = new StringBuilder();
		    s.append("Lignin ID,"+ligninid+"\n");
		    System.out.println(Arrays.asList(nodeLabels));
		    s.append("Adjacency matrix\n");
		    s.append(",");
		    for (int i = 0; i < numVertices; i++) {
		    	s.append((nodeLabels[i]!=""?nodeLabels[i]: i+1)+ ",");
		    }
		    s.append("\n");
		    
		    for (int i = 0; i < numVertices; i++) {
		      s.append((nodeLabels[i]!=""?nodeLabels[i]: i+1) + ",");
		      for (boolean j : adjMatrix[i]) {
		        s.append((j ? 1 : 0) + ",");
		      }
		      s.append("\n");
		    }	    
		    	
		    s.append("\n");
		    s.append("Connectivity matrix\n");
		    s.append(",");
		    for (int i = 0; i < numVertices; i++) {
		    	s.append((nodeLabels[i]!=""?nodeLabels[i]: i+1) + ",");
		    }
		    s.append("\n");
		    
		    for (int i = 0; i < numVertices; i++) {
		      s.append((nodeLabels[i]!=""?nodeLabels[i]: i+1) + ",");
		      for (int j = 0; j < numVertices; j++) {
		        s.append(connMatrix[i][j]=='\0'? '0'+", " : connMatrix[i][j] + ", ");
		      }
		      s.append("\n");
		    }	    
		    
		    return s.toString();
		  }
	  
	  
	  public void writetoCSV(String csvName)
	  {			 
		       
	        try {
	            PrintWriter output = new PrintWriter("matrix/"+csvName+".csv", "UTF-8");
	            output.write(getCSVformat(csvName));		           
	            output.close();
	            logger.info("Successfully written the CSV file"+csvName+".csv");
	        } catch (Exception e) {
	           
	           logger.error(e.getMessage());

	        }		    
		}	
	  
	  
	  public static void main(String args[]) {
		  
		  AdjacencyAndConnectivityMatrix g = new AdjacencyAndConnectivityMatrix(6);
		  logger.debug("test");
		  g.addEdge(0, 1,new char[]{'8','4'});		  
		  g.addEdge(1, 2, new char[]{'5','4'});		  
		  g.addEdge(2, 3, new char[]{'5','4'});
		  g.addEdge(4, 5, new char[]{'5','5'});		    
		  logger.debug(g.toString());
		  g.writetoCSV("1");
		  
	  }
	}