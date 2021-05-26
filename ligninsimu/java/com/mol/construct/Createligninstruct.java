package com.mol.construct;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.log4j.Logger;
import com.mol.mono.MonoUnit;
import com.mol.util.Constants;
import com.mol.util.AdjacencyAndConnectivityMatrix;
import com.mol.util.ConnectedGraphs;



public class Createligninstruct {
	
	 	//Inputs 
	  	static Logger logger = Logger.getLogger(Createligninstruct.class);

		static String[] bonds = { Constants.BO4, Constants.BB, Constants.B5}; // Bonds possible
		public static int[] bondper = {45,12,3};  // Ratio of the bonds
		
		static String[] branchingbnds = {Constants._4O5, Constants._55};
		static int[] branchingbndsper = {8, 25};

		static String[] mono = { "G", "S", "H" }; // Monomer units
		public static int[] monoPer = {100,0,0}; // Ratio of the units
		
		static int dp_units =100; // No of monomer units
		
		
		public static void main(String s[]) {		
			Date dt = new Date();
			getLiginList(dp_units, true);
			Date dt2 = new Date();
			System.out.println("Time taken="+ getDifference(dt,dt2));
			
		}	
		
				
		public static String getMtype()
		{
			if ( monoPer[1] > 0 && monoPer[0] > 0 && monoPer[2] > 0)
				return "sgh";
			else if (monoPer[1] == 0 &&  monoPer[0] > 0 && monoPer[2] > 0)
				return "gh";
			else if (monoPer[1] == 0 &&  monoPer[0] > 0 && monoPer[2] == 0)
				return "g";
			else if (monoPer[1] > 0 &&  monoPer[0] > 0 && monoPer[2] == 0)
				return "sg";
			return "";
				
		}
		
		public static String getMonoPer() {
			String monoStr = monoPer[0]+" / "+monoPer[1]+" / "+monoPer[2];		
			return monoStr;
		}
		public static List<Lignin> getLiginList(int dp, boolean branching)
		{
			String[] monoArr = encodercopy.getList(mono, monoPer, dp);
			long noOfPer = GeneratePolymer.factorial(dp);
			
			//System.out.println(Arrays.toString(monoArr));			
			//System.out.println("Number of permutation="+noOfPer);
						
			List<String[]> monoListStr = new ArrayList<String[]>();
			GeneratePolymer.permute(monoArr, 0, monoArr.length - 1, monoListStr, noOfPer);
			
			System.out.println("Before monomer list size="+monoListStr.size());
			
			/*for (String[] da : monoListStr)
			{
				System.out.println("Monomer permute="+Arrays.toString(da));
			}*/
						
			
			String[] bondArr = encodercopy.getList(bonds, bondper, dp);
			List<String[]> bondListStr = new ArrayList<String[]>();
			List<String[]> mArr;
			
			if (dp ==2)
			{
				bondListStr = generatepossiblebonds(bonds,bondper,dp, bondArr);			
				mArr = encodercopy.removeDuplicates(monoListStr);
			}
			else
			{
				GeneratePolymer.permute(bondArr, 0, bondArr.length-1, bondListStr, noOfPer);
				bondListStr =  encodercopy.removeDuplicates(encodercopy.possibleBonds(bondListStr));				
				mArr = encodercopy.removeDuplicates(monoListStr);
			}
			
			System.out.println("After monomer list size="+mArr.size());
			System.out.println("After bondListStr list size="+bondListStr.size());
			
			for (String[] da : mArr)
			{
				System.out.println(Arrays.toString(da));
			}
			//sop("bondListStr.size "+bondListStr.size());
			
			for (String[] da : bondListStr)
			{
				System.out.println(Arrays.toString(da));
			}
			
			
		    //sop("mArr.size "+mArr.size());			
		    
		    ArrayList<ArrayList<Integer>> monoseq = CombinationAlg.getCombination(dp,0);
		    //sop("combinations="+monoseq);
		    
		   
		    List<Lignin> ligninList = new ArrayList<Lignin>();		    
	        
			for(String[] monoList: mArr)
			{
				//System.out.println("MonoList ="+Arrays.toString(monoList));
				for (String[] bondList: bondListStr)
				{
					//System.out.println("BondList="+Arrays.toString(bondList));
					//System.out.println("BondList Length="+bondList.length);
					int mlen = 0;
				    Map<Integer,MonoUnit> monounitList = new HashMap<Integer,MonoUnit>();				    
				    for (int i=0;i<monoList.length; i++)
				    {
				    	monounitList.put(i+1, new MonoUnit(monoList[i],i+1));
				    }
					for (List<?> data: monoseq)
					{						
						String pBnd = ""+bondList[mlen].charAt(0);						
						String cBnd = ""+bondList[mlen].substring(1,bondList[mlen].length());
						
						if (pBnd.contentEquals("4"))
						{
							pBnd = ""+bondList[mlen].charAt(0)+bondList[mlen].charAt(1);
							cBnd = ""+bondList[mlen].charAt(bondList[mlen].length()-1);
						}
							
						//System.out.println("pBnd="+pBnd+",cBnd="+cBnd);
						
						int pID = (int)data.get(0);
						int cID = (int)data.get(1);
						MonoUnit pMono = monounitList.get(pID);
						MonoUnit cMono = monounitList.get(cID);						
						String pNode = pMono.getType();
						String cNode = cMono.getType();			
												
						pMono = checkBonds(pMono,pBnd, "p", pNode, pID);
						cMono = checkBonds(cMono,cBnd, "c",  pNode, pID);  //setpBnd(cBnd); // bond from Child C  - with P node
						if (pMono!= null && cMono != null)
						{
							if (pMono.getChild1()== null) pMono.setChild1(cMono);
							else if (pMono.getChild2()== null) pMono.setChild2(cMono);
							else if (pMono.getChild3()== null) pMono.setChild3(cMono);
							//System.out.println(pMono);						
							//System.out.println(cMono);
							monounitList.put((int)data.get(0), pMono);
							monounitList.put((int)data.get(1), cMono);
						}
						mlen++;
					}
					
					//System.out.println(monounitList);
					List<childNode> cNodeList = new ArrayList<childNode>();
					Lignin lignin=new Lignin();
					boolean ignore = false;
					for (int i:monounitList.keySet())
					{
						MonoUnit mn  = monounitList.get(i);
						//System.out.println(mn);									
						//System.out.println("INside"+mn);	
						if (mn.getChild1() != null)
						{
							MonoUnit cmn  = monounitList.get(mn.getChild1().getId());
							childNode cN = new childNode();						
							cN.setpType(cmn.getpID()+cmn.getpType());
							cN.setcType(cmn.getId()+cmn.getType());								
							cN.setBondType(mn.getBnd1()+cmn.getpBnd());	
							//System.out.println(cmn.getpType() +"-" + mn.getBnd1() +"-"+ cmn.getType() +"-"+cmn.getpBnd());							
							if ((cmn.getpType().equals("S") && (mn.getBnd1()== null || mn.getBnd1().equals("5"))) || (cmn.getType().equals("S") && (cmn.getpBnd()==null || cmn.getpBnd().equals("5"))))
							{
								ignore = true;								
							}
							if (!ignore) cNodeList.add(cN);
							//System.out.println("Ignore "+ignore);
							//System.out.println("Child1 "+cN);
						}
						
					}
					if (!ignore)
					{
						if (cNodeList!= null && cNodeList.size()>0)
						{
							AdjacencyAndConnectivityMatrix adjAndConnMx = new AdjacencyAndConnectivityMatrix(dp);
							adjAndConnMx.getAdjandConnMatrix(cNodeList);
							if (branching)   // This now only has two oligomers linking
							{
								cNodeList = ConnectedGraphs.addingBrancingUnits(cNodeList, monounitList, dp, branchingbndsper);
							}							
							adjAndConnMx.getAdjandConnMatrix(cNodeList);
							lignin.setPType(cNodeList.get(0).pType);
							lignin.setcNode(cNodeList);							
							lignin.setBndCnts(getBondRatiosforLignin(cNodeList));
							lignin.setAdjAndConnMtx(adjAndConnMx);
							ligninList.add(lignin);		
							
							//System.out.println(ligninList);
							//System.exit(0);
						}
					}
				}				
			}	
			System.out.println("Lignin Size="+ligninList.size());
			//for (Lignin lg : ligninList)
				//System.out.println(lg);
			return ligninList;
		}	
		
		
		
		static Map<String,Integer> getBondRatiosforLignin(List<childNode> cNodeList)
		{
			Map<String,Integer> bondpercent = new HashMap<String,Integer>();
			for (childNode cn : cNodeList)
			{	
				
				String bndType=cn.getBondType();
				if (cn.getBondType().equals(Constants._4OB))
					bndType = Constants.BO4;
				if (cn.getBondType().equals(Constants._5B))
					bndType = Constants.B5;
				if (cn.getBondType().equals(Constants._5O4))
					bndType = Constants._4O5;
				
				if (bondpercent.containsKey(bndType))
					bondpercent.replace(bndType,bondpercent.get(bndType)+1);
				else
					bondpercent.put(bndType, 1);
			}			
			return bondpercent;			
		}


		private static List<String[]> generatepossiblebonds(String[] bonds2, int[] bondper2, int dp2, String[] bondArr) 
		{
			List<String[]> lstBonds = new ArrayList<String[]>();
			if (dp2 == 2)
			{
				for (int i=0;i<bondper2.length;i++)
				{
					String [] possiblebonds = new String[1];
					if (bondper2[i]!=0)
						possiblebonds[0]=bonds2[i];
					if (possiblebonds.length>0 && possiblebonds[0] != null)
						lstBonds.add(possiblebonds);					
				}
			}
			else
			{
				GeneratePolymer.permute(bondArr, 0, bondArr.length-1, lstBonds, dp2);
			}
			//System.out.println("Possible bonds list size:"+lstBonds.size());
			return lstBonds;
		}

		static MonoUnit checkBonds(MonoUnit mono, String bnd, String type, String pNode, int pID)
		{
			boolean isupdated = false;
			if (bnd.equals("B") || bnd.equals("OB") && mono.isBetaC())
			{
				mono.setBetaC(false);
				mono = setBonds(mono,bnd,type, pNode, pID);
				if (mono != null) isupdated = true;
				
			}
			if (bnd.equals("5") && mono.isFivthC())
			{
				mono.setFivthC(false);
				mono = setBonds(mono,bnd, type, pNode, pID);
				if (mono != null) isupdated = true;
			}
			if (bnd.equals("O4") || bnd.equals("4") || bnd.equals("4O") && mono.isFourthC())
			{
				mono.setFourthC(false);
				mono = setBonds(mono,bnd,type, pNode, pID);
				if (mono != null) isupdated = true;
			}
			if (isupdated) return mono;
			return null;
		}
		static MonoUnit setBonds(MonoUnit mono, String bnd, String type, String pNode, int pID)
		{
			boolean isUpdated = false;
			if (type.equals("p"))
			{
				if (mono.getBnd1() == null)
				{
					mono.setBnd1(bnd);
					isUpdated = true;
				}
				else if (mono.getBnd2() == null)
				{
					mono.setBnd2(bnd);
					isUpdated = true;
				}
				else if (mono.getBnd3() == null)
				{
					mono.setBnd3(bnd);
					isUpdated = true;
				}
			}
			else if (type.equals("c"))
			{
				if (mono.getpBnd() == null)
				{
					mono.setpBnd(bnd);
					mono.setpType(pNode);
					mono.setpID(pID);
					isUpdated = true;
				}
			}
			if (isUpdated = true) return mono;
				
			return null;
		}
		
		static Map<String, Integer> setBondcount(List<String> bondList)
		{
			Map<String, Integer> boundcount = new HashMap<String, Integer>();
			boundcount.put(Constants.BO4, 0);
			boundcount.put(Constants.B5, 0);
			boundcount.put(Constants.BB, 0);
			boundcount.put(Constants._55, 0);
			boundcount.put(Constants.AO4, 0);
			for (String str : bondList) { 
		         if (boundcount.containsKey(str)) { 
		        	 boundcount.put(str, boundcount.get(str) + 1); 
		         } else { 
		        	 boundcount.put(str, 1); 
		         } 
		      } 
			return boundcount;
		}		
				
		
		 private static String getDifference(Date d1, Date d2) 
		 {
				String result = null;
				/** in milliseconds */
				long diff = d2.getTime() - d1.getTime();
				/**3 remove the milliseconds part */ 
				diff = diff / 1000;
				//long days = diff / (24 * 60 * 60);
				long hours = diff / (60 * 60) % 24;
				long minutes = diff / 60 % 60;
				long seconds = diff % 60;
				result = hours+":"+minutes+":"+seconds+"(Hrs:Mins:Secs)";
				return result;
		 }
		 
		
		
		static void sop(String s)
		{
			System.out.println(s);
		}
}
