package com.mol.construct;


import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class encodercopy {
	
	 //Inputs 
	 
	static String[] bonds = { "B04", "BB", "B5","405","55" }; // Bonds possible
	static int[] bondper = { 50, 30, 20, 0, 0};  // Ratio of the bonds

	static String[] mono = { "G", "S", "H" }; // Monomer units
	static int[] monoPer = {100,0,0}; // Ratio of the units
	
	static int dp =100; // No of monomer units
	
	public static String[] getBonds() {
		return bonds;
	}

	public static void setBonds(String[] bonds) {
		encodercopy.bonds = bonds;
	}

	public static int[] getBondper() {
		return bondper;
	}

	public static void setBondper(int[] bondper) {
		encodercopy.bondper = bondper;
	}

	public static String[] getMono() {
		return mono;
	}

	public static void setMono(String[] mono) {
		encodercopy.mono = mono;
	}

	public static String getMonoPer() {
		String monoStr = monoPer[0]+" / "+monoPer[1]+" / "+monoPer[2];		
		return monoStr;
	}

	public static void setMonoPer(int[] monoPer) {
		encodercopy.monoPer = monoPer;
	}

	public static int getDp() {
		return dp;
	}

	public static void setDp(int dp) {
		encodercopy.dp = dp;
	}

	public static void main(String s[])
	{
		 try {
			 //sop("here");
			for (Lignin ln :(generateEncodedFile(dp)))
			{
				sop(ln);
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Function to encode the monolist, bondlist to datafile
	 * @param dp
	 * @throws UnsupportedEncodingException 
	 * @throws FileNotFoundException 
	 */
	public static List<Lignin>  generateEncodedFile(int dp) throws FileNotFoundException, UnsupportedEncodingException
	{
		
		encodercopy en = new encodercopy();
		String[] monoArr = getList(mono, monoPer, dp);  // generating the list of monomers for the dp size
		
		GeneratePolymer gp = new GeneratePolymer();
		List<String[]> monoListStr = new ArrayList<String[]>();
		int n = monoArr.length;
		long noOfPer = GeneratePolymer.factorial(dp);
		GeneratePolymer.permute(monoArr, 0, n - 1, monoListStr, noOfPer);  // generating the permutation of monomers with the maximum permutation as 5000
		
		monoListStr=removeDuplicates(monoListStr);  // removing the repeated permutations in the list

		//sop(monoListStr.size());
		
		/*for(int i=0;i<monoListStr.size();i++)
		{
			sop(Arrays.toString(monoListStr.get(i)));
		}*/
		
		ArrayList<Integer> getSeq = Combination.printCombination(dp);  // possible combinations 
		//sop(getSeq);
		int nP = (getSeq.size()<=2)?getSeq.size():getSeq.size()/2;	
		String[] bondArr = getList(bonds, bondper, dp);
		//for (String S: bondArr) sop(S);
		List<String[]> bondListStr = new ArrayList<String[]>();

		n = bondArr.length;
		GeneratePolymer.permute(bondArr, 0, n - 1, bondListStr, noOfPer);
		bondListStr=removeDuplicates(possibleBonds(bondListStr));
		//sop(bondListStr);
		/* If dp is less than 3, tried to create bondpermunation list as bonds
		/*if (getSeq.size()<=2)
		{
			for(int i=0;i<bonds.length;i++)
			{
				String b[] = new String[1];
				b[0]=bonds[i];
				bondListStr.add(b);
			}
		}*/

		//sop(bondListStr.size());
		
		/* for(int i=0;i<bondListStr.size();i++)
		{
			sop(Arrays.toString(bondListStr.get(i)));
		} */
		
		PrintWriter writer;		
		writer = new PrintWriter("DataSet3.txt", "UTF-8");		
		List<String> data = new ArrayList<String>();
		List<Lignin> ligninList = new ArrayList<Lignin>();
		        
		for(int lp=0;lp<monoListStr.size();lp++)
		{
			monoArr = monoListStr.get(lp);
			//for (String s: monoArr) sop(s);
			for(int bp=0;bp<bondListStr.size();bp++)
			{
			bondArr = bondListStr.get(bp);
			//for (String s: bondArr) sop(s);
			List<childNode> cNodeList = new ArrayList<childNode>();
			ArrayList<monolignol> monoList = new ArrayList<monolignol>();
			ArrayList<String> bondList = new ArrayList<String>();
			String dataStr="";
			Lignin lignin=new Lignin();
			Edges edg = new Edges();
			
			for(int i=0;i<dp;i++)
			{
				monolignol m1 = new monolignol(monoArr[i], i+1);			
				monoList.add(m1);     
			}
			//sop (monoList);
			for(int i=0;i<bondArr.length;i++)
			{
				bondList.add(bondArr[i]);
			}
			edg.setEdWgt(bondList);
			// Mapping two nodes for the graph  (taking the combination values)
			for (int i=0;i<getSeq.size();i+=2)
			{
				//sop(getSeq.get(i)+"-"+getSeq.get(i+1));
				edg = (Edges) en.encode(edg, monoList.get(getSeq.get(i)-1),monoList.get(getSeq.get(i+1)-1));
				monoList.set(getSeq.get(i)-1,(monolignol) edg.getMonoL());						
			}			
			//sop (monoList);
			ArrayList<Integer> edgWtIn = new ArrayList<Integer>();
			for(monolignol mono: monoList)
			{			
				String[] ed= mono.getEdgeType();
				boolean setVal = false;
				for(int c=0;c<ed.length;c++)
				{
					setVal = false;
					for(int k=0;k<bonds.length;k++)
					{
						if (!setVal && ed[c]!=null && ed[c].equals(bonds[k]))
						{
							edgWtIn.add(k+1);
							setVal = true;
						}
					}
				}
			}
			edg.setEdgWgtIn(edgWtIn);
			//sop(edg);	
			//sop(monoList);
			int oldNode=0;
					
			for(int i=0;i<edg.getpNode().size();i++)
			{
				int nodeId = edg.getpNode().get(i);				
				monolignol monomer = monoList.get(nodeId-1);
				monolignol child1 = monomer.getchild1();
				monolignol child2 = monomer.getchild2();
				monolignol child3 = monomer.getchild3();				
				String nodeType = monomer.objId+monomer.getType();
				String connection1 = (child1 == null)?"":child1.objId+child1.getType();
				String connection2 = (child2 == null)?"":child2.objId+child2.getType();
				String connection3 = (child3 == null)?"":child3.objId+child3.getType();
				String[] connections = {connection1,connection2,connection3};
				String[] links = monomer.getEdgeType();
				
				String bondtype="";			
				int lbl = monomer.lbl.getIndex();
					
			    for (int b=0;b<monomer.getEdgeType().length;b++)
				{
					if (monomer.getEdgeType()!=null && monomer.getEdgeType()[b] != null)
					{
						bondtype = bondtype.concat(monomer.getEdgeType()[b]);
					}
				}
			    
				if (nodeId!= oldNode)
				{
					//writer.println(nodeId+"#\t\t"+nodeType+"\t\t"+connection1+"\t"+connection2+"\t"+connection3+"\t"+bondtype+"\t\t"+lbl);
					dataStr = dataStr.concat(nodeType+"\t("+bondtype+")"+"\t"+connection1+"\t"+connection2+"\t"+connection3+"\t|\t");								
					oldNode = nodeId;		
					
					for (int cnt=0;cnt<connections.length;cnt++)
					{
						if (connections[cnt]!=null && connections[cnt]!="" && links[cnt]!=null && links[cnt] !="" )
						{
							//sop("Link1="+nodeType+"--"+connections[cnt]+"--"+links[cnt]);
							String pn = nodeType.substring(1, nodeType.length());
							String cn = connections[cnt].substring(1, connections[cnt].length());
							String plnk = links[cnt].substring(0,1);
							String clnk = links[cnt].substring(1, links[cnt].length());
							//sop(pn+"-"+cn+"-"+plnk+"-"+clnk);
							//if (!((pn.equals("S") && plnk.equals("5")) || (cn.equals("S") && clnk.contentEquals("5"))))
							{								
								childNode cN = new childNode();						
								cN.setpType(nodeType);
								cN.setcType(connections[cnt]);
								cN.setBondType(links[cnt]);
								cNodeList.add(cN);
								//sop(cN);
							}
						}
					}	
					
				}
			}
			oldNode=0;
			for(int i=0;i<edg.getcNode().size();i++)		
			{
				int nodeId = edg.getcNode().get(i);
				monolignol monomer = monoList.get(nodeId-1);
				monolignol child1 = monomer.getchild1();
				monolignol child2 = monomer.getchild2();
				monolignol child3 = monomer.getchild3();
				String nodeType =monomer.objId+monomer.getType();
				String connection1 = (child1 == null)?"":child1.objId+child1.getType();
				String connection2 = (child2 == null)?"":child2.objId+child2.getType();
				String connection3 = (child3 == null)?"":child3.objId+child3.getType();
				String bondtype ="";
								
				String[] connections = {connection1,connection2,connection3};
				String[] links = monomer.getEdgeType();
				
				for (int b=0;b<monomer.getEdgeType().length;b++)
					if (monomer.getEdgeType()!=null && monomer.getEdgeType()[b] != null)
					{
						List<String> bndL = Arrays.asList(bonds);
						bondtype = bondtype.concat(" "+monomer.getEdgeType()[b]).concat(" ");										
					}				
				
				int lbl = monomer.lbl.getIndex();
				ArrayList<Integer> pNodeList = edg.getpNode();
				
				if (!pNodeList.contains(nodeId))
				{
					if (nodeId!= oldNode)
					{
						//writer.println(nodeId+"#\t\t"+nodeType+"\t\t"+connection1+"\t"+connection2+"\t"+connection3+"\t"+bondtype+"\t\t"+lbl);
						dataStr = dataStr.concat(nodeType+"\t("+bondtype+")"+"\t\t"+connection1+"\t"+connection2+"\t"+connection3);
						oldNode = nodeId;						
						for (int cnt=0;cnt<connections.length;cnt++)
						{		
							if (connections[cnt]!=null && connections[cnt]!="" && links[cnt]!=null && links[cnt] !="" )
							{
								//sop("Link2="+nodeType+"--"+connections[cnt]+"--"+links[cnt]);
								String pn = nodeType.substring(1, nodeType.length());
								String cn = connections[cnt].substring(1, connections[cnt].length());
								String plnk = links[cnt].substring(0,1);
								String clnk = links[cnt].substring(1, links[cnt].length());
								//sop(pn+"-"+cn+"-"+plnk+"-"+clnk);
								//if (!((pn.equals("S") && plnk.equals("5")) || (cn.equals("S") && clnk.contentEquals("5"))))
								{								
									childNode cN = new childNode();
									cN.setpType(nodeType);
									cN.setcType(connections[cnt]);
									cN.setBondType(links[cnt]);
									cNodeList.add(cN);
									//sop(cN);
								}
							}
						}	
					}
				}		
			}
			if (!data.contains(dataStr)) 
			{
				writer.println(dataStr);
				data.add(dataStr);
				lignin.setPType(cNodeList.get(0).pType);
				lignin.setcNode(cNodeList);
				ligninList.add(lignin);				
				lignin=new Lignin();
			}
			
		}
		}
		//sop(ligninList);
		writer.close();		
		return RemoveRepetitionfromLigninList(ligninList);
	}
	
	/**
	 * Function to map the child node to parent node.
	 * Mapping is set in the Monolignol object 
	 * @param ed
	 * @param l1
	 * @param l2
	 * @return
	 */
	public Object encode(Edges ed,monolignol l1, monolignol l2)
	{		
		if (checkID(l1.getchild1(),l2.objId) && checkID(l1.getchild2(),l2.objId) && checkID(l1.getchild3(),l2.objId))
		{	
			Random rn = new Random();
			//int random = rn.nextInt(2);
			//sop(random+1);
			sop("Edges"+ed.getEdWgt());
			int random=0;
			boolean isUpdated = false;
			switch (random+1)
			{
			case 3:				
					if (l1.getchild1() == null)
					{
						l1.setchild1(l2);
						isUpdated = true;
					}
					else
						if (l1.getchild2() == null)
						{
							l1.setchild2(l2);
							isUpdated = true;
						}
					else if (l1.getchild3() == null) 
					{
						l1.setchild3(l2);
						isUpdated = true;
					}
					break;
			case 2: 
				if (l1.getchild1() == null)
				{
					l1.setchild1(l2);
					isUpdated = true;
				}
				else if (l1.getchild2() == null)
				{
					l1.setchild2(l2);
					isUpdated = true;
				}
				break;
			case 1:
				if (l1.getchild1() == null)
				{
					l1.setchild1(l2);
					isUpdated = true;
				}
				break;
			}			
			if (!isFilled (l1.objId, ed) && isUpdated)
			{
				l1 = getEdgeType(ed,l1);		
				if (l1.isEdgUpdated)
				{
					ed.getpNode().add(l1.objId);
					ed.getcNode().add(l2.objId);				
					ed.setEdWgt(l1.getEdWgt());	
					ed.setMonoL(l1);		
				}
			}		
		}
		sop(ed.getEdWgt());
		return ed;			
	}
	
	static boolean checkID(monolignol obj, int id2)
	{	
		if (obj==null || obj.objId != id2)
		{
			return true;
		}
		return false;
	}
	
	static boolean isFilled(int id,Edges ed)
	{
		boolean isfill=false;
		int cnt=0;
		for(int pVer:ed.getpNode())
		{
			if (id == pVer)
			{
				cnt++;
			}
		}		
		if (cnt >= 3) isfill = true;
		return isfill;
	}
	
	static void sop(Object s)
	{
		System.out.println(s);
	}
	
	
	
	static String[] getList(String[] str, int[] bondper, int dp) {
		
		String[] strBonds = new String[dp];
		int[] totcnt = new int[str.length];
		// initializing the totcount[] for each bond with 1  
		for(int inc=0; inc < totcnt.length; inc++)
		{
			totcnt[inc]+=1;
		}
		// Loop to create the possible bond names (or) monomer names with the dp count
		int cnt = 0;
		boolean brkLoop = false;
		while (!brkLoop) {
			for (int j = 0; j < str.length; j++) {
				double cal = calPer(bondper[j], dp);	
				//sop(str[j]+"--"+cal);
				if (cal > 0 && totcnt[j] <= cal && cnt < dp) {						
					strBonds[cnt++] = str[j];	
					totcnt[j] = totcnt[j]+1;	
				}
			}
			brkLoop = true;
			for (int j = 0; j < str.length; j++) {
				double cal = calPer(bondper[j], dp);		
				if (cal > 0 && totcnt[j] <= cal && cnt < dp) {	
					brkLoop = false;
				}				
			}
		}
		//for (String s: strBonds) sop("getlist ->"+s);
		//Loop to add bonds at the end for chain length		
		while (cnt < dp)
		{
			for (int k=0;k<str.length;k++)
	      	 {
	            if (cnt >= dp)
	            {
	              break;
	            }
	            if (bondper[k]>0)
	            	strBonds[cnt++]=str[k];
	         }
		}				
		return strBonds;
	}
	
	static double calPer(int val, int dp) {
		if (val == 0)
			return val;
		double flrval = dp * (val / 100f);		
		return ((flrval > 0.5) ? Math.round(flrval) : 0);
	}
	
	static monolignol getEdgeType(Edges ed, monolignol m1)
	{	
		m1.isEdgUpdated = false;		
		ArrayList<String> edgArr = ed.getEdWgt();
		String c[]=new String[3];
	    List<String> newEdg = Arrays.asList(m1.getEdgeType());
	    int cnt=0;
	    for(String s: newEdg)
	    {
	    	if (s!=null)
	    	c[cnt++] = ""+s.charAt(0);
	    }	    
	    for(int k=0;k<m1.getEdgeType().length;k++)
		{	 
			  for (int j=0;j<edgArr.size();j++) 
			  {
				  String oldegd = edgArr.get(j);
				  String plnk = ""+oldegd.charAt(0);
				  String clnk = ""+oldegd.charAt(oldegd.length()-1);
				  String ctype = m1.getchild1().type;				 
				  if (k==1)
					  ctype= (m1.getchild2() != null?m1.getchild2().type:"");
				  else if (k==2)
					  ctype= (m1.getchild3() != null?m1.getchild3().type:"");
				  //sop(Arrays.toString(c)+"-"+plnk+"--"+clnk+"-"+m1.type+"-"+ctype);
				  //sop(Arrays.asList(m1.getEdgeType())+"-"+edgArr.get(j));
				  if (!((m1.edgeType.equals("S") && plnk.equals("5"))|| (ctype.equals("S") && clnk.equals("5"))))
				  {
					  if (!Arrays.asList(m1.getEdgeType()).contains(edgArr.get(j)) && !Arrays.asList(c).contains(plnk))
					  {			
						  //sop(m1.getEdgeType()[k]);
						  if (m1.getEdgeType()[k]==null || m1.getEdgeType()[k].isEmpty())
						  {
							  
							  m1.getEdgeType()[k]= edgArr.get(j);
							  edgArr.remove(edgArr.get(j));	
							  m1.setEdWgt(edgArr);		
							  m1.isEdgUpdated = true;
							  //sop(m1);
							  return m1;
						  }
					  }	
				  }
			  }			 
		  }	    
		return m1;				
	}
	
	static List<Lignin> RemoveRepetitionfromLigninList(List<Lignin> ligninlist)
	{
		ArrayList<Lignin> modifiedList = new ArrayList<Lignin>();
		//sop("Modifiedlist1="+ligninlist.size());
		for (Lignin lignin:ligninlist)
		{
			Lignin newLignin = new Lignin();
			newLignin.setPType(lignin.getPType());
			List<childNode> newChildList = new ArrayList<childNode>();
			List<String> cNode = new ArrayList<String>();
			for (childNode cN: lignin.getcNode())
			{	
				String cNodd_bond = cN.getcType()+"-"+cN.bondType.substring(1,cN.bondType.length());				
				//sop("cNodd_bond"+cNodd_bond +"ptype="+cN.getpType());
				boolean possibleBnd = true;
				String ctype = ""+cN.getcType().charAt(cN.getcType().length()-1);
				String ptype = ""+cN.getpType().charAt(cN.getpType().length()-1);
				String clink = ""+cN.bondType.charAt(cN.bondType.length()-1);
				String plink = ""+cN.bondType.charAt(0);
				
				if ((ctype.equals("S") && clink.equals("5"))||(ptype.equals("S") && plink.equals("5")))
				{
					possibleBnd = false;
				}
				if (!cNode.contains(cNodd_bond) && possibleBnd)
				{
					cNode.add(cNodd_bond);
					newChildList.add(cN);
				}
			}
			lignin.setcNode(newChildList);
			modifiedList.add(lignin);
			//sop(lignin);
		}		
		//sop("Modifiedlist="+modifiedList.size());
		return modifiedList;
	}
	
	static List<String[]> removeDuplicates(List<String[]> strData)
	{
		List<String[]> strList = new ArrayList<String[]>();
		for (int i=0;i<strData.size();i++){
			//sop("Entered loop");
			if (strList.isEmpty() || !isDuplicate(strList, strData.get(i)))
			{
				strList.add(strData.get(i));
			}
		}		
		return strList;
	}
	
	static boolean isDuplicate(List<String[]> strList, String[] strArr)
	{
		boolean isDup = false;		
		for (int i=0;i< strList.size();i++){
			String[] lclStrArr = strList.get(i);
			
			if (Arrays.toString(lclStrArr).equals(Arrays.toString(strArr)))				
					isDup = true;										
		}
		//sop (isDup);
		return isDup;
	}
	
	public static List<String[]> possibleBonds(List<String[]> bondListStr) {
		List<String[]> retArrr = new ArrayList<String[]>();
		for (String[] bnds : bondListStr)
		{
			char prevFirst = 0;
			char prevLast = 0;
			boolean ignore = false;
			//sop(Arrays.toString(bnds));
			List<String> bondArr = new ArrayList<String>();
			for (String bd: bnds)
			{
				char first = bd.charAt(0);
				char last = bd.charAt(bd.length()-1);
				char middle =  (bd.length()>2 ? bd.charAt(1): 0);
				if (prevLast == first)
				{
					if (first != last)
					{
						char temp = first;
						first = last;
						last = temp;	
						String bond = String.valueOf(first)+(middle!=0?String.valueOf(middle):"")+String.valueOf(last);
						bondArr.add(bond);
					}
					else
					{
						bondArr.add(" ");
						bondArr.add(bd);
					}
				}				
				else
				{
					bondArr.add(bd);
				}
				prevFirst = first;
				prevLast = last;
			}
			//sop(bondArr.toString());
			if (bondArr.size() > 0)
			{
				retArrr.add(bondArr.toArray(new String[0]));
			}
			//sop("Next>>");
			
		}
		return retArrr;
	}
}
