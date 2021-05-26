package com.mol.util;

import java.util.*;

import com.mol.construct.childNode;
import com.mol.mono.MonoUnit;

public class ConnectedGraphs {
	
	
	public static List<childNode> addingBrancingUnits(List<childNode> cNodeList,Map<Integer,MonoUnit> monounitList, int dp, int[] branchingbndsper)
	{
		int numberof405bnds= branchingbndsper[0];
		int numberof55bnds = branchingbndsper[1];
		
		boolean checkfornext = true;
		
		System.out.println("\n\n chain starts");
		while(checkfornext)
		{
			Map<String,Object> returnObj = connectedChains(cNodeList);
			List<List<String>> connectedgraphs = (List<List<String>>) returnObj.get("connectedchains");
			List<String> missingNodes = (List<String>) returnObj.get("missingnodes");
			int noofunitstoadd = (connectedgraphs !=null)?connectedgraphs.size():0;			
			System.out.println("\nnoofunitstoadd="+noofunitstoadd);
			System.out.println(cNodeList);
			System.out.println("connectedgraphs="+connectedgraphs);
			System.out.println("missingNodes="+missingNodes);		
			
			
			if (noofunitstoadd <=1 || (numberof405bnds<=0 && numberof55bnds<=0))
			{
				checkfornext=false;		
				 //System.exit(0);
			}
			if (checkfornext)
			{					
				boolean is4O5updated = false;
				boolean is55updated = false;
				String disconnectednode = connectedgraphs.get(connectedgraphs.size()-1).get(0);				
				int missingnodeid = Integer.parseInt(disconnectednode.substring(0,disconnectednode.length()-1));
				//System.out.println("disconnectednode"+disconnectednode+", missingnodeid"+missingnodeid);
						
			    if (missingnodeid !=0)
			    {
			    	MonoUnit newnode  = monounitList.get(missingnodeid);
					
					childNode missingnode = new childNode();							
					missingnode.setpType(newnode.getId()+newnode.getType());
					List<String> chainbefore= connectedgraphs.get(connectedgraphs.size()-2);
					conloop: for (int l1=0;l1<chainbefore.size(); l1++)
					{
						String nodeInMonoChain = chainbefore.get(l1);
						if (chainbefore.size()>l1+1)
							nodeInMonoChain = chainbefore.get(l1+1);		
						//System.out.println("nodeInMonoChain"+nodeInMonoChain);
						int elemInMonochain = Integer.parseInt(nodeInMonoChain.substring(0,nodeInMonoChain.length()-1));
						//System.out.println("nodeInMonoChain"+nodeInMonoChain+", elemInMonochain"+elemInMonochain);
						MonoUnit monoinchain  = monounitList.get(elemInMonochain);					
						
						if (newnode.isFivthC() && monoinchain != null && monoinchain.isFivthC() && numberof55bnds>0)
						{	
							newnode.setBnd1("5");	
							newnode.setFivthC(false);
							monoinchain.setpBnd("5");
							monoinchain.setFivthC(false);
							
							/* adding updated node to monolist - start*/
							monounitList.remove(missingnodeid);
							monounitList.put(missingnodeid, newnode);
							monounitList.remove(elemInMonochain);
							monounitList.put(elemInMonochain, monoinchain);
							/* adding updated node to monolist - end*/
							
							
							missingnode.setpType(monoinchain.getId()+monoinchain.getType());
							missingnode.setcType(newnode.getId()+newnode.getType());
							missingnode.setBondType(monoinchain.getpBnd()+newnode.getBnd1());
							cNodeList.add(missingnode);								
							is55updated=true;
						}	
						
						if (!is55updated & numberof405bnds>0)
						{
							if (newnode.isFivthC() && monoinchain != null && (monoinchain.isFourthC() && monoinchain.isFivthC()))
							{	
								newnode.setBnd1("5");	
								newnode.setFivthC(false);
								monoinchain.setpBnd("4O");
								monoinchain.setFourthC(false);
								
								/* adding updated node to monolist - start*/
								monounitList.remove(missingnodeid);
								monounitList.put(missingnodeid, newnode);
								monounitList.remove(elemInMonochain);
								monounitList.put(elemInMonochain, monoinchain);
								/* adding updated node to monolist - end*/
								
								missingnode.setpType(monoinchain.getId()+monoinchain.getType());
								missingnode.setcType(newnode.getId()+newnode.getType());
								missingnode.setBondType(monoinchain.getpBnd()+newnode.getBnd1());
								cNodeList.add(missingnode);
								//sop("After cNodeList 1="+cNodeList);
								is4O5updated=true;
							}
							else if (newnode.isFourthC() && monoinchain != null && monoinchain.isFivthC())
							{	
								newnode.setBnd1("O4");	
								newnode.setFourthC(false);
								monoinchain.setpBnd("5");
								monoinchain.setFivthC(false);
								
								/* adding updated node to monolist - start*/
								monounitList.remove(missingnodeid);
								monounitList.put(missingnodeid, newnode);
								monounitList.remove(elemInMonochain);
								monounitList.put(elemInMonochain, monoinchain);
								/* adding updated node to monolist - end*/
								
								
								missingnode.setpType(monoinchain.getId()+monoinchain.getType());
								missingnode.setcType(newnode.getId()+newnode.getType());
								missingnode.setBondType(monoinchain.getpBnd()+newnode.getBnd1());
								cNodeList.add(missingnode);
								//sop("After cNodeList 2="+cNodeList);
								is4O5updated=true;
							}
						}
												
						if (is4O5updated || is55updated) 
						{	
							if (is4O5updated && numberof405bnds > 0)	numberof405bnds-=1;		
							if (is55updated && numberof55bnds > 0)		numberof55bnds-=1;							
							break conloop;							
						}
					}							
			    }	
			    
			    //System.exit(0);
			   
			}	
			
		}
		
		System.out.println("\n chain ends \n\n");
		
		List<childNode> returncNodeList = new ArrayList<childNode>();
	    boolean isSorted = false;
	    int idx =0;
	    Map<Integer, Integer> nodespresent = new HashMap<Integer,Integer>();
	   
	    for (childNode cn : cNodeList)
		{
	      nodespresent.put(Integer.parseInt(cn.getcType().substring(0,cn.getcType().length()-1)), cNodeList.indexOf(cn));							
		}
	    //sop("nodespresent keys="+nodespresent.keySet());
	    List<Integer> sortedKeys = new ArrayList<Integer>(nodespresent.keySet()); 
	    Collections.sort(sortedKeys);
	    
	    //sop("sortedKeys="+sortedKeys);
	    
		for (int key : sortedKeys)
		{
			returncNodeList.add(cNodeList.get(nodespresent.get(key)));
		}
		
		//System.out.println("returncNodeList="+returncNodeList);
			
		
		return returncNodeList;			
		
	}
	
	private static Map<String,Object> connectedChains(List<childNode> cNodeList)
	{
	 	List<List<String>> connectedchains = new ArrayList<List<String>>();
	    List<String> connectednodes= new ArrayList<String>();
	    List<String> missingNodes = new ArrayList<String>();
	    Map<String,Object> returnMap = new HashMap<String, Object>();
	    boolean firstnode=true;		    
	    
	    List<childNode> sortedNodeList = new ArrayList<childNode>();		   
	    List<Float>  key = new ArrayList<Float>();
	    List<Integer>  value = new ArrayList<Integer>();
	    String[] chr = {"a","b","c"};
	    float ct= 0.1f;
	    for (childNode cn : cNodeList)
		{
	    	String nodeval = cn.getpType().substring(0,cn.getpType().length()-1);
	    	Float num = Float.parseFloat(nodeval);	    		    	
	    	if (key.contains(num))
	    	{
	    		ct += 0.1;
	    		num=num+ct;
	    	}
	    	key.add(num);
	    	value.add(cNodeList.indexOf(cn));
		}    
	    		    
	    System.out.println(key);
	    System.out.println(value);
	    List<Float> sortedKeys = new ArrayList<Float>(key);	    
	    Collections.sort(sortedKeys);
	    System.out.println(sortedKeys);
		for (Float sortval : sortedKeys)
		{
			sortedNodeList.add(cNodeList.get(value.get(key.indexOf(sortval))));
		}
		int len=0;
		System.out.println("conn="+sortedNodeList);
		for (childNode cn : sortedNodeList)
		{
			len++;
			String node1 = cn.getpType();
			String node2 = cn.getcType();	
			if (firstnode)
			{	
				connectednodes.add(node1);
				connectednodes.add(node2);					
				firstnode=false;
			}
			else if (connectednodes.contains(node1))
			{					
				connectednodes.add(node2);					
			}
			if (!connectednodes.contains(node1))
			{
				connectedchains.add(connectednodes);				
				if (len <= sortedNodeList.size()) missingNodes.add(node1);
				connectednodes = new ArrayList<String>();
				connectednodes.add(node1);
				connectednodes.add(node2);
			}
			if (len >= sortedNodeList.size())
			{
				connectedchains.add(connectednodes);				
			}
		}		
		if (connectedchains.size()>1)
		{   
			returnMap.put("connectedchains", connectedchains);
			returnMap.put("missingnodes", missingNodes);			
		}	
		
		return returnMap;
	}
}
