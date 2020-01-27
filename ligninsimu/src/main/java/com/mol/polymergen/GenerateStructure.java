package com.mol.polymergen;



import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Scanner;

import org.openscience.cdk.Bond;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import com.mol.mono.*;
import com.mol.construct.*;


public class GenerateStructure {

	
	public static void main(String s[]) throws IOException
	{
		Scanner in = new Scanner(System.in); 
		System.out.println("Enter the DP (no. of monomer units) : ");  
        String dp = in.nextLine(); 
        System.out.println("You entered dp = "+dp);
        System.out.println("G/H/S ratio : "+ encodercopy.getMonoPer());
        runData(Integer.parseInt(dp));
	}
	
	public static void runData(int dp) throws IOException
	{
		Date dt = new Date();
		List<Lignin> ligninList = encodercopy.generateEncodedFile(dp);
		System.out.println(ligninList.size());
		int i=0;
		PrintWriter writer;		
		writer = new PrintWriter("LigninSmiles.txt", "UTF-8");
		
		for (Lignin lignin:ligninList)
		{
			int parentNode = Integer.parseInt(""+lignin.getPType().substring(0,lignin.getPType().length()-1));
			String parentNodeType = ""+lignin.getPType().charAt(lignin.getPType().length()-1);
			//System.out.println(parentNode +"-"+parentNodeType);
			MonolignolBase g = MonolignolBase.getMonolignolUnit(parentNodeType, parentNode);
						
				//G_Conferyl g = new G_Conferyl(1);
				List<Integer> cc = new ArrayList<Integer>(); 
				List<childNode> cNList = lignin.getcNode();	
				for (childNode cN: cNList)
				{
					 //Child object initialization
					int gp = Integer.parseInt(""+cN.getpType().substring(0, cN.getpType().length()-1));
					String gptype = ""+cN.getpType().charAt(cN.getpType().length()-1);
					int gc = Integer.parseInt(""+cN.getcType().substring(0, cN.getcType().length()-1));
					String gctype = ""+cN.getcType().charAt(cN.getcType().length()-1);
					
					if (cN.getpType().equals(lignin.getPType()))
					{
						cc.add(gc);
						//System.out.println("Entered 1");
						//System.out.println(gc +"-"+gctype);
						MonolignolBase g2 = MonolignolBase.getMonolignolUnit(gctype, gc);							
						if (cN.getBondType().equals("B04"))
							g = generateStruct (g,g2,"bC_1","O2_"+gc);
						else if (cN.getBondType().equals("BB"))
							g = generateStruct (g,g2,"bC_1","bC_"+gc);
						else if (cN.getBondType().equals("55"))
							g = generateStruct (g,g2,"c5_1","c5_"+gc);
						else if (cN.getBondType().equals("A04"))
							g = generateStruct (g,g2,"aC_1","O2_"+gc);
						else if (cN.getBondType().equals("B5"))
							g = generateStruct (g,g2,"bC_1","c5_"+gc);										
					}
					else 
					{	
						//System.out.println("gc="+gp+"gc="+gc+"cc="+cc);
						if (cc.contains(gp) && cc.contains(gc))
						{	
							//System.out.println("Entered 2");
							if (cN.getBondType().equals("B04"))
								g = generateStruct (g,"bC_"+gp,"O2_"+gc);
							else if (cN.getBondType().equals("BB"))
								g = generateStruct (g,"bC_"+gp,"bC_"+gc);
							else if (cN.getBondType().equals("55"))
								g = generateStruct (g,"c5_"+gp,"c5_"+gc);
							else if (cN.getBondType().equals("A04"))
								g = generateStruct (g,"aC_"+gp,"O2_"+gc);
							else if (cN.getBondType().equals("B5"))
								g = generateStruct (g,"bC_"+gp,"c5_"+gc);							
						}	
						else if (!cc.contains(gc))
						{
							cc.add(gc);
							//System.out.println(gc +"-"+gctype);
							MonolignolBase g2 = MonolignolBase.getMonolignolUnit(gctype, gc);
							//System.out.println("Entered 3");
							if (cN.getBondType().equals("B04"))
								g = generateStruct (g,g2,"bC_"+gp,"O2_"+gc);
							else if (cN.getBondType().equals("BB"))
								g = generateStruct (g,g2,"bC_"+gp,"bC_"+gc);
							else if (cN.getBondType().equals("55"))
								g = generateStruct (g,g2,"c5_"+gp,"c5_"+gc);
							else if (cN.getBondType().equals("A04"))
								g = generateStruct (g,g2,"aC_"+gp,"O2_"+gc);
							else if (cN.getBondType().equals("B5"))
								g = generateStruct (g,g2,"bC_"+gp,"c5_"+gc);
						}	
					}
				}
				i++;
				g.toIMage(i, g.getMol());
				writer.println("Lignin="+i+",SMILEString="+g.getSmile(g.getMol()));
				System.out.println("Lignin="+i+",MolWt="+g.getMolWt(g.getMol())+"Struct="+lignin+"SMILEString="+g.getSmile(g.getMol()));
				
		}
		Date dt2 = new Date();		
		System.out.println("Time taken="+ getDifference(dt,dt2));
		writer.close();
		
	}
	
	static MonolignolBase generateStruct(MonolignolBase g, MonolignolBase g2, String b1, String b2)
	{
		IAtomContainer g1_mol = g.getMol();		
		IAtomContainer g2_mol = g2.getMol();
		IAtom atm1 = null, atm2 = null, atm3=null;
		
		for (int i=0;i<g1_mol.getAtomCount();i++)
		{
			if (g1_mol.getAtom(i).getAtomTypeName()!=null  && g1_mol.getAtom(i).getAtomTypeName().equals(b1))
			{
				atm1 = 	g1_mol.getAtom(i);				
				if (i>1)
				{
					atm3 = g1_mol.getAtom(i-1);					
					IBond bnd2 = g1_mol.getBond(atm3, atm1);
					if (bnd2!=null && bnd2.getOrder().equals(IBond.Order.DOUBLE))
					{
						bnd2.setOrder(IBond.Order.SINGLE);
					}
				}
			}
		}
		
		for (int i=0;i<g2_mol.getAtomCount();i++)
		{
			if (g2_mol.getAtom(i).getAtomTypeName()!= null && g2_mol.getAtom(i).getAtomTypeName().equals(b2))
				atm2 = 	g2_mol.getAtom(i);
			g1_mol.addAtom(g2_mol.getAtom(i));			
		}
		
		for (int i=0;i<g2_mol.getBondCount();i++)
		{
			g1_mol.addBond(g2_mol.getBond(i));			
		} 
		IBond bnd = new Bond(atm1,atm2);
		g1_mol.addBond(bnd);
		return g;
	}
	
	static MonolignolBase generateStruct(MonolignolBase g, String b1, String b2)
	{
		IAtomContainer g1_mol = g.getMol();		
		IAtom atm1 = null, atm2 = null, atm3=null;
		
		for (int i=0;i<g1_mol.getAtomCount();i++)
		{
			if (g1_mol.getAtom(i).getAtomTypeName()!=null  && g1_mol.getAtom(i).getAtomTypeName().equals(b1))
			{
				atm1 = 	g1_mol.getAtom(i);	
				//System.out.println(b1+"-"+atm1.getAtomTypeName());
				if (i>1)
				{
					atm3 = g1_mol.getAtom(i-1);					
					IBond bnd2 = g1_mol.getBond(atm3, atm1);
					if (bnd2!=null && bnd2.getOrder().equals(IBond.Order.DOUBLE))
					{
						bnd2.setOrder(IBond.Order.SINGLE);
					}
				}
			}			
		}
		
		for (int i=0;i<g1_mol.getAtomCount();i++)
		{
			if (g1_mol.getAtom(i).getAtomTypeName()!= null && g1_mol.getAtom(i).getAtomTypeName().equals(b2))
			{
				atm2 = 	g1_mol.getAtom(i);
				//System.out.println(b2+"-"+atm2.getAtomTypeName());				
			}						
		}	
		IBond bnd = new Bond(atm1,atm2);
		g.getMol().addBond(bnd);		
		return g;
	}
	
	 private static String getDifference(Date d1, Date d2) 
	 {
			String result = null;
			/** in milliseconds */
			long diff = d2.getTime() - d1.getTime();
			/** remove the milliseconds part */ 
			diff = diff / 1000;
			//long days = diff / (24 * 60 * 60);
			long hours = diff / (60 * 60) % 24;
			long minutes = diff / 60 % 60;
			long seconds = diff % 60;
			result = hours+":"+minutes+":"+seconds;
			return result;
	 }
}
