package com.mol.polymergen;



import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import org.json.simple.JSONObject;
import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

import com.mol.mono.*;
import com.mol.pojo.LigninsJson;
import com.mol.pojo.SingleStructDet;
import com.mol.util.Constants;
import com.mol.util.lcLogger;
import com.mol.construct.*;


public class GenerateStructure {

	
	static boolean addtojson = true;
	
	public static void main(String s[]) throws IOException, NumberFormatException, CDKException
	{
		
		Scanner in = new Scanner(System.in); 
		System.out.println("Enter the DP (no. of monomer units) : ");  
        String dp = in.nextLine(); 
        System.out.println("You entered dp = "+dp);
        System.out.println("G/H/S ratio : "+ Createligninstruct.getMonoPer());
        runData(Integer.parseInt(dp));    
	}
	
	public static void runData(int dp) throws IOException, CDKException
	{
		Date dt = new Date();
		//List<Lignin> ligninList = encodercopy.generateEncodedFile(dp);
		List<Lignin> ligninList = Createligninstruct.getLiginList(dp, true);	
				
		lcLogger.sop(""+ligninList.size(), "GenerateStructure");
		List<String> ligninSmiles = new ArrayList<String>();
		List<Double> molweight = new ArrayList<Double>();
		List<Map> bondkeys = new ArrayList<Map>();
						
		int i=0;
		int totlig=0;
		PrintWriter writer;		
		writer = new PrintWriter("LigninSmiles.txt", "UTF-8");
		double s_g_ratio = Double.valueOf(Createligninstruct.monoPer[1]) / Double.valueOf(Createligninstruct.monoPer[0]);
		String s_g = String.format("%.2f", s_g_ratio);
		
		
		LigninsJson jsonObjLignin = new LigninsJson();
		jsonObjLignin.setDp(dp);
		jsonObjLignin.setMonoType(Createligninstruct.getMtype());
		jsonObjLignin.setS_G_Ratio(new Double(s_g));
		jsonObjLignin.setName("polymer");
		
					   
		int temp =0;		
		
		for (Lignin lignin:ligninList)
		{			
		   
		    //temp++;
		    System.out.println(lignin);
		    System.out.println(lignin.getAdjAndConnMtx());
			int parentNode = Integer.parseInt(""+lignin.getPType().substring(0,lignin.getPType().length()-1));
			String parentNodeType = ""+lignin.getPType().charAt(lignin.getPType().length()-1);						
				
			List<Integer> cc = new ArrayList<Integer>(); 
			List<childNode> cNList = lignin.getcNode();	
			MonolignolBase g = null;
			System.out.println("dp="+dp+","+"childnodesize="+cNList.size());
				
			
			if (cNList.size() >= dp-4)
			{
			g = MonolignolBase.getMonolignolUnit(parentNodeType, parentNode);
			for (childNode cN: cNList)
			{
				 //Child object initialization
				int gp = Integer.parseInt(""+cN.getpType().substring(0, cN.getpType().length()-1));
				String gptype = ""+cN.getpType().charAt(cN.getpType().length()-1);
				int gc = Integer.parseInt(""+cN.getcType().substring(0, cN.getcType().length()-1));
				String gctype = ""+cN.getcType().charAt(cN.getcType().length()-1);			
				//System.out.println("gc="+gc+",gp="+gp);
				
				if (cN.getpType().equals(lignin.getPType()))
				{					
					cc.add(gc);
					//System.out.println("Entered 1");
					//System.out.println(gp +"-"+gc);
					MonolignolBase g2 = MonolignolBase.getMonolignolUnit(gctype, gc);		
					g = genStructTwo(g,g2,cN.getBondType(),gp,gc);													
				}
				else 
				{	
					//System.out.println("gc="+gp+"gc="+gc+"cc="+cc);
					if (cc.contains(gp) && cc.contains(gc))
					{	
						//System.out.println("Entered 2");							
						if (cN.getBondType().equals(Constants.BO4))
							g = generateStruct (g,"bC_"+gp,"O2_"+gc, cN.getBondType());
						else if (cN.getBondType().equals(Constants._4OB))
							g = generateStruct (g,"O2_"+gp,"bC_"+gc, cN.getBondType());
						else if (cN.getBondType().equals(Constants.BB))
							g = generateStruct (g,"bC_"+gp,"bC_"+gc, cN.getBondType());
						else if (cN.getBondType().equals(Constants._55))
							g = generateStruct (g,"c5_"+gp,"c5_"+gc, cN.getBondType());
						else if (cN.getBondType().equals(Constants.AO4))
							g = generateStruct (g,"aC_"+gp,"O2_"+gc, cN.getBondType());
						else if (cN.getBondType().equals(Constants.B5))
							g = generateStruct (g,"bC_"+gp,"c5_"+gc,cN.getBondType());
						else if (cN.getBondType().equals(Constants._5B))
							g = generateStruct (g,"c5_"+gp,"bC_"+gc, cN.getBondType());
						else if (cN.getBondType().equals(Constants._4O5))
							g = generateStruct (g,"O2_"+gp,"c5_"+gc,cN.getBondType());
						else if (cN.getBondType().equals(Constants._5O4))
							g = generateStruct (g,"c5_"+gp,"O2_"+gc,cN.getBondType());
					}	
					else if (!cc.contains(gc))
					{
						cc.add(gc);
						//System.out.println(gc +"-"+gctype);
						//System.out.println("Entered 3 "+ cN.getBondType());
						MonolignolBase g2 = MonolignolBase.getMonolignolUnit(gctype, gc);
						//System.out.println("Entered 3");						
						g = genStructTwo(g,g2,cN.getBondType(),gp,gc);								
						
					}	
				}
			 }
			}
			i++;
			//6System.out.println(lignin);
			if (g != null)
			{
				//System.out.println(lignin);
				
				g = updatedStruct(g);
				
				String smile = g.getSmile(g.getMol());
				double mwt = g.getMolWt(g.getMol());
				Map bndset = lignin.getBndCnts();
				//System.out.println("Bond "+ bndset +"/"+ bondkeys+"/"+bondkeys.contains(bndset));
				if (!ligninSmiles.contains(smile) ) 
				{
					ligninSmiles.add(smile);				
					molweight.add(mwt);
					bondkeys.add(lignin.getBndCnts());
					//writer.println("Lignin="+i+",SMILEString="+g.getSmile(g.getMol()));
					
					SingleStructDet oneChain = jsonObjLignin.getLignin();
					oneChain.setBndCnts(lignin.getBndCnts());
					oneChain.setMolWeight(g.getMolWt(g.getMol()));
					oneChain.setSmilestring(smile);		
					if (oneChain.getSmilestring() != null)
					{
						jsonObjLignin.getLigninchains().add(oneChain);
					}
					
					if (lignin.getBndCnts().containsKey(Constants._55) || lignin.getBndCnts().containsKey(Constants._4O5))
					{
						g.toIMage(dp+"_"+i, g.getMol());
						g.toCML(g.getMol(), dp+"_"+totlig);
						System.out.println(lignin);
						lignin.getAdjAndConnMtx().writetoCSV(dp+"_"+i);
					}
					//System.out.println("Molecular Weight :"+ g.getMolWt(g.getMol()));
					lignin.setMolWeight(g.getMolWt(g.getMol()));
					totlig++;	
					
					//System.out.println("Lignin="+i+",MolWt="+g.getMolWt(g.getMol())+"Struct="+lignin+"SMILEString="+g.getSmile(g.getMol()));
					//System.out.println(lignin.getAdjAndConnMtx().toString());
					
				}					
			}
							
		}
		if (ligninSmiles.size()>0)
		{
			jsonObjLignin.setNunber_of_structs(totlig);
			jsonObjLignin.createJSON(jsonObjLignin, writer);
		}
		
		Date dt2 = new Date();		
		System.out.println("Time taken="+ getDifference(dt,dt2));
		writer.close();
		System.out.println("Number of polymers generated="+totlig);
	}
	
	private static MonolignolBase genStructTwo(MonolignolBase g, MonolignolBase g2, String bondType, int gp, int gc) {
		
		if (bondType.equals(Constants.BO4))
			g = generateStruct (g,g2,"bC_"+gp,"O2_"+gc, bondType);
		else if (bondType.equals(Constants._4OB))
			g = generateStruct (g,g2,"O2_"+gp,"bC_"+gc, bondType);
		else if (bondType.equals(Constants.BB))
			g = generateStruct (g,g2,"bC_"+gp,"bC_"+gc, bondType);								
		else if (bondType.equals(Constants._55))
			g = generateStruct (g,g2,"c5_"+gp,"c5_"+gc, bondType);
		else if (bondType.equals(Constants.AO4))
			g = generateStruct (g,g2,"aC_"+gp,"O2_"+gc, bondType);
		else if (bondType.equals(Constants.B5))
			g = generateStruct (g,g2,"bC_"+gp,"c5_"+gc, bondType);
		else if (bondType.equals(Constants._5B))
			g = generateStruct (g,g2,"c5_"+gp,"bC_"+gc, bondType);
		else if (bondType.equals(Constants._4O5))
			g = generateStruct (g,g2,"O2_"+gp,"c5_"+gc, bondType);		
		else if (bondType.equals(Constants._5O4))
			g = generateStruct (g,g2,"c5_"+gp,"O2_"+gc,bondType);
		
		return g;
		
	}

	static MonolignolBase updatedStruct(MonolignolBase g) throws CDKException
	{
		IAtomContainer g1_mol = g.getMol();	
		for (int i=0;i<g1_mol.getAtomCount();i++)
		{			
			g1_mol.getAtom(i).setAtomTypeName(g1_mol.getAtom(i).getSymbol());			
		}		
		//g.toIMage(10, g.getMol());
		
		boolean isConnected = ConnectivityChecker.isConnected(g1_mol);
		System.out.println("isConnected="+isConnected);

	
		CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(g1_mol.getBuilder());
	    for (IAtom atom : g1_mol.atoms()) 
	    {
	     //System.out.println(atom);
	     IAtomType type = matcher.findMatchingAtomType(g1_mol, atom);
	     AtomTypeManipulator.configure(atom, type);
	    }		
		
		return g;
	}
	
	
	static MonolignolBase generateStruct(MonolignolBase g, MonolignolBase g2, String b1, String b2, String bondType)
	{
		IAtomContainer g1_mol = g.getMol();		
		IAtomContainer g2_mol = g2.getMol();
		IAtom atm1 = null, atm2 = null, alphaC1 = null, gammaCO1=null, alphaC2=null, gammaCO2=null, c4O=null, o5=null;
		IBond bondalpha_OH = null;
		//System.out.println("two="+b1+"-"+b2);
		
		String[] mol1= b1.split("_");
		String[] mol2= b2.split("_");
		//System.out.println(""+Arrays.toString(mol1)+","+Arrays.toString(mol2));
		for (int i=0;i<g1_mol.getAtomCount();i++)
		{
			// Changing double bond for B position bonding in Monomer1 or oligomer 
			if (g1_mol.getAtom(i).getAtomTypeName()!=null  && g1_mol.getAtom(i).getAtomTypeName().equals(b1) )
			{
				atm1 = 	g1_mol.getAtom(i);				
				if (i>1)
				{
					IAtom atm3 = g1_mol.getAtom(i-1);		
					if (mol1[0].equals("bC"))  // betaCarbon in the side chain
					{
						alphaC1=atm3;
						gammaCO1 = g1_mol.getAtom(i+2);
						//System.out.println("1Bond between (bC)="+atm1.getAtomTypeName()+","+alphaC1.getAtomTypeName()+","+gammaCO1.getAtomTypeName());
					}
					if (mol1[0].equals("c5"))  // 5th Carbon in the ring
					{
						c4O = g1_mol.getAtom(i+6);			
						//System.out.println("1Bond between (c5)="+atm1.getAtomTypeName()+","+c4O.getAtomTypeName());
					}
					
					//System.out.println("atm3="+atm3.getAtomTypeName()+",atm1="+atm1.getAtomTypeName());
					IBond bnd2 = g1_mol.getBond(atm3, atm1);					
					if (bnd2!=null && bnd2.getOrder().equals(IBond.Order.DOUBLE))
					{
						bnd2.setOrder(IBond.Order.SINGLE);
						o5 = new Atom(8);
						o5.setAtomTypeName("O5_"+i);
						bondalpha_OH = new Bond(atm3, o5);	
					}					
				}				
				
			}
		}
		
		for (int i=0;i<g2_mol.getAtomCount();i++)
		{
			// Changing double bond for B position in the attaching monomer
			if (g2_mol.getAtom(i).getAtomTypeName()!= null && g2_mol.getAtom(i).getAtomTypeName().equals(b2))
			{
				atm2 = 	g2_mol.getAtom(i);
				if (mol2[0].equals("bC")) // betaCarbon in the side chain
				{
					alphaC2 = g2_mol.getAtom(i-1);;
					gammaCO2 = g2_mol.getAtom(i+2);
					//System.out.println("Bond between="+atm2.getAtomTypeName()+","+alphaC2.getAtomTypeName()+","+gammaCO2.getAtomTypeName());
				}
				if (mol2[0].equals("c5"))  // 5th Carbon in the ring
				{
					c4O = g2_mol.getAtom(i+6);			
					//System.out.println("Bond between (c5)="+atm2.getAtomTypeName()+","+c4O.getAtomTypeName());
				}
				if (i>1)
				{
					IAtom atm3 = g2_mol.getAtom(i-1);								
					IBond bnd2 = g2_mol.getBond(atm3, atm2);
					//System.out.println("Bond between="+atm2.getAtomTypeName()+","+atm3.getAtomTypeName());
					if (bnd2!=null && bnd2.getOrder().equals(IBond.Order.DOUBLE))
					{
						bnd2.setOrder(IBond.Order.SINGLE);
						o5 = new Atom(8);
						o5.setAtomTypeName("O5_"+i);
						bondalpha_OH = new Bond(atm3, o5);	
					}
				}
			}
			g1_mol.addAtom(g2_mol.getAtom(i));		
			//System.out.println(b2+"-"+g2_mol.getAtom(i).getAtomTypeName());
		}
		
		for (int i=0;i<g2_mol.getBondCount();i++)
		{
			g1_mol.addBond(g2_mol.getBond(i));			
		} 
		IBond bnd = new Bond(atm1,atm2);
		if (bondType.equals(Constants.BB))
		{
			IBond bnd2 = new Bond(alphaC1,gammaCO2);
			IBond bnd3 = new Bond(alphaC2,gammaCO1);
			g1_mol.addBond(bnd2);
			g1_mol.addBond(bnd3);
		}
		if (bondType.equals(Constants.B5))
		{
			IBond bnd2 = new Bond(alphaC1,c4O);			
			g1_mol.addBond(bnd2);			
		}
		if (bondType.equals(Constants._5B))
		{
			IBond bnd2 = new Bond(c4O,alphaC2);			
			g1_mol.addBond(bnd2);			
		}
		g1_mol.addBond(bnd);
		if (bondalpha_OH != null && bondType.equals(Constants.BO4)) 
		{
			g.getMol().addAtom(o5);
			g.getMol().addBond(bondalpha_OH);
		}
		return g;
	}
	
	static MonolignolBase generateStruct(MonolignolBase g, String b1, String b2, String bondType)
	{
		IAtomContainer g1_mol = g.getMol();		
		IAtom atm1 = null, atm2 = null, atm3=null, alphaC1 = null, gammaCO1=null, alphaC2=null, gammaCO2=null, c4O=null, o5=null;
		IBond bondalpha_OH = null;
		//System.out.println("Single="+b1+"-"+b2);
		
		String[] mol1= b1.split("_");
		String[] mol2= b2.split("_");
		
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
					//System.out.println("2 Bond between="+atm1.getAtomTypeName()+","+atm3.getAtomTypeName());			
					if (bnd2!=null && bnd2.getOrder().equals(IBond.Order.DOUBLE))
					{
						bnd2.setOrder(IBond.Order.SINGLE);
						o5 = new Atom(8);
						o5.setAtomTypeName("O5_"+i);
						bondalpha_OH = new Bond(atm3, o5);						
					}
				}
			}			
		}
		
		for (int i=0;i<g1_mol.getAtomCount();i++)
		{
			if (g1_mol.getAtom(i).getAtomTypeName()!= null && g1_mol.getAtom(i).getAtomTypeName().equals(b2))
			{
				atm2 = 	g1_mol.getAtom(i);
				//System.out.println("Bond between="+atm1.getAtomTypeName()+","+atm2.getAtomTypeName());
				//System.out.println(b2+"-"+atm2.getAtomTypeName());				
			}						
		}	
		IBond bnd = new Bond(atm1,atm2);
		g.getMol().addBond(bnd);		
		if (bondalpha_OH != null && bondType.equals(Constants.BO4)) 
		{
			g.getMol().addAtom(o5);
			g.getMol().addBond(bondalpha_OH);
		}
		return g;
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
	 
	
}
