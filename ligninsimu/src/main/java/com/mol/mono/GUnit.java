package com.mol.mono;

import java.io.IOException;
import java.util.ArrayList;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

public class GUnit extends MonolignolBase {
	
	IAtomContainer mol = null;
	ArrayList<IAtom> bondingAtom = new ArrayList<IAtom>();
	public static void main(String s[]) throws IOException, CDKException
	{
		GUnit gunit = new GUnit(1);
		gunit.toIMage(101, gunit.mol);		
	}
	
	
	public IAtomContainer getMol() {
		return mol;
	}


	public void setMol(IAtomContainer mol) {
		this.mol = mol;
	}

	public GUnit(int i)  {
		super();
		mol = generateBaseMol(i);
		IAtom o3,c7,c3 = null;
		IBond b13,b14;
		
		
		c7 = new Atom("C");	
		
		c7.setID("c7_"+i);
		
		o3 = new Atom("O");
		o3.setID("O3_"+i);
		
		
		for (int i1=0;i1<mol.getAtomCount();i1++)
		{
			if (mol.getAtom(i1).getID()!=null  && mol.getAtom(i1).getID().equals("c5_"+i))
			{
				c5 = mol.getAtom(i1);
			}
			if (mol.getAtom(i1).getID()!=null  && mol.getAtom(i1).getID().equals("c3_"+i))
			{
				c3 = mol.getAtom(i1);
			}
		}
				
		b13 = new Bond(c3, o3);
		b14 = new Bond(o3, c7);
		mol.addAtom(c7);		
		mol.addAtom(o3);		
		mol.addBond(b13);
		mol.addBond(b14);
		//mol = configureMol(mol);
		setBondingAtom(bondingAtom);
		
	}
	
	public ArrayList<IAtom> getBondingAtom() {
		return bondingAtom;
	}
	
	private void setBondingAtom(ArrayList<IAtom> bondingAtom) {
		bondingAtom = new ArrayList<IAtom>();
		bondingAtom.add(c4);		
		bondingAtom.add(c1);
		bondingAtom.add(betaC);
		bondingAtom.add(alphaC);		
		this.bondingAtom = bondingAtom;
	}

}
