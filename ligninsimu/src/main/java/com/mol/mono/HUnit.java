package com.mol.mono;


import java.io.IOException;
import java.util.ArrayList;


import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;


public class HUnit extends MonolignolBase {
	
	IAtomContainer mol = null;
	ArrayList<IAtom> bondingAtom = new ArrayList<IAtom>();
	
	public static void main(String s[]) throws IOException
	{
		HUnit hunit = new HUnit(1);
		hunit.toIMage(102, hunit.mol);		
	}

	public HUnit(int i) {
		super();
		mol = generateBaseMol(i);		
		setBondingAtom(bondingAtom);
		
	}
	
	public IAtomContainer getMol() {
		return mol;
	}


	public void setMol(IAtomContainer mol) {
		this.mol = mol;
	}
	
	public ArrayList<IAtom> getBondingAtom() {
		return bondingAtom;
	}
	
	private void setBondingAtom(ArrayList<IAtom> bondingAtom) {
		bondingAtom = new ArrayList<IAtom>();
		bondingAtom.add(c4);		
		bondingAtom.add(c1);
		bondingAtom.add(c5);
		bondingAtom.add(c3);
		bondingAtom.add(betaC);
		bondingAtom.add(alphaC);		
		this.bondingAtom = bondingAtom;
	}

}

