package com.mol.mono;
import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;

public class Atom_Local {
	
	protected String lcl_atomName = "";
	
	protected IAtom lcl_Atom;

	public String getLcl_atomName(int atomNumber) {
		return lcl_atomName;
	}

	public void setLcl_atomName(String lcl_atomName) {
		this.lcl_atomName = lcl_atomName;
	}

	public IAtom getLcl_Atom(int i) {
		lcl_Atom = new Atom(i);
		return lcl_Atom;
	}

	public Atom_Local(String lcl_atomName, int i) {
		super();
		this.lcl_atomName = lcl_atomName;
		this.lcl_Atom = new Atom(i);
	}

	public void setLcl_Atom(IAtom lcl_Atom) {
		this.lcl_Atom = lcl_Atom;
	}
	
	

}
