package com.mol.mono;

import org.openscience.cdk.Atom;

public class LcAtom extends Atom   {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	protected String atomId;

	public LcAtom(int i) {
		super.atomicNumber = i;
	}

	public String getAtomId() {
		return atomId;
	}

	public void setAtomId(String atomId) {
		this.atomId = atomId;
	}

}
