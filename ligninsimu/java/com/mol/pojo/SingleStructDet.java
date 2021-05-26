package com.mol.pojo;

import java.util.Map;

public class SingleStructDet {

	String smilestring;
	double molWeight;
	Map<String,Integer> bndCnts;
	public String getSmilestring() {
		return smilestring;
	}
	public void setSmilestring(String smilestring) {
		this.smilestring = smilestring;
	}
	public double getMolWeight() {
		return molWeight;
	}
	public void setMolWeight(double molWeight) {
		this.molWeight = molWeight;
	}
	public Map<String, Integer> getBndCnts() {
		return bndCnts;
	}
	public void setBndCnts(Map<String, Integer> bndCnts) {
		this.bndCnts = bndCnts;
	}
}
