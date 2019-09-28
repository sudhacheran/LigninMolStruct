package com.mol.construct;


import java.util.ArrayList;
import java.util.Arrays;

public class monolignol {
	static int NoOfBond = 3;
	static int id=0;
	String type;
	monolignol bondl;
	monolignol bond2;
	monolignol bond3;
	int objId;	
	String[] edgeType = new String[5];
	ArrayList<String> edWgt= new ArrayList<String>();
	Labels lbl;
	boolean isEdgUpdated = false;	
		
	public monolignol(String type) {
		super();
		this.type = type;
		this.objId = ++id;	
		lbl = new Labels(this.objId,this.type);
	}
	
	
	public monolignol(String type, int idV) {
		super();
		this.type = type;
		this.objId = idV;	
		lbl = new Labels(this.objId,this.type);
	}
	
	public String getType() {
		return type;
	}
	public void setType(String type) {
		this.type = type;
	}
	public monolignol getBondl() {
		return bondl;
	}
	public void setBondl(monolignol bondl) {
		if (bondl != null) 
		this.bondl = bondl;
	}
	public monolignol getBond2() {
		
		return bond2;
	}
	public void setBond2(monolignol bond2) {
		if (bond2 != null) 
		this.bond2 = bond2;
	}
	public monolignol getBond3() {
		return bond3;
	}
	public void setBond3(monolignol bond3) {
		if (bond3 != null) 
		this.bond3 = bond3;
	}
	
	public ArrayList<String> getEdWgt() {
		return edWgt;
	}
	public void setEdWgt(ArrayList<String> edWgt) {
		this.edWgt = edWgt;
	}
	public String[] getEdgeType() {
		return edgeType;
	}
	public void setEdgeType(String[] edgeType) {
		this.edgeType = edgeType;
	}
	@Override
	public String toString() {
		return "monolignol [type=" + type + ", bondl=" + ((bondl==null)?"":bondl.objId+bondl.getType()) + ", bond2=" + ((bond2==null)?"":bond2.objId+bond2.getType())  + ", bond3=" + ((bond3==null)?"":bond3.objId+bond3.getType())  + ", objId="
				+ objId + ", edgeType=" + Arrays.toString(edgeType) + ", lbl=" + lbl + "]\n";
	}
		 

}


