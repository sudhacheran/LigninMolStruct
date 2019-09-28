package com.mol.construct;
import java.util.List;

public class Lignin {

	
	String PType;
	List<childNode> cNode;
	public String getPType() {
		return PType;
	}
	public void setPType(String pType) {
		PType = pType;
	}
	public List<childNode> getcNode() {
		return cNode;
	}
	public void setcNode(List<childNode> cNode) {
		this.cNode = cNode;
	}
	@Override
	public String toString() {
		return "Lignin [PType=" + PType + ", cNode=" + cNode + "]";
	}
}

