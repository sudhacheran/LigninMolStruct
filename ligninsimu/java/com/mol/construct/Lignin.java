package com.mol.construct;
import java.util.List;
import java.util.Map;

import com.mol.util.AdjacencyAndConnectivityMatrix;

public class Lignin {

	
	String PType;
	List<childNode> cNode;
	Map<String,Integer> bndCnts;
	AdjacencyAndConnectivityMatrix adjAndConnMtx;	
	double molWeight;
	
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
	
	
	public AdjacencyAndConnectivityMatrix getAdjAndConnMtx() {
		return adjAndConnMtx;
	}
	public void setAdjAndConnMtx(AdjacencyAndConnectivityMatrix adjAndConnMtx) {
		this.adjAndConnMtx = adjAndConnMtx;
	}
	@Override
	public String toString() {
		return "Lignin [PType=" + PType + ", cNode=" + cNode + "] \n Bond counts: "+bndCnts;
	}
}

