package com.mol.pojo;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;

public class LigninsJson {
	
	String monoType;
	double S_G_Ratio;
	String name;
	int dp;
	int nunber_of_structs;
	SingleStructDet lignin;	
	List<SingleStructDet> ligninchains;
	
	
	public String getMonoType() {
		return monoType;
	}
	public void setMonoType(String monoType) {
		this.monoType = monoType;
	}
	public double getS_G_Ratio() {
		return S_G_Ratio;
	}
	public void setS_G_Ratio(double s_G_Ratio) {
		S_G_Ratio = s_G_Ratio;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public int getDp() {
		return dp;
	}
	public void setDp(int dp) {
		this.dp = dp;
	}
	public int getNunber_of_structs() {
		return nunber_of_structs;
	}
	public void setNunber_of_structs(int nunber_of_structs) {
		this.nunber_of_structs = nunber_of_structs;
	}
	public List<SingleStructDet> getLigninchains() {
		if (ligninchains == null) ligninchains = new ArrayList<SingleStructDet>();
		return ligninchains;
	}
	public void setLigninchains(List<SingleStructDet> ligninchains) {
		this.ligninchains = ligninchains;
	}
	public SingleStructDet getLignin() {
		return new SingleStructDet();
	}
	public void setLignin(SingleStructDet lignin) {
		this.lignin = lignin;
	}
	
	public void createJSON(Object ligninDet, PrintWriter jsonfile)
	{
		  ObjectMapper mapper = new ObjectMapper();
	      //Converting the Object to JSONString
	      String jsonString;
		try {
			jsonString = mapper.writeValueAsString(ligninDet);
			jsonfile.write(jsonString);
			System.out.println("Successfully Copied JSON Object to File...");
			System.out.println(jsonString);
		
		} catch (JsonProcessingException e) {
			
			e.printStackTrace();
		}
	}
	      
	
}





