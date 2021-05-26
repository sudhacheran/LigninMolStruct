package com.mol.construct;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class GeneratePolymer {
		static final int noOfPer_limit = 10000;
	

public static void main(String[] args) {
		
		int dp = 6;
		
		long noOfPer = factorial(dp);
		System.out.println("No of permutation="+noOfPer);

		String[] str = { "B04", "BB", "55", "B5", "B1", "405" };
		int[] bondper = { 50, 30, 20, 0, 0, 0 };

		String[] mono = { "G", "S", "H" };
		int[] monoPer = { 100, 0, 0 };
		
		GeneratePolymer gp = new GeneratePolymer();
		List<String> polymentString = gp.getPossiblePolymer(str, bondper, mono, monoPer, dp,noOfPer);
		
		for(int i=0;i<polymentString.size();i++)
		{
			System.out.println(polymentString.get(i));
		}
		
}

public List<String> getCorrectedPolymers(String[] str, int[] bondper, String[] mono, int[] monoPer, int dp, long noOfPer)
{
	List<String> polymerList = getPossiblePolymer (str, bondper, mono, monoPer, dp, noOfPer);	
	
	for(int i=0;i<polymerList.size();i++)
	{
		System.out.println(polymerList.get(i));
	}
			
	return polymerList;
}


	
/**
 * Method returns the possible polymer chain permutations (within the specified limit parameter)
 * @param str - Bond array
 * @param bondper - Bond percentage
 * @param mono - Monomer unit array
 * @param monoPer - Monomer percentage
 * @param dp  - Chain length (no. of monomer units)
 * @param noOfPer - Limiting the permutation
 * @return
 */
public List<String> getPossiblePolymer(String[] str, int[] bondper, String[] mono, int[] monoPer, int dp, long noOfPer)
{

		String[] strBonds = getList(str, bondper, dp);
		String[] monoBonds = getList(mono, monoPer, dp);

		System.out.println("No of Bonds=" + strBonds.length);
		System.out.println("No of Monomers=" + monoBonds.length);

		List<String[]> bondList = new ArrayList<String[]>();

		int n = strBonds.length;
		permute(strBonds, 0, n - 1, bondList, noOfPer);

		List<String[]> monoList = new ArrayList<String[]>();
		n = monoBonds.length;
		permute(monoBonds, 0, n - 1, monoList, noOfPer);

		System.out.println("Bond Permutations" + bondList.size());
		System.out.println("Monomer Permutations" + monoList.size());
		List<String> polymer = new ArrayList<String>();

		//int d=0;
		for (String[] s4 : monoList) {
			for (String[] s2 : bondList) {
				//String[] s4 = monoList.get(d++);
				String polymerStr = "";
				int e = 0;
				for (String tt : s2) {
					polymerStr += s4[e++] + "-" + tt + "-";
				}		
				if (!polymer.contains(polymerStr))					
				polymer.add(polymerStr);
			}
		}
		//sop(""+polymer.size());
		return polymer;
	}

	/**
	 * Method returns the possible bond names and monomer names with in the dp 
	 * @param str - Bond or monomer array
	 * @param bondper - Bond or monomer percentile list
	 * @param dp  - length of the polymer chain
	 * @return
	 */
	public String[] getList(String[] str, int[] bondper, int dp) {
		
		String[] strBonds = new String[dp];
		int[] totcnt = new int[str.length];
		// initializing the totcount[] for each bond with 1  
		for(int inc=0; inc < totcnt.length; inc++)
		{
			totcnt[inc]+=1;
		}
		// Loop to create the possible bond names (or) monomer names with the dp count
		int cnt = 0;
		boolean brkLoop = false;
		while (!brkLoop) {
			for (int j = 0; j < str.length; j++) {
				double cal = calPer(bondper[j], dp);		
				if (cal > 0 && totcnt[j] <= cal & cnt < dp) {	
					sop(cal+"<"+totcnt[j]+"="+str[j]+" ="+cnt);
					strBonds[cnt++] = str[j];	
					totcnt[j] = totcnt[j]+1;
				}
			}
			sop("entered");
			brkLoop = true;
			for (int j = 0; j < str.length; j++) {
				double cal = calPer(bondper[j], dp);		
				if (cal > 0 && totcnt[j] <= cal && cnt < dp) {	
					brkLoop = false;
				}				
			}
		}
		//Loop to add bonds at the end for chain length
		//sop("cnt="+cnt);
		while (cnt < dp)
		{
			for (int k=0;k<str.length;k++)
	      	 {
	            if (cnt >= dp)
	            {
	              break;
	            }
	            strBonds[cnt++]=str[k];
	         }
		}
				
		return strBonds;
	}

	/**
	 * This is a recursive method to create the permutations of bond (or) monomer array
	 * @param str -  Bond array / Monomer array
	 * @param l   - Initial index of the array to permute
	 * @param r   - end index of the array to permute
	 * @param bondList  - List to store the possible combinations of Bond / Monomers
	 * @param noOfPermutation
	 */
	public static void permute(String[] str, int l, int r, List<String[]> bondList, long noOfPermutation) {
		if (bondList.size() >= noOfPermutation) {
			//System.out.print("Enteret empty return");
			return;
		}

		if (l == r) {
			String[] st3 = str.clone(); // cloning the string array to add it to ArrayList
			//if (!bondList.contains(st3)) {
				bondList.add(st3);
			//}
		} else
			for (int i = l; i <= r; i++) {
				if (r%2 == 0)	
					str = swap(str, i, l);
				else
					str = swap(str, l, i);
				permute(str, l + 1, r, bondList, noOfPermutation); // Recursive call
				str = swap(str, l, 0);
			}

	}

	/**
	 * Method to calculate the ratio
	 * @param val - Percentage value
	 * @param dp   - polymer chain length
	 * @return
	 */
	public double calPer(int val, int dp) {
		if (val == 0)
			return val;
		double flrval = Math.round(dp * (val / 100f));
		
		return ((flrval > 0) ? flrval : 1);
	}

	/**
	 * Method to swap tow values in a string array
	 * @param a
	 * @param i
	 * @param j
	 * @return
	 */
	public static String[] swap(String[] a, int i, int j) {
		String temp;
		temp = a[i];
		a[i] = a[j];
		a[j] = temp;
		return a;
	}
	
	public static long factorial(int n)
	{
		if (n > 25) return noOfPer_limit;
		long fact=1;
		for (int i=1;i<=n;i++)
		{
			fact=fact*i;
		}
		if (fact > noOfPer_limit) return noOfPer_limit;
		return fact;
	}

	public static void sop(String s) {
		System.out.println(s);
	}
}

