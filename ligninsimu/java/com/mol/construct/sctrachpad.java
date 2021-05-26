package com.mol.construct;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class sctrachpad {




private static List<String[]> possibleBonds(List<String[]> bondListStr) {
	List<String[]> retArrr = new ArrayList<String[]>();
	for (String[] bnds : bondListStr)
	{
		char prevFirst = 0;
		char prevLast = 0;
		boolean ignore = false;
		sop(Arrays.toString(bnds));
		List<String> bondArr = new ArrayList<String>();
		for (String bd: bnds)
		{
			char first = bd.charAt(0);
			char last = bd.charAt(bd.length()-1);
			char middle =  (bd.length()>2 ? bd.charAt(1): 0);
			if (prevLast == first)
			{
				if (first != last)
				{
					char temp = first;
					first = last;
					last = temp;	
					String bond = String.valueOf(first)+(middle!=0?String.valueOf(middle):"")+String.valueOf(last);
					bondArr.add(bond);
				}						
			}
			else
			{
				bondArr.add(bd);
			}
			prevFirst = first;
			prevLast = last;
		}
		sop(bondArr.toString());
		if (bondArr.size() > 0)
		{
			retArrr.add(bondArr.toArray(new String[0]));
		}
		sop("Next>>");
		
	}
	return retArrr;
}

private static List<String[]> removeduplicates(List<String[]> monoListStr) {
	List<String> str = new ArrayList<String>();
	List<String[]> retArrr = new ArrayList<String[]>();
	for(String[] marr: monoListStr)
	{
		String s1="";
		for (String sr: marr)
		{ 
			s1 = s1.concat(sr);
		}
		
		if (!str.contains(s1))
		{
			retArrr.add(marr);
			str.add(s1);
		}
		
	}
	return retArrr;
}

public static void sop(String s)
{
	System.out.print(s);
}


}