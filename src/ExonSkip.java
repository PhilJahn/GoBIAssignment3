import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import AugmentedTree.Interval;
import AugmentedTree.IntervalTree;

public class ExonSkip implements Interval{

	private HashSet<String> sv_prot;
	private HashSet<String> wt_prot;
	
	private int start;
	private int stop;
	private HashSet<Region> wt_introns;
	private HashSet<Region> wt_exons;
	
	private int minEx;
	private int maxEx;
	
	private int minBase;
	private int maxBase;
	
	
	public ExonSkip(HashSet<String> sv_prot, HashSet<String> wt_prot, HashSet<Region> wt_introns, HashSet<Region> wt_exons, int minEx, int maxEx, int minBase, int maxBase, int start, int stop){
		this.sv_prot = sv_prot;
		this.wt_prot = wt_prot;
		this.wt_introns = wt_introns;
		this.wt_exons = wt_exons;
		
		this.start = start;
		this.stop = stop;
		
		this.minEx = minEx;
		this.maxEx = maxEx;
		
		this.minBase = minBase;
		this.maxBase = maxBase;
	}
	

		
	public HashSet<String> getSVProt(){
		return sv_prot;
	}
	
	public HashSet<String> getWTProt(){
		return wt_prot;
	}
	
	
	public int getStart(){
		return start;
	}
	
	public int getStop(){
		return stop;
	}
	
	public HashSet<Region> getWTIntrons(){
		return wt_introns;
		
	}
	
	public HashSet<Region> getWTExons(){
		return wt_exons;
	}	
	
	public int getMinEx(){
		return minEx;
	}
	
	public int getMaxEx(){
		return maxEx;
	}
	
	public int getMinBase(){
		return minBase;
	}
	
	public int getMaxBase(){
		return maxBase;
	}
	
	public String toString(){
		return start + ":" + stop + " sv: " + sv_prot.toString() + " wt: " + wt_prot.toString();
	}
	
}
