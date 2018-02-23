import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import AugmentedTree.IntervalTree;

public class Gene extends Region{

	
//	public static void main(String[] args){
//		Annotation a = new Annotation("id", "name", "str", '-', "test");
//		IntervalTree<Region> test = new IntervalTree<Region>();
//		Region r1 = new Region(1,20,a);
//		Region r2 = new Region(1,8,a);
//		Region r3 = new Region(0,10,a);
//		Region r4 = new Region(7,10,a);
//		Region r5 = new Region(6,10,a);
//		test.add(r1);
//		test.add(r2);
//		test.add(r3);
//		test.add(r4);
//		test.add(r5);
//		ArrayList<Region> posCheck = new ArrayList<Region>();
//		posCheck = test.getIntervalsEndAt(10, posCheck);
//		System.out.println(test.toTreeString());
//		System.out.println(posCheck.toString());
//	}
	
	private HashMap<String,Transcript> transcripts;
	
	private HashSet<String> transReadIds;
	
	private String gene_biotype;
	
	public Gene(int x1, int x2, Annotation annotation, String gene_biotype){
		super(x1,x2,annotation);
		transcripts = new HashMap <String,Transcript>();
		transReadIds = new  HashSet<String>();
		this.gene_biotype = gene_biotype;
	}

//	public Gene(int x1, int x2, Annotation annotation, Annotation subannotation){
//		super(x1,x2,annotation);
//		transcripts = new HashMap <String,Transcript>();
//		transcripts.put(subannotation.getId(),new Transcript(x1,x2,subannotation));
//	}

	public int getTranscriptNumber(){
		int n = 0;
		for( String k : transcripts.keySet() ){
			if(transcripts.get(k).getRegionsTree().size() > 0){
				n ++;
			}
		}
		return n;
	}
	
	public int getProteinNumber(){
		int n = 0;
		for( String k : transcripts.keySet() ){
			n += transcripts.get(k).getProtNum();
		}
		return n;
	}
	
	public void removeEmpty(){
		for( String k : transcripts.keySet() ){
			if(transcripts.get(k).getRegionsTree().size() == 0){
				transcripts.remove(k);
			}
		}
	}
	
	public int hashCode(){
		return this.getAnnotation().hashCode();
	}
	
	public String toString(){
		String output = "Gene: " + this.getAnnotation().getId() + " " + super.toString() + "\n";
		for( String k : transcripts.keySet() ){
			output += transcripts.get(k).toString();
		}
		return output;
	}

	public void add(Transcript trans) {
		this.transcripts.put(trans.getAnnotation().getId(), trans);
	}
	
	class StartRegionComparator implements Comparator<Region>
	{
	    public int compare(Region x1, Region x2)
	    {
	        return x1.getStart() - x2.getStart();
	    }
	}
	
	public HashSet<String> getTransReadIds(){
		return transReadIds;
	}
	
	public String getBiotype(){
		return gene_biotype;
	}

	public boolean inTranscript(Read curRead) {
		boolean transcriptomic = false;
//		if(curRead.getReadName().equals("711")){
//			System.out.println("711");
//		}
		for( String k : transcripts.keySet() ){
			if(!transcriptomic){
				transcriptomic |= transcripts.get(k).inTranscript(curRead);
			}
		}
//		if(curRead.getReadName().equals("7075")){
//			System.out.println(transcriptomic);
//		}
		if(transcriptomic){
			this.transReadIds.add(curRead.getReadName());
		}
		
		return transcriptomic;
	}
}

