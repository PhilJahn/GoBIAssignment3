import java.util.ArrayList;
import java.util.Arrays;
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
	
	private IntervalTree<RegionBlock> mergedTrans;
	
	public Gene(int x1, int x2, Annotation annotation, String gene_biotype){
		super(x1,x2,annotation);
		transcripts = new HashMap <String,Transcript>();
		transReadIds = new  HashSet<String>();
		this.gene_biotype = gene_biotype;
		mergedTrans = new IntervalTree<RegionBlock>();
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

	public ArrayList<String> inTranscript(Read curRead) {
		ArrayList<String> transcriptList = new ArrayList<String>();
//		if(curRead.getReadName().equals("711")){
//			System.out.println("711");
//		}
		for( String k : transcripts.keySet() ){
			if(transcripts.get(k).inTranscript(curRead)){
				transcriptList.add(k);
			}
		}
		
		return transcriptList;
	}
	
	public boolean inMerged(Read curRead) {
		if(!(mergedTrans.size() > 0)){
			this.mergeTrans();
//			if(this.getAnnotation().getId().equals("ENSG00000227232") && curRead.getReadName().equals("4995990")){
//				System.out.println("I was here");
//			}
		}
//		if(this.getAnnotation().getId().equals("YBL087C") && curRead.getReadName().equals("4995990")){
//			System.out.println("YBL087C");
//			
//			for( String k : transcripts.keySet() ){
//				System.out.println(k);
//				System.out.println(transcripts.get(k).getExons());
//			}
//			
//			System.out.println("Merged");
//			System.out.println(this.mergedTrans.toString());
//			System.out.println(this.mergedTrans.size());
//			System.out.println("4995990");
//			System.out.println(curRead.getAlignmentBlocks().toString());
//		}
		ArrayList<RegionBlock> readRV =  new ArrayList<RegionBlock>(curRead.getAlignmentBlocks());
		
		HashSet<RegionBlock> cont;
		for(RegionBlock rb : readRV){
			cont = new HashSet<RegionBlock>();
			cont = mergedTrans.getIntervalsSpanning(rb.getStart(), rb.getStop()-1, cont);
			if(cont.size() == 0){
				return false;
			}
		}
		
		return true;
	}
	
	private void mergeTrans(){
		IntervalTree<RegionBlock> mergeRegions = new IntervalTree<RegionBlock>();
		for( String k : transcripts.keySet() ){
			mergeRegions.addAll(transcripts.get(k).getExons());
		}
		
//		if(this.getAnnotation().getId().equals("ENSG00000227232")){
//			System.out.println("Merge of ENSG00000227232");
//			System.out.println(mergeRegions.toString());	
//			
//		}
		
		Iterator<Set<RegionBlock>> mergereg_it = mergeRegions.groupIterator();
		while(mergereg_it.hasNext()){
			Set<RegionBlock> overlap = mergereg_it.next();
			RegionBlock[] overlapArray = new RegionBlock[overlap.size()];
			overlapArray = overlap.toArray(overlapArray);
			Arrays.sort(overlapArray,new StartRegionBlockComparator());
			int start = overlapArray[0].getStart();
			Arrays.sort(overlapArray,new StopRegionBlockComparator());
			int stop = overlapArray[overlapArray.length-1].getStop();
			
			RegionBlock newblock = new RegionBlock(start,stop);
			mergedTrans.add(newblock);
		}
	}
}

class StartRegionBlockComparator implements Comparator<RegionBlock>
{
    public int compare(RegionBlock x1, RegionBlock x2)
    {
        return x1.getStart() - x2.getStart();
    }
}

class StopRegionBlockComparator implements Comparator<RegionBlock>
{
    public int compare(RegionBlock x1, RegionBlock x2)
    {
        return x1.getStop() - x2.getStop();
    }
}
