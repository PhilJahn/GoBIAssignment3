import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

import AugmentedTree.IntervalTree;

public class Transcript extends RegionVector{
	
	private IntervalTree<Region> introns; 
	
	private IntervalTree<RegionBlock> exons; 

	public Transcript(int start, int stop, Annotation annotation)  {
		super(start, stop, annotation);
		exons= new IntervalTree<RegionBlock>();
	}
	
	public int getProtNum(){
		ArrayList<Annotation> annos = new ArrayList<Annotation>();
		for(Region r: this.getRegionsTree()){
			Annotation anno = r.getAnnotation();
			if(!annos.contains(anno)){
				annos.add(anno);
			}
		}
		return annos.size();
		
	}
	
	public boolean setIntrons(){
		if(this.getRegionsTree().size() > 1){
			introns = this.invert().getRegionsTree();
			return true;
		}
		else{
			return false;
		}
	}
	
	public IntervalTree<Region> getIntrons(){
		return introns;
	}
	
	public boolean equals(Object o){
		if(o.getClass() == this.getClass()){
			return equals((Transcript) o);
		}
		return false;
	}
	
	public boolean equals(Transcript t){
		return this.getAnnotation().equals(t.getAnnotation());
	}
	
	public boolean add(RegionBlock exon){
		return exons.add(exon);
	}
	
	public int hashCode(){
		return this.getAnnotation().hashCode();
	}
	
	public String toString(){
		return "Transcript: " + super.toString();
	}

	public boolean inTranscript(Read curRead) {
		
//		System.out.println(curRead.toString());
		
		
		boolean in = true;
		ArrayList<RegionBlock> exonsfop = curRead.getAlignmentBlocksFoP();
		ArrayList<RegionBlock> exonssop = curRead.getAlignmentBlocksSoP();
		ArrayList<RegionBlock> intronsR = curRead.getIntronBlocks();
		
//		if(curRead.getReadName().equals("7075")){
//			System.out.println("7075");
//			System.out.println(this.exons.toString());
//			System.out.println("FOP " + exonsfop.toString());
//			System.out.println("SOP " +exonssop.toString());
//			System.out.println("Intron " +intronsR.toString());
//		}
//		
		in = inTranscriptExon(exonsfop);
//		if(curRead.getReadName().equals("7075")){
//			System.out.println("FOP " + in);
//		}
		if(in){
			in = inTranscriptExon(exonssop);
//			if(curRead.getReadName().equals("7075")){
//				System.out.println("SOP " + in);
//			}
		}
		
		if(in && intronsR.size() > 0){
			HashSet<RegionBlock> cont;
			for(RegionBlock intron: intronsR){
				cont = exons.getIntervalsSpannedBy(intron.getStart()-1,intron.getStop()-1, new HashSet<RegionBlock>());
				if(cont.size() != 0){
//					if(curRead.getReadName().equals("7075")){
//						System.out.println("Intron " + false);
//					}
					return false;
				}
			}
//			if(curRead.getReadName().equals("7075")){
//				System.out.println("Intron " + in);
//			}
		}
		
		
		return in;
	}
	
	private boolean inTranscriptExon(ArrayList<RegionBlock> exons){
		boolean in = true;
		ArrayList<RegionBlock> cont;
		
		
		RegionBlock curExon = exons.get(0);
		
		if(exons.size() > 1){
			cont = this.exons.getIntervalsEndAt(curExon.getStop()-1, new ArrayList<RegionBlock>());
			if(cont.size() > 0){
				in &= cont.get(0).getStart() <= curExon.getStart();
			}
			else{
				return false;
			}
			for(int i = 1; i < exons.size()-1; i++){
				if(in){
					curExon = exons.get(i);
					cont = this.exons.getIntervalsEqual(curExon.getStart(),curExon.getStop()-1, new ArrayList<RegionBlock>());
					in &= cont.size() > 0;
				}
			}
			if(in){
				curExon = exons.get(exons.size()-1);
				cont = this.exons.getIntervalsBeginAt(curExon.getStart(), new ArrayList<RegionBlock>());
				if(cont.size() > 0){
					in &= cont.get(0).getStop() >= curExon.getStop()-1;
				}
				else{
					return false;
				}
			}
		}
		else{
			cont = this.exons.getIntervalsSpanning(curExon.getStart(),curExon.getStop()-1, new ArrayList<RegionBlock>());
			in &= cont.size() > 0;
		}
		return in;
	}
	
}
