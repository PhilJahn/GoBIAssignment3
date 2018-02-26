import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.SortedSet;

import AugmentedTree.Interval;
import AugmentedTree.IntervalTree;
import net.sf.samtools.*;

public class Read implements Interval{

	private int start;
	private int stop;
	private String readname;
	//false =forward, true = reverse
	private char strand;
	private String chr;
	private boolean incons;
	private int mm;
	private int clip;
	private int split;
	
	private IntervalTree<RegionBlock> alignment_blocks;
	
	private ArrayList<RegionBlock> exons_fop;
	private ArrayList<RegionBlock> exons_sop;
	
	
	
	private ArrayList<RegionBlock> introns;
	
	public Read(SAMRecord samRecord1, SAMRecord samRecord2) {
		SAMRecord fop;
		SAMRecord sop;
		if(samRecord1.getFirstOfPairFlag()){
			fop = samRecord1;
			sop = samRecord2;
		}
		else{
			sop = samRecord1;
			fop = samRecord2;
		}
		
		
//		mm = 0;
//		if( fop.getAttribute("NM") != null ){
//			mm += (Integer) fop.getAttribute("NM");
//		}
//		if( fop.getAttribute("nM") != null ){
//			mm += (Integer) fop.getAttribute("nM");
//		}
//		if( fop.getAttribute("XM") != null ){
//			mm += (Integer) fop.getAttribute("XM");
//		}
//		if( sop.getAttribute("NM") != null ){
//			mm += (Integer) sop.getAttribute("NM");
//		}
//		if( sop.getAttribute("nM") != null ){
//			mm += (Integer) sop.getAttribute("nM");
//		}
//		if( sop.getAttribute("XM") != null ){
//			mm += (Integer) sop.getAttribute("XM");
//		}
		
		Integer mmfop = (Integer) fop.getAttribute("NM");
		mmfop = (mmfop != null) ? mmfop : (Integer) fop.getAttribute("nM");
		mmfop = (mmfop != null) ? mmfop : (Integer) fop.getAttribute("XM");
		
		Integer mmsop = (Integer) sop.getAttribute("NM");
		mmsop = (mmsop != null) ? mmsop : (Integer) sop.getAttribute("nM");
		mmsop = (mmsop != null) ? mmsop : (Integer) sop.getAttribute("XM");
		
		this.mm = mmfop + mmsop;
		
		this.start = Math.min(fop.getAlignmentStart(), sop.getAlignmentStart());
		int o_start = Math.max(fop.getAlignmentStart(), sop.getAlignmentStart());
		this.stop = Math.max(fop.getAlignmentEnd(), sop.getAlignmentEnd());
		int o_stop = Math.min(fop.getAlignmentEnd(), sop.getAlignmentEnd());
		
		readname = fop.getReadName();
		
		chr = fop.getReferenceName();
		
		if(fop.getReadNegativeStrandFlag()){
			strand = '-';
		}
		else{
			strand = '+';
		}
		
		ArrayList<RegionBlock> exons_fop = new ArrayList<RegionBlock>();
		for (AlignmentBlock block : fop.getAlignmentBlocks()){
			int start = block.getReferenceStart();
			int stop = start + block.getLength();
			exons_fop.add(new RegionBlock(start,stop));
		}
		
		
//		exons_fop.clear();
//		
//		exons_fop.add(new RegionBlock(114357090, 114347830));
//		exons_fop.add(new RegionBlock(114347890, 114347902));
		
		exons_fop.sort(new StartRegionBlockComparator());
		
		IntervalTree<RegionBlock> exons_fop_tree = new IntervalTree<RegionBlock>(exons_fop);
		exons_fop.clear();
		Iterator<Set<RegionBlock>> exon_fop_it = exons_fop_tree.groupIterator();
		while(exon_fop_it.hasNext()){
			Set<RegionBlock> overlap = exon_fop_it.next();
			RegionBlock[] overlapArray = new RegionBlock[overlap.size()];
			overlapArray = overlap.toArray(overlapArray);
			Arrays.sort(overlapArray,new StartRegionBlockComparator());
			int start = overlapArray[0].getStart();
			Arrays.sort(overlapArray,new StopRegionBlockComparator());
			int stop = overlapArray[overlapArray.length-1].getStop();
			
			RegionBlock newblock = new RegionBlock(start,stop);
			exons_fop.add(newblock);
		}
		
		exons_fop.sort(new StartRegionBlockComparator());
		
		ArrayList<RegionBlock> introns = new ArrayList<RegionBlock>();
		for(int i =1 ; i < exons_fop.size(); i++){
			int start = exons_fop.get(i-1).getStop();
			int stop = exons_fop.get(i).getStart();
			introns.add(new RegionBlock(start+1,stop));
		}

		ArrayList<RegionBlock> exons_sop = new ArrayList<RegionBlock>();
		for (AlignmentBlock block : sop.getAlignmentBlocks()){
			int start = block.getReferenceStart();
			int stop = start + block.getLength();
			exons_sop.add(new RegionBlock(start,stop));
		}
		
//		exons_sop.clear();
//		
//		exons_sop.add(new RegionBlock(114357290, 114357336));
//		exons_sop.add(new RegionBlock(114357090, 114357168));
//		exons_sop.add(new RegionBlock(114347890, 114347902));
		
		exons_sop.sort(new StartRegionBlockComparator());
		
		IntervalTree<RegionBlock> exons_sop_tree = new IntervalTree<RegionBlock>(exons_sop);
		exons_sop.clear();
		Iterator<Set<RegionBlock>> exon_sop_it = exons_sop_tree.groupIterator();
		while(exon_sop_it.hasNext()){
			Set<RegionBlock> overlap = exon_sop_it.next();
			RegionBlock[] overlapArray = new RegionBlock[overlap.size()];
			overlapArray = overlap.toArray(overlapArray);
			Arrays.sort(overlapArray,new StartRegionBlockComparator());
			int start = overlapArray[0].getStart();
			Arrays.sort(overlapArray,new StopRegionBlockComparator());
			int stop = overlapArray[overlapArray.length-1].getStop();
			
			RegionBlock newblock = new RegionBlock(start,stop);
			exons_sop.add(newblock);
		}
		
		exons_sop.sort(new StartRegionBlockComparator());
		
		for(int i =1 ; i < exons_sop.size(); i++){
			int start = exons_sop.get(i-1).getStop();
			int stop = exons_sop.get(i).getStart();
			if(start < stop){
				introns.add(new RegionBlock(start+1,stop));
			}
		}
		
		this.exons_fop = exons_fop;
		this.exons_sop = exons_sop;
		this.introns = introns;
		
		IntervalTree<RegionBlock> exons = new IntervalTree<RegionBlock>();
		exons.addAll(exons_sop);
		exons.addAll(exons_fop);
		
		incons = false;
		
//		if(readname.equals("4995990")){
//			System.out.println("Exons_FoP: " + exons_fop.toString());
//			System.out.println("Exons_SoP: "  + exons_sop.toString());
//			System.out.println("Introns: " + introns.toString());
//			System.out.println("Start: " + start);
//			System.out.println("Stop: " + stop);
//			System.out.println("o_start: " + o_start);
//			System.out.println("o_stop: " + o_stop);
//			
//		}
		
		for( RegionBlock intron : introns){
			if(intron.getStart() <= o_stop+1 && intron.getStop() > o_start){
				
//				System.out.println(intron);
				
				HashSet<RegionBlock> cont = new HashSet<RegionBlock>();
				cont = exons.getIntervalsIntersecting(intron.getStart(), intron.getStop()-1, cont);
				ArrayList<RegionBlock> bord = new ArrayList<RegionBlock>();
				bord =exons.getIntervalsIntersecting(intron.getStart()-1, intron.getStop(), bord);
				
//				if(readname.equals("9887855")){
//					System.out.println(intron);
//					System.out.println("cont " + cont.size());
//					System.out.println("bord " + cont.size());
//				}
//				
				
				if(cont.size() != 0){
					incons = true;
				}
				
				if(bord.size() != 4){
					incons = true;
				}
			}
		}
		
		alignment_blocks = new IntervalTree<RegionBlock>();
		
		Iterator<Set<RegionBlock>> exon_it = exons.groupIterator();
		while(exon_it.hasNext()){
			Set<RegionBlock> overlap = exon_it.next();
			RegionBlock[] overlapArray = new RegionBlock[overlap.size()];
			overlapArray = overlap.toArray(overlapArray);
			Arrays.sort(overlapArray,new StartRegionBlockComparator());
			int start = overlapArray[0].getStart();
			Arrays.sort(overlapArray,new StopRegionBlockComparator());
			int stop = overlapArray[overlapArray.length-1].getStop();
			
			RegionBlock newblock = new RegionBlock(start,stop);
			alignment_blocks.add(newblock);
		}
		

		this.introns.sort(new StartRegionBlockComparator());
		
		IntervalTree<RegionBlock> introns_tree = new IntervalTree<RegionBlock>(this.introns);
		this.introns.clear();
		Iterator<Set<RegionBlock>> introns_it = introns_tree.groupIterator();
		while(introns_it.hasNext()){
			Set<RegionBlock> overlap = introns_it.next();
			RegionBlock[] overlapArray = new RegionBlock[overlap.size()];
			overlapArray = overlap.toArray(overlapArray);
			Arrays.sort(overlapArray,new StartRegionBlockComparator());
			int start = overlapArray[0].getStart();
			Arrays.sort(overlapArray,new StopRegionBlockComparator());
			int stop = overlapArray[overlapArray.length-1].getStop();
			
			RegionBlock newblock = new RegionBlock(start,stop);
			this.introns.add(newblock);
		}
		
		
		
		int fop_start_diff = fop.getAlignmentStart() - fop.getUnclippedStart();
		int fop_stop_diff = fop.getUnclippedEnd() - fop.getAlignmentEnd();
		
		int sop_start_diff = sop.getAlignmentStart() - sop.getUnclippedStart();
		int sop_stop_diff = sop.getUnclippedEnd() - sop.getAlignmentEnd();
		
		clip = fop_start_diff + fop_stop_diff + sop_start_diff + sop_stop_diff;
		
//		if(incons){
//			System.out.println(readname);
//			System.out.println("FoP: " + exons_fop.toString());
//			System.out.println("SoP: " + exons_sop.toString());
//			
//			System.out.println("Intron: " + introns.toString());
//			
//			System.out.println("Ges: "+ this.start + ":" + this.stop);
//			System.out.println("Other: " + o_start + ":" + o_stop);
//			System.out.println(alignment_blocks.toTreeString());
//		}
	}
	
	public String toString(){
		String tab = "\t";
		String brk = "\n";
		char sep = ':';
		char sip = '|';
		char lin = '-';
		StringBuilder resultBuilder = new StringBuilder();
		resultBuilder.append(readname);
		resultBuilder.append(tab);
		resultBuilder.append(start);
		resultBuilder.append(sep);
		resultBuilder.append(stop);
		resultBuilder.append(tab);
		resultBuilder.append(incons);
		resultBuilder.append(tab);
		resultBuilder.append(strand);
		resultBuilder.append(tab);
		ArrayList<RegionBlock> blocks= new ArrayList<RegionBlock>(alignment_blocks);
		resultBuilder.append(blocks.toString());
		resultBuilder.append(tab);
		resultBuilder.append(exons_fop.toString());
		resultBuilder.append(tab);
		resultBuilder.append(exons_sop.toString());
		
		return resultBuilder.toString();
	}

	@Override
	public int getStart() {
		return start;
	}

	@Override
	public int getStop() {
		return stop;
	}
	
	public String getReadName(){
		return readname;
	}
	
	public IntervalTree<RegionBlock> getAlignmentBlocks(){
		return alignment_blocks;
	}
	
	public ArrayList<RegionBlock> getAlignmentBlocksFoP(){
		return exons_fop;
	}
	
	public ArrayList<RegionBlock> getAlignmentBlocksSoP(){
		return exons_sop;
	}
	
	public ArrayList<RegionBlock> getIntronBlocks(){
		return introns;
	}
	
	public boolean isConsistent(){
		return !incons;
	}
	
	public String getChromosome(){
		return chr;
	}
	
	public char getStrand(){
		return strand;
	}
	
	public int getMM(){
		return mm;
	}
	
	public int getClip(){
		return clip;
	}

	public int getSplit(){
		return introns.size();
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
	
}
