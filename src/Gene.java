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
	
	public Gene(int x1, int x2, Annotation annotation){
		super(x1,x2,annotation);
		transcripts = new HashMap <String,Transcript>();
		transReadIds = new  HashSet<String>();
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
	
	public ArrayList<ExonSkip> getExonSkips() {
//		long startTime = System.currentTimeMillis();
		IntervalTree<Region> introns = new IntervalTree<Region>();
		
		for( String k : transcripts.keySet() ){
			Transcript curTran = transcripts.get(k);
			if(curTran.setIntrons()){
				introns.addAll(curTran.getIntrons());
			}
		}
		//introns holds all introns
		
//		long stopTime = System.currentTimeMillis();
//		System.out.println("All Introns: " + (stopTime-startTime));
		
//		startTime = System.currentTimeMillis();
		IntervalTree<Region> mergeIntrons = new IntervalTree<Region>();
		Iterator<Region> iterator = introns.iterator();
		int lastStart = -1;
		int lastStop = -1;
		int curStart = -2;
		int curStop = -2;
		while(iterator.hasNext()){
			Region curIntron = iterator.next();
			curStart = curIntron.getStart();
			curStop = curIntron.getStop();
			if(lastStart != curStart || lastStop != curStop){
				ArrayList<Region> sameIntron = new ArrayList<Region>();
				sameIntron = introns.getIntervalsEqual(curStart, curStop, sameIntron);
				Region mergeIntron = curIntron;
				for( Region r : sameIntron){
					mergeIntron = mergeIntron.merge(r);
				}
				mergeIntrons.add(mergeIntron);
				lastStart = curStart;
				lastStop = curStop;
			}
		}
		//mergeIntrons holds all distinct introns
		
//		stopTime = System.currentTimeMillis();
//		System.out.println("merge Introns: " + (stopTime-startTime));
		
//		System.out.println(mergeIntrons.toTreeString());
		
		ArrayList<ExonSkip> skips = new ArrayList<ExonSkip>();
		
		Iterator<Region> mergeIterator = mergeIntrons.iterator();
//		startTime = System.currentTimeMillis();
		while(mergeIterator.hasNext()){
			Region sv = mergeIterator.next();
			int svstart = sv.getStart();
			int svstop = sv.getStop();
			
			HashSet<String> sv_transcript_ids = sv.getAnnotation().getSuperIds();
			
			HashSet<String> wt_start_transcript_ids = new HashSet<String>();
			
			HashSet<Region> wt_candidatesStart = new HashSet<Region>();
			wt_candidatesStart = mergeIntrons.getIntervalsBeginAt(svstart, wt_candidatesStart);
			
			for(Region wt_candidate : wt_candidatesStart){
				wt_start_transcript_ids.addAll(wt_candidate.getAnnotation().getSuperIds());
			}

			HashSet<String> wt_stop_transcript_ids = new HashSet<String>();
			
			HashSet<Region> wt_candidatesStop = new HashSet<Region>();
			wt_candidatesStop = mergeIntrons.getIntervalsEndAt(svstop, wt_candidatesStop);
			
			for(Region wt_candidate : wt_candidatesStop){
				wt_stop_transcript_ids.addAll(wt_candidate.getAnnotation().getSuperIds());
			}
			
			HashSet<String> wt_transcript_ids = new HashSet<String>();
			
			wt_transcript_ids.addAll(wt_stop_transcript_ids);
			
			wt_transcript_ids.retainAll(wt_start_transcript_ids);
			wt_transcript_ids.removeAll(sv_transcript_ids);
			
			if(wt_transcript_ids.size() > 0){
			
			int minBase = svstop -svstart;
			int maxBase = 0;
			
			int minEx = svstop -svstart;
			int maxEx = 0;
			
			HashSet<String> wt_prot_ids = new HashSet<String>();
			
			HashSet<Region> wt_skips = new HashSet<Region>();
			
			HashSet<Region> wt_exons = new HashSet<Region>();
			
//			System.out.println("Transcript IDs " + wt_transcript_ids);
			
			for( String id: wt_transcript_ids){
				Transcript wt_transcript = transcripts.get(id);
				HashSet<Region> wt_id_skips = new HashSet<Region>();
				wt_id_skips = wt_transcript.getIntrons().getIntervalsSpannedBy(svstart, svstop, wt_id_skips);
				
				HashSet<Region> wt_id_exons = new HashSet<Region>();
				wt_id_exons = wt_transcript.getRegionsTree().getIntervalsSpannedBy(svstart, svstop, wt_id_exons);
				int ex = wt_id_skips.size()-1;
				maxEx = Math.max(maxEx, ex);
				minEx = Math.min(minEx, ex);
				
				int base = svstop-svstart;
				for(Region wt_exon: wt_id_exons){	
					wt_exons.add(wt_exon);
					
				}
				for(Region wt_skip: wt_id_skips){
					HashSet<String> wt_skip_prot = new HashSet<String>();
					base -= wt_skip.length();
	
					wt_skips.add(wt_skip);
					
					wt_skip_prot = wt_skip.getAnnotation().getIds();
					
					for(String skip_prot: wt_skip_prot){
						wt_prot_ids.add(skip_prot);
					}
				}
				
				
				maxBase = Math.max(maxBase, base);
				minBase = Math.min(minBase, base);
			}
			

			ExonSkip skip = new ExonSkip(sv.getAnnotation().getIds(),wt_prot_ids,wt_skips,wt_exons,minEx,maxEx,minBase,maxBase, svstart, svstop);
			skips.add(skip);
			}
		}
//		stopTime = System.currentTimeMillis();
//		System.out.println("Skip Identifizierung: " + (stopTime-startTime));
		
		
		return skips;
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

