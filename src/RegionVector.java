import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import AugmentedTree.IntervalTree;

public class RegionVector extends Region{
	
//Start inclusive, Stop exclusive


	private IntervalTree<Region> regions;
	
	public RegionVector(int x1, int x2, Annotation annotation){
		super(x1,x2,annotation);
		regions = new IntervalTree<Region>();
	}

	public RegionVector(int x1, int x2, Annotation annotation, Annotation subannotation){
		super(x1,x2,annotation);
		regions = new IntervalTree<Region>();
		regions.add(new Region(x1,x2,subannotation));
	}
	
	public RegionVector(IntervalTree<Region> region, Annotation annotation){
		super(region.getStart(),region.getStop(),annotation);
		this.regions = region;
	}
	
	public RegionVector(Vector<Region> region, Annotation annotation){
		super(region.get(0).getStart(),region.get(region.size()-1).getStop(),annotation);
		regions = new IntervalTree<Region>();
		regions.addAll(region);
	}
	
	public RegionVector(Region[] region, Annotation annotation){
		super(region[0].getStart(),region[region.length-1].getStop(), annotation);
		regions = new IntervalTree<Region>();
		for(int i = 0; i < region.length; i++){
			regions.add(region[i]);
		}
	}
	
	public Region[] getRegionsArray(){
		return regions.toArray(new Region[0]);
	}
	
	public IntervalTree<Region> getRegionsTree(){
		return regions;
	}
	
	public ArrayList<Region> getRegions(){
		ArrayList <Region> r = new ArrayList<Region>();
		Region[] regionArray = regions.toArray(new Region[0]);
		for(int i = 0; i < regionArray.length; i++){
			r.add(regionArray[i]);
		}
		return r;
	}
	
	public RegionVector merge(){
		
		Vector<Region> resultV = new Vector<Region>();

		Iterator<Set<Region>> iterator = regions.groupIterator();
		while(iterator.hasNext()){
			Collection<Region> overlap = (Collection<Region>) iterator.next();
			Vector<Region> overlapVector = new Vector<Region>(overlap);
			int start = overlapVector.get(0).getStart();
			overlapVector.sort(new StopRegionComparator());
			int stop = overlapVector.lastElement().getStop();
			resultV.add(new Region(start,stop,overlapVector.firstElement().getAnnotation()));
		}
		
		return new RegionVector(resultV,this.getAnnotation());
   
	}
	
	public void merge(RegionVector rv){
		
		IntervalTree<Region> regions = this.regions.clone(); 

		regions.addAll(rv.getRegionsTree());
	}
	
	public RegionVector subtract(RegionVector rv){
		IntervalTree<Region> regions = this.regions.clone();
		rv = rv.merge();
		Region[] rvarray = rv.getRegionsArray();
		for(int i = 0; i < rvarray.length; i++){
			Vector<Region> rvvector = new Vector<Region>();
			int rvstart = rvarray[i].getStart();
			int rvstop = rvarray[i].getStop();
			rvvector = regions.getIntervalsIntersecting(rvstart, rvstop, rvvector );
			
			for(int j = 0; j < rvvector.size(); j++){
				Region curreg = rvvector.get(j);
				int start = curreg.getStart();
				int stop = curreg.getStop();
				
				if(stop <= rvstop){
					if(start < rvstart){
						regions.add(new Region(start,rvstart,curreg.getAnnotation()));
						regions.remove(curreg);
					}
					else{
						regions.remove(curreg);
					}
				}
				else{
					if(start < rvstart){
						regions.add(new Region(start,rvstart,curreg.getAnnotation()));
						regions.add(new Region(rvstop,stop,curreg.getAnnotation()));
						regions.remove(curreg);
					}
					else{
						regions.add(new Region(rvstop,stop,curreg.getAnnotation()));
						regions.remove(curreg);
					}
				}
			}
		}
		return new RegionVector(regions,this.getAnnotation());
		
	}
	
	public RegionVector invert(){
		IntervalTree<Region> result = new IntervalTree<Region>();
		ArrayList<Region> regionsArray = new ArrayList<Region>(regions);
		regionsArray.sort(new StartRegionComparator());
		HashSet<Region> resultSet = new HashSet<Region>();
		int start;
		int stop;
		for(int i = 1; i < regionsArray.size(); i++){
			Region last = regionsArray.get(i-1);
			Region cur = regionsArray.get(i);
			start = last.getStop();
			stop = cur.getStart();
			Annotation anno;
			if(last.getAnnotation().equals(cur.getAnnotation())){
				anno = cur.getAnnotation();
			}
			else{
				anno = new Annotation(last.getAnnotation(),cur.getAnnotation());
			}
			resultSet.add(new Region(start+1,stop,anno));
		}
		result.addAll(resultSet);
		return new RegionVector(result,this.getAnnotation());
	}
	
	public int coveredLength(){
		int l = 0;
		
		Iterator<Set<Region>> iterator = regions.groupIterator();
		while(iterator.hasNext()){
			Collection<Region> overlap = (Collection<Region>) iterator.next();
			Vector<Region> overlapVector = new Vector<Region>(overlap);
			int start = overlapVector.get(0).getStart();
			overlapVector.sort(new StopRegionComparator());
			int stop = overlapVector.lastElement().getStop();
			l += stop-start;
		}
		
		return l;
	}
	
	public ArrayList<Region> getCoveredRegion(RegionVector rv){
		ArrayList<Region> results = new ArrayList<Region>();
		rv = rv.merge();
		Region[] rvarray = rv.getRegionsArray();
		for(int i = 0; i < rvarray.length; i++){
			ArrayList<Region> rvlist = new ArrayList<Region>();
			int rvstart = rvarray[i].getStart();
			int rvstop = rvarray[i].getStop();
			rvlist = regions.getIntervalsIntersecting(rvstart+1, rvstop-1, rvlist );
			results.addAll(rvlist);
		}
		
		if(results.size() > 0){
			return results;
		}else{
			return null;
		}

	}
	
	public ArrayList<Region> getCoveredRegion(Region r){
		ArrayList<Region> results = new ArrayList<Region>();
		int rvstart = r.getStart();
		int rvstop = r.getStop();
		results = regions.getIntervalsIntersecting(rvstart+1, rvstop-1, results );
		
		if(results.size() > 0){
			return results;
		}else{
			return null;
		}

	}
	
	public RegionVector getRegion(Region r){
		IntervalTree<Region> results = new IntervalTree<Region>();
		Vector<Region> rvvector = new Vector<Region>();
		int rvstart = r.getStart();
		int rvstop = r.getStop();
		rvvector = regions.getIntervalsIntersecting(rvstart+1, rvstop-1, rvvector );
		results.addAll(rvvector);
		if(results.size() > 0){
			return new RegionVector(results,this.getAnnotation());
		}else{
			return null;
		}
	}

	
	public boolean isSub(RegionVector rv){
		return regions.containsAll(rv.getRegions());
	}
	
	public boolean add(Region region){
		if(region.getStart() < this.getStart()){
			this.setStart(region.getStart());
		}
		if(region.getStop() > this.getStop()){
			this.setStop(region.getStop());
		}
		return regions.add(region);
	}
	
	public boolean remove(Region region){
		return regions.remove(region);
	}	
	
	public int hashCode(){
		return this.getAnnotation().hashCode();
	}
	
	public String toString(){
		return super.toString() + "\n" + regions.toTreeString();
	}
	
	class StartRegionComparator implements Comparator<Region>
	{
	    public int compare(Region x1, Region x2)
	    {
	        return x1.getStart() - x2.getStart();
	    }
	}
	
	class StopRegionComparator implements Comparator<Region>
	{
	    public int compare(Region x1, Region x2)
	    {
	        return x1.getStop() - x2.getStop();
	    }
	}
	

}
