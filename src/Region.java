import java.util.Comparator;

import AugmentedTree.Interval;

public class Region implements Interval {
    private int start;
    private int stop;
    
    private Annotation annotation;

//    public Region(int start, int stop) 
//    {
//        this.start = start;
//        this.stop = stop;
//        this.annotation = null;
//    }
    
    public Region(int start, int stop, Annotation annotation) 
    {
        this.start = start;
        this.stop = stop;
        this.annotation = annotation;
    }

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }
    
    public int length(){
    	return stop - start;
    }
    
    public Annotation getAnnotation(){
    	return annotation;
    }
    
    public void setStart(int start){
    	this.start = start;
    }
    
    public void setStop(int stop){
    	this.stop = stop;
    }
    
    public Region merge(Region r){
    	int nstart = Math.min(start,r.getStart());
    	int nstop = Math.max(stop, r.getStop());
    	Annotation nannotation = new Annotation(annotation,r.getAnnotation());
		return new Region(nstart,nstop,nannotation);
    }
    
    @Override
    public String toString(){
		return start + ":" + stop;
    	
    }
    
    public boolean equals(Object o){
    	if(this.getClass() == o.getClass()){
    		return this.equals((Region) o);
    	}
		return false;
    }
    
    public boolean equals(Interval r){
    	return r.getStart() == this.getStart() && r.getStop() == this.getStop();
    }
    
    public int hashCode(){
    	String hashString = this.getStart() + "" + this.getStop();
    	return hashString.hashCode();
    }
}
