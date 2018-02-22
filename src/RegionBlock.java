import AugmentedTree.Interval;

public class RegionBlock implements Interval{

	private int start;
	private int stop;
	
	public RegionBlock(int start,int stop){
		this.start = start;
		this.stop = stop;
	}
	
	@Override
	public String toString(){
		return (start + ":" +stop);
	}
	
	@Override
	public int getStart() {
		return start;
	}

	@Override
	public int getStop() {
		return stop;
	}

}
