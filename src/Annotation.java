import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

public class Annotation {
	
	private HashSet<String> id;
	private String chr;
	private char strand;
	private HashSet<String> super_id;
	private String super_super_id;
	private String name;
	private String type;
	private String gene_name;
	
	public Annotation(String id, String name, String chr, char strand, String type){
		this.id = new HashSet<String>();
		this.id.add(id);
		this.name = name;
		this.chr = chr;
		this.strand = strand;
		this.super_id = new HashSet<String>();
		this.super_id.add("");
		this.super_super_id = "";
		this.type = type;
		this.gene_name = name;
	}
	
	public Annotation(String id, String name, String chr, char strand, String super_id, String type, String gene_name){
		this.id = new HashSet<String>();
		this.id.add(id);
		this.name = name;
		this.chr = chr;
		this.strand = strand;
		this.super_id = new HashSet<String>();
		this.super_id.add(super_id);
		this.super_super_id = "";
		this.type = type;
		this.gene_name = gene_name;
	}
	
	public Annotation(String id, String chr, char strand, String super_id, String super_super_id, String type, String gene_name){
		this.id = new HashSet<String>();
		this.id.add(id);
		this.name = "";
		this.chr = chr;
		this.strand = strand;
		this.super_id = new HashSet<String>();
		this.super_id.add(super_id);
		this.super_super_id = super_super_id;
		this.type = type;
		this.gene_name = gene_name;
	}
	
	public Annotation(Annotation a, Annotation b) {
		id = new HashSet<String>();
		id.addAll(new ArrayList<String>(a.getIds()));
		id.addAll(new ArrayList<String>(b.getIds()));
		super_id = new HashSet<String>();
		super_id.addAll(new ArrayList<String>(a.getSuperIds()));
		super_id.addAll(new ArrayList<String>(b.getSuperIds()));
		if(a.getName().equals("")){
			this.name = b.getName();
		}
		else{
			this.name = a.getName() + "|" + b.getName();
		}
		this.chr = a.getChromosome();
		this.strand = a.getStrand();
		this.super_super_id = a.getSuperSuperId();
		this.gene_name = a.getGeneName();
		String typea = a.getType();
		String typeb = b.getType();
		if(typea.equals(typeb)){
			this.type = typea;
		}
		else{
			this.type = typea + "|" + typeb;
		}
	}

	public String getId(){
		Iterator<String> idIterator = id.iterator();
		return idIterator.next();
	}
	
	public HashSet<String> getIds(){
		return id;
	}
	
	public String getName(){
		return name;
	}
	
	public String getChromosome(){
		return chr;
	}
	
	public char getStrand(){
		return strand;
	}
	
	public String getSuperId(){
		Iterator<String> sidIterator = super_id.iterator();
		return sidIterator.next();
	}
	
	public HashSet<String> getSuperIds(){
		return super_id;
	}
	
	public String getSuperSuperId(){
		return super_super_id;
	}
	
	public String getType(){
		return type;
	}
	
	public String getGeneName(){
		return gene_name;
	}
	
	public boolean isSup(Annotation a){
		if(a != null){
			return a.getSuperId().equals(id);
		}
		else{
			return false;
		}
	}
	
	public boolean isSub(Annotation a){
		if(a != null){
			boolean b = true;
			for( String aid : a.getIds()){
				b &= super_id.contains(aid);
			}
			return b;
		}
		else{
			return false;
		}
	}
	
	public String toString(){
		if(name.equals("")){
			return id.toString() + "; " + super_id.toString();
		}
		else{
			return id.toString() + "; " + name.toString() + "; " + super_id.toString();
		}
	}
	
	public boolean equals(Annotation a){
		return a.getIds().equals(id) && a.getSuperIds().equals(super_id) && a.getSuperSuperId().equals(super_super_id);
	}
	
	
	@Override
	public int hashCode(){
		return id.hashCode();
	}
	
}
