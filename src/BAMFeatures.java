import java.awt.List;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import AugmentedTree.IntervalTree;
import net.sf.samtools.*;

public class BAMFeatures {

	public static void main(String[] args) {
//		long startTime = System.currentTimeMillis();		
		String gtfPath ="";
		String bamPath ="";
		String outputPath ="";
		String frstrand = "";
		for(int i =0; i < args.length-1; i++){
			if(args[i].equals("-gtf")){
				gtfPath = args[i+1];
				i++;
			}
			else if(args[i].equals("-bam")){
				bamPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-o")){
				outputPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-frstrand")){
				frstrand = args[i+1];
				i++; 
			}
		}
		
		if(gtfPath.equals("") || bamPath.equals("") || outputPath.equals("")){
			System.out.println("Usage Info:\n-gtf <filepath for GTF>\n-bam <filepath for BAM>\n-o <filepath for ouput>\n[-frstrand true/false <strandness>]");
		}
		else{
			int strandness = 0;
			if(frstrand.equals("true")){
				strandness = 1;
			}
			else if(frstrand.equals("false")){
				strandness = 2;	
			}
			
			Path gtfFilePath = Paths.get(gtfPath);
			Path bamFilePath = Paths.get(bamPath);
			Path outputFilePath = Paths.get(outputPath);
			
			File gtfFile = gtfFilePath.toFile();
			File bamFile = bamFilePath.toFile();
			
			BAMFeatures bamf = new BAMFeatures(gtfFile, bamFile, strandness);
		
			String output = bamf.getOutput();
			
			ArrayList<String> outputAsList = new ArrayList<String>();
			outputAsList.add(output);
			try {
				Files.write(outputFilePath, outputAsList, Charset.forName("UTF-8"));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
//		long stopTime = System.currentTimeMillis();
//		System.out.println("Input:" + (stopTime-startTime));	
	}

	private HashMap<Integer,Gene> geneSet;
	private HashMap<String,HashMap<Character,IntervalTree<Gene>>> geneTree;
	private HashMap<String,HashMap<Character,IntervalTree<Read>>> readTree;
	
	private SAMFileReader bam_reader;
	
	public BAMFeatures (File gtfFile, File bamFile, int strandness) {
//		long startTime = System.currentTimeMillis();
		
		bam_reader = new SAMFileReader(bamFile);
		bam_reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		
		geneSet = new HashMap<Integer,Gene>();
		
	    try {
	        BufferedReader br = new BufferedReader (new FileReader(gtfFile));
	        String line;
	        Gene curGene = null;
	        Annotation curGAnno = null;
	        Transcript curTrans = null;
	        Annotation curTAnno = null;
	        Annotation curCAnno = null;
	        ArrayList<Transcript> unknownGene = new ArrayList<Transcript>();
	        ArrayList<Region> unknownTranscript = new ArrayList<Region>();
	        HashMap<String,Annotation> mapAnno = new HashMap<String,Annotation>(); 
	        while ((line = br.readLine()) != null){
	        	String[] lineSplit = line.split("\t");
	        	if(lineSplit.length >= 8){
	        	String[] attrSplit = lineSplit[8].split(";");
	        	HashMap<String,String> attr;
	        	if(lineSplit[2].equals("gene")){
	        		attr = getAttributes(attrSplit);
	        		String gene_id = attr.get("gene_id");
	        		String gene_name = attr.get("gene_name");
	        		String gene_biotype = attr.get("gene_biotype");
	        		String type = lineSplit[1];
	        		char strand = lineSplit[6].charAt(0);
	        		String chr = lineSplit[0];
	        		int start = Integer.parseInt(lineSplit[3]);
	        		int stop = Integer.parseInt(lineSplit[4]);
	        		
	        		Annotation geneAnno = new Annotation(gene_id,gene_name,chr,strand,type);
	        		Gene gene = new Gene(start, stop, geneAnno, gene_biotype);
	        		geneSet.put(gene.hashCode(),gene);
	        		curGene = gene;
	        		curGAnno = geneAnno;
	        	}
	        	if(lineSplit[2].equals("CDS")){
	        		attr = getAttributes(attrSplit);
	        		String super_super_id = attr.get("gene_id");
	        		String gene_name = attr.get("gene_name");
	        		String super_id = attr.get("transcript_id");
	        		String id = attr.get("protein_id");
	        		String type = lineSplit[1];
	        		char strand = lineSplit[6].charAt(0);
	        		String chr = lineSplit[0];
	        		int start = Integer.parseInt(lineSplit[3]);
	        		int stop = Integer.parseInt(lineSplit[4]);
	        		Region cds;
	        		Annotation cdsAnno;
	        		if(curCAnno != null && curCAnno.getId().equals(id) && curCAnno.getSuperId().equals(super_id) && curCAnno.getSuperSuperId().equals(super_super_id)){
	        			cdsAnno = curCAnno;
	        		}
	        		else{
	        			cdsAnno = mapAnno.get(id+super_id+super_super_id);
	        			if(!(cdsAnno != null)){
	        				cdsAnno = new Annotation(id,chr,strand,super_id,super_super_id,type,gene_name);
	        				mapAnno.put(id+super_id+super_super_id,cdsAnno);
	        			}
	        			curCAnno = cdsAnno;
	        		}
	        		cds = new Region(start,stop,cdsAnno);
	        		
	        		
	        		if(cdsAnno.isSub(curTAnno)){
	        			curTrans.add(cds);
	        		}
	        		else {
		        		String transcript_name = attr.get("transcript_name");	        			
		        		Annotation transAnno = new Annotation(super_id,transcript_name,chr,strand,super_super_id,type,gene_name);
		        		Transcript trans = new Transcript(start,stop,transAnno);

		        		if(transAnno.isSub(curGAnno)){
		        			curGene.add(trans);
		        		}
		        		else{
		        			System.out.println(transAnno.getSuperId() + " is missing");
//			        		Annotation geneAnno = new Annotation(super_super_id,gene_name,chr,strand,type);
//			        		Gene gene = new Gene(start, stop, geneAnno);
//			        		geneSet.put(gene.hashCode(),gene);
//			        		curGene = gene;
//			        		curGAnno = geneAnno;
//			        		curGene.add(trans);
		        		}
		        		
		        		curTAnno = transAnno;
		        		curTrans = trans;
		        		curTrans.add(cds);
	        		}
	        		
	        	}
	        	if(lineSplit[2].equals("exon")){
	        		attr = getAttributes(attrSplit);
	        		String super_super_id = attr.get("gene_id");
	        		String gene_name = attr.get("gene_name");
	        		String super_id = attr.get("transcript_id");
	        		String id = attr.get("protein_id");
	        		String type = lineSplit[1];
	        		char strand = lineSplit[6].charAt(0);
	        		String chr = lineSplit[0];
	        		int start = Integer.parseInt(lineSplit[3]);
	        		int stop = Integer.parseInt(lineSplit[4]);
	        		RegionBlock exon;
	        		Annotation exonAnno;
	        		if(curCAnno != null && curCAnno.getId().equals(id) && curCAnno.getSuperId().equals(super_id) && curCAnno.getSuperSuperId().equals(super_super_id)){
	        			exonAnno = curCAnno;
	        		}
	        		else{
	        			exonAnno = mapAnno.get(id+super_id+super_super_id);
	        			if(!(exonAnno != null)){
	        				exonAnno = new Annotation(id,chr,strand,super_id,super_super_id,type,gene_name);
	        				mapAnno.put(id+super_id+super_super_id,exonAnno);
	        			}
	        		}
	        		exon = new RegionBlock(start,stop);
	        		
	        		
	        		if(exonAnno.isSub(curTAnno)){
	        			curTrans.add(exon);
	        		}
	        		else {
		        		String transcript_name = attr.get("transcript_name");	        			
		        		Annotation transAnno = new Annotation(super_id,transcript_name,chr,strand,super_super_id,type,gene_name);
		        		Transcript trans = new Transcript(start,stop,transAnno);

		        		if(transAnno.isSub(curGAnno)){
		        			curGene.add(trans);
		        		}
		        		else{
		        			System.out.println(transAnno.getSuperId() + " is missing");
//			        		Annotation geneAnno = new Annotation(super_super_id,gene_name,chr,strand,type);
//			        		Gene gene = new Gene(start, stop, geneAnno);
//			        		geneSet.put(gene.hashCode(),gene);
//			        		curGene = gene;
//			        		curGAnno = geneAnno;
//			        		curGene.add(trans);
		        		}
		        		
		        		curTAnno = transAnno;
		        		curTrans = trans;
		        		curTrans.add(exon);
	        		}
	        	}
	        	
	        	else{}
	        	}
	        }
	        br.close();	
	    	
	    } catch (Exception e) {
			e.printStackTrace();
		}
	    
	    geneTree = new HashMap<String,HashMap<Character,IntervalTree<Gene>>>();
	    
	    for(int key: geneSet.keySet()){
	    	Gene curGene = geneSet.get(key);
	    	Annotation curAnno = curGene.getAnnotation();
	    	char str = curAnno.getStrand();
	    	String chr = curAnno.getChromosome();
	    	if(geneTree.containsKey(chr)){
	    		HashMap<Character,IntervalTree<Gene>> chrTree = geneTree.get(chr);
	    		if(chrTree.containsKey(str)){
	    			IntervalTree<Gene> strTree = chrTree.get(str);
	    			strTree.add(curGene);
	    		}
	    		else{
	    			IntervalTree<Gene> strTree = new IntervalTree<Gene>();
	    			strTree.add(curGene);
	    			chrTree.put(str, strTree);
	    		}
	    	}
	    	else{
    			IntervalTree<Gene> strTree = new IntervalTree<Gene>();
    			strTree.add(curGene);
    			HashMap<Character,IntervalTree<Gene>> chrTree = new HashMap<Character,IntervalTree<Gene>>();
    			chrTree.put(str, strTree);
    			geneTree.put(chr, chrTree);
	    	}
	    }
	    
	    
//	    System.out.println(geneTree.keySet().toString());
//	    
//	    for(String chr : geneTree.keySet()){
//	    	System.out.println(geneTree.get(chr).keySet().toString());
//	    	for(char str: geneTree.get(chr).keySet()){
//	    		if(readTree.containsKey(chr)){
//		    		if(readTree.get(chr).containsKey(str)){
//		    			System.out.println(geneTree.get(chr).get(str).size() + " : " + readTree.get(chr).get(str).size());
//		    		}
//	    		}
//	    	}
//	    }
	    
//		long stopTime = System.currentTimeMillis();
//		System.out.println("Input:" + (stopTime-startTime));
	}
	
	public static HashMap<String,String> getAttributes (String[] attrs){
		HashMap <String,String> attrMap = new HashMap <String,String>();
		for(int i = 0; i < attrs.length; i++){
			String curattr = attrs[i];
			curattr = curattr.trim();
			String[] attr = curattr.split("\\s+");
			attr[0] = attr[0].trim();
			attr[1] = attr[1].trim();
			attr[1] = attr[1].replaceAll("\"", "");
			attrMap.put(attr[0],attr[1]);
		}
		return attrMap;
	}
	
	public ArrayList<Gene> getGenes(){
		ArrayList<Gene> result = new ArrayList<Gene>();
		for (Integer key: geneSet.keySet()) {
		   result.add(geneSet.get(key));
		}
		return result;
	}
	
	public String getOutput(){
		String tab = "\t";
		String brk = "\n";
		char dop = ':';
		char sip = '|';
		char lin = '-';
		char com = ',';
		
		String incStr = "split-inconsistent:true";
		String mmStr = "mm:";
		String clipStr = "clipping:";
		String splStr = "nsplit:";
		String couStr = "gcount:";
		String distStr = "gdist:";
		String ansStr = "antisense:";
		String pcrStr = "pcrindex:";
		String merStr = "MERGED";
		String intStr = "INTRON";

		StringBuilder resultBuilder = new StringBuilder("");

		SAMRecordIterator bam_it = bam_reader.iterator();

		HashMap<String,SAMRecord> store = new HashMap<String,SAMRecord>();
		
		readTree = new HashMap<String,HashMap<Character,IntervalTree<Read>>>();
		
		while(bam_it.hasNext()){

			SAMRecord samr = bam_it.next();
			
			boolean ignore = samr.getNotPrimaryAlignmentFlag() || samr.getReadUnmappedFlag() || samr.getMateUnmappedFlag() || (samr.getMateNegativeStrandFlag() == samr.getReadNegativeStrandFlag());

			boolean inPair = samr.getFirstOfPairFlag()||samr.getSecondOfPairFlag();

			if(inPair && !ignore){
				String readname = samr.getReadName();
//				&& readname.equals("87178")	
//				Integer x= count.put(readname,1);
//				if(x != null){
//					count.put(readname, x+1);
//				}
				if(store.containsKey(readname)){
//					System.out.println(readname);
					Read curRead = new Read(samr, store.get(readname));
					
					resultBuilder.append(readname);
					
					if(curRead.isConsistent()){
						
						resultBuilder.append(tab);
						resultBuilder.append(mmStr);
						int mm = curRead.getMM();
						resultBuilder.append(mm);
						
						resultBuilder.append(tab);
						resultBuilder.append(clipStr);
						int clip = curRead.getClip();
						resultBuilder.append(clip);
						
						resultBuilder.append(tab);
						resultBuilder.append(splStr);
						int spl = curRead.getSplit();
						resultBuilder.append(spl);
						
					    char str = curRead.getStrand();
					    String chr = curRead.getChromosome();
					    
					    int start = curRead.getStart();
					    int stop = curRead.getStop();
					    
					    boolean transcriptomic = false;
					   
					    HashSet<Gene> cgenes = new HashSet<Gene>();
					    cgenes = geneTree.get(chr).get(str).getIntervalsSpanning(start, stop, cgenes);
					    HashSet<Gene> igenes = new HashSet<Gene>();
					    igenes = geneTree.get(chr).get(str).getIntervalsSpannedBy(start, stop, igenes);
					    
					    
					    if(cgenes.size() > 0){
					    	for(Gene g : cgenes){
//								if(curRead.getReadName().equals("7075")){
//									System.out.println(g.toString());
//								}
					    		
					    		transcriptomic |= g.inTranscript(curRead);
					    	}
					    }
					    else{
					    	resultBuilder.append(tab);
				    		resultBuilder.append(couStr);
				    		resultBuilder.append(0);
				    		
					    	if(igenes.size() > 0){
					    		resultBuilder.append(tab);
					    		resultBuilder.append("skip");
					    	}
					    	else{
					    		
					    		resultBuilder.append(tab);
					    		resultBuilder.append(distStr);
					    		
					    		ArrayList<Gene> rneigh = new ArrayList<Gene>();
					    		rneigh = geneTree.get(chr).get(str).getIntervalsRightNeighbor(start, stop, rneigh);
					    		long rdist = -1;
					    		if(rneigh.size() > 0){
					    			rdist = rneigh.get(0).getStart() - stop -1;
					    		}
					    		
					    		ArrayList<Gene> lneigh = new ArrayList<Gene>();
					    		lneigh = geneTree.get(chr).get(str).getIntervalsLeftNeighbor(start, stop, lneigh);
					    		long ldist = -1;
					    		if(lneigh.size() > 0){
					    			ldist = start - lneigh.get(0).getStop() -1;
					    		}
					    		
					    		if(ldist < 0){
					    			if(rdist < 0){
					    				resultBuilder.append(0);
					    			}
					    			else{
					    				resultBuilder.append(rdist);
					    			}
					    		}
					    		else{
					    			if(rdist < 0){
					    				resultBuilder.append(ldist);
					    			}
					    			else{
					    				resultBuilder.append(Math.min(ldist, rdist));
					    			}
					    		}
					    		
					    	}
					    }
					    
					    if(transcriptomic){
						    if(readTree.containsKey(chr)){
						    	HashMap<Character,IntervalTree<Read>> chrTree = readTree.get(chr);
						    	if(chrTree.containsKey(str)){
						    		IntervalTree<Read> strTree = chrTree.get(str);
						    		strTree.add(curRead);
						    	}
						    	else{
						    		IntervalTree<Read> strTree = new IntervalTree<Read>();
						    		strTree.add(curRead);
						    		chrTree.put(str, strTree);
						    	}
						   	}
					    	else{
					    		IntervalTree<Read> strTree = new IntervalTree<Read>();
					    		strTree.add(curRead);
					    		HashMap<Character,IntervalTree<Read>> chrTree = new HashMap<Character,IntervalTree<Read>>();
					    		chrTree.put(str, strTree);
					   			readTree.put(chr, chrTree);
					    	}
					    }
					}
					else{
						resultBuilder.append(tab);
						resultBuilder.append(incStr);
					}
					store.remove(readname);
					resultBuilder.append(brk);
				}
				else{
					store.put(readname, samr);
				}
			}
		}
		bam_it.close();
		
		return resultBuilder.toString();
		
	}
	
	class StartRegionComparator implements Comparator<Region>
	{
	    public int compare(Region x1, Region x2)
	    {
	        return x1.getStart() - x2.getStart();
	    }
	}
	
	class StartReadComparator implements Comparator<Read>
	{
	    public int compare(Read x1, Read x2)
	    {
	        return x1.getStart() - x2.getStart();
	    }
	}

}
