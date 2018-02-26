import java.awt.List;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
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
import java.util.TreeMap;
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
			
			File gtfFile = gtfFilePath.toFile();
			File bamFile = bamFilePath.toFile();
			
			BAMFeatures bamf = new BAMFeatures(gtfFile, bamFile, strandness);
			
			try {
				bamf.getOutput(outputPath);
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
	
	private int strandness;
	
	private SAMFileReader bam_reader;
	
	public BAMFeatures (File gtfFile, File bamFile, int strandness) {
//		long startTime = System.currentTimeMillis();
		this.strandness = strandness;
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
//	        	if(lineSplit[2].equals("CDS")){
//	        		attr = getAttributes(attrSplit);
//	        		String super_super_id = attr.get("gene_id");
//	        		String gene_name = attr.get("gene_name");
//	        		String super_id = attr.get("transcript_id");
//	        		String id = attr.get("protein_id");
//	        		String type = lineSplit[1];
//	        		char strand = lineSplit[6].charAt(0);
//	        		String chr = lineSplit[0];
//	        		int start = Integer.parseInt(lineSplit[3]);
//	        		int stop = Integer.parseInt(lineSplit[4]);
//	        		Region cds;
//	        		Annotation cdsAnno;
//	        		if(curCAnno != null && curCAnno.getId().equals(id) && curCAnno.getSuperId().equals(super_id) && curCAnno.getSuperSuperId().equals(super_super_id)){
//	        			cdsAnno = curCAnno;
//	        		}
//	        		else{
//	        			cdsAnno = mapAnno.get(id+super_id+super_super_id);
//	        			if(!(cdsAnno != null)){
//	        				cdsAnno = new Annotation(id,chr,strand,super_id,super_super_id,type,gene_name);
//	        				mapAnno.put(id+super_id+super_super_id,cdsAnno);
//	        			}
//	        			curCAnno = cdsAnno;
//	        		}
//	        		cds = new Region(start,stop,cdsAnno);
//	        		
//	        		
//	        		if(cdsAnno.isSub(curTAnno)){
//	        			curTrans.add(cds);
//	        		}
//	        		else {
//		        		String transcript_name = attr.get("transcript_name");	        			
//		        		Annotation transAnno = new Annotation(super_id,transcript_name,chr,strand,super_super_id,type,gene_name);
//		        		Transcript trans = new Transcript(start,stop,transAnno);
//
//		        		if(transAnno.isSub(curGAnno)){
//		        			curGene.add(trans);
//		        		}
//		        		else{
//		        			System.out.println(transAnno.getSuperId() + " is missing");
////			        		Annotation geneAnno = new Annotation(super_super_id,gene_name,chr,strand,type);
////			        		Gene gene = new Gene(start, stop, geneAnno);
////			        		geneSet.put(gene.hashCode(),gene);
////			        		curGene = gene;
////			        		curGAnno = geneAnno;
////			        		curGene.add(trans);
//		        		}
//		        		
//		        		curTAnno = transAnno;
//		        		curTrans = trans;
//		        		curTrans.add(cds);
//	        		}
//	        		
//	        	}
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
	
	public void getOutput(String outputPath) throws IOException{
		String tab = "\t";
		String brk = "\n";
		char dop = ':';
		char sip = '|';
		char lin = '-';
		char com = ',';
		char spa = ' ';
		
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
		
		// <Stop <Start, Count>>
		TreeMap<Integer,HashMap<Integer,HashMap<String,Pair<Integer>>>> readMap = new TreeMap<Integer,HashMap<Integer,HashMap<String,Pair<Integer>>>>();
		String curChr = "";
		FileWriter outputWriter = new FileWriter(outputPath,false);
		outputWriter.write(resultBuilder.toString());
		
		boolean nostrand = strandness == 0;
		
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
					
//					if(readname.equals("9887855")){
//						System.out.println("Combined: " +curRead.getAlignmentBlocks().toString());
//						System.out.println("FoP: " + curRead.getAlignmentBlocksFoP().toString());
//						System.out.println("SoP: " +curRead.getAlignmentBlocksSoP().toString());
//						System.out.println("Consistent? " +curRead.isConsistent());
//					}
					char str = curRead.getStrand();
					int mm = curRead.getMM();
					int clip = curRead.getClip();
					int spl = curRead.getSplit();

				    char revstr;
				    
				    if(strandness == 2){
					    if(str == '-'){
					    	str = '+';
					    	revstr = '-';
					    }
					    else{
					    	str = '-';
					    	revstr = '+';
					    }
				    }
				    else{
					    if(str == '-'){
					    	revstr = '+';
					    }
					    else{
					    	revstr = '-';
					    }
				    }
				    
				    String chr = curRead.getChromosome();
				    
				    int start = curRead.getStart();
				    int stop = curRead.getStop();
				    
				    boolean skipped = false;
				    
					if(curRead.isConsistent()){

					    HashSet<Gene> cgenes = new HashSet<Gene>();
					    cgenes = geneTree.get(chr).get(str).getIntervalsSpanning(start, stop, cgenes);
					    HashSet<Gene> igenes = new HashSet<Gene>();
					    igenes = geneTree.get(chr).get(str).getIntervalsSpannedBy(start, stop, igenes);
					    
					    HashSet<Gene> cgenesRev = new HashSet<Gene>();;
					    HashSet<Gene> igenesRev= new HashSet<Gene>();;
					    
					    if(nostrand){

						    cgenesRev = geneTree.get(chr).get(revstr).getIntervalsSpanning(start, stop, cgenesRev);

						    igenesRev = geneTree.get(chr).get(revstr).getIntervalsSpannedBy(start, stop, igenesRev);
						    
					    }
					    
					    if(cgenes.size() > 0 || (cgenesRev.size() > 0  && nostrand)){
					    	resultBuilder.append(readname);
					    	
							resultBuilder.append(tab);
							resultBuilder.append(mmStr);
							resultBuilder.append(mm);
							
							resultBuilder.append(tab);
							resultBuilder.append(clipStr);
							resultBuilder.append(clip);
							
							resultBuilder.append(tab);
							resultBuilder.append(splStr);
							resultBuilder.append(spl);
					    	
					    	ArrayList<String> gList = new ArrayList<String>();
					    	ArrayList<String> tgList = new ArrayList<String>();
					    	StringBuilder tgBuilder = new StringBuilder();
					    	for(Gene g : cgenes){
					    		tgList = g.inTranscript(curRead);
					    		if(tgList.size() > 0){
					    			tgBuilder.append(g.getAnnotation().getId());
					    			tgBuilder.append(com);
					    			tgBuilder.append(g.getBiotype());
					    			tgBuilder.append(dop);
					    			tgBuilder.append(tgList.get(0));
					    			for( int i = 1; i < tgList.size(); i++){
					    				tgBuilder.append(com);
						    			tgBuilder.append(tgList.get(i));
					    			}
					    			
					    			gList.add(tgBuilder.toString());
					    			tgBuilder.setLength(0);
					    		}
					    	}
					    	for(Gene g : cgenesRev){
					    		tgList = g.inTranscript(curRead);
					    		if(tgList.size() > 0){
					    			tgBuilder.append(g.getAnnotation().getId());
					    			tgBuilder.append(com);
					    			tgBuilder.append(g.getBiotype());
					    			tgBuilder.append(dop);
					    			tgBuilder.append(tgList.get(0));
					    			for( int i = 1; i < tgList.size(); i++){
					    				tgBuilder.append(com);
						    			tgBuilder.append(tgList.get(i));
					    			}
					    			
					    			gList.add(tgBuilder.toString());
					    			tgBuilder.setLength(0);
					    		}
					    	}
					    	
					    	if(gList.size() > 0){
						    	resultBuilder.append(tab);
					    		resultBuilder.append(couStr);
					    		resultBuilder.append(gList.size());
					    		
					    		resultBuilder.append(tab);
					    		resultBuilder.append(gList.get(0));
					    		for( int i = 1; i < gList.size(); i++){
				    				resultBuilder.append(sip);
					    			resultBuilder.append(gList.get(i));
				    			}
					    	}
					    	else{
					    		for(Gene g : cgenes){
						    		if(g.inMerged(curRead)){
						    			tgBuilder.append(g.getAnnotation().getId());
						    			tgBuilder.append(com);
						    			tgBuilder.append(g.getBiotype());
						    			tgBuilder.append(dop);
						    			tgBuilder.append(merStr);
						    			gList.add(tgBuilder.toString());
						    			tgBuilder.setLength(0);
						    		}
					    		}
					    		for(Gene g : cgenesRev){
						    		if(g.inMerged(curRead)){
						    			tgBuilder.append(g.getAnnotation().getId());
						    			tgBuilder.append(com);
						    			tgBuilder.append(g.getBiotype());
						    			tgBuilder.append(dop);
						    			tgBuilder.append(merStr);
						    			gList.add(tgBuilder.toString());
						    			tgBuilder.setLength(0);
						    		}
					    		}
					    		if(gList.size() > 0){
							    	resultBuilder.append(tab);
						    		resultBuilder.append(couStr);
						    		resultBuilder.append(gList.size());
						    		
						    		resultBuilder.append(tab);
						    		resultBuilder.append(gList.get(0));
						    		for( int i = 1; i < gList.size(); i++){
					    				resultBuilder.append(sip);
						    			resultBuilder.append(gList.get(i));
					    			}
						    	}
					    		else{
						    		for(Gene g : cgenes){
							    		tgBuilder.append(g.getAnnotation().getId());
							    		tgBuilder.append(com);
							    		tgBuilder.append(g.getBiotype());
							    		tgBuilder.append(dop);
							    		tgBuilder.append(intStr);
							    		gList.add(tgBuilder.toString());
							    		tgBuilder.setLength(0);
						    		}
						    		for(Gene g : cgenesRev){
							    		tgBuilder.append(g.getAnnotation().getId());
							    		tgBuilder.append(com);
							    		tgBuilder.append(g.getBiotype());
							    		tgBuilder.append(dop);
							    		tgBuilder.append(intStr);
							    		gList.add(tgBuilder.toString());
							    		tgBuilder.setLength(0);
						    		}
							    	resultBuilder.append(tab);
						    		resultBuilder.append(couStr);
						    		resultBuilder.append(gList.size());
						    		
						    		resultBuilder.append(tab);
						    		resultBuilder.append(gList.get(0));
						    		for( int i = 1; i < gList.size(); i++){
					    				resultBuilder.append(sip);
						    			resultBuilder.append(gList.get(i));
					    			}
					    		}
					    	}
					    }
					    else{
				    		
					    	if(igenes.size() > 0){
					    		skipped = true;
					    	}
					    	else{
					    		resultBuilder.append(readname);
								resultBuilder.append(tab);
								resultBuilder.append(mmStr);
								resultBuilder.append(mm);
								
								resultBuilder.append(tab);
								resultBuilder.append(clipStr);
								resultBuilder.append(clip);
								
								resultBuilder.append(tab);
								resultBuilder.append(splStr);
								resultBuilder.append(spl);
								
						    	resultBuilder.append(tab);
					    		resultBuilder.append(couStr);
					    		resultBuilder.append(0);
					    		
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
					    		
					    		if(nostrand){
						    		ArrayList<Gene> rneighRev = new ArrayList<Gene>();
						    		rneighRev = geneTree.get(chr).get(revstr).getIntervalsRightNeighbor(start, stop, rneighRev);
						    		long rdistRev = -1;
						    		if(rneighRev.size() > 0){
						    			rdistRev = rneighRev.get(0).getStart() - stop -1;
						    		}
						 
						    		rdist = Math.min(rdist, rdistRev);
						    		
						    		ArrayList<Gene> lneighRev = new ArrayList<Gene>();
						    		lneighRev = geneTree.get(chr).get(revstr).getIntervalsLeftNeighbor(start, stop, lneighRev);
						    		long ldistRev = -1;
						    		if(lneighRev.size() > 0){
						    			ldistRev = start - lneighRev.get(0).getStop() -1;
						    		}
						    		

						    		ldist = Math.min(ldist, ldistRev);

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
					    		
					    		boolean antisense = false;
					    		ArrayList<Gene> antigenes = new ArrayList<Gene>();
							    if(geneTree.get(chr).containsKey(revstr)){
							    	
							    	antigenes = geneTree.get(chr).get(revstr).getIntervalsSpanning(start, stop, antigenes);
							    	if(antigenes.size() > 0){
							    		antisense = true;
							    	}
							    }
							    
							    resultBuilder.append(tab);
					    		resultBuilder.append(ansStr);
					    		resultBuilder.append(false);
					    	}
					    }
					    
					    if(!skipped){
					    	
					    	if(!(curChr.equals(chr))){
					    		readMap.clear();
					    	}
					    	
					    	readMap.headMap(start).clear();
					    	int c = 0;
					    	
					    	//Pair ( + , -)
					    	
					    	String readString = curRead.getAlignmentBlocks().toString();
					    	
					    	if(readMap.containsKey(stop)){
					    		HashMap<Integer,HashMap<String,Pair<Integer>>> stopMap = readMap.get(stop);
					    		if(stopMap.containsKey(start)){
					    			HashMap<String,Pair<Integer>> startMap = stopMap.get(start);
					    			
					    			if(startMap.containsKey(readString)){
					    				Pair<Integer> p = startMap.get(readString);
					    				if(str == '+'){
					    					c = p.p1;
					    					p = new Pair<Integer>(p.p1+1, p.p2);
					    				
					    				}
					    				else{
					    					c = p.p2;
					    					p = new Pair<Integer>(p.p1, p.p2+1);
					    				}
					    				startMap.put(readString, p);
					    			}
					    			else{
					    				if(str == '+'){
						    				startMap.put(readString, new Pair<Integer>(1,0));
						    			}
						    			else{
						    				startMap.put(readString, new Pair<Integer>(0,1));
						    			}
					    			}
					    		}
					    		else{
					    			HashMap<String,Pair<Integer>> startMap = new HashMap<String,Pair<Integer>>();
					    			if(str == '+'){
					    				startMap.put(readString, new Pair<Integer>(1,0));
					    			}
					    			else{
					    				startMap.put(readString, new Pair<Integer>(0,1));
					    			}
					    			stopMap.put(start,startMap);
					    			
					    		}
					    	}
					    	else{
					    		HashMap<String,Pair<Integer>> startMap = new HashMap<String,Pair<Integer>>();
					    		HashMap<Integer,HashMap<String,Pair<Integer>>> stopMap = new HashMap<Integer,HashMap<String,Pair<Integer>>>();
				    			if(str == '+'){
				    				startMap.put(readString, new Pair<Integer>(1,0));
				    			}
				    			else{
				    				startMap.put(readString, new Pair<Integer>(0,1));
				    			}
				    			stopMap.put(start,startMap);
					    		readMap.put(stop, stopMap);
					    	}
					    	
						    resultBuilder.append(tab);
				    		resultBuilder.append(pcrStr);
				    		resultBuilder.append(spa);
				    		resultBuilder.append(c);
				    		
//				    		resultBuilder.append(tab);
//				    		resultBuilder.append(readMap.size());
				    		
				    		resultBuilder.append(brk);
					    }
					    
					    curChr = chr;
					}
					else{
						
					    HashSet<Gene> cgenes = new HashSet<Gene>();
					    cgenes = geneTree.get(chr).get(str).getIntervalsSpanning(start, stop, cgenes);
					    HashSet<Gene> igenes = new HashSet<Gene>();
					    igenes = geneTree.get(chr).get(str).getIntervalsSpannedBy(start, stop, igenes);
					    
					    
					    skipped = cgenes.size() == 0 && igenes.size() > 0;
					    
					    if(!skipped){
							resultBuilder.append(readname);
							resultBuilder.append(tab);
							resultBuilder.append(incStr);
							resultBuilder.append(brk);
					    }
					}
					store.remove(readname);
					outputWriter.write(resultBuilder.toString());
					resultBuilder.setLength(0);
				}
				else{
					store.put(readname, samr);
				}
			}
		}
		bam_it.close();
		
		outputWriter.close();
		
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
	
	class Pair<T> {
	    T p1;
	    T p2;
	    Pair(T p1, T p2) {
	        this.p1 = p1;
	        this.p2 = p2;
	    }
	}

}
