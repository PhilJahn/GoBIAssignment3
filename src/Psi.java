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

public class Psi {

	public static void main(String[] args) {
//		long startTime = System.currentTimeMillis();		
		String gtfPath ="";
		String bamPath ="";
		String outputPath ="";
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
		}
		
		if(gtfPath.equals("") || bamPath.equals("") || outputPath.equals("")){
			System.out.println("Usage Info:\n-gtf <filepath for GTF>\n-bam <filepath for BAM>\n-o <filepath for ouput>");
		}
		else{
			
			String baiPath = bamPath.concat(".bai");
			
			Path gtfFilePath = Paths.get(gtfPath);
			Path bamFilePath = Paths.get(bamPath);
			Path baiFilePath = Paths.get(baiPath);
			Path outputFilePath = Paths.get(outputPath);
			
			File gtfFile = gtfFilePath.toFile();
			File bamFile = bamFilePath.toFile();
			File baiFile = baiFilePath.toFile();
			
			Psi psi = new Psi(gtfFile, bamFile, baiFile);
		
			String output = psi.getOutput();
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
	
	public Psi (File gtfFile, File bamFile, File baiFile) {
//		long startTime = System.currentTimeMillis();
		
		SAMFileReader bam_reader = new SAMFileReader(bamFile,baiFile);
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
	        		String type = lineSplit[1];
	        		char strand = lineSplit[6].charAt(0);
	        		String chr = lineSplit[0];
	        		int start = Integer.parseInt(lineSplit[3]);
	        		int stop = Integer.parseInt(lineSplit[4]);
	        		
	        		Annotation geneAnno = new Annotation(gene_id,gene_name,chr,strand,type);
	        		Gene gene = new Gene(start, stop, geneAnno);
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
	    
	    this.setReadTree(bam_reader);
	    
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
	
	
	public String getOutput(){
//		long startTime = System.currentTimeMillis();
		String tab = "\t";
		String brk = "\n";
		char sep = ':';
		char sip = '|';
		char lin = '-';
//		StringBuilder resultBuilder = new StringBuilder("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tWT_exons\tWT_prots\tSV_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases\n");

		StringBuilder resultBuilder = new StringBuilder("gene\texon\tnum_incl_reads\tnum_excl_reads\tnum_total_reads\tpsi\n");
//
//		for(int key: geneSet.keySet()){
//			Gene gene = geneSet.get(key);
//			Annotation geneAnno = gene.getAnnotation();
//			resultBuilder.append(geneAnno.getId());
//			resultBuilder.append(tab);
//			resultBuilder.append(geneAnno.getName());
//			resultBuilder.append(tab);
//			resultBuilder.append(gene.getStart());
//			resultBuilder.append(sep);
//			resultBuilder.append(gene.getStop());
//			resultBuilder.append(brk);
//		}
		
//		String keyString = "ENSG00000143839.12";
//		int keyInt = keyString.hashCode();
//		Gene gene = geneSet.get(keyInt);
//		geneSet.clear();
//		geneSet.put(keyInt, gene);

		
		for (Integer key: geneSet.keySet()) {
			Gene curGene = geneSet.get(key);
			Annotation curGAnno = curGene.getAnnotation();
//			System.out.println("getExonSkips");
			ArrayList<ExonSkip> skips = curGene.getExonSkips();
			
//			printExonSkip(resultBuilder, curGene, curGAnno, skips, tab, brk, sep, sip);
			
			String geneid = curGAnno.getId();
			char str = curGAnno.getStrand();
			String chr = curGAnno.getChromosome();
			
//			resultBuilder.append((curGene.getStart() + " " + curGene.getStop()));
//			resultBuilder.append(brk);
			
			HashSet<Read> readSet = new HashSet<Read>();
			if(readTree.containsKey(chr)){
	    		if(readTree.get(chr).containsKey(str)){
	    			readSet = readTree.get(chr).get(str).getIntervalsSpannedBy(curGene.getStart(), curGene.getStop(), readSet);
	    		}
			}
			HashSet<String> readIds = curGene.getTransReadIds();
	    	IntervalTree<Read> reads = new IntervalTree<Read>(readSet);
			
//			System.out.println(reads.size());
			
//			System.out.println(reads.toTreeString());
			
			
			
			for(ExonSkip skip: skips){
				ArrayList<Region> exons = new ArrayList<Region>(skip.getWTExons());
				exons.sort(new StartRegionComparator());
				for( Region exon : exons){
					int start = exon.getStart();
					int stop = exon.getStop()+1;
					
					int c = 0;
					int d = 0;
					int ges = 0;
					
					ArrayList<Read> exonReads = new ArrayList<Read>();
					exonReads = reads.getIntervalsIntersecting(start, stop, exonReads);
					
//					System.out.println(exonReads.size());					
//					resultBuilder.append(brk);
					
					exonReads.sort(new StartReadComparator());
					
//					int[] incl = new int[]{67102,67144,67100,67145,67106,67105,67108,67109,67187,67184,67114,67111,67078,67079,67117,67116,67190,67072,67074,67124,67201,67168,67089,67122,67123,67128,67129,67206,67080,67083,67084,67081,67164,67087,67121,67088,67162,67086,67135,67133,67134,67139,67137,67138,67090,67095,67093,67131,67098,67132,67099,67130};
//					int[] excl = new int[]{7022,7043,6991,6971,6970,6896,6929,6909,7019,6941,6988,6987,6924,6926,6947,7010,7153,6980,7034,7015,6984,7075,7030,6916,6974,6996,6930,6977,6976,7229,6978,6937,6914};
//					int[] almexcl = new int[]{7043,6980,6977,6987,6991,7034,7010,6971,7019,6924,6976,7015,6916,6974,6978,6970,6926,6947,7229,6896,6929,7153,6909,7022,6914,6988,6984,6941,7030,6996,6937,6930};
//
//					HashSet<Integer> exclSet = new HashSet<Integer>();
//					for(int i = 0; i < excl.length; i++){
//						exclSet.add(excl[i]);
//					}
//					for(int i = 0; i < almexcl.length; i++){
//						exclSet.remove(almexcl[i]);
//					}
//					
//					System.out.println(exclSet.toString());
					
//					System.out.println(incl.length + " " + excl.length + " " + (incl.length+excl.length));
//					
					for(Read er : exonReads){
						if(readIds.contains(er.getReadName())){
						boolean inincl = false;
						boolean inexcl = false;
//						resultBuilder.append(er.toString());
						
						ges++;
						
						ArrayList<RegionBlock> in = new ArrayList<RegionBlock>();
						in = er.getAlignmentBlocks().getIntervalsSpannedBy(start, stop, in);
						
						ArrayList<RegionBlock> overl = new ArrayList<RegionBlock>();
						overl = er.getAlignmentBlocks().getIntervalsSpannedBy(start, stop, overl);
						
//						resultBuilder.append(tab);						
						if(in.size() > 0){
							if(overl.size() == in.size()){
								c++;
//								resultBuilder.append("incl");
							}
//							else{
//								resultBuilder.append("overl");
//							}
						}
						else{
							d++;
//							resultBuilder.append("excl");
						}
						
//						for(int id : incl){
//							if(er.getReadName().equals("" + id)){
//								inincl = true;
////								System.out.println("id" + id);
//							}
//						}
//						for(int id : excl){
//							if(er.getReadName().equals("" + id)){
//								inexcl = true;
//								System.out.print("," + id);
//							}
//						}
//						
//						if(inincl){
//							resultBuilder.append(tab);
//							resultBuilder.append("incl");
//						}
//						if(inexcl){
//							resultBuilder.append(tab);
//							resultBuilder.append("excl");
//						}
//						else{
//							resultBuilder.append(tab);
//							resultBuilder.append("miss");
//						}
//						
//						resultBuilder.append(brk);
						}
					}
					
					if(ges > 0){
						resultBuilder.append(geneid);
						resultBuilder.append(tab);
						resultBuilder.append(start);
						resultBuilder.append(lin);
						resultBuilder.append(stop);
					
						resultBuilder.append(tab);
						resultBuilder.append(c);
					
						resultBuilder.append(tab);
						resultBuilder.append(d);
	
						resultBuilder.append(tab);
						resultBuilder.append(ges);
					
						double psi = (c+0.0)/ges;
					
						resultBuilder.append(tab);
						resultBuilder.append(psi);
						resultBuilder.append(brk);
					}
					
				}
				
			}
			
		}
//		long stopTime = System.currentTimeMillis();
//		System.out.println("Skips:" + (stopTime-startTime));
		return resultBuilder.toString();
	}
	
	public ArrayList<Gene> getGenes(){
		ArrayList<Gene> result = new ArrayList<Gene>();
		for (Integer key: geneSet.keySet()) {
		   result.add(geneSet.get(key));
		}
		return result;
	}
	
	public void setReadTree(SAMFileReader bam_reader){
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
					if(curRead.isConsistent()){
					    char str = curRead.getStrand();
					    String chr = curRead.getChromosome();
					    
					    int start = curRead.getStart();
					    int stop = curRead.getStop();
					    
					    boolean transcriptomic = false;
					   
					    HashSet<Gene> contin = new HashSet<Gene>();
					    contin = geneTree.get(chr).get(str).getIntervalsSpanning(start, stop, contin);
					    
					    if(contin.size() > 0){
					    	for(Gene g : contin){
//								if(curRead.getReadName().equals("7075")){
//									System.out.println(g.toString());
//								}
					    		
					    		transcriptomic |= g.inTranscript(curRead);
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
					store.remove(readname);
				}
				else{
					store.put(readname, samr);
				}
			}
		}
		bam_it.close();
		
		
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
