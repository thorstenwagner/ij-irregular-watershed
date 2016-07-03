package de.biomedical_imaging.ij.irregular_watershed;


import ij.*;
import ij.process.*;
import ij.blob.Blob;
import ij.blob.ManyBlobs;
import ij.gui.*;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;

import ij.plugin.filter.*;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.frame.PlugInFrame;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.plugin.filter.EDM;
import ij.plugin.ImageCalculator;
import ij.plugin.filter.ParticleAnalyzer;

/**
 * MIT License
 * 
 * Copyright (c) 2016 Thorsten Wagner
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * Please cite BioVoxxel and Jan Brocher when you publish results 
 * obtained by usage of this plugin or a modified version of it
 * 
 * Thank you
 * 
 * May 04, 2015: Bug fix: Now correctly works with stacks (Thorsten Wagner, wagner@biomedical-imaging.de)
 * May 21, 2015: Now uses Jan's EDM based erosion and adds convexity based watershed which is scale invariant 
 *
 */


public class Watershed_Irregular_Features implements ExtendedPlugInFilter, DialogListener {
	ImagePlus imp;
	private double erosions = 1;
	private double convexityThreshold = 0;
	private PlugInFilterRunner pfr;
	private int nPasses = 1;
	private int pass;
	private int flags = DOES_8G|KEEP_PREVIEW|SNAPSHOT;
	public int setup(String arg, ImagePlus imp) {
		this.imp = imp; 
		if(!imp.getProcessor().isBinary()) {
			IJ.error("works only on 8-bit binary images");
			return DONE;
		} else {
			return DOES_8G | DOES_STACKS;
		}
		
	
	}
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
		GenericDialog gd = new GenericDialog("Watershed Irregular Features");
		gd.addNumericField("erosion cycle number:", erosions, 0, 5, "");
		gd.addNumericField("Convexity_Threshsold", 0, 2);
		gd.addPreviewCheckbox(pfr); // passing pfr makes the filter ready for preview
		gd.addDialogListener(this); // the DialogItemChanged method will be called on user input
		gd.addHelp("http://fiji.sc/BioVoxxel_Toolbox#Watershed_Irregular_Structures");
		gd.showDialog(); // display the dialog; preview runs in the background now
		if (gd.wasCanceled()) {
			return DONE;
		}
		IJ.register(this.getClass()); // protect static class variables (filter parameters) from garbage collection
		this.pfr = pfr;
		return IJ.setupDialog(imp, flags); // ask whether to process all slices of stack (if a stack)
	}
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		erosions = (int) gd.getNextNumber();
		convexityThreshold = gd.getNextNumber();
		if (gd.invalidNumber() || erosions<1) {
			Checkbox previewCheckbox = (Checkbox) gd.getCheckboxes().firstElement();
			if (previewCheckbox.getState()) {
				previewCheckbox.setSize(130, 20);
				previewCheckbox.setLabel("Invalid number");
			}
			return false;
		} else {
			return true;
		}
	}
	public void run(ImageProcessor ip) {
		// ip.snapshot();
		Prefs.blackBackground = true;
		boolean invertedLut = ip.isInvertedLut();
		if(invertedLut) {
			ip.invertLut();
		}
		ImagePlus origImp = new ImagePlus("",ip.duplicate());
		ImageProcessor erosionIP = ip.duplicate();
		ImagePlus erosionImp = new ImagePlus("", erosionIP);
		ImageProcessor watershedIP = ip.duplicate();
		ImagePlus watershedImp = new ImagePlus("", watershedIP);
		
		
		EDM edm = new EDM();
		
		
		if(convexityThreshold==0){
			//If the convexity threshold is set to 0 simple erode the image n times.
			edm.toEDM(erosionIP);
			erosionIP.threshold((int)erosions);
		}
		else{
			//Do a seperate erosion for each connectec component (CC). If
			//a CC has a convexity larger than the convexity threshold stop the erosion
			//for this object
			ArrayList<Blob> objects = new ArrayList<Blob>();
			
			
			erosionImp = new ImagePlus("", erosionIP);
			ManyBlobs mb = new ManyBlobs(erosionImp);
			mb.setBackground(0);
			mb.findConnectedComponents();
			while(mb.size()>0) {
				for (Blob blob : mb) {
					if(blob.getConvexity() >convexityThreshold){
						objects.add(blob);
						erosionIP.setColor(Color.black);
						
						erosionIP.fillPolygon(blob.getOuterContour());
					}
				}
				edm.toEDM(erosionIP);
				erosionIP.threshold(2);
				
				erosionImp = new ImagePlus("", erosionIP);
	
				mb = new ManyBlobs(erosionImp);
				mb.setBackground(0);
				mb.findConnectedComponents();

				
			}
			erosionIP.set(0);
			for (Blob blob : objects) {
				blob.setDefaultColor(Color.white);
				blob.draw(erosionIP);
			}
		}
		
		//separate original objects with the normal watershed algorithm from IJ
		
		edm.toWatershed(watershedIP);
		
		
		//first, the watershed separation lines are extracted
		//then, the second image calculation keeps only those separators
		//which overlap with an eroded particle
		ImageCalculator calculateImages = new ImageCalculator();
		ImagePlus extractedSeparatorsImp = calculateImages.run("XOR create", watershedImp, origImp);
		ImagePlus remainingSeparatorsImp = calculateImages.run("AND create", extractedSeparatorsImp, erosionImp);
		
		//the remaining separator lines are analyzed to get their starting position
		//for later selection
		int options = ParticleAnalyzer.RECORD_STARTS;
		int measurements = Measurements.CENTROID;
		ResultsTable resultsTable = new ResultsTable();
		ParticleAnalyzer pa = new ParticleAnalyzer(options, measurements, resultsTable, 0.0, 999999999.9);
		pa.analyze(remainingSeparatorsImp);
		int xStart, yStart;
		watershedIP.setValue(255.0);
		
		//the remaining separation lines in their original size and orientation
		//are copied into the watersheded image to close undesired separations
		for(int r=0; r<resultsTable.getCounter(); r++) {
			xStart = (int)resultsTable.getValue("XStart", r);
			yStart = (int)resultsTable.getValue("YStart", r);
			IJ.doWand(extractedSeparatorsImp, xStart, yStart, 0.0, "8-connected");
			Roi selectedSeparator = extractedSeparatorsImp.getRoi();
			watershedIP.fill(selectedSeparator);
		}
		//watershedImp.show();
		//the corrected watersheded image is copied into the original to be able to
		//undo the watershed if undesired
		//origImp.getProcessor().setPixels(watershedImp.getProcessor().getPixels());
		
		if(invertedLut) {
			ip.invertLut();
		}
		for(int i = 0; i < ip.getWidth(); i++){
			for(int j =0; j < ip.getHeight(); j++){
				ip.putPixel(i, j, watershedImp.getProcessor().getPixel(i, j));
			}
		}

	}
	public void setNPasses (int nPasses) {
		this.nPasses = nPasses;
		pass = 0;
	}
	void showProgress(double percent) {
		percent = (double)(pass-1)/nPasses + percent/nPasses;
		IJ.showProgress(percent);
	}
}
