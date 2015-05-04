import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.frame.PlugInFrame;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.plugin.filter.EDM;
import ij.plugin.ImageCalculator;
import ij.plugin.filter.ParticleAnalyzer;

/**
 * Watershed Irregular Features
 * 
 * The plugin might serve to reduce artifacts created by the normal
 * watershed algorithm implicated in ImageJ when applied to objects 
 * of irregular shape containing relatively small connections between them.
 * 
 * Copyright (C), 2014, Jan Brocher / BioVoxxel
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * IN NO EVENT WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES 
 * AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, 
 * INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING 
 * OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO 
 * LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR 
 * THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), 
 * EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF 
 * SUCH DAMAGES.
 * 
 * Please cite BioVoxxel and Jan Brocher when you publish results 
 * obtained by usage of this plugin or a modified version of it
 * 
 * Thank you
 *
 */


public class Watershed_Irregular_Features implements ExtendedPlugInFilter, DialogListener {
	ImagePlus imp;
	private double erosions = 1;
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
			return DOES_8G;
		}
		
	}

	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
		GenericDialog gd = new GenericDialog("Watershed Irregular Features");
			gd.addNumericField("erosion cycle number:", erosions, 1, 5, "");
	                gd.addPreviewCheckbox(pfr);	// passing pfr makes the filter ready for preview
	                gd.addDialogListener(this);	// the DialogItemChanged method will be called on user input
	                gd.addHelp("http://fiji.sc/BioVoxxel_Toolbox#Watershed_Irregular_Structures");
	                gd.showDialog();		// display the dialog; preview runs in the background now
	                if (gd.wasCanceled()) {
	                	return DONE;
	                }
	               	IJ.register(this.getClass());	// protect static class variables (filter parameters) from garbage collection
	        this.pfr = pfr;
	        return IJ.setupDialog(imp, flags); // ask whether to process all slices of stack (if a stack)
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		erosions = (int) gd.getNextNumber();
		if (gd.invalidNumber() || erosions<1) {
			IJ.error("invalid number");
			return false;
		} else {
			return true;
		}
	}

	public void run(ImageProcessor ip) {
		ip.snapshot();
		//boolean blackBackground = Prefs.blackBackground;
		
		boolean invertedLut = ip.isInvertedLut();
		if(invertedLut) {
			ip.invertLut();
		}

		ImagePlus erosionImp = imp.duplicate();
		ImageProcessor erosionIP = erosionImp.getProcessor();

		ImagePlus watershedImp = erosionImp.duplicate();
		ImageProcessor watershedIP = watershedImp.getProcessor();

		//dilate objects according to erosion cycle number entered by user
		for(int n=0; n<=erosions; n++) {
			erosionIP.dilate();
		}

		//separate original objects with the normal watershed algorithm from IJ
		EDM watershedEDM = new EDM();
		watershedEDM.toWatershed(watershedIP);

		//first, the watershed separation lines are extracted
		//then, the second image calculation keeps only those separators 
		//which overlap with an eroded particle
		ImageCalculator calculateImages = new ImageCalculator();
		ImagePlus extractedSeparatorsImp = calculateImages.run("XOR create", watershedImp, imp);	
		
		ImagePlus remainingSeparatorsImp = calculateImages.run("AND create", extractedSeparatorsImp, erosionImp);

		//the remaining separator lines are analyzed to get their starting position
		//for later selection
		int options = ParticleAnalyzer.CLEAR_WORKSHEET|ParticleAnalyzer.RECORD_STARTS;
		int measurements = Measurements.CENTROID;
		ResultsTable rt = new ResultsTable();
		ParticleAnalyzer pa = new ParticleAnalyzer(options, measurements, rt, 0.0, 999999999.9);

		pa.analyze(remainingSeparatorsImp);
		int xStart, yStart;

		watershedIP.setValue(255.0);
		
		//the remaining separation lines in their original size and orientation
		//are copied into the watersheded image to close undesired separations
		for(int r=0; r<rt.getCounter(); r++) {
			xStart = (int) rt.getValue("XStart", r);
			yStart = (int) rt.getValue("YStart", r);
			IJ.doWand(extractedSeparatorsImp, xStart, yStart, 0.0, "8-connected");
			Roi selectedSeparator = extractedSeparatorsImp.getRoi();
			watershedIP.fill(selectedSeparator);
			
		}

		//the corrected watersheded image is copied into the original to be able to
		//undo the watershed if undesired
		calculateImages.run("Copy", imp, watershedImp);
		
		if(invertedLut) {
			ip.invertLut();
		}
		
		imp.updateAndDraw();
				
		//prepare all unused ImagePlus, ImageProcessor and Constructors for carbage collection
		watershedImp.close();
		watershedEDM = null;
		erosionImp.close();
		extractedSeparatorsImp.close();
		remainingSeparatorsImp.close();
		rt = null;
		pa = null;
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
