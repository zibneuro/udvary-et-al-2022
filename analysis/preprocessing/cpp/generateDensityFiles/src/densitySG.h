/****************************************************************************/
/*                                                                          */
/* File:    densitySG.h                                                   	*/
/*                                                                          */
/* Purpose: Compute Density of SpatialGraph						        	*/
/*                                                                          */
/* Author: 	Daniel Udvary													*/
/* 			In Silico Brain Sciences										*/
/*			Max Planck Institute for Neurobiology of Behavior – caesar		*/
/*			Ludwig-Erhard-Allee 2											*/
/*			53175 Bonn, Germany												*/
/*                                                                          */
/* e-mail: 	daniel.udvary@mpinb.mpg.de										*/
/*                                                                          */
/* Date: 	16.10.2018                                                  	*/
/*                                                                          */
/****************************************************************************/

#pragma once
#include "typedefs.h"
#include "basics.h"
#include "amiraReader.h"
#include "helper.h"

#ifndef DENSITYSG
#define DENSITYSG

class densitySG
{
	public:
		/* Constructor either with given SpatialGraph or with path/to/spatialGraph*/
		densitySG(AmiraSpatialGraph * sg, int label); 
		densitySG(const char * inputfilename, int label); 
		
		/* Length Functions */
		ImageDataPointerType computeLength();
		double getLengthInBoxSimpleTotal(double input_bounds[6]);
		double getLengthInBoxSimple(double input_bounds[6]);
		double getLengthInBoxSimple(double input_bounds[6], bool showErrorMessage);

		void addIntersectionPts(double input_bounds[6]);
		void addIntersectionPts(ImageDataPointerType volume);
		bool addIntersectionPtsInVoxel(double bounds[6]);

		// Old FUNCTIONS
		//bool arePtsInSameBox(double pt1[3], double pt2[3], double minBox[3]);
		//bool arePtsInSameRange(double x1, double x2, double minX, double maxX);
		//bool onEdge(double x, double minX, double maxX);
		//std::list< double > getIntersectAllT(double xBox[2], double yBox[2], double zBox[2], double x1, double y1, double z1, double x2, double y2, double z2);
		//double computeIntersectT(double Dst1, double Dst2);
		//void getIntersectFromT(double t, double pt1[3], double pt2[3], double intersectPt[3]);
		
		/* BP Functions */
		ImageDataPointerType computeBP();
		std::list< int > getBranchPointIDs();
		bool isSamePoint(double x1, double y1, double z1, double x2, double y2, double z2); 
		int getBranchingPtInBox(double input_bounds[6]);
		
		/* TP Functions */
		ImageDataPointerType computeTP();
		double getTerminalPtInBox(double input_bounds[6]);
		
		/* Additional Functions */
		void computeParas(); 
		void extendBounds(double inputbounds[6]);
		void printProperties();
		void getVoxelBoundingBox(double boundsVoxel[6],double origin[3],
									int x,int y,int z, double InputVoxelSize);
		int ptInBox(double x, double minX, double maxX);
		bool inBox(double input_bounds[6], double x, double y, double z);

		/* Getter and Setter Methods */
		void setVoxelSize(double VoxelSize); 
		double getVoxelSize() {return VoxelSize;};
		void setLabel(int label); 
		int getLabel() {return label;};
		void setOutputpath(std::string outputpath); 
		std::string getOutputpath() {return outputpath;};
		bool getBoolAmirafile() {return boolAmirafile;}
		void setBoolAmirafile(bool boolAmirafile);
		void setEpsilon(double epsilon); 
		double getEpsilon() {return epsilon;};
		void setBounds(double inputbounds[6]); 
		void getBounds(double inputbounds[6]);
		AmiraSpatialGraph * getSpatialGraph() {return sg;};
	
	private:
		AmiraSpatialGraph * sg; 
		std::string outputpath; 
		double VoxelSize;
		int label; 
		double epsilon; 
		bool boolAmirafile;
		double bounds[6]; 
};

#endif
