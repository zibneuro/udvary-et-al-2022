/****************************************************************************/
/*                                                                          */
/* Program:                                                                 */
/*                                                                          */
/* File:      amiraReader.h                                                 */
/*                                                                          */
/* Purpose:   interface class between Amira SpatialGraph, Surface and       */
/*            other files and the internal spatial graph data structure     */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*            5353 Parkside Drive                                           */
/*            Jupiter, Florida 33458                                        */
/*            USA                                                           */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   22.12.2010                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/* Revised by: 	Daniel Udvary												*/
/* 				In Silico Brain Sciences									*/
/*				Max Planck Institute for Neurobiology of Behavior – caesar	*/
/*				Ludwig-Erhard-Allee 2										*/
/*				53175 Bonn, Germany											*/
/*                                                                          */
/* e-mail: 		daniel.udvary@mpinb.mpg.de									*/
/*                                                                          */
/* Date: 		04.09.2018                                                  */
/*                                                                          */
/****************************************************************************/

#pragma once
#include "typedefs.h"
#include "basics.h"

#ifndef AMIRAREADER
#define AMIRAREADER


class Reader
{
	public:
		Reader(const char * filename);
		Reader(const char * filename, const char * outputFilename);
		~Reader();
		
		// default method for reading Amira SpatialGraph files
		void readSpatialGraphFile(bool applyTransform/*, const char * fname=NULL*/);
		static void readSpatialGraphFile(bool applyTransform, const char * fname, AmiraSpatialGraph * sg);
		// default method for writing Amira SpatialGraph files
		void writeSpatialGraphFile();
		// alternate method for writing Amira SpatialGraph files
		void writeSpatialGraphFileFromEdges();
		
		// default method for loading SpatialGraphSets
		static void readSpatialGraphSetFile(const char * fname, std::vector< unsigned int >& originalGraphIndices,
										std::vector< unsigned int >& spatialGraphSetLabels, std::vector< double * >& spatialGraphTransforms,
										std::vector< std::string >& originalGraphFiles, std::map< unsigned int, std::string >& cellTypeIDLabels);
// 		SpatialGraphSet * loadSpatialGraphSet();
		
		// default method for reading NeuroConv hoc files
		void readHocFile();
		// method for writing NeuroConv hoc files
		void writeHocFile();
		// method for writing NeuroConv hoc files used during
		// registration. writes separate files for neuron
		// morphology and anatomical landmarks
		void writeSeparateHocFiles();
		
		// Method to read hoc file with all contours
		// Pia and alpha -> Pia
		void readHocFileAllPiaContours();

		// default method for writing Amira Surface files
		void writeAmiraSurfaceFile(PolyDataPointerType triangleData);
		// default method for reading Amira Surface files
		PolyDataPointerType readAmiraSurfaceFile();
		
		// default method for reading Amira Landmark files
		PointsPointerType readLandmarkFile();
		// method for reading Amira Landmark files representing
		// closed contours and converting them into closed contour
		// graphs (i.e., PolyData format)
		void readLandmarkFile(bool applyTransform);
		// default method for writing Amira Landmark files
		void writeLandmarkFile(PointsPointerType pts);
		
		// default method for reading Amira scalar field files
		ImageDataPointerType readScalarField();
		// default method for writing Amira scalar field files
		void writeScalarField(ImageDataPointerType field);
		
		// default method for reading Amira Vectorfield files
		ImageDataPointerType readVectorField();
		ImageDataPointerType readNVectorField();
		// default method for writing Amira Vectorfield files
		void writeVectorField(ImageDataPointerType field);
		void writeNVectorField(ImageDataPointerType field);
		
		// read NeuroNet connection matrix tables
		void readConnectionMatrix(ConnectionMatrix * connectome);
		
		// read NeuroNet connection matrix tables
		void readSynapsesPerCellTable(CellTable * table);
		
		// read NeuroNet convergence matrix tables
		void readConvergenceTableCSV(CellTable * table);

		AmiraSpatialGraph * getSpatialGraph() { return inputSpatialGraph; }
		void setSpatialGraph(AmiraSpatialGraph * outputSpatialGraph) { inputSpatialGraph = outputSpatialGraph; }
		
	private:
		const char * inputFilename;
		const char * outputFilename;
		std::string letters;
		std::string numbers;
		std::string signs;
		std::string otherChars;
		std::string whitespace;
		// cell type/column names
		std::list< int > barrelLabels;
		std::map< int, const char * > int2Labels;
		std::map< std::string, int > labels2Int;
		std::list< const char * > hocLabels;
		std::list< const char * > hocNeuronLabels;
		std::list< const char * > hocLandmarkLabels;
		std::list< unsigned int > postTypes;
		std::list< unsigned int > preTypes;
		std::map< std::string, unsigned int > celltypeLabels2Int;
		std::map< unsigned int, const char * > int2CelltypeLabels;
		
		AmiraSpatialGraph * inputSpatialGraph;
		
		void initializeConstants();
};



#endif
