/****************************************************************************/
/*                                                                          */
/* File:      barrel_field.h                                                */
/*                                                                          */
/* Purpose:   class providing interface and common computations for         */
/*            standard barrel field                                         */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*            5353 Parkside Drive                                           */
/*            Jupiter, Florida 33458                                        */
/*            USA                                                           */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   26.09.2011                                                    */
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
/* Date: 		04.12.2017                                                  */
/*                                                                          */
/****************************************************************************/

#pragma once
#include "typedefs.h"
#include "basics.h"
#include "amiraReader.h"

#ifndef BARREL_FIELD
#define BARREL_FIELD


class BarrelField
{
	public:
		BarrelField();
		BarrelField(const char * filepath);
		~BarrelField();
		
		std::list< int > barrelLabels;
		std::list< int > borderBarrels;
		std::map< int, int > neighborBarrel;
		std::map< int, std::list< int > > barrelGrid;
		std::map< int, double > avgTopDist;
		std::map< int, double > avgPiaWMDist;
		std::map< int, double > avgBarrelArea;
		std::map< int, const char * > int2Labels;
		std::map< std::string, int > labels2Int;
		
		std::map< int, Column * > avgColumns;
		std::map< int, Column * > avgBarrels;
		std::map< int, double * > avgAxes;
		std::map< int, double * > avgCenters;
		Surface * avgPiaSurface;
		Surface * avgWMSurface;
		Surface * avgL4UpperSurface;
		Surface * avgL4LowerSurface;
		
		// returns local z axis at point x
		void localZAxis(double x[3], double zAxis[3]);
		// returns local z axis at and closest point on column axis to point x
		void localZAxis(double x[3], double closestPt[3], double zAxis[3]);
		// returns local z axis at and closest point on column axis to point x
		// (verbose version, not used)
		void localZAxis(double x[3], double closestPt[3], bool verbose, double zAxis[3]);
		// returns (column-) z axis of barrel HBID
		void localZAxis(int HBID, double zAxis[3]);
		// creates unit vectors of coordinate system at point x
		// and stores them in user-supplied double[3]
		void localCoordinateSystem(double x[3], double xAxis[3], double yAxis[3], double zAxis[3]);
		// computes coordinates of point x in cylindrical coordinate system
		// centered on closest column and stores them in user-supplied double[3]
		void localCylindricalCoordinates(double x[3], double coordinates[3]);
		// computes coordinates of point x in cylindrical coordinate system
		// centered on column 'HBID' and stores them in user-supplied double[3]
		void localCylindricalCoordinates(double x[3], int HBID, double coordinates[3]);
		
		// returns ID of barrel (column)
		// closest to point x
		int closestBarrel(double x[3]);
		// checks whether point x actually
		// lies inside S1, taking both convex hull
		// and Pia and WM into account.
		// preferred method for checking this.
		bool isInsideS1(double x[3]);
		// returns the laminar position
		// of point x:
		// 1 - SUPRA; 2 - GRAN; 3 - INFRA
		// 0 - Other (e.g. below WM)
		int laminarPosition(double x[3]);
		// returns column ID of closest
		// barrel column if inside
		// otherwise, returns 0 (septum)
		int insideColumn(double x[3]);
		// returns distance to Pia along local
		// z-axis; returns -1 in case of errors
		double piaDistance(double x[3]);
	
	private:
		AmiraSpatialGraph * avgBarrelField;
		ImageDataPointerType axisVecField;
		PolyDataPointerType avgAxesField;
		PolyDataPointerType avgPia;
		PolyDataPointerType avgWM;
		PolyDataPointerType avgL4Upper;
		PolyDataPointerType avgL4Lower;
		PolyDataPointerType S1HullData;
		ClosedSurface * S1ConvexHull;
		
		void initializeConstants();
		void readStandardBarrelField();
		void readStandardBarrelField(const char * filepath);
		
};

#endif
