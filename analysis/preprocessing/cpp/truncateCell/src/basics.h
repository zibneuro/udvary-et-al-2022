/****************************************************************************/
/*                                                                          */
/* Program:   CortexCoordinates                                             */
/*                                                                          */
/* File:      basics.h                                                      */
/*                                                                          */
/* Purpose:   basic classes for the internal data structure                 */
/*            SpatialGraph, Edge, Vertex(deprecated) etc.                   */
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
/* Date: 		21.11.2019                                                  */
/*                                                                          */
/****************************************************************************/

#pragma once
#include "typedefs.h"

#ifndef BASICS
#define BASICS

class Vertex
{
	public:
		Vertex(float * _coordinates, int _label);
		Vertex(double * _coordinates, int _label);
		Vertex(Vertex * otherVertex);
		~Vertex();
// 	private:
		double coordinates[3];
		int label;
};

class Edge
{
	public:
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates, float _radius);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates, double _radius);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates, std::list< double > radList );
		Edge(Edge * otherEdge);
		~Edge();
// 	private:
		int edgeConnectivity[2];
		int numEdgePoints;
		int label;
		std::list< double * > edgePointCoordinates;
		double radius;
		std::list< double > radiusList;
		// for hoc file structure
		int fatherID;
		double segmentLength();
};

class AmiraSpatialGraph
{
	public:
		AmiraSpatialGraph();
		~AmiraSpatialGraph();
		
		// appends otherSpatialGraph to this SpatialGraph
		void mergeSpatialGraph(AmiraSpatialGraph * otherSpatialGraph);
		// clears all data from internal data storage
		// use with caution!
		void clear();
		
		void addVertex( Vertex * newVertex );
		void addEdge( Edge * newEdge );
		void addPolyDataObject(PolyDataPointerType object, int label);
		void addLine(double start[3], double end[3]);
		void addLine(double start[3], double end[3], int ID);
		// clears all structures of type 'label' from internal data storage
		// use with caution!
		void removeLabel(int label);
		void removeLandmarkLabels();
		
		// sets internal transformation matrix
		// from double[4][4]
		void setTransformation(double ** newTransform);
		// sets internal transformation matrix
		// from vtkTransform
		void setTransformation(TransformPointerType newTransform);
		// applies transformation matrix to all structures.
		// can only be applied once; to apply more than once,
		// set transformation before each use of this method
		void applyTransformation();
		// applies transformation matrix to all structures of type 'label'.
		// can only be applied once; to apply more than once,
		// set transformation before each use of this method
		void applyTransformation(int label);
		
		void printTransformation();
		
		std::vector< Vertex * >::iterator verticesBegin();
		std::vector< Edge * >::iterator edgesBegin();
		std::vector< Vertex * >::iterator verticesEnd();
		std::vector< Edge * >::iterator edgesEnd();
		std::vector< Vertex * > * verticesPointer(){return &vertices;}
		std::vector< Edge * > * edgesPointer(){return &edges;}
		
		bool isLabelInSpatialGraph(int checkLabel);
		// DEPRECATED: extract all points of landmark 'label' planewise;
		// return their plane indices in zIndexList;
		// sanity check: returns false if no points found
		bool extractLandmark(int label, std::list< std::list< double * > >& planewisePointList, std::list< int >& zIndexList);
		// stores structures of type 'label' in polyData and their plane indices in zIndexList
		// sanity check: returns false if no points found
		bool extractLandmark(int label, PolyDataPointerType polyData, std::list< int >& zIndexList);
		// stores structures of type 'label' in polyData
		// preferred method to interface SpatialGraph structures with VTK
		// sanity check: returns false if no points found
		bool extractLandmark(int label, PolyDataPointerType polyData);
		
		unsigned int getNumberOfVertices() { return vertices.size(); }
		unsigned int getNumberOfEdges() { return edges.size(); }
		unsigned int getNumberOfPoints();
		
		void getBoundingBox(double bounds[6]);
		void getBoundingBox(int label, double bounds[6]);
		
		void setHomeBarrel(int ID){this->homeBarrel = ID;}
		int getHomeBarrel(){return this->homeBarrel;}
		// calculate total segment length from edge startID to soma
		// calculate total segment angle from edge startID to soma
		double totalSegmentLength(int startID);
		double totalSegmentLengthPrecise(int startID);
		double cumulatedSegmentAngle(int startID);
		
		// Manipulates/Truncates Spatial Graph
		// bounds determine slicing bounds of slice [-150 150]
		// cut_COORD determines cutting dimension (X_COORD, Y_COORD, Z_COORD)
		void truncateSpatialGraph(double bounds[2], int cut_COORD);

		void vesselsToPoints();

		void setFatherID();
	
	private:
// 		unsigned int numberOfVertices, numberOfEdges, numberOfPoints;
		
		std::vector< Vertex * > vertices;
		std::vector< Edge * > edges;
		
		int homeBarrel;
		
		bool isIdentity;	//avoid going through all points if transformation == 1
		bool transformationApplied;	//flag to ensure transformation is at most applied once
		double transformation[4][4];
		
		void removeDuplicatePoints();
};

// contains top and bottom center pts as double *
// top and bottom contours in vtkPolyData in order top, bottom
class Column
{
	public:
		Column();
		Column(Column * otherColumn);
		Column(PolyDataPointerType contours, double * top, double * bottom);
		~Column();
		
		double getHeight(){return sqrt((top[0] - bottom[0])*(top[0] - bottom[0]) + (top[1] - bottom[1])*(top[1] - bottom[1]) + (top[2] - bottom[2])*(top[2] - bottom[2]));}
		void getCenter ( double center[3] );
		void translateColumn(const double * shift);
		void rotateColumn(gsl_matrix * rot);
		void rotateColumn(HomogeneousMatrixPointerType mat);
		
		PolyDataPointerType contours;
		double * top;
		double * bottom;
};

class Surface
{
	public:
		Surface(PolyDataPointerType mesh);
		virtual ~Surface() {};
		
		void intersectLine(double * axis, double * center);
		void intersectLineInDirection(double * axis, double * center);
		
		double * getLastIntersectPoint();
		void getLastIntersectPoint(double pt[3]);
		
		vtkIdType getLastIntersectCellID();
		bool isValid(){return dataValid;}
		bool isIntersectionFound(){return intersectionFound;}
		
		PolyDataPointerType ptr(){return surfaceMesh;}
// 		void setSurfaceMesh(PolyDataPointerType mesh){surfaceMesh = mesh;}
		
	protected:
	       PolyDataPointerType surfaceMesh;
	       CellLocatorPointerType locator;
	       
	       bool dataValid, intersectionFound;
	       double * intersectPt;
	       vtkIdType intersectID;
};

class ClosedSurface : public Surface
{
	public:
		ClosedSurface ( PolyDataPointerType mesh );
		~ClosedSurface();
		
		bool isPointInsideSurface ( double pt[3] );
		
	private:
		SelectEnclosedPointsFilterType insideSurfaceFilter;
};

class ConnectionMatrix
{
	public:
		ConnectionMatrix();
		~ConnectionMatrix();
		
		std::map< unsigned int, std::pair< unsigned int, unsigned int > > IDColumnCelltypeMap;
		
		std::list< unsigned int > preTypes;
		std::list< unsigned int > postTypes;
		std::map< unsigned int, unsigned int > preTypeNumbers;
		std::map< unsigned int, unsigned int > postTypeNumbers;
		std::map< unsigned int, std::vector< unsigned int > > preTypeIDs;
		std::map< unsigned int, std::vector< unsigned int > > postTypeIDs;
		
		std::list< unsigned int > preColumns;
		std::list< unsigned int > postColumns;
		std::map< unsigned int, std::vector< unsigned int > > preColumnIDs;
		std::map< unsigned int, std::vector< unsigned int > > postColumnIDs;
		
		std::map< unsigned int, unsigned int > preColumnNumbers;
		std::map< unsigned int, unsigned int > postColumnNumbers;
		std::map< unsigned int, std::map<unsigned int, unsigned int > > columnCelltypePreNumbers;
		std::map< unsigned int, std::map<unsigned int, unsigned int > > columnCelltypePostNumbers;
		
// 		std::map< unsigned int, std::map<unsigned int, float > > prePostMatrix;
// 		std::map< unsigned int, std::map<unsigned int, float > > postPreMatrix;
		std::map< MatrixIndexType, float > matrix;
		
		SelectionType getPreColumnCelltypeSelection(unsigned int column, unsigned int cellType);
		SelectionType getUniqueSelection(SelectionType preSelection);
		SelectionType getPostColumnCelltypeSelection(unsigned int column, unsigned int cellType);
		
		bool simpleMatrix; // if no duplicated innervation values are stored
		std::map< unsigned int, unsigned int > preDuplicatedIDtoSingleID;

		double getAverageInnervation(SelectionType preSelection, SelectionType postSelection);
		double getAverageConvergence(SelectionType preSelection, SelectionType postSelection);
		double getInnervationSum(SelectionType preSelection, SelectionType postSelection);

		/* For mean and SD */
		std::list<double> getInnervationValuesPerPair(SelectionType preSelection, SelectionType postSelection);
		std::list<double> getInnervationValuesUniquePerPair(SelectionType preSelection, SelectionType postSelection);
		std::list<double> getConvergenceValuesPerPair(SelectionType preSelection, SelectionType postSelection);
		std::list<double> getConvergenceValuesUniquePerPair(SelectionType preSelection, SelectionType postSelection);

		std::list<double> getInnervationValuesPerPostsynapse(SelectionType preSelection, SelectionType postSelection);
		std::list<double> getConvergenceValuesPerPostsynapse(SelectionType preSelection, SelectionType postSelection);

		std::list<double> getInnervationValuesPerPresynapse(SelectionType preSelection, SelectionType postSelection);
		std::list<double> getInnervationValuesUniquePerPresynapse(SelectionType preSelection, SelectionType postSelection);
		std::list<double> getConvergenceValuesPerPresynapse(SelectionType preSelection, SelectionType postSelection);
		std::list<double> getConvergenceValuesUniquePerPresynapse(SelectionType preSelection, SelectionType postSelection);

		double getAverageConvergenceOptimized(SelectionType preSelection, SelectionType postSelection);
		double getAverageConvergenceIterateSelection(SelectionType preSelection, SelectionType postSelection);
		double getAverageConvergenceIterateMatrix(SelectionType preSelection, SelectionType postSelection);

		void writeConnectionMatrixAsImage(const char * outputFilename, SelectionType preSelection, SelectionType postSelection);
		void writeSimpleConnectionMatrixAsImage(const char* outputFilename, SelectionType preSelection, SelectionType postSelection,
													bool CP, unsigned int downSamplingFactor);

		
// 	private:
		
		
};

struct CellTableRow
{
	unsigned int cellID;
	unsigned int cellType;
	unsigned int column;
	bool insideColumn;
	std::vector< float > somaLocation;
	float totalSynapses;
	std::vector< float > synapsesPerPreTypeColumn;
};

class CellTable
{
	public:
		CellTable();
		~CellTable();
		// header contains cell type/column pairs pointing
		// to index in synapsesPerPreTypeColumn vector for each row
		std::map< ColumnCellTypePair, unsigned int > header;
		std::vector< CellTableRow * > rows;
		
		// returns vector of all rows that match any cell type/
		// column combination supplied. If any of the column/
		// cell type input lists is empty, it uses only the other
		// list for filtering rows (i.e. leaving the columns list
		// empty will return rows with matching cell type in all columns)
		std::vector< CellTableRow * > getPostColumnCelltypeRows(std::list< unsigned int > columns, std::list< unsigned int > cellTypes);
		// returns list with indices for CellTableRow data vectors
		// of all column/cell type combinations as supplied by the input.
		// If any of the column/cell type input lists is empty, it uses only the other
		// list for filtering rows (i.e. leaving the columns list
		// empty will return rows with matching cell type in all columns)
		std::list< unsigned int > getPreColumnCelltypeColumns(std::list< unsigned int > columns, std::list< unsigned int > cellTypes);

		// Summarize SymLocal Subtyeps
		// Convergence Values have to be weighted according to occurrence of SymLocal Cell Types
		// Requires path to NumberOfCellsTable
		// e.g., fname = "/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_3x3_final/nrCells.csv"
		void mergeLocal(const char* fname);
		// Returns cell count of local subtypes to weight them properly
		std::map< ColumnCellTypePair, unsigned int > getCellNumbers(const char* fname);
};


#endif
