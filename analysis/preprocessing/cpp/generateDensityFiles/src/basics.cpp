/****************************************************************************/
/*                                                                          */
/* Program:   CortexCoordinates                                             */
/*                                                                          */
/* File:      basics.cpp                                                    */
/*                                                                          */
/* Purpose:   basic classes for the internal data structure                 */
/*            SpatialGraph, Edge, Vertex(deprecated)                        */
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
/* Date: 		16.02.2020                                                  */
/*                                                                          */
/****************************************************************************/

#include "basics.h"

Vertex::Vertex(float * _coordinates, int _label)
{
	coordinates[0] = (double)_coordinates[0];
	coordinates[1] = (double)_coordinates[1];
	coordinates[2] = (double)_coordinates[2];
	label = _label;
};

Vertex::Vertex(double * _coordinates, int _label)
{
	coordinates[0] = _coordinates[0];
	coordinates[1] = _coordinates[1];
	coordinates[2] = _coordinates[2];
	label = _label;
};

Vertex::Vertex(Vertex * otherVertex)
{
	this->coordinates[0] = otherVertex->coordinates[0];
	this->coordinates[1] = otherVertex->coordinates[1];
	this->coordinates[2] = otherVertex->coordinates[2];
	this->label = otherVertex->label;
}

Vertex::~Vertex()
{
	//tbd
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	edgePointCoordinates = _edgePointCoordinates;
	radius = 0.0;
	fatherID = -1;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates, double _radius)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	edgePointCoordinates = _edgePointCoordinates;
	radius = _radius;
	fatherID = -1;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates, std::list< double > radList )
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	edgePointCoordinates = _edgePointCoordinates;
	radiusList = radList;
	fatherID = -1;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	std::list< float * >::iterator iter;
	for(iter = _edgePointCoordinates.begin(); iter != _edgePointCoordinates.end(); ++iter)
		edgePointCoordinates.push_back((double *)(*iter));
	radius = 0.0;
	fatherID = -1;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates, float _radius)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	std::list< float * >::iterator iter;
	for(iter = _edgePointCoordinates.begin(); iter != _edgePointCoordinates.end(); ++iter)
		edgePointCoordinates.push_back((double *)(*iter));
	radius = _radius;
	fatherID = -1;
};

Edge::Edge(Edge * otherEdge)
{
	this->edgeConnectivity[0] = otherEdge->edgeConnectivity[0];
	this->edgeConnectivity[1] = otherEdge->edgeConnectivity[1];
	this->label = otherEdge->label;
	this->numEdgePoints = otherEdge->numEdgePoints;
	this->radius = otherEdge->radius;
	std::list< double >::const_iterator radiusListIt;
	for(radiusListIt = otherEdge->radiusList.begin(); radiusListIt != otherEdge->radiusList.end(); ++radiusListIt)
		this->radiusList.push_back(*radiusListIt);
	std::list< double *  >::const_iterator edgePtListIt;
	for(edgePtListIt = otherEdge->edgePointCoordinates.begin(); edgePtListIt != otherEdge->edgePointCoordinates.end(); ++edgePtListIt)
	{
		double * pt = new double [3];
		pt[0] = (*edgePtListIt)[0], pt[1] = (*edgePtListIt)[1], pt[2] = (*edgePtListIt)[2];
		this->edgePointCoordinates.push_back(pt);
	}
	this->fatherID = otherEdge->fatherID;
}

Edge::~Edge()
{
	edgePointCoordinates.clear();
	radiusList.clear();
}

double Edge::segmentLength()
{
	if(!edgePointCoordinates.size())
		return 0;
	
	double length = 0, * lastPt = edgePointCoordinates.front();
	std::list< double * >::const_iterator ptIt = edgePointCoordinates.begin();
	++ptIt;
	while(ptIt != edgePointCoordinates.end())
	{
		double * nextPt = *ptIt;
		length += sqrt((nextPt[0] - lastPt[0])*(nextPt[0] - lastPt[0]) + (nextPt[1] - lastPt[1])*(nextPt[1] - lastPt[1]) + (nextPt[2] - lastPt[2])*(nextPt[2] - lastPt[2]));
		lastPt = nextPt;
		++ptIt;
	}
	return length;
};

AmiraSpatialGraph::AmiraSpatialGraph()
{
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
		{
			if(ii != jj)
				transformation[ii][jj] = 0;
			else
				transformation[ii][jj] = 1;
		}
	isIdentity = 1;
	transformationApplied = 0;
	homeBarrel = 0;
};

AmiraSpatialGraph::~AmiraSpatialGraph()
{
	std::vector< Edge * >::iterator it1;
	std::vector< Vertex * >::iterator it2;
	for(it1 = edges.begin(); it1 != edges.end(); ++it1)
		delete *it1;
	for(it2 = vertices.begin(); it2 != vertices.end(); ++it2)
		delete *it2;
	
	vertices.clear();
	edges.clear();
};

unsigned int AmiraSpatialGraph::getNumberOfPoints()
{
	unsigned int number = 0;
	std::vector< Edge * >::iterator it1;
	for(it1 = edges.begin(); it1 != edges.end(); ++it1)
	{
		number += (*it1)->numEdgePoints;
	}
	
	return number;
}

void AmiraSpatialGraph::getBoundingBox(double bounds[6])
{
	double xMin = 1E09, xMax = -1E09, yMin = 1E09, yMax = -1E09, zMin = 1E09, zMax = -1E09;
	bool hasPoints = 0;
	for(int ii = 0; ii < edges.size(); ++ii)
	{
		std::list< double * >::const_iterator edgePtListIt;
		for(edgePtListIt = edges[ii]->edgePointCoordinates.begin(); edgePtListIt != edges[ii]->edgePointCoordinates.end(); ++edgePtListIt)
		{
			hasPoints = 1;
			double * pt = *edgePtListIt;
			if(pt[0] < xMin)
				xMin = pt[0];
			if(pt[0] > xMax)
				xMax = pt[0];
			if(pt[1] < yMin)
				yMin = pt[1];
			if(pt[1] > yMax)
				yMax = pt[1];
			if(pt[2] < zMin)
				zMin = pt[2];
			if(pt[2] > zMax)
				zMax = pt[2];
		}
	}
	
	if(hasPoints)
	{
		bounds[0] = xMin;
		bounds[1] = xMax;
		bounds[2] = yMin;
		bounds[3] = yMax;
		bounds[4] = zMin;
		bounds[5] = zMax;
	}
	else
	{
		bounds[0] = 0;
		bounds[1] = 0;
		bounds[2] = 0;
		bounds[3] = 0;
		bounds[4] = 0;
		bounds[5] = 0;
	}
}

void AmiraSpatialGraph::getBoundingBox(int label, double bounds[6])
{
	if(isLabelInSpatialGraph(label))
	{
		double xMin = 1E09, xMax = -1E09, yMin = 1E09, yMax = -1E09, zMin = 1E09, zMax = -1E09;
		bool hasPoints = 0;
		for(int ii = 0; ii < edges.size(); ++ii)
		{
			if(edges[ii]->label == label)
			{
				std::list< double * >::const_iterator edgePtListIt;
				for(edgePtListIt = edges[ii]->edgePointCoordinates.begin(); edgePtListIt != edges[ii]->edgePointCoordinates.end(); ++edgePtListIt)
				{
					hasPoints = 1;
					double * pt = *edgePtListIt;
					if(pt[0] < xMin)
						xMin = pt[0];
					if(pt[0] > xMax)
						xMax = pt[0];
					if(pt[1] < yMin)
						yMin = pt[1];
					if(pt[1] > yMax)
						yMax = pt[1];
					if(pt[2] < zMin)
						zMin = pt[2];
					if(pt[2] > zMax)
						zMax = pt[2];
				}
			}
		}
		
		if(hasPoints)
		{
			bounds[0] = xMin;
			bounds[1] = xMax;
			bounds[2] = yMin;
			bounds[3] = yMax;
			bounds[4] = zMin;
			bounds[5] = zMax;
		}
		else
		{
			bounds[0] = 0;
			bounds[1] = 0;
			bounds[2] = 0;
			bounds[3] = 0;
			bounds[4] = 0;
			bounds[5] = 0;
		}
	}
	else
	{
		bounds[0] = 0;
		bounds[1] = 0;
		bounds[2] = 0;
		bounds[3] = 0;
		bounds[4] = 0;
		bounds[5] = 0;
	}
}

void AmiraSpatialGraph::addVertex(Vertex * newVertex)
{
	vertices.push_back(newVertex);
};

void AmiraSpatialGraph::addEdge(Edge * newEdge)
{
	edges.push_back(newEdge);
};

void AmiraSpatialGraph::setTransformation(double ** newTransform)
{
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
			transformation[ii][jj] = newTransform[ii][jj];
	
	bool isId = 1;
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
		{
			if(ii != jj)
				if(fabs(transformation[ii][jj]) > 1E-6)
				{
					isId = 0;
// 					std::cout << "T[" << ii << "][" << jj << "] = " << transformation[ii][jj] << std::endl;
				}
			
			else
				if(fabs(transformation[ii][jj] - 1) > 1E-6)
				{
					isId = 0;
// 					std::cout << "T[" << ii << "][" << jj << "] = " << transformation[ii][jj] << std::endl;
				}
		}
	
	isIdentity = isId;
	transformationApplied = 0;
};

void AmiraSpatialGraph::setTransformation(TransformPointerType newTransform)
{
	HomogeneousMatrixPointerType mat = newTransform->GetMatrix();
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
			transformation[ii][jj] = mat->GetElement(ii, jj);
	
	bool isId = 1;
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
		{
			if(ii != jj)
				if(fabs(transformation[ii][jj]) > 1E-6)
				{
					isId = 0;
// 					std::cout << "T[" << ii << "][" << jj << "] = " << transformation[ii][jj] << std::endl;
				}
			
			else
				if(fabs(transformation[ii][jj] - 1) > 1E-6)
				{
					isId = 0;
// 					std::cout << "T[" << ii << "][" << jj << "] = " << transformation[ii][jj] << std::endl;
				}
		}
	
	isIdentity = isId;
	transformationApplied = 0;
};

void AmiraSpatialGraph::applyTransformation()
{
	if(!isIdentity && !transformationApplied)
	{
// 		printTransformation();
		transformationApplied = 1;
// 		std::flush(std::cout << "Transforming " << vertices.size() << " vertices..." << std::endl);
		std::vector< Vertex * >::iterator vertexIt;
// 		int vertexCount = 1;
		for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
		{
// 			std::cout << "Transforming Vertex " << vertexCount << " of " << vertices.size() << std::endl;
// 			++vertexCount;
			double oldCoords[4], newCoords[4];
			for(int ii = 0; ii < 3; ++ii)
			{
				oldCoords[ii] = (*vertexIt)->coordinates[ii];
				newCoords[ii] = 0;
			}
			oldCoords[3] = 1;
			newCoords[3] = 1;
			for(int ii = 0; ii < 3; ++ii)
				for(int jj = 0; jj < 4; ++jj)
					newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
			
			for(int ii = 0; ii < 3; ++ii)
				(*vertexIt)->coordinates[ii] = newCoords[ii];
		}
		
// 		std::flush(std::cout << "Transforming " << edges.size() << " edges..." << std::endl);
		std::vector< Edge * >::iterator edgeIt;
// 		int edgeCount = 1;
		for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
		{
// 			std::cout << "Transforming Edge " << edgeCount << " of " << edges.size() << std::endl;
// 			std::cout << (*edgeIt)->edgePointCoordinates.size() << " points" << std::endl;
// 			++edgeCount;
			std::list< double * >::iterator edgeListIt;
			for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{
				double oldCoords[4], newCoords[4];
				for(int ii = 0; ii < 3; ++ii)
				{
					oldCoords[ii] = (*edgeListIt)[ii];
					newCoords[ii] = 0;
				}
				oldCoords[3] = 1;
				newCoords[3] = 1;
				for(int ii = 0; ii < 3; ++ii)
					for(int jj = 0; jj < 4; ++jj)
						newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
				
				for(int ii = 0; ii < 3; ++ii)
					(*edgeListIt)[ii] = newCoords[ii];
			}
		}
// 		std::flush(std::cout << "done!" << std::endl);
	}
};

// apply transformation only to one specific label
void AmiraSpatialGraph::applyTransformation(int label)
{
	if(!isIdentity && !transformationApplied)
	{
// 		printTransformation();
		transformationApplied = 1;
		std::vector< Vertex * >::iterator vertexIt;
		int vertexCount = 0;
		for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
		{
			if((*vertexIt)->label == label)
			{
				++vertexCount;
				double oldCoords[4], newCoords[4];
				for(int ii = 0; ii < 3; ++ii)
				{
					oldCoords[ii] = (*vertexIt)->coordinates[ii];
					newCoords[ii] = 0;
				}
				oldCoords[3] = 1;
				newCoords[3] = 1;
				for(int ii = 0; ii < 3; ++ii)
					for(int jj = 0; jj < 4; ++jj)
						newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
				
				for(int ii = 0; ii < 3; ++ii)
					(*vertexIt)->coordinates[ii] = newCoords[ii];
			}
		}
// 		std::cout << "Transforming " << vertexCount << " vertices..." << std::endl;
		
		std::vector< Edge * >::iterator edgeIt;
		int edgeCount = 0;
		for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
		{
			if((*edgeIt)->label == label)
			{
				++edgeCount;
				std::list< double * >::iterator edgeListIt;
				for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
				{
					double oldCoords[4], newCoords[4];
					for(int ii = 0; ii < 3; ++ii)
					{
						oldCoords[ii] = (*edgeListIt)[ii];
						newCoords[ii] = 0;
					}
					oldCoords[3] = 1;
					newCoords[3] = 1;
					for(int ii = 0; ii < 3; ++ii)
						for(int jj = 0; jj < 4; ++jj)
							newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
					
					for(int ii = 0; ii < 3; ++ii)
						(*edgeListIt)[ii] = newCoords[ii];
				}
			}
		}
// 		std::cout << "Transforming " << edgeCount << " edges..." << std::endl;
	}
};

void AmiraSpatialGraph::printTransformation()
{
	std::cout << "SpatialGraph transformation matrix:" << std::endl << std::endl;
	for(int ii = 0; ii < 4; ++ii)
	{
		std:: cout << "[";
		for(int jj = 0; jj < 3; ++jj)
			std::cout << transformation[ii][jj] << ", ";
		std::cout << transformation[ii][3] << "]" << std::endl;
	}
	std::cout << std::endl;
	std::cout << "isIdentity = " << isIdentity << std::endl;
	std::cout << "transformationApplied = " << transformationApplied << std::endl;
};

std::vector< Vertex * >::iterator AmiraSpatialGraph::verticesBegin()
{
	std::vector< Vertex * >::iterator it = vertices.begin();
	return it;
};

std::vector< Edge * >::iterator AmiraSpatialGraph::edgesBegin()
{
	std::vector< Edge * >::iterator it = edges.begin();
	return it;
};

std::vector< Vertex * >::iterator AmiraSpatialGraph::verticesEnd()
{
	std::vector< Vertex * >::iterator it = vertices.end();
	return it;
};

std::vector< Edge * >::iterator AmiraSpatialGraph::edgesEnd()
{
	std::vector< Edge * >::iterator it = edges.end();
	return it;
};

void AmiraSpatialGraph::mergeSpatialGraph(AmiraSpatialGraph * otherSpatialGraph)
{
	unsigned int oldVertexNr = vertices.size();
	std::vector< Vertex * >::iterator vertexIt;
	for(vertexIt = otherSpatialGraph->verticesBegin(); vertexIt != otherSpatialGraph->verticesEnd(); ++vertexIt)
	{
		Vertex * newVertex = new Vertex(*vertexIt);
		this->vertices.push_back(newVertex);
	}
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = otherSpatialGraph->edgesBegin(); edgeIt != otherSpatialGraph->edgesEnd(); ++edgeIt)
	{
		Edge * newEdge = new Edge(*edgeIt);
		newEdge->edgeConnectivity[0] += oldVertexNr;
		newEdge->edgeConnectivity[1] += oldVertexNr;
		this->edges.push_back(newEdge);
	}
};

// resets SpatialGraph. use with extreme caution!!! all data will be gone
void AmiraSpatialGraph::clear()
{
	std::vector< Edge * >::iterator edgeIt;
	std::vector< Vertex * >::iterator vertexIt;
	for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
		delete *edgeIt;
	for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
		delete *vertexIt;
	edges.clear(), vertices.clear();
	homeBarrel = 0;
	isIdentity = 1;
	transformationApplied = 0;
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
			transformation[ii][jj] = (ii == jj) ? 1 : 0;
};

void AmiraSpatialGraph::vesselsToPoints()
{
	std::vector< Edge * >::iterator edgeIt;
	std::vector< Vertex * >::iterator vertexIt;
	std::vector< Vertex * > verticesVec;
	
// 	std::flush(std::cout << "converting blood vessels to points..." << std::endl);
	
	for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
	{
		verticesVec.push_back(*vertexIt);
	}
	
	for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
	{
		if((*edgeIt)->label == Vessel)
		{
			double circ = 0;
			double * center = new double[3];
			center[0] = 0, center[1] = 0, center[2] = 0;
			double * curr;
			double * last;
			std::list< double * >::iterator it = (*edgeIt)->edgePointCoordinates.begin();
			last = *it;
			++it;
			while(it != (*edgeIt)->edgePointCoordinates.end())
			{
				curr = *it;
				circ += std::sqrt((curr[0] - last[0])*(curr[0] - last[0]) + (curr[1] - last[1])*(curr[1] - last[1]) + (curr[2] - last[2])*(curr[2] - last[2]));
				center[0] += curr[0];
				center[1] += curr[1];
				center[2] += curr[2];
				last = curr;
				++it;
			}
			center[0] /= (double)((*edgeIt)->numEdgePoints - 1);
			center[1] /= (double)((*edgeIt)->numEdgePoints - 1);
			center[2] /= (double)((*edgeIt)->numEdgePoints - 1);
			(*edgeIt)->numEdgePoints = 1;
			(*edgeIt)->edgePointCoordinates.clear();
			(*edgeIt)->edgePointCoordinates.push_back(center);
			(*edgeIt)->radius = circ/(2*PI);
			
			for(int ii = 0; ii < 3; ++ii)
			{
				verticesVec[(*edgeIt)->edgeConnectivity[0]]->coordinates[ii] = center[ii];
				verticesVec[(*edgeIt)->edgeConnectivity[1]]->coordinates[ii] = center[ii];
			}
		}
	}
	
	int ii = 0;
	for(vertexIt = vertices.begin(); vertexIt != vertices.end() && ii < verticesVec.size(); ++vertexIt, ++ii)
	{
		for(int jj = 0; jj < 3; ++jj)
			(*vertexIt)->coordinates[jj] = verticesVec[ii]->coordinates[jj];
	}
	
	removeDuplicatePoints();
};

void AmiraSpatialGraph::removeDuplicatePoints()
{
	std::vector< Edge * >::iterator edgeIt = edges.begin();
	std::vector< Vertex * >::iterator vertexIt = vertices.begin();
	bool duplicate = 0;
	while(edgeIt != edges.end() && vertexIt != vertices.end())
	{
		bool duplicateVertex = (bool)((*edgeIt)->edgeConnectivity[0] != (*edgeIt)->edgeConnectivity[1]);
		duplicate = 0;
		
		if((*edgeIt)->label == Vessel)
		{
			std::vector< Edge * >::iterator checkIt = edgeIt;
			std::vector< Vertex * >::iterator checkIt2 = vertexIt;
			++checkIt;
			++checkIt2;
			if(duplicateVertex)
				++checkIt2;
			
			while(checkIt != edges.end() && checkIt2 != vertices.end())
			{
				if((*edgeIt)->label == Vessel)
				{
					double * tmp1 = (*edgeIt)->edgePointCoordinates.front();
					double * tmp2 = (*checkIt)->edgePointCoordinates.front();
					
					if(std::abs(tmp1[2] - tmp2[2]) > 1.0)
					{
						++checkIt;
						++checkIt2;
						if(duplicateVertex)
							++checkIt2;
						continue;
					}
					
					float dist = std::sqrt((tmp1[0] - tmp2[0])*(tmp1[0] - tmp2[0]) + (tmp1[1] - tmp2[1])*(tmp1[1] - tmp2[1]) + (tmp1[2] - tmp2[2])*(tmp1[2] - tmp2[2]));
					if(dist < (*edgeIt)->radius || dist < (*checkIt)->radius)
					{
						if((*edgeIt)->radius < (*checkIt)->radius)
						{
							duplicate = 1;
							edgeIt = edges.erase(edgeIt);
							vertexIt = vertices.erase(vertexIt);
							if(duplicateVertex)
								vertexIt = vertices.erase(vertexIt);
							
							std::vector< Edge * >::iterator updateIt;
							for(updateIt = edgeIt; updateIt != edges.end(); ++updateIt)
							{
								if(duplicateVertex)
								{
									(*updateIt)->edgeConnectivity[0] -= 2;
									(*updateIt)->edgeConnectivity[1] -= 2;
								}
								else
								{
									(*updateIt)->edgeConnectivity[0] -= 1;
									(*updateIt)->edgeConnectivity[1] -= 1;
								}
							}
							++checkIt;
							++checkIt2;
							if(duplicateVertex)
								++checkIt2;
						}
						else
						{
							duplicate = 0;
							checkIt = edges.erase(checkIt);
							checkIt2 = vertices.erase(checkIt2);
							if(duplicateVertex)
								checkIt2 = vertices.erase(checkIt2);
							
							std::vector< Edge * >::iterator updateIt;
							for(updateIt = checkIt; updateIt != edges.end(); ++updateIt)
							{
								if(duplicateVertex)
								{
									(*updateIt)->edgeConnectivity[0] -= 2;
									(*updateIt)->edgeConnectivity[1] -= 2;
								}
								else
								{
									(*updateIt)->edgeConnectivity[0] -= 1;
									(*updateIt)->edgeConnectivity[1] -= 1;
								}
							}
						}
					}
					else
					{
						++checkIt;
						++checkIt2;
						if(duplicateVertex)
							++checkIt2;
					}
				}
				else
				{
					++checkIt;
					++checkIt2;
					if(duplicateVertex)
						++checkIt2;
				}
			}
		}
		
		if(!duplicate)
		{
			++edgeIt;
			++vertexIt;
			if(duplicateVertex)
				++vertexIt;
		}
	}
};

// calculate total segment length from edge startID to soma
// NOTE: Only computes euclidean distance between Nodes/Vertex Points
double AmiraSpatialGraph::totalSegmentLength(int startID)
{
	if(startID >= edges.size() || startID < 0)
		return 0;
// 	double totalLength = edges[startID]->segmentLength();
	double totalLength = 0;
	double thisVec[3], thisNorm = 0;
	for(int ii = 0; ii < 3; ++ii)
	{
		thisVec[ii] = edges[startID]->edgePointCoordinates.back()[ii] - edges[startID]->edgePointCoordinates.front()[ii];
		thisNorm += thisVec[ii]*thisVec[ii];
	}
	totalLength = sqrt(thisNorm);

	int fatherID = edges[startID]->fatherID;
// 	int segmentCnt = 1;
	while(fatherID != -1)
	{
		int currID = fatherID;
// 		totalLength += edges[currID]->segmentLength();
		thisNorm = 0;
		for(int ii = 0; ii < 3; ++ii)
		{
			thisVec[ii] = edges[currID]->edgePointCoordinates.back()[ii] - edges[currID]->edgePointCoordinates.front()[ii];
			thisNorm += thisVec[ii]*thisVec[ii];
		}
		totalLength += sqrt(thisNorm);
		fatherID = edges[currID]->fatherID;
// 		++segmentCnt;
	}
// 	std::cout << "nr of segments to soma = " << segmentCnt << std::endl;
	return totalLength;
};

// calculate total segment length from edge startID to soma
// NOTE: Consider points on edges. More precise than totalSegmentLength
double AmiraSpatialGraph::totalSegmentLengthPrecise(int startID)
{
	if(startID >= edges.size() || startID < 0)
		return 0;
	double totalLength = edges[startID]->segmentLength();

	int fatherID = edges[startID]->fatherID;
	while(fatherID != -1)
	{
		int currID = fatherID;
		totalLength += edges[currID]->segmentLength();
		fatherID = edges[currID]->fatherID;
	}
	return totalLength;
};

double AmiraSpatialGraph::cumulatedSegmentAngle(int startID)
{
	if(startID >= edges.size() || startID < 0)
		return 0;
	double totalAngle = 0;
	double thisVec[3], thisNorm = 0;
	for(int ii = 0; ii < 3; ++ii)
	{
		thisVec[ii] = edges[startID]->edgePointCoordinates.back()[ii] - edges[startID]->edgePointCoordinates.front()[ii];
		thisNorm += thisVec[ii]*thisVec[ii];
	}
	thisNorm = sqrt(thisNorm);
	if(thisNorm)
		for(int ii = 0; ii < 3; ++ii)
			thisVec[ii] /= thisNorm;
	int fatherID = edges[startID]->fatherID;
	int segmentCnt = 1;
	while(fatherID != -1)
	{
		++segmentCnt;
		int currID = fatherID;
		double nextVec[3], nextNorm = 0;
		for(int ii = 0; ii < 3; ++ii)
		{
			nextVec[ii] = edges[currID]->edgePointCoordinates.back()[ii] - edges[currID]->edgePointCoordinates.front()[ii];
			nextNorm += nextVec[ii]*nextVec[ii];
		}
		nextNorm = sqrt(nextNorm);
		if(nextNorm)
			for(int ii = 0; ii < 3; ++ii)
				nextVec[ii] /= nextNorm;
		double angle = 0;
		for(int ii = 0; ii < 3; ++ii)
		{
			angle += thisVec[ii]*nextVec[ii];
			thisVec[ii] = nextVec[ii];
		}
		angle = acos(angle);
		totalAngle += angle;
		fatherID = edges[currID]->fatherID;
	}
	return totalAngle;
// 	return totalAngle/(double)segmentCnt;
};

bool AmiraSpatialGraph::isLabelInSpatialGraph(int checkLabel)
{
	for(int ii = 0; ii < edges.size(); ++ii)
		if(edges[ii]->label == checkLabel)
		{
			return 1;
		}
	return 0;
};

//extract all points of landmark 'label' planewise; return their plane indices in zIndexList; return false if empty
bool AmiraSpatialGraph::extractLandmark(int label, std::list< std::list< double * > >& planewisePointList, std::list< int >& zIndexList)
{
	std::vector< Edge * >::iterator edgeIt;
	edgeIt = this->edges.begin();
	if(edgeIt != this->edges.end())
	{
		int zIndex = 0;
		while(edgeIt != this->edges.end() && (*edgeIt)->label != label)
			++edgeIt;
		
		zIndex = lround((*edgeIt)->edgePointCoordinates.back()[2]);
		zIndexList.push_back(zIndex);
		std::list< double * > planeCoordinates;
		std::list< double * >::iterator pointListIt;
		while(edgeIt != this->edges.end())
		{
			if((*edgeIt)->label != label)
			{
				++edgeIt;
				continue;
			}
			int tmpZ = lround((*edgeIt)->edgePointCoordinates.back()[2]);
			// 			std::flush(std::cout << "in plane " << tmpZ << std::endl);
			if(tmpZ == zIndex)
			{
				for(pointListIt = (*edgeIt)->edgePointCoordinates.begin(); pointListIt != (*edgeIt)->edgePointCoordinates.end(); ++pointListIt)
				{
					double * tmpPoint = *pointListIt;
// 					tmpPoint[2] = (double)(int)(tmpPoint[2] + 0.5);
					tmpPoint[2] = round(tmpPoint[2]);
					planeCoordinates.push_back(*pointListIt);
					planeCoordinates.back()[2] = tmpPoint[2];
				}
			}
			else
			{
				std::list< double * > tmpList(planeCoordinates);
				planewisePointList.push_back(tmpList);
				planeCoordinates.clear();
				
				for(pointListIt = (*edgeIt)->edgePointCoordinates.begin(); pointListIt != (*edgeIt)->edgePointCoordinates.end(); ++pointListIt)
				{
					double * tmpPoint = *pointListIt;
// 					tmpPoint[2] = (double)(int)(tmpPoint[2] + 0.5);
					tmpPoint[2] = round(tmpPoint[2]);
					planeCoordinates.push_back(*pointListIt);
					planeCoordinates.back()[2] = tmpPoint[2];
				}
				zIndex = tmpZ;
				zIndexList.push_back(zIndex);
			}
			++edgeIt;
		}
		zIndexList.sort();
		planewisePointList.push_back(planeCoordinates);
		std::list< std::list< double * > >::const_iterator controlit;
		for(controlit = planewisePointList.begin(); controlit != planewisePointList.end(); ++controlit)
			if(controlit->size())
				return 1;
		
		return 0;
	}
	
	return 0;
};

//extract all points of landmark 'label' as PolyData; return false if empty
bool AmiraSpatialGraph::extractLandmark(int label, PolyDataPointerType polyData)
{
	if(!polyData)
		polyData = PolyDataPointerType::New();
	polyData->Allocate(1);
	PointsPointerType points = PointsPointerType::New();
	points->SetDataTypeToFloat();
	int lastID = 0;
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = this->edges.begin(); edgeIt != this->edges.end(); ++edgeIt)
	{
		if((*edgeIt)->label != label)
			continue;
		
		int end = (*edgeIt)->edgePointCoordinates.size();
		if(label >= Landmark && label < ZAxis)	// vtkPolygon does NOT use the same point twice on a contour as a SpatialGraph does
			--end;
		PolygonPointerType poly = PolygonPointerType::New();
		poly->GetPointIds()->SetNumberOfIds(end);
		
		std::list< double * >::iterator pointListIt;
		pointListIt = (*edgeIt)->edgePointCoordinates.begin();
		for(int ii = 0; ii < end; ++pointListIt, ++ii)
		{
			double * tmp = new double[3];
			tmp[0] = (*pointListIt)[0], tmp[1] = (*pointListIt)[1], tmp[2] = (*pointListIt)[2];
			points->InsertNextPoint(tmp);
			poly->GetPointIds()->SetId(ii, ii + lastID);
		}
		lastID += end;
		polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
	}
	polyData->SetPoints(points);
	polyData->Update();
	
	if(polyData->GetNumberOfPoints())
		return 1;
	
	return 0;
};

//extract all points of landmark 'label' as PolyData; return their plane indices in zIndexList; return false if empty
bool AmiraSpatialGraph::extractLandmark(int label, PolyDataPointerType polyData, std::list< int >& zIndexList)
{
	std::vector< Edge * >::iterator edgeIt;
	edgeIt = this->edges.begin();
	if(edgeIt != this->edges.end())
	{
		polyData->Allocate(1);
		PointsPointerType points = PointsPointerType::New();
		points->SetDataTypeToFloat();
		int lastID = 0;
		while(edgeIt != this->edges.end() && (*edgeIt)->label != label)
			++edgeIt;
		
		std::list< double * >::iterator pointListIt;
		while(edgeIt != this->edges.end())
		{
			if((*edgeIt)->label != label)
			{
				++edgeIt;
				continue;
			}
			int end = (*edgeIt)->edgePointCoordinates.size();
			if(label >= Landmark)	// vtkPolygon does NOT use the same point twice on a contour as a SpatialGraph does
				--end;
			PolygonPointerType poly = PolygonPointerType::New();
			poly->GetPointIds()->SetNumberOfIds(end);
			pointListIt = (*edgeIt)->edgePointCoordinates.begin();
			for(int ii = 0; ii < end; ++pointListIt, ++ii)
			{
				double * tmp = *pointListIt;
				points->InsertNextPoint(tmp);
				poly->GetPointIds()->SetId(ii, ii + lastID);
			}
			lastID += end;
			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
			++edgeIt;
		}
		polyData->SetPoints(points);
		polyData->Update();
		
		for(int ii = 0; ii < polyData->GetNumberOfCells(); ++ii)
		{
			double * bounds = polyData->GetCell(ii)->GetBounds();
			if(lround(bounds[4]) == lround(bounds[5]))
				zIndexList.push_back(lround(bounds[4]));
		}
		zIndexList.sort();
		zIndexList.unique();
		
		if(polyData->GetNumberOfPoints())
			return 1;
			
		return 0;
	}
	
	return 0;
};

//writes all cells in object as closed loops in SpatialGraph
void AmiraSpatialGraph::addPolyDataObject(PolyDataPointerType object, int label)
{
	if(object->GetNumberOfCells())
	{
		for(int ii = 0; ii < object->GetNumberOfCells(); ++ii)
		{
			PointsPointerType cellPoints = object->GetCell(ii)->GetPoints();
			if(cellPoints->GetNumberOfPoints())
			{
				double * start = new double[3];
				double * end = new double[3];
				cellPoints->GetPoint(0, start);
				int endIndex = (label >= Landmark) ? 0 : cellPoints->GetNumberOfPoints()-1;
				cellPoints->GetPoint(0, end);
				std::list< double * > edgePts;
				for(int jj = 0; jj < cellPoints->GetNumberOfPoints(); ++jj)
				{
					double * pt = new double[3];
					cellPoints->GetPoint(jj, pt);
					edgePts.push_back(pt);
				}
				if(!endIndex)
					edgePts.push_back(end);
				int connectivity[2];
				int nrEdgePts = edgePts.size();
				connectivity[0] = getNumberOfVertices();
				connectivity[1] = getNumberOfVertices() + 1;
				Vertex * startVertex = new Vertex(start, label);
				Vertex * endVertex = new Vertex(end, label);
				Edge * newEdge = new Edge(connectivity, nrEdgePts, label, edgePts);
				addVertex(startVertex);
				addVertex(endVertex);
				addEdge(newEdge);
			}
		}
	}
};

void AmiraSpatialGraph::addLine(double start[3], double end[3])
{
	double * endPoint = new double[3];
	double * bottomPoint = new double[3];
	for(int ii = 0; ii < 3; ++ii)
	{
		endPoint[ii] = end[ii];
		bottomPoint[ii] = start[ii];
	}
	Vertex * newVert1 = new Vertex(endPoint, ZAxis);
	Vertex * newVert2 = new Vertex(bottomPoint, ZAxis);
	this->addVertex(newVert2);
	this->addVertex(newVert1);
	int connectionIndex[2];
	if(!this->getNumberOfVertices())
	{
		connectionIndex[0] = 0;
		connectionIndex[1] = 1;
	}
	else
	{
		connectionIndex[0] = this->getNumberOfVertices() - 2;
		connectionIndex[1] = this->getNumberOfVertices() - 1;
	}
	std::list< double * > axisCoords;
	axisCoords.push_back(bottomPoint);
	axisCoords.push_back(endPoint);
	Edge * newAxis = new Edge(connectionIndex, 2, ZAxis, axisCoords);
	this->addEdge(newAxis);
};

void AmiraSpatialGraph::addLine(double start[3], double end[3], int ID)
{
	double * endPoint = new double[3];
	double * bottomPoint = new double[3];
	for(int ii = 0; ii < 3; ++ii)
	{
		endPoint[ii] = end[ii];
		bottomPoint[ii] = start[ii];
	}
	Vertex * newVert1 = new Vertex(endPoint, ID);
	Vertex * newVert2 = new Vertex(bottomPoint, ID);
	this->addVertex(newVert2);
	this->addVertex(newVert1);
	int connectionIndex[2];
	if(!this->getNumberOfVertices())
	{
		connectionIndex[0] = 0;
		connectionIndex[1] = 1;
	}
	else
	{
		connectionIndex[0] = this->getNumberOfVertices() - 2;
		connectionIndex[1] = this->getNumberOfVertices() - 1;
	}
	std::list< double * > axisCoords;
	axisCoords.push_back(bottomPoint);
	axisCoords.push_back(endPoint);
	Edge * newAxis = new Edge(connectionIndex, 2, ID, axisCoords);
	this->addEdge(newAxis);
};

//careful! impossible to undo
// Removes all landmark labels
// only keeps labels between Neuron and below Landmark
// keeps Dendrite, ApicalDendrite, BasalDendrite, Axon and Soma
void AmiraSpatialGraph::removeLandmarkLabels()
{
	std::map< int, int > vertexIDLUT;
	std::vector< Vertex * > newVertices;
	std::vector< Edge * > newEdges;
	int newID = 0;
	for(int ii = 0; ii < this->vertices.size(); ++ii)
	{
		if(this->vertices[ii]->label < Neuron || this->vertices[ii]->label>= Landmark)
		{
			continue;
		}
		else
		{
			Vertex * newVertex = new Vertex(this->vertices[ii]);
			newVertices.push_back(newVertex);
			vertexIDLUT.insert(std::pair< int, int >(ii, newID));
			++newID;
		}
	}
	for(int ii = 0; ii < this->edges.size(); ++ii)
	{
		if(this->edges[ii]->label < Neuron || this->edges[ii]->label>= Landmark)
		{
			continue;
		}
		else
		{
			Edge * newEdge = new Edge(this->edges[ii]);
			int oldID0 = newEdge->edgeConnectivity[0];
			int oldID1 = newEdge->edgeConnectivity[1];
			newEdge->edgeConnectivity[0] = vertexIDLUT[oldID0];
			newEdge->edgeConnectivity[1] = vertexIDLUT[oldID1];
			newEdges.push_back(newEdge);
		}
	}
	
	std::vector< Vertex * >::iterator vertexIt;
	for(vertexIt = this->vertices.begin(); vertexIt != this->vertices.end(); ++vertexIt)
	{
		delete *vertexIt;
	}
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = this->edges.begin(); edgeIt != this->edges.end(); ++edgeIt)
	{
		delete *edgeIt;
	}
	this->vertices.clear();
	this->vertices = newVertices;
	this->edges.clear();
	this->edges = newEdges;
};

//careful! impossible to undo
void AmiraSpatialGraph::removeLabel(int label)
{
	std::map< int, int > vertexIDLUT;
	std::vector< Vertex * > newVertices;
	std::vector< Edge * > newEdges;
	int newID = 0;
	for(int ii = 0; ii < this->vertices.size(); ++ii)
	{
		if(this->vertices[ii]->label == label)
		{
			continue;
		}
		else
		{
			Vertex * newVertex = new Vertex(this->vertices[ii]);
			newVertices.push_back(newVertex);
			vertexIDLUT.insert(std::pair< int, int >(ii, newID));
			++newID;
		}
	}
	for(int ii = 0; ii < this->edges.size(); ++ii)
	{
		if(this->edges[ii]->label == label)
		{
			continue;
		}
		else
		{
			Edge * newEdge = new Edge(this->edges[ii]);
			int oldID0 = newEdge->edgeConnectivity[0];
			int oldID1 = newEdge->edgeConnectivity[1];
			newEdge->edgeConnectivity[0] = vertexIDLUT[oldID0];
			newEdge->edgeConnectivity[1] = vertexIDLUT[oldID1];
			newEdges.push_back(newEdge);
		}
	}

	std::vector< Vertex * >::iterator vertexIt;
	for(vertexIt = this->vertices.begin(); vertexIt != this->vertices.end(); ++vertexIt)
	{
		delete *vertexIt;
	}
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = this->edges.begin(); edgeIt != this->edges.end(); ++edgeIt)
	{
		delete *edgeIt;
	}
	this->vertices.clear();
	this->vertices = newVertices;
	this->edges.clear();
	this->edges = newEdges;

// 	std::vector< Vertex * >::iterator vertexIt;
// 	vertexIt = this->vertices.begin();
// 	while(vertexIt != this->vertices.end())
// 	{
// 		if((*vertexIt)->label == label)
// 		{
// 			delete *vertexIt;
// 			vertexIt = this->vertices.erase(vertexIt);
// 		}
// 		else
// 			++vertexIt;
// 	}
// 	std::vector< Edge * >::iterator edgeIt;
// 	edgeIt = this->edges.begin();
// 	while(edgeIt != this->edges.end())
// 	{
// 		if((*edgeIt)->label == label)
// 		{
// 			delete *edgeIt;
// 			edgeIt = this->edges.erase(edgeIt);
// 		}
// 		else
// 			++edgeIt;
// 	}
};

/* Sets the father ID in case it is not set previously (for SpatialGraphs from Amira format)
 * 	fatherID is set if SpatialGraph is loaded from .HOC. Otherwise fatherID == -1 for all edges
 * 	otherwise fatherID==-1 if father node is soma
 * Use EdgeConnectivity to set fatherIDs properly
 * Tested for one example. Also works for .HOC files, fatherIDs stay the same.
 */
void AmiraSpatialGraph::setFatherID()
{
	// Find IDs of Edges that are labeled as soma
	std::vector<int> somaEdgeIDs; // For amira SpatialGraph this can be several short edges, thus use a vector

	// Check whether FatherID has been set already, i.e. -1 for all edges
	bool fatherIDOK = true;

	// EdgeConnectivity -> EdgeID
	// EdgeConnectivity[0] -> VertexID of father edge
	// EdgeConnectivity[1] -> VertexID of current edge
	// to get fatherID:
	// 		edgeConnectivity1toEdgeID[ edgeIDtoEdgeConnectivity0[ edgeID ] ]
	int edgeID = 0;
	std::map<int,int> edgeIDtoEdgeConnectivity0;
	std::map<int,int> edgeConnectivity1toEdgeID;

	for(std::vector< Edge * >::iterator edgeIt = this->edgesBegin(); edgeIt != this->edgesEnd(); ++edgeIt)
	{
		// Check and find soma label
		if ((*edgeIt)->label == Soma)
		{
			somaEdgeIDs.push_back(edgeID);
		}

		// Check fatherID
		if (-1 != (*edgeIt)->fatherID)
		{
			fatherIDOK = false;
		}

		// edgeIDtoEdgeConnectivity0
		edgeIDtoEdgeConnectivity0.insert(std::pair<int,int>(edgeID,(*edgeIt)->edgeConnectivity[0]));
		// edgeConnectivity1toEdgeID
		if (edgeConnectivity1toEdgeID.count((*edgeIt)->edgeConnectivity[1])>0)
		{
			std::cout << "ERROR (setFatherID)! EdgeConnectivity[1] was already added to map!" << std::endl;
			return;
		}
		else
		{
			edgeConnectivity1toEdgeID.insert(std::pair<int,int>((*edgeIt)->edgeConnectivity[1],edgeID));
		}
		edgeID++;
	}

	if (somaEdgeIDs.empty())
	{
		std::cout << "ERROR (setFatherID)! No Soma Found! Cannot set fatherID!" << std::endl;
		return;
	}

	if (!fatherIDOK)
	{
		std::cout << "WARNING! FatherID has already been set previously (not all of them are -1)!" << std::endl;
		std::cout << "  FatherIDs are set again and might be modified!" << std::endl;
		std::cout << "  Proceed with caution!" << std::endl;
	}

	// Set FatherID
	edgeID = 0;
	for(std::vector< Edge * >::iterator edgeIt = this->edgesBegin(); edgeIt != this->edgesEnd(); ++edgeIt)
	{
		if ((*edgeIt)->label != Soma)
		{
			// fatherID = edgeConnectivity1toEdgeID[ edgeIDtoEdgeConnectivity0[ edgeID ] ]
			if (edgeIDtoEdgeConnectivity0.count(edgeID)==0)
			{
				std::cout << "ERROR (setFatherID)! edgeID not found in Map!" << std::endl;
				return;
			}

			int edgeConnectivity0_tmp = edgeIDtoEdgeConnectivity0[edgeID];

			if (edgeConnectivity1toEdgeID.count(edgeConnectivity0_tmp)==0)
			{
				std::cout << "ERROR (setFatherID)! edgeConnectivity[1] not found in Map!" << std::endl;
				return;
			}

			int edgeID_tmp = edgeConnectivity1toEdgeID[edgeConnectivity0_tmp];

			// Check if EdgeID is ID of Soma, then set to -1
			std::vector<int>::iterator it = find (somaEdgeIDs.begin(), somaEdgeIDs.end(), edgeID_tmp);
			if (it != somaEdgeIDs.end())
			{
				(*edgeIt)->fatherID = -1;
			}
			else
			{
				(*edgeIt)->fatherID = edgeID_tmp;
			}
		}
		else
		{
			(*edgeIt)->fatherID = -1;
		}
		edgeID++;
	}
};

void AmiraSpatialGraph::truncateSpatialGraph(double bounds[2], int cut_COORD)
{
	std::vector< Edge * > newEdges;
	std::vector< int > VertexIDtoDel;
	std::vector< int > VertexIDEnd;
	int counter = 0;

	// Go through all edges
	for(int ii = 0; ii < this->edges.size(); ++ii)
	{
		bool edgeOut = false;
		bool modifyEdge = false;
		//flush(std::cout << "--- Edge " << ii << "/" << this->edges.size() << " ---" << std::endl);

		// If connected vertex is on deletion list, if so, delete whole edge as well
		int vertexID1 = this->edges[ii]->edgeConnectivity[0];
		if (std::find(VertexIDEnd.begin(), VertexIDEnd.end(), vertexID1) != VertexIDEnd.end())
		{
			//flush(std::cout << "   Edge is daughter of cut Edge"<< std::endl);
			edgeOut = true;
		}
		else
		{	// If first point of Edge is already outside of slicing box, delete edge
			std::list< double >::iterator edgeRadiusListIt = this->edges[ii]->radiusList.begin();
			//flush(std::cout << "Size of Radius List " << this->edges[ii]->radiusList.size() << std::endl);
			double previousRad = *edgeRadiusListIt;
			std::list< double * >::iterator edgeListIt = this->edges[ii]->edgePointCoordinates.begin();
			double * previousPt = *edgeListIt;

			if (previousPt[cut_COORD]<bounds[0] || previousPt[cut_COORD]>bounds[1])
			{
				edgeOut = true;
			}
			else
			{ 	// If first point of Edge is not outside of slicing box, check whether remaining points are outside
				// if outside -> modify edge
				// if inside -> keep edge
				std::list< double * > newEdgePointCoordinates;
				std::list< double > newRadList;

				++edgeRadiusListIt;
				++edgeListIt;
				double currentRad;
				double * currentPt;
				int jj = 0;

				for(; edgeListIt != this->edges[ii]->edgePointCoordinates.end(); ++edgeListIt, ++edgeRadiusListIt, ++jj)
				{
					currentPt = *edgeListIt;
					currentRad = *edgeRadiusListIt;

					// Add radius and EdgePoint List to new Edge
					newEdgePointCoordinates.push_back(previousPt);
					newRadList.push_back(previousRad);

					if (currentPt[cut_COORD]<bounds[0] || currentPt[cut_COORD]>bounds[1])
					{
						modifyEdge = true;
						break;
					}
					else
					{
						previousPt = currentPt;
						previousRad = currentRad;
					}
					//flush(std::cout << "Pt0: " << previousPt[0] << "," << previousPt[1] << "," << previousPt[2] << std::endl);
					//flush(std::cout << "Pt1: " << currentPt[0] << "," << currentPt[1] << "," << currentPt[2] << std::endl);
				}

				// If current point is outside of slicing box
				if (modifyEdge)
				{
					//flush(std::cout << "Modify Edge!" << std::endl);
					double coord = 0.0;

					if (currentPt[cut_COORD]<bounds[0])
					{
						coord = bounds[0];
					}
					else if (currentPt[cut_COORD]>bounds[1])
					{
						coord = bounds[1];
					}
					else
					{
						std::cout << "Something is wrong! Point not outside of bounding box" << std::endl;
						return;
					}

					// Compute intersecting Point via standing vector plus weighted direction vector
					double Dst1 = previousPt[cut_COORD]-coord;
					double Dst2 = currentPt[cut_COORD]-coord;
					double t = 0.0;

					if ((Dst1*Dst2) >= 0.0)
					{
						t = 0.0;
					}
					else if (Dst1 == Dst2)
					{
						t = 0.0;
					}
					else
					{
						t = (-Dst1/(Dst2-Dst1 /*+ 1E-8 */));
					}
					if (t<0.0 || t>1.0)
					{
						std::cout << "Something is wrong! t is equal or below 0.0 or above 1.0!" << std::endl;
						std::cout << "Dist1: " << Dst1 << " Dist2: " << Dst2 << " t: " << t << std::endl;
						std::cout << "coord: " << coord << " previousPt: " << previousPt[cut_COORD] << " currentPt: " << currentPt[cut_COORD] << std::endl;
						t = 0.0;
						return;
					}

					// Assign new radius
					double newRad = previousRad + (currentRad - previousRad) * t;

					// New Point
					double * newPt = new double[3];
					for (int i = 0; i<3; i++)
					{
						newPt[i] = previousPt[i] + (currentPt[i]-previousPt[i]) * t;
						this->vertices[this->edges[ii]->edgeConnectivity[1]]->coordinates[i] = newPt[i];
					}
					VertexIDEnd.push_back(this->edges[ii]->edgeConnectivity[1]);

					//flush(std::cout << "--- Interpolate Edge " << ii << "/" << this->edges.size() << " newEdgeID " << counter << std::endl);
					//flush(std::cout << "    Previous: " << previousPt[0] << " " << previousPt[1] << " " << previousPt[2] << " r = " << previousRad << std::endl);
					//flush(std::cout << "    Current: " << currentPt[0] << " " << currentPt[1] << " " << currentPt[2] << " r = " << currentRad << std::endl);
					//flush(std::cout << "    NewPt: " << newPt[0] << " " << newPt[1] << " " << newPt[2] << " r = " << newRad << std::endl);
					//flush(std::cout << this->edges[ii]->segmentLength() << " " << this->edges[ii]->numEdgePoints << std::endl);

					// Add new point and new radius to list of EdgePoints and Radius
					newEdgePointCoordinates.push_back(newPt);
					newRadList.push_back(newRad);

					// Create new modified Edge
					int newEdgeConnectivity[2];
					newEdgeConnectivity[0] = this->edges[ii]->edgeConnectivity[0];
					newEdgeConnectivity[1] = this->edges[ii]->edgeConnectivity[1];
					int newNumEdgePoints = newEdgePointCoordinates.size();

					//for(std::list< double * >::iterator ptIt = newEdgePointCoordinates.begin(); ptIt != newEdgePointCoordinates.end(); ++ptIt)
					//{
					//	flush(std::cout << "newEdgePts " << (*ptIt)[0] << " " << (*ptIt)[1] << " " << (*ptIt)[2] << std::endl);
					//}

					Edge * newEdge = new Edge(newEdgeConnectivity, newNumEdgePoints, this->edges[ii]->label, newEdgePointCoordinates, newRadList);
					newEdges.push_back(newEdge);
				}
				else
				{
					//flush(std::cout << "--- Keep Edge " << ii << "/" << this->edges.size() << " newEdgeID " << counter << std::endl);
					Edge * newEdge = new Edge(this->edges[ii]);
					newEdges.push_back(newEdge);
				}
				counter++;
			}
		}

		// Change Vertex point! (don't delete this point)
		// if Edge deleted, add vertex id
		if(edgeOut)
		{
			//flush(std::cout << "--- Delete Edge " << ii << "/" << this->edges.size() << std::endl);
			int vertexID2 = this->edges[ii]->edgeConnectivity[1];
			VertexIDtoDel.push_back(vertexID2);
			VertexIDEnd.push_back(vertexID2);
		}
	}

	// Delete corresponding Vertices
	std::map< int, int > vertexIDLUT;
	std::vector< Vertex * > newVertices;
	int newID = 0;
	for(int ii = 0; ii < this->vertices.size(); ++ii)
	{
		if (std::find(VertexIDtoDel.begin(), VertexIDtoDel.end(), ii) != VertexIDtoDel.end())
		{
			continue;
		}
		else
		{
			Vertex * newVertex = new Vertex(this->vertices[ii]);
			newVertices.push_back(newVertex);
			vertexIDLUT.insert(std::pair< int, int >(ii, newID));
			++newID;
		}
	}

	//std::cout << "VERTEX SIZE " << newVertices.size() << " Delete: " << VertexIDtoDel.size() << "/" << this->vertices.size() << std::endl;
	// Update EdgeConnectivity
	int counter2 = 0;
	for(std::vector< Edge * >::iterator edgeIt = newEdges.begin(); edgeIt != newEdges.end(); ++edgeIt, ++counter2)
	{
		// Update Vertex ID
		int oldID0 = (*edgeIt)->edgeConnectivity[0];
		int oldID1 = (*edgeIt)->edgeConnectivity[1];
		(*edgeIt)->edgeConnectivity[0] = vertexIDLUT[oldID0];
		(*edgeIt)->edgeConnectivity[1] = vertexIDLUT[oldID1];

		// Check whether Vertices and EdgePoints match
		// If not, display Error Message
		double * ptE1 = (*edgeIt)->edgePointCoordinates.front();
		double * ptV1 = newVertices[(*edgeIt)->edgeConnectivity[0]]->coordinates;
		double diff1 = sqrt(pow((ptE1[X_COORD]-ptV1[X_COORD]),2.0)+pow((ptE1[Y_COORD]-ptV1[Y_COORD]),2.0)+pow((ptE1[Z_COORD]-ptV1[Z_COORD]),2.0));
		if (diff1>0)
		{
			flush(std::cout << "Error: Vertex Coordinates do not equal Point Coordinates (begin) Bounds: [" << bounds[0] << " " << bounds[1] << "]" << std::endl);
			flush(std::cout << "  EdgeConnectivity: [" << (*edgeIt)->edgeConnectivity[0] << " " << (*edgeIt)->edgeConnectivity[1] << "] " << std::endl);
			flush(std::cout << "  EdgeID: " << counter2 << " NumEdges: " << (*edgeIt)->numEdgePoints << " EdgeLabel: " << (*edgeIt)->label << std::endl);
			flush(std::cout << "  EdgePointBegin: " << ptE1[0] << " " << ptE1[1] << " " << ptE1[2] << std::endl);
			flush(std::cout << "  VertexPoint: " << ptV1[0] << " " << ptV1[1] << " " << ptV1[2] << " VertexLabel: " << newVertices[(*edgeIt)->edgeConnectivity[0]]->label << std::endl);
		}
		double * ptE2 = (*edgeIt)->edgePointCoordinates.back();
		double * ptV2 = newVertices[(*edgeIt)->edgeConnectivity[1]]->coordinates;
		double diff2 = sqrt(pow((ptE2[X_COORD]-ptV2[X_COORD]),2.0)+pow((ptE2[Y_COORD]-ptV2[Y_COORD]),2.0)+pow((ptE2[Z_COORD]-ptV2[Z_COORD]),2.0));
		if (diff2>0)
		{
			flush(std::cout << "Error: Vertex Coordinates do not equal Point Coordinates (end) Bounds: [" << bounds[0] << " " << bounds[1] << "]" << std::endl);
			flush(std::cout << "  EdgeConnectivity: [" << (*edgeIt)->edgeConnectivity[0] << " " << (*edgeIt)->edgeConnectivity[1] << "] " << std::endl);
			flush(std::cout << "  EdgeID: " << counter2 << " NumEdges: " << (*edgeIt)->numEdgePoints << " EdgeLabel: " << (*edgeIt)->label << std::endl);
			flush(std::cout << "  EdgePointEnd: " << ptE2[0] << " " << ptE2[1] << " " << ptE2[2] << std::endl);
			flush(std::cout << "  VertexPoint: " << ptV2[0] << " " << ptV2[1] << " " << ptV2[2] << " VertexLabel: " << newVertices[(*edgeIt)->edgeConnectivity[1]]->label << std::endl);
		}
	}

	// Delete Vertices and Edges that are outside of slicing box
	for(std::vector< Vertex * >::iterator vertexIt = this->vertices.begin(); vertexIt != this->vertices.end(); ++vertexIt)
	{
		delete *vertexIt;
	}
	this->vertices.clear();
	this->vertices = newVertices;

	for(std::vector< Edge * >::iterator edgeIt = this->edges.begin(); edgeIt != this->edges.end(); ++edgeIt)
	{
		delete *edgeIt;
	}
	this->edges.clear();
	this->edges = newEdges;

	// Display
//	std::cout << "NumVertices: " << this->vertices.size() << " NumEdges: " << this->edges.size() << std::endl;
//	for(int ii = 0; ii < this->vertices.size(); ++ii)
//	{
//		flush(std::cout << "VertexID " << ii << ": " << this->vertices[ii]->coordinates[0] << " " << this->vertices[ii]->coordinates[1] << " " << this->vertices[ii]->coordinates[2] << " " << std::endl);
//	}
//
//	for(int ii = 0; ii < this->edges.size(); ++ii)
//	{
//		//std::list< double * >::iterator edgeListIt1 = this->edges[ii]->edgePointCoordinates.begin();
//		double * pt1 = this->edges[ii]->edgePointCoordinates.front();
//		//std::list< double * >::iterator edgeListIt2 = this->edges[ii]->edgePointCoordinates.end();
//		double * pt2 = this->edges[ii]->edgePointCoordinates.back();
//
//		flush(std::cout << "EdgeID " << ii << ": EdgeConnectivity: [" << this->edges[ii]->edgeConnectivity[0] << " " << this->edges[ii]->edgeConnectivity[1] << "] " << std::endl);
//		flush(std::cout << "EdgeID " << ii << ": NumEdges: " << this->edges[ii]->numEdgePoints << std::endl);
//		flush(std::cout << "EdgeID " << ii << ": Pt1: " << pt1[0] << " " << pt1[1] << " " << pt1[2] << std::endl);
//		flush(std::cout << "EdgeID " << ii << ": Pt2: " << pt2[0] << " " << pt2[1] << " " << pt2[2] << std::endl);
//		int jjPt = 0;
//		for (std::list< double * >::iterator edgeListIt = this->edges[ii]->edgePointCoordinates.begin(); edgeListIt != this->edges[ii]->edgePointCoordinates.end(); ++edgeListIt, ++jjPt)
//		{
//			double * pt1 = *edgeListIt;
//			flush(std::cout << "EdgeID " << ii << ": Pt" << jjPt << ": " << pt1[0] << " " << pt1[1] << " " << pt1[2] << std::endl);
//		}
//	}
};

Column::Column()
{
	contours = NULL;
	this->top = new double[3];
	this->bottom =  new double[3];
};

Column::Column ( Column* otherColumn )
{
	this->contours = PolyDataPointerType::New();
	this->contours->DeepCopy(otherColumn->contours);
	this->top = new double[3];
	this->top[0] = otherColumn->top[0], this->top[1] = otherColumn->top[1], this->top[2] = otherColumn->top[2];
	this->bottom =  new double[3];
	this->bottom[0] = otherColumn->bottom[0], this->bottom[1] = otherColumn->bottom[1], this->bottom[2] = otherColumn->bottom[2];
}

Column::Column(PolyDataPointerType contours, double * top, double * bottom)
{
	this->contours = contours;
	this->top = new double[3];
	this->top[0] = top[0], this->top[1] = top[1], this->top[2] = top[2];
	this->bottom =  new double[3];
	this->bottom[0] = bottom[0], this->bottom[1] = bottom[1], this->bottom[2] = bottom[2];
};

Column::~Column()
{
	if(this->top)
		delete top;
	if(this->bottom)
		delete bottom;
};

void Column::getCenter ( double center[3] )
{
	center[0] = 0.5*(top[0] + bottom[0]);
	center[1] = 0.5*(top[1] + bottom[1]);
	center[2] = 0.5*(top[2] + bottom[2]);
}

void Column::translateColumn(const double * shift)
{
	for(int ii = 0; ii < 3; ++ii)
	{
		this->top[ii] += shift[ii];
		this->bottom[ii] += shift[ii];
	}
	
	if(this->contours)
	{
		TransformFilterType transform = TransformFilterType::New();
		TransformPointerType translation = TransformPointerType::New();
		translation->Translate(shift);
		transform->SetTransform(translation);
		transform->SetInput(this->contours);
		transform->Update();
		this->contours->DeepCopy(transform->GetOutput());
		this->contours->Update();
	}
}

void Column::rotateColumn(gsl_matrix * rot)
{
	TransformFilterType transform = TransformFilterType::New();
	TransformPointerType rotation = TransformPointerType::New();
	HomogeneousMatrixPointerType mat = HomogeneousMatrixPointerType::New();
	for(int ii = 0; ii < 3; ++ii)
		for(int jj = 0; jj < 3; ++jj)
			mat->SetElement(ii, jj, gsl_matrix_get(rot, ii, jj));
	for(int ii = 0; ii < 3; ++ii)
	{
		mat->SetElement(ii, 3, 0);
		mat->SetElement(3, ii, 0);
	}
	mat->SetElement(3, 3, 1);
	
	double hTop[4], hBottom[4];
	for(int ii = 0; ii < 3; ++ii)
	{
		hTop[ii] = this->top[ii];
		hBottom[ii] = this->bottom[ii];
	}
	hTop[3] = hBottom[3] = 1;
	mat->MultiplyPoint(hTop, hTop);
	mat->MultiplyPoint(hBottom, hBottom);
	for(int ii = 0; ii < 3; ++ii)
	{
		this->top[ii] = hTop[ii];
		this->bottom[ii] = hBottom[ii];
	}
	
	if(this->contours)
	{
		rotation->SetMatrix(mat);
		transform->SetTransform(rotation);
		transform->SetInput(this->contours);
		transform->Update();
		this->contours->DeepCopy(transform->GetOutput());
		this->contours->Update();
	}
};

void Column::rotateColumn(HomogeneousMatrixPointerType mat)
{
	TransformFilterType transform = TransformFilterType::New();
	TransformPointerType rotation = TransformPointerType::New();
	double hTop[4], hBottom[4];
	for(int ii = 0; ii < 3; ++ii)
	{
		hTop[ii] = this->top[ii];
		hBottom[ii] = this->bottom[ii];
	}
	hTop[3] = hBottom[3] = 1;
	mat->MultiplyPoint(hTop, hTop);
	mat->MultiplyPoint(hBottom, hBottom);
	for(int ii = 0; ii < 3; ++ii)
	{
		this->top[ii] = hTop[ii];
		this->bottom[ii] = hBottom[ii];
	}
	
	if(this->contours)
	{
		rotation->SetMatrix(mat);
		transform->SetTransform(rotation);
		transform->SetInput(this->contours);
		transform->Update();
		this->contours->DeepCopy(transform->GetOutput());
		this->contours->Update();
	}
};

Surface::Surface(PolyDataPointerType mesh)
{
	dataValid = 0;
	intersectionFound = 0;
	surfaceMesh = PolyDataPointerType::New();
	surfaceMesh->Allocate(1);
	surfaceMesh->DeepCopy(mesh);
	locator = CellLocatorPointerType::New();
	locator->AutomaticOn();
	locator->CacheCellBoundsOn();
	locator->SetDataSet(surfaceMesh);
	locator->BuildLocator();
};

void Surface::intersectLine(double * axis, double * center)
{
	double * pt = new double[3];
	double a0[3], a1[3], tol = 1E-03, t, pcoords[3];
	int subId;
	vtkIdType cellID;
	GenericCellPointerType intersectCell = GenericCellPointerType::New();
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] = center[jj];
		a1[jj] = center[jj];
	}
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] += axis[jj]*4000;
		a1[jj] -= axis[jj]*4000;
	}
	int intersection = locator->IntersectWithLine(a0, a1, tol, t, pt, pcoords, subId, cellID, intersectCell);
	dataValid = 1;
	if(intersection)
	{
		intersectPt = pt;
		intersectID = cellID;
		intersectionFound = 1;
	}
	else
	{
		delete [] pt;
		intersectPt = NULL;
		intersectID = -1;
		intersectionFound = 0;
	}
};

void Surface::intersectLineInDirection ( double* axis, double* center )
{
	double * pt = new double[3];
	double a0[3], a1[3], tol = 1E-03, t, pcoords[3];
	int subId;
	vtkIdType cellID;
	GenericCellPointerType intersectCell = GenericCellPointerType::New();
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] = center[jj];
		a1[jj] = center[jj];
	}
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] += axis[jj]*4000;
	}
	int intersection = locator->IntersectWithLine(a0, a1, tol, t, pt, pcoords, subId, cellID, intersectCell);
	dataValid = 1;
	if(intersection)
	{
		intersectPt = pt;
		intersectID = cellID;
		intersectionFound = 1;
	}
	else
	{
		delete [] pt;
		intersectPt = NULL;
		intersectID = -1;
		intersectionFound = 0;
	}
}

double * Surface::getLastIntersectPoint()
{
	if(dataValid)
	{
		if(intersectPt)
		{
			double * pt = new double[3];
			pt[0] = intersectPt[0], pt[1] = intersectPt[1], pt[2] = intersectPt[2];
			return pt;
		}
		return NULL;
	}
	else
		return NULL;
};

void Surface::getLastIntersectPoint(double pt[3])
{
	if(dataValid && intersectPt)
		pt[0] = intersectPt[0], pt[1] = intersectPt[1], pt[2] = intersectPt[2];
};

vtkIdType Surface::getLastIntersectCellID()
{
	if(dataValid)
	{
		return intersectID;
	}
	return 0;
};

ClosedSurface::ClosedSurface ( PolyDataPointerType mesh ) : Surface ( mesh )
{
	insideSurfaceFilter = SelectEnclosedPointsFilterType::New();
	insideSurfaceFilter->Initialize(mesh);
	insideSurfaceFilter->CheckSurfaceOn();
// 	insideSurfaceFilter->Print(std::cout);
}

ClosedSurface::~ClosedSurface()
{
	if(insideSurfaceFilter) insideSurfaceFilter->Complete();
}

// use vtkSelectEnclosedPoints filter
bool ClosedSurface::isPointInsideSurface ( double pt[3] )
{
	return insideSurfaceFilter->IsInsideSurface(pt);
}

ConnectionMatrix::ConnectionMatrix()
{
	simpleMatrix = false; // set simpleMatrix to false by default!
}

ConnectionMatrix::~ConnectionMatrix()
{
// TBI
}

SelectionType ConnectionMatrix::getPreColumnCelltypeSelection(unsigned int column, unsigned int cellType)
{
	SelectionType intersection;
	SelectionType colIDs = preColumnIDs[column];
	SelectionType celltypeIDs = preTypeIDs[cellType];
	std::set_intersection(colIDs.begin(), colIDs.end(), celltypeIDs.begin(), celltypeIDs.end(), std::back_inserter(intersection));
	
	return intersection;
}

SelectionType ConnectionMatrix::getPostColumnCelltypeSelection(unsigned int column, unsigned int cellType)
{
	SelectionType intersection;
	SelectionType colIDs = postColumnIDs[column];
	SelectionType celltypeIDs = postTypeIDs[cellType];
	std::set_intersection(colIDs.begin(), colIDs.end(), celltypeIDs.begin(), celltypeIDs.end(), std::back_inserter(intersection));
	
	return intersection;
}

SelectionType ConnectionMatrix::getUniqueSelection(SelectionType preSelection)
{
	if (!simpleMatrix)
	{
		return preSelection;
	}

	SelectionType preSelectionUnique;

	for(int i = 0; i < preSelection.size(); ++i)
	{
		std::map<unsigned int, unsigned int>::iterator it = preDuplicatedIDtoSingleID.find(preSelection[i]);
		if (it == preDuplicatedIDtoSingleID.end())
		{
			std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
			std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
			std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
			return preSelection;
		}
		else
		{
			if (std::find(preSelectionUnique.begin(), preSelectionUnique.end(),(*it).second) != preSelectionUnique.end())
				continue;
			else
				preSelectionUnique.push_back((*it).second);
		}
	}

	return preSelectionUnique;
}

double ConnectionMatrix::getInnervationSum(SelectionType preSelection, SelectionType postSelection)
{	// Sum of Innervation
	double innervation = 0;
	std::map<unsigned int, unsigned int>::iterator it;

	for(int i = 0; i < preSelection.size(); ++i)
	{
		if (simpleMatrix)
		{
			it = preDuplicatedIDtoSingleID.find(preSelection[i]);
			if (it == preDuplicatedIDtoSingleID.end())
			{
				std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
				std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
				std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
				return -1;
			}
		}

		for(int j = 0; j < postSelection.size(); ++j)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				innervation = innervation + matrixIt->second;
			}
		}

	}

	return innervation;
	// To get average per Presynapse: I_sum / selectionPre.size();
	// To get average per Postsynapse: I_sum / selectionPost.size();
	// To get average per Pair: I_sum / (selectionPre.size() * selectionPost.size() );
}

double ConnectionMatrix::getAverageInnervation(SelectionType preSelection, SelectionType postSelection)
{
	double innervation = ConnectionMatrix::getInnervationSum(preSelection, postSelection);
	innervation = innervation / (preSelection.size()*postSelection.size());
	return innervation;
}

double ConnectionMatrix::getAverageConvergenceIterateMatrix(SelectionType preSelection, SelectionType postSelection)
{
	double convergence = 0;
	std::map<unsigned int, unsigned int> PreCellIDCounter;

	if (simpleMatrix)
	{
		for(int i = 0; i < preSelection.size(); ++i)
		{
			std::map<unsigned int, unsigned int>::iterator it = preDuplicatedIDtoSingleID.find(preSelection[i]);
			if (it == preDuplicatedIDtoSingleID.end())
			{
				std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
				std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
				std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
			}
			else
			{
				if (PreCellIDCounter.count(it->second)>0)
				{
					PreCellIDCounter[it->second]++;
				}
				else
				{
					PreCellIDCounter.insert(std::pair< unsigned int, unsigned int >(it->second,1));
				}
			}
		}
	}

	for (std::map< MatrixIndexType, float >::iterator it = matrix.begin(); it != matrix.end(); ++it)
	{
		unsigned int preCellID = it->first.first;
		unsigned int postCellID = it->first.second;

		if (std::find(postSelection.begin(), postSelection.end(), postCellID) != postSelection.end())
		{
			if (simpleMatrix)
			{
				if (PreCellIDCounter.count(preCellID)>0)
				{
					float innervation = it->second;
					convergence = convergence + (1 - exp(-1*innervation))*PreCellIDCounter[preCellID];
				}
			}
			else if (std::find(preSelection.begin(), preSelection.end(), preCellID) != preSelection.end())
			{
				float innervation = it->second;
				convergence = convergence + (1 - exp(-1*innervation));
			}
		}
	}

	convergence = convergence / (preSelection.size()*postSelection.size());
	return convergence;
}

double ConnectionMatrix::getAverageConvergenceOptimized(SelectionType preSelection, SelectionType postSelection)
{
	unsigned int loopSelection = preSelection.size()*postSelection.size();
	unsigned int loopMatrix = matrix.size() + preSelection.size();

	if (loopMatrix<loopSelection)
	{
		return getAverageConvergenceIterateMatrix(preSelection, postSelection);
	}
	else
	{
		 /*if (loopSelection>(6000^2))
		 {
			 std::cout << "WARNING! THIS HAS CAUSED ISSUES BEFORE! NOW FIXED!" << std::endl;
		 }*/
		// might not have worked
		// getAverageConvergenceIterate(SelectionpreSelection, postSelection);

		return getAverageConvergence(preSelection, postSelection);
	}
}

double ConnectionMatrix::getAverageConvergenceIterateSelection(SelectionType preSelection, SelectionType postSelection)
{
	return getAverageConvergence(preSelection, postSelection);

	/* HERES A BUG, PARRALLIZATION DID NOT WORK, therefore don't use parralllization
	// Parallized!
	double convergence = 0;
	std::map<unsigned int, unsigned int>::iterator it;

	for(int i = 0; i < preSelection.size(); ++i)
	{
		if (simpleMatrix)
		{
			it = preDuplicatedIDtoSingleID.find(preSelection[i]);
			if (it == preDuplicatedIDtoSingleID.end())
			{
				std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
				std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
				std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
				return -1;
			}
		}

		//#pragma omp parallel for
		//HERES A BUG, PARRALLIZATION DID NOT WORK, therefore don't use parralllization
		for(int j = 0; j < postSelection.size(); ++j)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				convergence = convergence + (1 - exp(-1*matrixIt->second));
			}
		}
	}

	convergence = convergence / (preSelection.size()*postSelection.size());
	return convergence;*/
}

double ConnectionMatrix::getAverageConvergence(SelectionType preSelection, SelectionType postSelection)
{
	double convergence = 0;
	std::map<unsigned int, unsigned int>::iterator it;

	for(int i = 0; i < preSelection.size(); ++i)
	{
		if (simpleMatrix)
		{
			it = preDuplicatedIDtoSingleID.find(preSelection[i]);
			if (it == preDuplicatedIDtoSingleID.end())
			{
				std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
				std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
				std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
				return -1;
			}
		}

		for(int j = 0; j < postSelection.size(); ++j)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				convergence = convergence + (1 - exp(-1*matrixIt->second));
			}
		}
	}

	convergence = convergence / (preSelection.size()*postSelection.size());
	return convergence;
}

std::list<double> ConnectionMatrix::getInnervationValuesPerPair(SelectionType preSelection, SelectionType postSelection)
{
	std::list<double> innervationValues;
	std::map<unsigned int, unsigned int>::iterator it;

	for(int i = 0; i < preSelection.size(); ++i)
	{
		if (simpleMatrix)
		{
			it = preDuplicatedIDtoSingleID.find(preSelection[i]);
			if (it == preDuplicatedIDtoSingleID.end())
			{
				std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
				std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
				std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
				return innervationValues;
			}
		}

		for(int j = 0; j < postSelection.size(); ++j)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				innervationValues.push_back(matrixIt->second);
			}
			else
			{
				innervationValues.push_back(0);
			}
		}
	}

	return innervationValues;
}

std::list<double> ConnectionMatrix::getInnervationValuesPerPostsynapse(SelectionType preSelection, SelectionType postSelection)
{
	std::list<double> innervationValues;

	for(int j = 0; j < postSelection.size(); ++j)
	{
		double I = 0;

		for(int i = 0; i < preSelection.size(); ++i)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				std::map<unsigned int, unsigned int>::iterator it = preDuplicatedIDtoSingleID.find(preSelection[i]);
				if (it == preDuplicatedIDtoSingleID.end())
				{
					std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
					std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
					std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
					return innervationValues;
				}
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				I = I + matrixIt->second;
			}
		}

		innervationValues.push_back(I);
	}

	return innervationValues;
}

std::list<double> ConnectionMatrix::getConvergenceValuesPerPostsynapse(SelectionType preSelection, SelectionType postSelection)
{
	std::list<double> convergenceValues;

	for(int j = 0; j < postSelection.size(); ++j)
	{
		double convergence = 0;

		for(int i = 0; i < preSelection.size(); ++i)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				std::map<unsigned int, unsigned int>::iterator it = preDuplicatedIDtoSingleID.find(preSelection[i]);
				if (it == preDuplicatedIDtoSingleID.end())
				{
					std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
					std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
					std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
					return convergenceValues;
				}
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				convergence = convergence + (1 - exp(-1*matrixIt->second));
			}
		}

		convergence = convergence/preSelection.size();
		convergenceValues.push_back(convergence);
	}

	return convergenceValues;
}

std::list<double> ConnectionMatrix::getInnervationValuesPerPresynapse(SelectionType preSelection, SelectionType postSelection)
{
	std::list<double> innervationValues;
	std::map<unsigned int, unsigned int>::iterator it;

	for(int i = 0; i < preSelection.size(); ++i)
	{
		if (simpleMatrix)
		{
			it = preDuplicatedIDtoSingleID.find(preSelection[i]);
			if (it == preDuplicatedIDtoSingleID.end())
			{
				std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
				std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
				std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
				return innervationValues;
			}
		}

		double innervation = 0;

		for(int j = 0; j < postSelection.size(); ++j)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				innervation = innervation + matrixIt->second;
			}
		}
		innervationValues.push_back(innervation);
	}
	return innervationValues;
}

std::list<double> ConnectionMatrix::getInnervationValuesUniquePerPresynapse(SelectionType preSelection, SelectionType postSelection)
{
	std::list<double> innervationValues;
	std::map<unsigned int, unsigned int>::iterator it;
	std::vector<unsigned int> singleIDlist;

	for(int i = 0; i < preSelection.size(); ++i)
	{
		if (simpleMatrix)
		{
			it = preDuplicatedIDtoSingleID.find(preSelection[i]);
			if (it == preDuplicatedIDtoSingleID.end())
			{
				std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
				std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
				std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
				return innervationValues;
			}

			// Check if ID has already been used for computation, if so skip
			if (std::find(singleIDlist.begin(), singleIDlist.end(), (*it).second) != singleIDlist.end())
				continue;
			else
				singleIDlist.push_back((*it).second);
		}

		double innervation = 0;

		for(int j = 0; j < postSelection.size(); ++j)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				innervation = innervation + matrixIt->second;
			}
		}
		innervationValues.push_back(innervation);
	}
	return innervationValues;
}

std::list<double> ConnectionMatrix::getConvergenceValuesPerPresynapse(SelectionType preSelection, SelectionType postSelection)
{
	std::list<double> convergenceValues;
	std::map<unsigned int, unsigned int>::iterator it;

	for(int i = 0; i < preSelection.size(); ++i)
	{
		if (simpleMatrix)
		{
			it = preDuplicatedIDtoSingleID.find(preSelection[i]);
			if (it == preDuplicatedIDtoSingleID.end())
			{
				std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
				std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
				std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
				return convergenceValues;
			}
		}

		double convergence = 0;

		for(int j = 0; j < postSelection.size(); ++j)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				convergence = convergence + (1 - exp(-1*matrixIt->second));
			}
		}

		convergence = convergence/postSelection.size();
		convergenceValues.push_back(convergence);
	}
	return convergenceValues;
}

std::list<double> ConnectionMatrix::getConvergenceValuesUniquePerPresynapse(SelectionType preSelection, SelectionType postSelection)
{
	std::list<double> convergenceValues;

	std::map<unsigned int, unsigned int>::iterator it;
	std::vector<unsigned int> singleIDlist;

	for(int i = 0; i < preSelection.size(); ++i)
	{
		if (simpleMatrix)
		{
			it = preDuplicatedIDtoSingleID.find(preSelection[i]);
			if (it == preDuplicatedIDtoSingleID.end())
			{
				std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
				std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
				std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
				return convergenceValues;
			}

			// Check if ID has already been used for computation, if so skip
			if (std::find(singleIDlist.begin(), singleIDlist.end(), (*it).second) != singleIDlist.end())
				continue;
			else
				singleIDlist.push_back((*it).second);
		}

		double convergence = 0;

		for(int j = 0; j < postSelection.size(); ++j)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				convergence = convergence + (1 - exp(-1*matrixIt->second));
			}
		}
		convergence = convergence/postSelection.size();
		convergenceValues.push_back(convergence);
	}

	return convergenceValues;
}

std::list<double> ConnectionMatrix::getInnervationValuesUniquePerPair(SelectionType preSelection, SelectionType postSelection)
{
	std::list<double> innervationValues;
	std::map<unsigned int, unsigned int>::iterator it;
	std::vector<unsigned int> singleIDlist;

	for(int i = 0; i < preSelection.size(); ++i)
	{
		if (simpleMatrix)
		{
			it = preDuplicatedIDtoSingleID.find(preSelection[i]);
			if (it == preDuplicatedIDtoSingleID.end())
			{
				std::cout << "ERROR! PreCellID " << preSelection[i] << " was not found!" << std::endl;
				std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[i]].first;
				std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[i]].second << std::endl;
				return innervationValues;
			}

			// Check if ID has already been used for computation, if so skip
			if (std::find(singleIDlist.begin(), singleIDlist.end(), (*it).second) != singleIDlist.end())
				continue;
			else
				singleIDlist.push_back((*it).second);
		}

		for(int j = 0; j < postSelection.size(); ++j)
		{
			MatrixIndexType index;

			if (simpleMatrix)
			{
				index = std::make_pair((*it).second, postSelection[j]);
			}
			else
			{
				index = std::make_pair(preSelection[i], postSelection[j]);
			}

			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);

			if(matrixIt != matrix.end())
			{
				innervationValues.push_back(matrixIt->second);
			}
			else
			{
				innervationValues.push_back(0);
			}
		}
	}

	return innervationValues;
}

std::list<double> ConnectionMatrix::getConvergenceValuesPerPair(SelectionType preSelection, SelectionType postSelection)
{
	std::list<double> innervationValues = getInnervationValuesPerPair(preSelection, postSelection);
	std::list<double> convergenceValues;

	for (std::list<double>::iterator it=innervationValues.begin(); it != innervationValues.end(); ++it)
	{
		convergenceValues.push_back(1 - exp(-1*(*it)));
	}
	return convergenceValues;
}

std::list<double> ConnectionMatrix::getConvergenceValuesUniquePerPair(SelectionType preSelection, SelectionType postSelection)
{
	std::list<double> innervationValues = getInnervationValuesUniquePerPair(preSelection, postSelection);
	std::list<double> convergenceValues;

	for (std::list<double>::iterator it=innervationValues.begin(); it != innervationValues.end(); ++it)
	{
		convergenceValues.push_back(1 - exp(-1*(*it)));
	}
	return convergenceValues;
}

// bool CP: true: connection probability; false: innervation value
void ConnectionMatrix::writeSimpleConnectionMatrixAsImage(const char * outputFilename,
															SelectionType preSelection, SelectionType postSelection,
															bool CP, unsigned int downSamplingFactor)
{
	if (!simpleMatrix)
	{
		std::cout << "ERROR! writeSimpleConnectionMatrixAsImage only works for simpleMatrix!" << std::endl;
		return;
	}

	CalcImage2DType::RegionType matrixRegion;
	CalcImage2DType::SizeType matrixSize;
	matrixSize[0] = postSelection.size();
	matrixSize[1] = preSelection.size();
	matrixRegion.SetSize(matrixSize);

	CalcImage2DType::Pointer matrixImage = CalcImage2DType::New();
	matrixImage->SetRegions(matrixRegion);
	matrixImage->Allocate();
	matrixImage->FillBuffer(0);

	float minValue = 1e6;
	float maxValue = 0.0;

	std::map<unsigned int, unsigned int>::iterator it;

	//std::flush(std::cout << "CREATED IMAGE" << std::endl);

	for(int ii = 0; ii < preSelection.size(); ++ii)
	{
		it = preDuplicatedIDtoSingleID.find(preSelection[ii]);

		if (it == preDuplicatedIDtoSingleID.end())
		{
			std::cout << "ERROR! PreCellID " << preSelection[ii] << " was not found!" << std::endl;
			std::cout << "	ColumnID: " << IDColumnCelltypeMap[preSelection[ii]].first;
			std::cout << "	CellTypeID: " << IDColumnCelltypeMap[preSelection[ii]].second << std::endl;
			return;
		}

		for(int jj = 0; jj < postSelection.size(); ++jj)
		{

			MatrixIndexType index = std::make_pair((*it).second, postSelection[jj]);
			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);
			if(matrixIt != matrix.end())
			{
				CalcImage2DType::IndexType imageIndex;
				imageIndex[0] = jj;
				imageIndex[1] = ii;

				float value = matrixIt->second;
				if (CP)
					value = 1-exp(-value);

				matrixImage->SetPixel(imageIndex, value);
				if(value < minValue)
				{
					minValue = value;
				}
				if(value > maxValue)
				{
					maxValue = value;
				}
			}
			else
			{
				minValue = 0.0;
			}
		}
	}

	//std::flush(std::cout << "DONE FILLING IMAGE" << std::endl);

	Image2DType::RegionType scaledRegion;
	Image2DType::SizeType scaledSize;
	scaledSize[0] = postSelection.size();
	scaledSize[1] = preSelection.size();
	scaledRegion.SetSize(scaledSize);

	Image2DType::Pointer scaledMatrixImage = Image2DType::New();
	scaledMatrixImage->SetRegions(scaledRegion);
	scaledMatrixImage->Allocate();
	CalcImage2DToImage2DFilterType::Pointer rescaledMatrixFilter = CalcImage2DToImage2DFilterType::New();
	rescaledMatrixFilter->SetOutputMinimum(0);
	rescaledMatrixFilter->SetOutputMaximum(255);
	rescaledMatrixFilter->SetInput(matrixImage);
	rescaledMatrixFilter->Update();
	scaledMatrixImage = rescaledMatrixFilter->GetOutput();
	scaledMatrixImage->Update();

	Image2DWriterType::Pointer planeWriter = Image2DWriterType::New();
	planeWriter->SetFileName(outputFilename);

	if (downSamplingFactor>1)
	{
		// Rescale / Downsample
		ShrinkImageFilterType::Pointer shrinkFilter = ShrinkImageFilterType::New();
		shrinkFilter->SetInput(scaledMatrixImage);
		shrinkFilter->SetShrinkFactor(0, downSamplingFactor); // shrink the first dimension by a factor of 2 (i.e. 100 gets changed to 50)
		shrinkFilter->SetShrinkFactor(1, downSamplingFactor); // shrink the second dimension by a factor of 2 (i.e. 100 gets changed to 50)
		shrinkFilter->Update();
		planeWriter->SetInput(shrinkFilter->GetOutput());
	}
	else
	{
		planeWriter->SetInput(scaledMatrixImage);
	}

	try
	{
		planeWriter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "WriterExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}

	std::string imageScaleName(outputFilename);
	imageScaleName += "_scale.txt";
	std::ofstream imageScaleFile(imageScaleName.c_str());
	imageScaleFile << "Minimum value = " << minValue << std::endl;
	imageScaleFile << "Maximum value = " << maxValue << std::endl;
	imageScaleFile.close();
}

void ConnectionMatrix::writeConnectionMatrixAsImage(const char* outputFilename, SelectionType preSelection, SelectionType postSelection)
{
	if (simpleMatrix)
	{
		std::cout << "ERROR! function writeConnectionMatrixAsImage does not support simpleMatrix!" << std::endl;
		std::cout << "       Use function writeSimpleConnectionMatrixAsImage instead!" << std::endl;
		return;
	}

	CalcImage2DType::RegionType matrixRegion;
	CalcImage2DType::SizeType matrixSize;
	matrixSize[0] = postSelection.size();
	matrixSize[1] = preSelection.size();
	matrixRegion.SetSize(matrixSize);
	
	CalcImage2DType::Pointer matrixImage = CalcImage2DType::New();
	matrixImage->SetRegions(matrixRegion);
	matrixImage->Allocate();
	matrixImage->FillBuffer(0);
	
	float minInnervation = 1e6;
	float maxInnervation = 0.0;
	for(int ii = 0; ii < preSelection.size(); ++ii)
	{
		for(int jj = 0; jj < postSelection.size(); ++jj)
		{
			MatrixIndexType index(preSelection[ii], postSelection[jj]);
			std::map< MatrixIndexType, float >::const_iterator matrixIt = matrix.find(index);
			if(matrixIt != matrix.end())
			{
				CalcImage2DType::IndexType imageIndex;
				imageIndex[0] = jj;
				imageIndex[1] = ii;
				matrixImage->SetPixel(imageIndex, matrixIt->second);
				if(matrixIt->second < minInnervation)
				{
					minInnervation = matrixIt->second;
				}
				if(matrixIt->second > maxInnervation)
				{
					maxInnervation = matrixIt->second;
				}
			}
			else
			{
				minInnervation = 0.0;
			}
		}
	}
	
	Image2DType::RegionType scaledRegion;
	Image2DType::SizeType scaledSize;
	scaledSize[0] = postSelection.size();
	scaledSize[1] = preSelection.size();
	scaledRegion.SetSize(scaledSize);
	
	Image2DType::Pointer scaledMatrixImage = Image2DType::New();
	scaledMatrixImage->SetRegions(scaledRegion);
	scaledMatrixImage->Allocate();
	CalcImage2DToImage2DFilterType::Pointer rescaledMatrixFilter = CalcImage2DToImage2DFilterType::New();
	rescaledMatrixFilter->SetOutputMinimum(0);
	rescaledMatrixFilter->SetOutputMaximum(255);
	rescaledMatrixFilter->SetInput(matrixImage);
	rescaledMatrixFilter->Update();
	scaledMatrixImage = rescaledMatrixFilter->GetOutput();
	scaledMatrixImage->Update();
	
	Image2DWriterType::Pointer planeWriter = Image2DWriterType::New();
	planeWriter->SetFileName(outputFilename);
	planeWriter->SetInput(scaledMatrixImage);
	try
	{
		planeWriter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "WriterExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}
	
	std::string imageScaleName(outputFilename);
	imageScaleName += "_scale.txt";
	std::ofstream imageScaleFile(imageScaleName.c_str());
	imageScaleFile << "Minimum innervation = " << minInnervation << std::endl;
	imageScaleFile << "Maximum innervation = " << maxInnervation << std::endl;
	imageScaleFile.close();
}


CellTable::CellTable()
{
	//TBI
}

CellTable::~CellTable()
{
	for(int i = 0; i < rows.size(); ++i)
	{
		delete rows[i];
	}
	rows.clear();
}

std::list< unsigned int > CellTable::getPreColumnCelltypeColumns(std::list< unsigned int > columns, std::list< unsigned int > cellTypes)
{
	std::list< unsigned int > selection;
	std::map< ColumnCellTypePair, unsigned int >::const_iterator headerIt;
	
	if(!columns.size() && !cellTypes.size())
	{
		for(headerIt = header.begin(); headerIt != header.end(); ++headerIt)
		{
			selection.push_back(headerIt->second);
		}
	}
	if(columns.size() && !cellTypes.size())
	{
		for(headerIt = header.begin(); headerIt != header.end(); ++headerIt)
		{
			unsigned int column = headerIt->first.first;
			if(std::find(columns.begin(), columns.end(), column) != columns.end())
			{
				selection.push_back(headerIt->second);
			}
		}
	}
	if(!columns.size() && cellTypes.size())
	{
		for(headerIt = header.begin(); headerIt != header.end(); ++headerIt)
		{
			unsigned int cellType = headerIt->first.second;
			if(std::find(cellTypes.begin(), cellTypes.end(), cellType) != cellTypes.end())
			{
				selection.push_back(headerIt->second);
			}
		}
	}
	if(columns.size() && cellTypes.size())
	{
		for(headerIt = header.begin(); headerIt != header.end(); ++headerIt)
		{
			unsigned int column = headerIt->first.first;
			unsigned int cellType = headerIt->first.second;
			if(std::find(columns.begin(), columns.end(), column) != columns.end() 
				&& std::find(cellTypes.begin(), cellTypes.end(), cellType) != cellTypes.end())
			{
				selection.push_back(headerIt->second);
			}
		}
	}
	
	return selection;
}

std::vector< CellTableRow * > CellTable::getPostColumnCelltypeRows(std::list< unsigned int > columns, std::list< unsigned int > cellTypes)
{
	std::vector< CellTableRow * > selection;
	
	if(!columns.size() && !cellTypes.size())
	{
		selection = rows;
	}
	if(columns.size() && !cellTypes.size())
	{
		for(int i = 0; i < rows.size(); ++i)
		{
			unsigned int column = rows[i]->column;
			if(std::find(columns.begin(), columns.end(), column) != columns.end())
			{
				selection.push_back(rows[i]);
			}
		}
	}
	if(!columns.size() && cellTypes.size())
	{
		for(int i = 0; i < rows.size(); ++i)
		{
			unsigned int cellType = rows[i]->cellType;
			if(std::find(cellTypes.begin(), cellTypes.end(), cellType) != cellTypes.end())
			{
				selection.push_back(rows[i]);
			}
		}
	}
	if(columns.size() && cellTypes.size())
	{
		for(int i = 0; i < rows.size(); ++i)
		{
			unsigned int column = rows[i]->column;
			unsigned int cellType = rows[i]->cellType;
			if(std::find(columns.begin(), columns.end(), column) != columns.end() 
				&& std::find(cellTypes.begin(), cellTypes.end(), cellType) != cellTypes.end())
			{
				selection.push_back(rows[i]);
			}
		}
	}
	
	return selection;
}

void CellTable::mergeLocal(const char* fname)
{
	std::map< ColumnCellTypePair, unsigned int > numCells = getCellNumbers(fname);
	if (numCells.empty())
		return;

	std::vector <int> PreColumnList;
	PreColumnList.push_back(A1);
	PreColumnList.push_back(A2);
	PreColumnList.push_back(A3);
	PreColumnList.push_back(A4);
	PreColumnList.push_back(B1);
	PreColumnList.push_back(B2);
	PreColumnList.push_back(B3);
	PreColumnList.push_back(B4);
	PreColumnList.push_back(C1);
	PreColumnList.push_back(C2);
	PreColumnList.push_back(C3);
	PreColumnList.push_back(C4);
	PreColumnList.push_back(D1);
	PreColumnList.push_back(D2);
	PreColumnList.push_back(D3);
	PreColumnList.push_back(D4);
	PreColumnList.push_back(E1);
	PreColumnList.push_back(E2);
	PreColumnList.push_back(E3);
	PreColumnList.push_back(E4);
	PreColumnList.push_back(Alpha);
	PreColumnList.push_back(Beta);
	PreColumnList.push_back(Gamma);
	PreColumnList.push_back(Delta);

	// Get Indices of SymLocal Types for synapsesPerPreColumn
	std::vector< int > indicesLocal;
	std::vector< int > indicesLocalDel;
	std::vector< int > freqLocal;

	for (int i = 0; i!=PreColumnList.size(); ++i)
	{
		for (std::map< ColumnCellTypePair, unsigned int>::iterator headerit=header.begin(); headerit!=header.end(); ++headerit)
		{
			if (PreColumnList.at(i)==headerit->first.first && headerit->first.second >= SymLocal1axon && headerit->first.second <= SymLocal6axon)
			{
				indicesLocal.push_back(headerit->second);

				// Number of Cells of respective Column/CellType Pair
				unsigned int n = numCells.find(ColumnCellTypePair(headerit->first.first,headerit->first.second))->second;
				freqLocal.push_back(n);

				if (headerit->first.second != SymLocal1axon)
				{
					indicesLocalDel.push_back(headerit->second);
				}
			}
		}
	}

	std::sort (indicesLocalDel.begin(),indicesLocalDel.end());

	// Go through rows (individual cells)
	for (std::vector< CellTableRow * >::iterator row = rows.begin(); row != rows.end(); ++row)
	{
		// Change CellType ID of Postsynaptic Type
		if ((*row)->cellType >= SymLocal1 && (*row)->cellType <= SymLocal6)
		{
			(*row)->cellType = SymLocal;
		}

		// Add up Presynaptic SymLocal Types (all to SymLocal1) in synapsesPerPreTypeColumn
		for (int i=0; i<indicesLocal.size()/6; i++)
		{
			float sum = 0;
			for (int ii=0; ii<6; ii++)
			{
				sum = sum + float(freqLocal.at(i*6+ii));
			}

			// Weight SymLocals according to their occurrence
			float tmp = 0;
			for (int ii=0; ii<6; ii++)
			{
				float p = (*row)->synapsesPerPreTypeColumn.at(indicesLocal.at(i*6+ii));
				float w = (freqLocal.at(i*6+ii))/sum;
				tmp = tmp + w * p;
			}
			(*row)->synapsesPerPreTypeColumn.at(indicesLocal.at(i*6)) = tmp;
		}

		// Delete SymLocal2-6 from synapsesPerPreTypeColumn
		for (int i = indicesLocalDel.size(); i-- > 0 ; )
		{
			(*row)->synapsesPerPreTypeColumn.erase((*row)->synapsesPerPreTypeColumn.begin()+indicesLocalDel.at(i));
		}
	}

	// Create a new header, insert former values, skip SymLocal2-6, rename SymLocal1 to SymLocal, adjust columnID
	std::map< ColumnCellTypePair, unsigned int > header_backup = header;
	header.clear();

	for (std::map< ColumnCellTypePair, unsigned int >::iterator headerit=header_backup.begin(); headerit!=header_backup.end(); ++headerit)
	{
		unsigned int column = headerit->first.first;
		unsigned int cellType = headerit->first.second;
		unsigned int columnID = headerit->second;
		int subtmp = 0;

		// Subtract Indices of SymLocal2-6 from ColumnID
		for (std::vector<int>::iterator it = indicesLocalDel.begin() ; it != indicesLocalDel.end(); ++it)
		{
			if ((columnID) >= (*it))
			{
				subtmp++;
			}
		}
		columnID = columnID-subtmp;

		// Rename SymLocal1 to SymLocal in header
		if (cellType == SymLocal1axon)
		{
			cellType = SymLocalaxon;
		} // Skip this
		else if (cellType >= SymLocal2axon && cellType <= SymLocal6axon)
		{
			continue;
		}

		// Insert Values into header
		header.insert(std::pair< ColumnCellTypePair, unsigned int >(ColumnCellTypePair(column, cellType), columnID));
	}
}

std::map< ColumnCellTypePair, unsigned int > CellTable::getCellNumbers(const char* fname)
{
	std::ifstream inputStream(fname);
	std::map< ColumnCellTypePair, unsigned int > numCells;

	if(!inputStream.fail())
	{
		std::map< std::string, unsigned int> ColumnLabels2Int;
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("Alpha"), Alpha));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("A1"), A1));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("A2"), A2));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("A3"), A3));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("A4"), A4));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("Beta"), Beta));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("B1"), B1));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("B2"), B2));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("B3"), B3));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("B4"), B4));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("Gamma"), Gamma));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C1"), C1));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C2"), C2));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C3"), C3));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C4"), C4));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C5"), C5));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("C6"), C6));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("Delta"), Delta));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D1"), D1));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D2"), D2));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D3"), D3));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D4"), D4));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D5"), D5));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("D6"), D6));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E1"), E1));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E2"), E2));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E3"), E3));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E4"), E4));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E5"), E5));
		ColumnLabels2Int.insert(std::pair< std::string, int >(std::string("E6"), E6));

		std::string currentLine;
		std::vector< int > colVector;

		bool headerFound = false;
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				if (currentLine.find("CELL TYPE/COLUMN TOTAL") != std::string::npos)
				{
					std::getline(inputStream, currentLine);

					std::size_t delim = currentLine.find("\t");

					while(delim != std::string::npos)
					{
						currentLine = currentLine.substr(delim+1);
						delim = currentLine.find("\t");
						colVector.push_back(ColumnLabels2Int[currentLine.substr(0,delim)]);
					}
					headerFound = true;
				}

				int celltypeID = -1;

				if (headerFound)
				{
					if (currentLine.find("SymLocal1axon") != std::string::npos)
					{
						celltypeID = SymLocal1axon;
					}
					else if (currentLine.find("SymLocal2axon") != std::string::npos)
					{
						celltypeID = SymLocal2axon;
					}
					else if (currentLine.find("SymLocal3axon") != std::string::npos)
					{
						celltypeID = SymLocal3axon;
					}
					else if (currentLine.find("SymLocal4axon") != std::string::npos)
					{
						celltypeID = SymLocal4axon;
					}
					else if (currentLine.find("SymLocal5axon") != std::string::npos)
					{
						celltypeID = SymLocal5axon;
					}
					else if (currentLine.find("SymLocal6axon") != std::string::npos)
					{
						celltypeID = SymLocal6axon;
					}
				}

				if (celltypeID != -1)
				{
					std::size_t delim = currentLine.find("\t");

					int i = 0;

					while(delim != std::string::npos)
					{
						if (i>=colVector.size())
							break;

						currentLine = currentLine.substr(delim+1);
						delim = currentLine.find("\t");
						std::string tmp = currentLine.substr(0,delim);
						unsigned int n = atoi(tmp.c_str());
						numCells.insert(std::pair< ColumnCellTypePair, unsigned int >(ColumnCellTypePair(colVector.at(i), celltypeID), n));
						i++;
					}

					if ((i)!=colVector.size())
					{
						flush(std::cout << "Header size (" << i << ") and Column Vector size (" << colVector.size() << ") do not match!" << std::endl);
					}
				}

				if (currentLine.find("CELL TYPE/COLUMN INSIDE COLUMN") != std::string::npos)
				{
					break;
				}
			}
		}
	}
	else
	{
		std::cout << "ERROR! Reading nrCells table " << fname << " failed!" << std::endl;
	}

	return numCells;
}
