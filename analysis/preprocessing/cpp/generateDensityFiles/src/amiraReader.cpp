/****************************************************************************/
/*                                                                          */
/* Program:   					                                            */
/*                                                                          */
/* File:      amiraReader.cpp                                               */
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
/* Date: 		19.08.2018                                                  */
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "typedefs.h"
#include "basics.h"
#include "amiraReader.h"
// #include <boost/config/posix_features.hpp>

// #define DEBUG

const char * letters_const = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char * numbers_const = "0123456789";
const char * signs_const = "+-";
const char * otherChars_const = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
const char * whitespace_const = "\t ";
std::list< int > barrelLabels_const;
std::map< int, const char * > int2Labels_const;
std::list< const char * > hocLabels_const;
std::list< const char * > hocNeuronLabels_const;
std::list< const char * > hocLandmarkLabels_const;
void initializeConstants_const();

Reader::Reader(const char * filename)
{
	inputSpatialGraph = NULL;
	this->inputFilename = filename;
	this->letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	this->numbers = "0123456789";
	this->signs = "+-";
	this->otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
	this->whitespace = "\t ";
	this->initializeConstants();
};

Reader::Reader(const char * filename, const char * outputFilename)
{
	inputSpatialGraph = NULL;
	this->inputFilename = filename;
	this->outputFilename = outputFilename;
	this->letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	this->numbers = "0123456789";
	this->signs = "+-";
	this->otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
	this->whitespace = "\t ";
	this->initializeConstants();
};

Reader::~Reader()
{
	//tbd
};

void Reader::readSpatialGraphFile(bool applyTransform)
{
	std::ifstream inputStream(inputFilename);
	
	if(!inputStream.fail())
	{
// 		std::cout << "Reading SpatialGraph file " << inputFilename << std::endl;
		std::string currentLine;
		unsigned int line = 0;
		
		double ** transformation = new double *[4];
		for(int ii = 0; ii < 4; ++ii)
		{
			transformation[ii] = new double[4];
			for(int jj = 0; jj < 4; ++jj)
			{
				if(ii != jj)
					transformation[ii][jj] = 0;
				else
					transformation[ii][jj] = 1;
			}
		}
		
		bool parameters = 0;
		bool transform = 0;
		unsigned int brackets = 0;
		unsigned int vertexTransformIndex = 1000000, edgeTransformIndex = 1000000;
		unsigned int vertexCoordIndex = 1000000, vertexLabelIndex = 1000000, edgeConnectivityIndex = 1000000, edgePointIndex = 1000000, edgeLabelIndex = 1000000, edgePointCoordIndex = 1000000, edgeRadiusIndex = 1000000;
		unsigned int currentIndex = 0;
		unsigned int vertex = 0, edge = 0, point = 0;
		
		std::vector< Vertex * > inputVertices;
		std::vector< Edge * > inputEdges;
		std::list< double * > tmpVertices;
		std::list< int > tmpVertexLabels;
		std::list< int * > tmpEdgeConnections;
		std::list< int > tmpNoEdgePoints;
		std::list< int > tmpEdgeLabels;
		std::list< double * > edgePoints;
		std::list< double > edgePointRadius;
		
		while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
		{
// 			if(!parameters)
// 			{
// 				std::cout << currentLine << std::endl;
// 				++line;
// 			}
			
			if(currentLine.size())
			{
				if(currentLine.find("@", 0) == 0)
				{
					char * tmp = new char[currentLine.size() - 1];
					currentLine.copy(tmp, currentLine.size() - 1, 1);
					currentIndex = atoi(tmp);
// 					std::cout << "Reading data section " << currentIndex << std::endl;
					delete [] tmp;
					continue;
				}
				
				if(currentIndex == 0)
				{
					std::string::size_type loc = currentLine.find("define", 0);
					if(loc != std::string::npos)
					{
						if(currentLine.find("VERTEX", 7) != std::string::npos)
						{
							char * tmp = new char[currentLine.size() - 14];
							currentLine.copy(tmp, currentLine.size() - 14, 14);
							vertex = atoi(tmp);
// 							std::cout << "vertex = " << vertex << std::endl;
// 							inputVertices.resize(vertex);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("EDGE", 7) != std::string::npos)
						{
							char * tmp = new char[currentLine.size() - 12];
							currentLine.copy(tmp, currentLine.size() - 12, 12);
							edge = atoi(tmp);
// 							std::cout << "edges = " << edge << std::endl;
// 							inputEdges.resize(edge);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("POINT", 7) != std::string::npos)
						{
							char * tmp = new char[currentLine.size() - 13];
							currentLine.copy(tmp, currentLine.size() - 13, 13);
							point = atoi(tmp);
// 							std::cout << "points = " << point << std::endl;
							delete [] tmp;
							continue;
						}
					}
					
					loc = currentLine.find("Parameters", 0);
					if(loc == 0)
					{
						parameters = 1;
						if(currentLine.find("{", 0) != std::string::npos)
							brackets = 1;
						continue;
					}
// 					if(parameters && currentLine.find("{", 0) == std::string::npos && currentLine.find("}", 0) == std::string::npos)
// 						continue;
					if(parameters)
					{
						std::string::size_type startPos = 0;
						for(std::string::size_type bPos = currentLine.find("{", startPos); bPos != std::string::npos; )
						{
							++brackets;
							if(bPos == currentLine.size() - 1)
								break;
							bPos = currentLine.find("{", bPos+1);
						}
						for(std::string::size_type bPos = currentLine.find("}", startPos); bPos != std::string::npos; )
						{
							--brackets;
							if(bPos == currentLine.size() - 1)
								break;
							bPos = currentLine.find("}", bPos+1);
						}
						if(!brackets)
							parameters = 0;
					}
// 					if(parameters && currentLine.find("{", 0) != std::string::npos)
// 					{
// 						++brackets;
// 						continue;
// 					}
// 					if(parameters && currentLine.find("}", 0) != std::string::npos)
// 					{
// 						--brackets;
// 						if(!brackets)
// 							parameters = 0;
// 						continue;
// 					}
					
					if(parameters && currentLine.find("TransformationMatrix ", 0) != std::string::npos)
					{
// 						std::cout << "found correct section transform parameters!" << std::endl;
						unsigned int count = 0;
						std::string::size_type loc1, loc2, loc3;
						loc1 = currentLine.find_first_of(numbers, 0);
						loc2 = currentLine.find_first_of(signs, 0);
						if(loc2 != std::string::npos)
							if(loc2 < loc1)
								loc1 = loc2;
						loc2 = currentLine.find_first_of(whitespace, loc1 + 1);	//ignores last value: is always 1 anyways
						while(loc2 != std::string::npos && count < 16)
						{
							char * tmp1 = new char[loc2 - loc1];
							currentLine.copy(tmp1, loc2 - loc1, loc1);
							double ftmp1 = atof(tmp1);
							transformation[count%4][count/4]= ftmp1;	// amira files are columns after each other
							loc3 = loc2;
							loc1 = currentLine.find_first_of(numbers, loc3);
							loc2 = currentLine.find_first_of(signs, loc3);
							if(loc2 != std::string::npos)
								if(loc2 < loc1)
									loc1 = loc2;
							loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
							++count;
							delete [] tmp1;
						}
// 						std::cout << "transformation matrix:" << std::endl;
// 						for(int ii = 0; ii < 4; ++ii)
// 						{
// 							std::cout << "[";
// 							for(int jj = 0; jj < 4; ++jj)
// 							{
// 								if(jj < 3)
// 									std::cout << transformation[ii][jj] << ",\t";
// 								else
// 									std::cout << transformation[ii][jj];
// 							}
// 							std::cout << "]" << std::endl;
// 						}
						//remove numeric artifacts from z-axis:
						for(int ii = 0; ii < 2; ++ii)
						{
							transformation[2][ii] = 0;
							transformation[ii][2] = 0;
						}
						transformation[2][2] = 1;
					}
					loc = currentLine.find("VERTEX", 0);
					if(loc == 0)
					{
						if(currentLine.find("VertexCoordinates", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							vertexCoordIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("GraphLabels", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							vertexLabelIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("TransformInfo", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							vertexTransformIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
					}
					
					loc = currentLine.find("EDGE", 0);
					if(loc == 0)
					{
						if(currentLine.find("EdgeConnectivity", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeConnectivityIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("NumEdgePoints", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgePointIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("GraphLabels", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeLabelIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("TransformInfo", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeTransformIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
					}
					
					loc = currentLine.find("POINT", 0);
					if(loc == 0)
					{
						if(currentLine.find("EdgePointCoordinates", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgePointCoordIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("Radius", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeRadiusIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
					}
				}
				
				if(currentIndex == vertexTransformIndex || currentIndex == edgeTransformIndex)
					continue;
				
				if(currentIndex == vertexCoordIndex)
				{
					const char * thisLine = currentLine.c_str();
					double * tmpCoords = new double[3];
					char ** endptr = new char*;
					tmpCoords[0] = strtod(thisLine, endptr);
					tmpCoords[1] = strtod(*endptr, endptr);
					tmpCoords[2] = strtod(*endptr, endptr);
					
// 					std::string::size_type loc1 = currentLine.find(" ", 0);
// 					std::string::size_type loc2 = currentLine.find(" ", loc1 + 1);
// 					char * tmp1 = new char[loc1];
// 					char * tmp2 = new char[loc2 - loc1];
// 					char * tmp3 = new char[currentLine.size() - loc2 - 1];
// 					currentLine.copy(tmp1, loc1, 0);
// 					currentLine.copy(tmp2, loc2 - loc1, loc1 + 1);
// 					currentLine.copy(tmp3, currentLine.size() - loc2 - 1, loc2 + 1);
// 					double * tmpCoords = new double[3];
// // 					tmpCoords[0] = atof(tmp1);
// // 					tmpCoords[1] = atof(tmp2);
// // 					tmpCoords[2] = atof(tmp3);
// 					char ** endptr = new char*;
// 					tmpCoords[0] = strtod(tmp1, endptr);
// 					tmpCoords[1] = strtod(tmp2, endptr);
// 					tmpCoords[2] = strtod(tmp3, endptr);
					
					tmpVertices.push_back(tmpCoords);
// 					delete [] tmp1, delete [] tmp2, delete [] tmp3;
					delete endptr;
				}
				
				if(currentIndex == vertexLabelIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					int tmplabel = atoi(tmp);
					
					tmpVertexLabels.push_back(tmplabel);
					delete [] tmp;
				}
				
				if(currentIndex == edgeConnectivityIndex)
				{
					const char * thisLine = currentLine.c_str();
					int * tmpCoords = new int[2];
					char ** endptr = new char*;
					tmpCoords[0] = static_cast< int >(strtol(thisLine, endptr, 10));
					tmpCoords[1] = static_cast< int >(strtol(*endptr, endptr, 10));
					
// 					std::string::size_type loc = currentLine.find(" ", 0);
// 					char * tmp1 = new char[loc];
// 					char * tmp2 = new char[currentLine.size() - loc - 1];
// 					currentLine.copy(tmp1, loc, 0);
// 					currentLine.copy(tmp2, currentLine.size() - loc - 1, loc + 1);
// 					int * tmpCoords = new int[2];
// 					tmpCoords[0] = atoi(tmp1);
// 					tmpCoords[1] = atoi(tmp2);
					
					tmpEdgeConnections.push_back(tmpCoords);
// 					delete [] tmp1, delete [] tmp2;
					delete endptr;
				}
				
				if(currentIndex == edgePointIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					int tmplabel = atoi(tmp);
					
					tmpNoEdgePoints.push_back(tmplabel);
					delete [] tmp;
				}
				
				if(currentIndex == edgeLabelIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					int tmplabel = atoi(tmp);
					
					tmpEdgeLabels.push_back(tmplabel);
					delete [] tmp;
				}
				
				if(currentIndex == edgePointCoordIndex)
				{
					const char * thisLine = currentLine.c_str();
					double * tmpCoords = new double[3];
					char ** endptr = new char*;
					tmpCoords[0] = strtod(thisLine, endptr);
					tmpCoords[1] = strtod(*endptr, endptr);
					tmpCoords[2] = strtod(*endptr, endptr);
					
// 					std::string::size_type loc1 = currentLine.find(" ", 0);
// 					std::string::size_type loc2 = currentLine.find(" ", loc1 + 1);
// 					char * tmp1 = new char[loc1];
// 					char * tmp2 = new char[loc2 - loc1];
// 					char * tmp3 = new char[currentLine.size() - loc2 - 1];
// 					currentLine.copy(tmp1, loc1, 0);
// 					currentLine.copy(tmp2, loc2 - loc1, loc1 + 1);
// 					currentLine.copy(tmp3, currentLine.size() - loc2 - 1, loc2 + 1);
// 					double * tmpCoords = new double[3];
// // 					tmpCoords[0] = atof(tmp1);
// // 					tmpCoords[1] = atof(tmp2);
// // 					tmpCoords[2] = atof(tmp3);
// 					char ** endptr = new char*;
// 					tmpCoords[0] = strtod(tmp1, endptr);
// 					tmpCoords[1] = strtod(tmp2, endptr);
// 					tmpCoords[2] = strtod(tmp3, endptr);
					
					edgePoints.push_back(tmpCoords);
// 					delete [] tmp1, delete [] tmp2, delete [] tmp3;
					delete endptr;
				}
				
				if(currentIndex == edgeRadiusIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					double tmpradius = atof(tmp);
					edgePointRadius.push_back(tmpradius);
					delete [] tmp;
				}
			}
		}
		
		std::list< double * >::iterator tmpvertexit;
		std::list< int >::iterator tmpvertexlabelit;
		std::list< int * >::iterator tmpedgeconnectivityit;
		std::list< int >::iterator tmpnumberedgepointsit;
		std::list< int >::iterator tmpedgelabelit;
		std::list< double * >::iterator edgepointit;
		std::list< double >::iterator edgeradiusit;
		
		for(tmpvertexit = tmpVertices.begin(), tmpvertexlabelit = tmpVertexLabels.begin(); tmpvertexit != tmpVertices.end() && tmpvertexlabelit != tmpVertexLabels.end(); ++tmpvertexit, ++tmpvertexlabelit)
		{
			Vertex * tmpVertex = new Vertex(*tmpvertexit, *tmpvertexlabelit);
			inputVertices.push_back(tmpVertex);
		}
		
		edgepointit = edgePoints.begin();
		edgeradiusit = edgePointRadius.begin();
		for(tmpedgeconnectivityit = tmpEdgeConnections.begin(), tmpnumberedgepointsit = tmpNoEdgePoints.begin(), tmpedgelabelit = tmpEdgeLabels.begin();
		tmpedgeconnectivityit != tmpEdgeConnections.end() && tmpnumberedgepointsit != tmpNoEdgePoints.end() && tmpedgelabelit != tmpEdgeLabels.end();
		++tmpedgeconnectivityit, ++tmpnumberedgepointsit, ++tmpedgelabelit)
		{
			std::list< double * >::iterator tmpit = edgepointit;
			std::list< double >::iterator tmpit2 = edgeradiusit;
			for(int ii = 0; ii < *tmpnumberedgepointsit; ++ii)
			{
				++tmpit;
				++tmpit2;
			}
			std::list< double * > tmpPoints(edgepointit, tmpit);
			std::list< double > tmpRadius(edgeradiusit, tmpit2);
			Edge * tmpEdge = new Edge(*tmpedgeconnectivityit, *tmpnumberedgepointsit, *tmpedgelabelit, tmpPoints, tmpRadius);
			
			inputEdges.push_back(tmpEdge);
			edgepointit = tmpit;
			edgeradiusit = tmpit2;
		}
		
		inputSpatialGraph = new AmiraSpatialGraph;
		
		std::vector< Vertex * >::iterator vertexIter;
		std::vector< Edge * >::iterator edgeIter;
		
		for(vertexIter = inputVertices.begin(); vertexIter != inputVertices.end(); ++ vertexIter)
			inputSpatialGraph->addVertex(*vertexIter);
		
		for(edgeIter = inputEdges.begin(); edgeIter != inputEdges.end(); ++edgeIter)
			inputSpatialGraph->addEdge(*edgeIter);
// 		bool isId = 1;
// 		for(int ii = 0; ii < 4; ++ii)
// 			for(int jj = 0; jj < 4; ++jj)
// 			{
// 				if(ii != jj)
// 					if(transformation[ii][jj] != 0)
// 						isId = 0;
// 				
// 				else
// 					if(transformation[ii][jj] != 1)
// 						isId = 0;
// 			}
// 		std::cout << "transformation matrix:" << std::endl;
// 		for(int ii = 0; ii < 4; ++ii)
// 		{
// 			std::cout << "[";
// 			for(int jj = 0; jj < 4; ++jj)
// 			{
// 				if(jj < 3)
// 					std::cout << transformation[ii][jj] << ",\t";
// 				else
// 					std::cout << transformation[ii][jj];
// 			}
// 			std::cout << "]" << std::endl;
// 		}
		if(applyTransform)
		{
// 			inputSpatialGraph->printTransformation();
			inputSpatialGraph->setTransformation(transformation);
			inputSpatialGraph->applyTransformation();
		}
		for(int ii = 0; ii < 4; ++ii)
			delete [] transformation[ii];
		delete [] transformation;
		
// 		std::cout << "Vertex number = " << inputSpatialGraph->getNumberOfVertices() << std::endl;
// 		std::cout << "Edge number = " << inputSpatialGraph->getNumberOfEdges() << std::endl;
// 		std::cout << "Point number = " << inputSpatialGraph->getNumberOfPoints() << std::endl;
		
// 		std::cout << "VertexCoordinates @" << vertexCoordIndex << std::endl;
// 		std::cout << "Vertex GraphLabels @" << vertexLabelIndex << std::endl;
// 		std::cout << "Vertex TransformInfo @" << vertexTransformIndex << std::endl;
// 		std::cout << "EdgeConnectivity @" << edgeConnectivityIndex << std::endl;
// 		std::cout << "NumEdgePoints @" << edgePointIndex << std::endl;
// 		std::cout << "Edge GraphLabels @" << edgeLabelIndex << std::endl;
// 		std::cout << "Edge TransformInfo @" << edgeTransformIndex << std::endl;
// 		std::cout << "EdgePointCoordinates @" << edgePointCoordIndex << std::endl;
// 		std::cout << "Edge Radius @" << edgeRadiusIndex << std::endl;
	}
	
	inputStream.close();
};

void Reader::readSpatialGraphFile(bool applyTransform, const char* fname, AmiraSpatialGraph* sg)
{
	std::ifstream inputStream(fname);
	
	if(!inputStream.fail())
	{
// 		std::cout << "Reading SpatialGraph file " << inputFilename << std::endl;
		std::string currentLine;
		unsigned int line = 0;
		
		double ** transformation = new double *[4];
		for(int ii = 0; ii < 4; ++ii)
		{
			transformation[ii] = new double[4];
			for(int jj = 0; jj < 4; ++jj)
			{
				if(ii != jj)
					transformation[ii][jj] = 0;
				else
					transformation[ii][jj] = 1;
			}
		}
		
		bool parameters = 0;
		bool transform = 0;
		unsigned int brackets = 0;
		unsigned int vertexTransformIndex = 1000000, edgeTransformIndex = 1000000;
		unsigned int vertexCoordIndex = 1000000, vertexLabelIndex = 1000000, edgeConnectivityIndex = 1000000, edgePointIndex = 1000000, edgeLabelIndex = 1000000, edgePointCoordIndex = 1000000, edgeRadiusIndex = 1000000;
		unsigned int currentIndex = 0;
		unsigned int vertex = 0, edge = 0, point = 0;
		
		std::vector< Vertex * > inputVertices;
		std::vector< Edge * > inputEdges;
		std::list< double * > tmpVertices;
		std::list< int > tmpVertexLabels;
		std::list< int * > tmpEdgeConnections;
		std::list< int > tmpNoEdgePoints;
		std::list< int > tmpEdgeLabels;
		std::list< double * > edgePoints;
		std::list< double > edgePointRadius;
		
		while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
		{
// 			if(!parameters)
// 			{
// 				std::cout << currentLine << std::endl;
// 				++line;
// 			}
			
			if(currentLine.size())
			{
				if(currentLine.find("@", 0) == 0)
				{
					char * tmp = new char[currentLine.size() - 1];
					currentLine.copy(tmp, currentLine.size() - 1, 1);
					currentIndex = atoi(tmp);
// 					std::cout << "Reading data section " << currentIndex << std::endl;
					delete [] tmp;
					continue;
				}
				
				if(currentIndex == 0)
				{
					std::string::size_type loc = currentLine.find("define", 0);
					if(loc != std::string::npos)
					{
						if(currentLine.find("VERTEX", 7) != std::string::npos)
						{
							char * tmp = new char[currentLine.size() - 14];
							currentLine.copy(tmp, currentLine.size() - 14, 14);
							vertex = atoi(tmp);
// 							std::cout << "vertex = " << vertex << std::endl;
// 							inputVertices.resize(vertex);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("EDGE", 7) != std::string::npos)
						{
							char * tmp = new char[currentLine.size() - 12];
							currentLine.copy(tmp, currentLine.size() - 12, 12);
							edge = atoi(tmp);
// 							std::cout << "edges = " << edge << std::endl;
// 							inputEdges.resize(edge);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("POINT", 7) != std::string::npos)
						{
							char * tmp = new char[currentLine.size() - 13];
							currentLine.copy(tmp, currentLine.size() - 13, 13);
							point = atoi(tmp);
// 							std::cout << "points = " << point << std::endl;
							delete [] tmp;
							continue;
						}
					}
					
					loc = currentLine.find("Parameters", 0);
					if(loc == 0)
					{
						parameters = 1;
						if(currentLine.find("{", 0) != std::string::npos)
							brackets = 1;
						continue;
					}
// 					if(parameters && currentLine.find("{", 0) == std::string::npos && currentLine.find("}", 0) == std::string::npos)
// 						continue;
					if(parameters)
					{
						std::string::size_type startPos = 0;
						for(std::string::size_type bPos = currentLine.find("{", startPos); bPos != std::string::npos; )
						{
							++brackets;
							if(bPos == currentLine.size() - 1)
								break;
							bPos = currentLine.find("{", bPos+1);
						}
						for(std::string::size_type bPos = currentLine.find("}", startPos); bPos != std::string::npos; )
						{
							--brackets;
							if(bPos == currentLine.size() - 1)
								break;
							bPos = currentLine.find("}", bPos+1);
						}
						if(!brackets)
							parameters = 0;
					}
// 					if(parameters && currentLine.find("{", 0) != std::string::npos)
// 					{
// 						++brackets;
// 						continue;
// 					}
// 					if(parameters && currentLine.find("}", 0) != std::string::npos)
// 					{
// 						--brackets;
// 						if(!brackets)
// 							parameters = 0;
// 						continue;
// 					}
					
					if(parameters && currentLine.find("TransformationMatrix ", 0) != std::string::npos)
					{
// 						std::cout << "found correct section transform parameters!" << std::endl;
						unsigned int count = 0;
						std::string::size_type loc1, loc2, loc3;
						loc1 = currentLine.find_first_of(numbers_const, 0);
						loc2 = currentLine.find_first_of(signs_const, 0);
						if(loc2 != std::string::npos)
							if(loc2 < loc1)
								loc1 = loc2;
						loc2 = currentLine.find_first_of(whitespace_const, loc1 + 1);	//ignores last value: is always 1 anyways
						while(loc2 != std::string::npos && count < 16)
						{
							char * tmp1 = new char[loc2 - loc1];
							currentLine.copy(tmp1, loc2 - loc1, loc1);
							double ftmp1 = atof(tmp1);
							transformation[count%4][count/4]= ftmp1;	// amira files are columns after each other
							loc3 = loc2;
							loc1 = currentLine.find_first_of(numbers_const, loc3);
							loc2 = currentLine.find_first_of(signs_const, loc3);
							if(loc2 != std::string::npos)
								if(loc2 < loc1)
									loc1 = loc2;
							loc2 = currentLine.find_first_of(whitespace_const, loc1 + 1);
							++count;
							delete [] tmp1;
						}
// 						std::cout << "transformation matrix:" << std::endl;
// 						for(int ii = 0; ii < 4; ++ii)
// 						{
// 							std::cout << "[";
// 							for(int jj = 0; jj < 4; ++jj)
// 							{
// 								if(jj < 3)
// 									std::cout << transformation[ii][jj] << ",\t";
// 								else
// 									std::cout << transformation[ii][jj];
// 							}
// 							std::cout << "]" << std::endl;
// 						}
						//remove numeric artifacts from z-axis:
						for(int ii = 0; ii < 2; ++ii)
						{
							transformation[2][ii] = 0;
							transformation[ii][2] = 0;
						}
						transformation[2][2] = 1;
					}
					loc = currentLine.find("VERTEX", 0);
					if(loc == 0)
					{
						if(currentLine.find("VertexCoordinates", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							vertexCoordIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("GraphLabels", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							vertexLabelIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("TransformInfo", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							vertexTransformIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
					}
					
					loc = currentLine.find("EDGE", 0);
					if(loc == 0)
					{
						if(currentLine.find("EdgeConnectivity", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeConnectivityIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("NumEdgePoints", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgePointIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("GraphLabels", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeLabelIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("TransformInfo", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeTransformIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
					}
					
					loc = currentLine.find("POINT", 0);
					if(loc == 0)
					{
						if(currentLine.find("EdgePointCoordinates", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgePointCoordIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("Radius", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							edgeRadiusIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
					}
				}
				
				if(currentIndex == vertexTransformIndex || currentIndex == edgeTransformIndex)
					continue;
				
				if(currentIndex == vertexCoordIndex)
				{
					const char * thisLine = currentLine.c_str();
					double * tmpCoords = new double[3];
					char ** endptr = new char*;
					tmpCoords[0] = strtod(thisLine, endptr);
					tmpCoords[1] = strtod(*endptr, endptr);
					tmpCoords[2] = strtod(*endptr, endptr);
					
// 					std::string::size_type loc1 = currentLine.find(" ", 0);
// 					std::string::size_type loc2 = currentLine.find(" ", loc1 + 1);
// 					char * tmp1 = new char[loc1];
// 					char * tmp2 = new char[loc2 - loc1];
// 					char * tmp3 = new char[currentLine.size() - loc2 - 1];
// 					currentLine.copy(tmp1, loc1, 0);
// 					currentLine.copy(tmp2, loc2 - loc1, loc1 + 1);
// 					currentLine.copy(tmp3, currentLine.size() - loc2 - 1, loc2 + 1);
// 					double * tmpCoords = new double[3];
// // 					tmpCoords[0] = atof(tmp1);
// // 					tmpCoords[1] = atof(tmp2);
// // 					tmpCoords[2] = atof(tmp3);
// 					char ** endptr = new char*;
// 					tmpCoords[0] = strtod(tmp1, endptr);
// 					tmpCoords[1] = strtod(tmp2, endptr);
// 					tmpCoords[2] = strtod(tmp3, endptr);
					
					tmpVertices.push_back(tmpCoords);
// 					delete [] tmp1, delete [] tmp2, delete [] tmp3;
					delete endptr;
				}
				
				if(currentIndex == vertexLabelIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					int tmplabel = atoi(tmp);
					
					tmpVertexLabels.push_back(tmplabel);
					delete [] tmp;
				}
				
				if(currentIndex == edgeConnectivityIndex)
				{
					const char * thisLine = currentLine.c_str();
					int * tmpCoords = new int[2];
					char ** endptr = new char*;
					tmpCoords[0] = static_cast< int >(strtol(thisLine, endptr, 10));
					tmpCoords[1] = static_cast< int >(strtol(*endptr, endptr, 10));
					
// 					std::string::size_type loc = currentLine.find(" ", 0);
// 					char * tmp1 = new char[loc];
// 					char * tmp2 = new char[currentLine.size() - loc - 1];
// 					currentLine.copy(tmp1, loc, 0);
// 					currentLine.copy(tmp2, currentLine.size() - loc - 1, loc + 1);
// 					int * tmpCoords = new int[2];
// 					tmpCoords[0] = atoi(tmp1);
// 					tmpCoords[1] = atoi(tmp2);
					
					tmpEdgeConnections.push_back(tmpCoords);
// 					delete [] tmp1, delete [] tmp2;
					delete endptr;
				}
				
				if(currentIndex == edgePointIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					int tmplabel = atoi(tmp);
					
					tmpNoEdgePoints.push_back(tmplabel);
					delete [] tmp;
				}
				
				if(currentIndex == edgeLabelIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					int tmplabel = atoi(tmp);
					
					tmpEdgeLabels.push_back(tmplabel);
					delete [] tmp;
				}
				
				if(currentIndex == edgePointCoordIndex)
				{
					const char * thisLine = currentLine.c_str();
					double * tmpCoords = new double[3];
					char ** endptr = new char*;
					tmpCoords[0] = strtod(thisLine, endptr);
					tmpCoords[1] = strtod(*endptr, endptr);
					tmpCoords[2] = strtod(*endptr, endptr);
					
// 					std::string::size_type loc1 = currentLine.find(" ", 0);
// 					std::string::size_type loc2 = currentLine.find(" ", loc1 + 1);
// 					char * tmp1 = new char[loc1];
// 					char * tmp2 = new char[loc2 - loc1];
// 					char * tmp3 = new char[currentLine.size() - loc2 - 1];
// 					currentLine.copy(tmp1, loc1, 0);
// 					currentLine.copy(tmp2, loc2 - loc1, loc1 + 1);
// 					currentLine.copy(tmp3, currentLine.size() - loc2 - 1, loc2 + 1);
// 					double * tmpCoords = new double[3];
// // 					tmpCoords[0] = atof(tmp1);
// // 					tmpCoords[1] = atof(tmp2);
// // 					tmpCoords[2] = atof(tmp3);
// 					char ** endptr = new char*;
// 					tmpCoords[0] = strtod(tmp1, endptr);
// 					tmpCoords[1] = strtod(tmp2, endptr);
// 					tmpCoords[2] = strtod(tmp3, endptr);
					
					edgePoints.push_back(tmpCoords);
// 					delete [] tmp1, delete [] tmp2, delete [] tmp3;
					delete endptr;
				}
				
				if(currentIndex == edgeRadiusIndex)
				{
					char * tmp = new char[currentLine.size()];
					currentLine.copy(tmp, currentLine.size(), 0);
					double tmpradius = atof(tmp);
					edgePointRadius.push_back(tmpradius);
					delete [] tmp;
				}
			}
		}
		
		std::list< double * >::iterator tmpvertexit;
		std::list< int >::iterator tmpvertexlabelit;
		std::list< int * >::iterator tmpedgeconnectivityit;
		std::list< int >::iterator tmpnumberedgepointsit;
		std::list< int >::iterator tmpedgelabelit;
		std::list< double * >::iterator edgepointit;
		std::list< double >::iterator edgeradiusit;
		
		for(tmpvertexit = tmpVertices.begin(), tmpvertexlabelit = tmpVertexLabels.begin(); tmpvertexit != tmpVertices.end() && tmpvertexlabelit != tmpVertexLabels.end(); ++tmpvertexit, ++tmpvertexlabelit)
		{
			Vertex * tmpVertex = new Vertex(*tmpvertexit, *tmpvertexlabelit);
			inputVertices.push_back(tmpVertex);
		}
		
		edgepointit = edgePoints.begin();
		edgeradiusit = edgePointRadius.begin();
		for(tmpedgeconnectivityit = tmpEdgeConnections.begin(), tmpnumberedgepointsit = tmpNoEdgePoints.begin(), tmpedgelabelit = tmpEdgeLabels.begin();
		tmpedgeconnectivityit != tmpEdgeConnections.end() && tmpnumberedgepointsit != tmpNoEdgePoints.end() && tmpedgelabelit != tmpEdgeLabels.end();
		++tmpedgeconnectivityit, ++tmpnumberedgepointsit, ++tmpedgelabelit)
		{
			std::list< double * >::iterator tmpit = edgepointit;
			std::list< double >::iterator tmpit2 = edgeradiusit;
			for(int ii = 0; ii < *tmpnumberedgepointsit; ++ii)
			{
				++tmpit;
				++tmpit2;
			}
			std::list< double * > tmpPoints(edgepointit, tmpit);
			std::list< double > tmpRadius(edgeradiusit, tmpit2);
			Edge * tmpEdge = new Edge(*tmpedgeconnectivityit, *tmpnumberedgepointsit, *tmpedgelabelit, tmpPoints, tmpRadius);
			
			inputEdges.push_back(tmpEdge);
			edgepointit = tmpit;
			edgeradiusit = tmpit2;
		}
		
// 		AmiraSpatialGraph * sg = new AmiraSpatialGraph;
		
		std::vector< Vertex * >::iterator vertexIter;
		std::vector< Edge * >::iterator edgeIter;
		
		for(vertexIter = inputVertices.begin(); vertexIter != inputVertices.end(); ++ vertexIter)
			sg->addVertex(*vertexIter);
		
		for(edgeIter = inputEdges.begin(); edgeIter != inputEdges.end(); ++edgeIter)
			sg->addEdge(*edgeIter);
// 		bool isId = 1;
// 		for(int ii = 0; ii < 4; ++ii)
// 			for(int jj = 0; jj < 4; ++jj)
// 			{
// 				if(ii != jj)
// 					if(transformation[ii][jj] != 0)
// 						isId = 0;
// 				
// 				else
// 					if(transformation[ii][jj] != 1)
// 						isId = 0;
// 			}
// 		std::cout << "transformation matrix:" << std::endl;
// 		for(int ii = 0; ii < 4; ++ii)
// 		{
// 			std::cout << "[";
// 			for(int jj = 0; jj < 4; ++jj)
// 			{
// 				if(jj < 3)
// 					std::cout << transformation[ii][jj] << ",\t";
// 				else
// 					std::cout << transformation[ii][jj];
// 			}
// 			std::cout << "]" << std::endl;
// 		}
		if(applyTransform)
		{
// 			sg->printTransformation();
			sg->setTransformation(transformation);
			sg->applyTransformation();
		}
		for(int ii = 0; ii < 4; ++ii)
			delete [] transformation[ii];
		delete [] transformation;
		
// 		std::cout << "Vertex number = " << sg->getNumberOfVertices() << std::endl;
// 		std::cout << "Edge number = " << sg->getNumberOfEdges() << std::endl;
// 		std::cout << "Point number = " << sg->getNumberOfPoints() << std::endl;
		
// 		std::cout << "VertexCoordinates @" << vertexCoordIndex << std::endl;
// 		std::cout << "Vertex GraphLabels @" << vertexLabelIndex << std::endl;
// 		std::cout << "Vertex TransformInfo @" << vertexTransformIndex << std::endl;
// 		std::cout << "EdgeConnectivity @" << edgeConnectivityIndex << std::endl;
// 		std::cout << "NumEdgePoints @" << edgePointIndex << std::endl;
// 		std::cout << "Edge GraphLabels @" << edgeLabelIndex << std::endl;
// 		std::cout << "Edge TransformInfo @" << edgeTransformIndex << std::endl;
// 		std::cout << "EdgePointCoordinates @" << edgePointCoordIndex << std::endl;
// 		std::cout << "Edge Radius @" << edgeRadiusIndex << std::endl;
	
// 		return sg;
	}
	
	inputStream.close();
	return;
}

void Reader::writeSpatialGraphFile()
{
// 	std::list<std::list<Compartment * > >::iterator edge_list_it;
// 	std::list<Compartment * >::iterator edge_it;
	std::vector< Vertex * >::iterator vertexIt;
	std::vector< Edge * >::iterator edgeIt;
	
	int number_of_edge_points = inputSpatialGraph->getNumberOfPoints();
	
// 	for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
// 	{
// 		number_of_edge_points += (*edge_list_it).size();
// 	}
// 	for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
// 	{
// 		number_of_edge_points += (*edge_list_contour_it).size();
// 	}
// 	for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
// 	{
// 		number_of_edge_points += (*edge_list_contour_it).size();
// 	}
	
	std::string format = outputFilename;
	format += ".am";
	
	#ifdef DEBUG
	std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
	//std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
	#endif
	
	std::ofstream NeuroMorphData( format.c_str() );
	
	NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	NeuroMorphData << "# This SpatialGraph file was created by the Neuron Registration Tool NeuroMap " << std::endl;
	NeuroMorphData << "# NeuroMap was programmed by Robert Egger," << std::endl;
	NeuroMorphData << "# Max-Planck-Florida Institute, Jupiter, Florida " << std::endl;
	
	NeuroMorphData << "define VERTEX " << inputSpatialGraph->getNumberOfVertices() << std::endl;
	NeuroMorphData << "define EDGE " << inputSpatialGraph->getNumberOfEdges()  << std::endl;
// 	NeuroMorphData << "define GRAPH " << /*amira_spatial_graph->vertice_list.size() +*/ /*amira_contour_graph->vertice_list.size() +*/ inputSpatialGraph->getNumberOfVertices() + /*amira_spatial_graph->edge_list.size() +*/ /*amira_contour_graph->edge_list.size() +*/ inputSpatialGraph->getNumberOfEdges() << std::endl;
	NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
	
	NeuroMorphData << "Parameters {GraphLabels {"                                	<<std::endl;
	NeuroMorphData << "        Neuron { "                                	<<std::endl;
	NeuroMorphData << "            Dendrite {"                           	<<std::endl;
	NeuroMorphData << "                ApicalDendrite {"                 	<<std::endl;
	NeuroMorphData << "                    Color 1 0.5 0.5,"          	<<std::endl;
	NeuroMorphData << "                    Id 4 }"                     	<<std::endl;
	NeuroMorphData << "                BasalDendrite {"         		<<std::endl;
	NeuroMorphData << "                    Color 0.8 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                    Id 5 }"				<<std::endl;
	NeuroMorphData << "                Color 1 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 3 }"				<<std::endl;
	NeuroMorphData << "            Axon {"					<<std::endl;
	NeuroMorphData << "                Color 0 0 1,"			<<std::endl;
	NeuroMorphData << "                Id 6 }"				<<std::endl;
	NeuroMorphData << "            Soma {"					<<std::endl;
	NeuroMorphData << "                Color 1 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 7 }"				<<std::endl;
	NeuroMorphData << "            Color 1 0 0,"				<<std::endl;
	NeuroMorphData << "            Id 2 }"					<<std::endl;
	NeuroMorphData << "        Landmark {"					<<std::endl;
	NeuroMorphData << "            Pia {"					<<std::endl;
	NeuroMorphData << "                Color 0 1 0.5,"			<<std::endl;
	NeuroMorphData << "                Id 9 }"				<<std::endl;
	NeuroMorphData << "            Vessel {"				<<std::endl;
	NeuroMorphData << "                Color 1 0.5 0,"			<<std::endl;
	NeuroMorphData << "                Id 10 }"				<<std::endl;
	NeuroMorphData << "            Barrel {"				<<std::endl;
	NeuroMorphData << "                aRow {"				<<std::endl;
	NeuroMorphData << "                    A1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 13 }"			<<std::endl;
	NeuroMorphData << "                    A2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 14 }"			<<std::endl;
	NeuroMorphData << "                    A3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 15 }"			<<std::endl;
	NeuroMorphData << "                    A4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 16 }"			<<std::endl;
	NeuroMorphData << "                Color 1 0.2 0.2,"			<<std::endl;
	NeuroMorphData << "                Id 12 }"				<<std::endl;
	NeuroMorphData << "                bRow {"				<<std::endl;
	NeuroMorphData << "                    B1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 18 }"			<<std::endl;
	NeuroMorphData << "                    B2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 19 }"			<<std::endl;
	NeuroMorphData << "                    B3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 20 }"			<<std::endl;
	NeuroMorphData << "                    B4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 21 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                    Id 17 }"				<<std::endl;
	NeuroMorphData << "                cRow {"				<<std::endl;
	NeuroMorphData << "                    C1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 23 }"			<<std::endl;
	NeuroMorphData << "                    C2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 24 }"			<<std::endl;
	NeuroMorphData << "                    C3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 25 }"			<<std::endl;
	NeuroMorphData << "                    C4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 26 }"			<<std::endl;
	NeuroMorphData << "                    C5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 27 }"			<<std::endl;
	NeuroMorphData << "                    C6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 28 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                    Id 22 }"				<<std::endl;
	NeuroMorphData << "                dRow {"				<<std::endl;
	NeuroMorphData << "                    D1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 30 }"			<<std::endl;
	NeuroMorphData << "                    D2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 31 }"			<<std::endl;
	NeuroMorphData << "                    D3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 32 }"			<<std::endl;
	NeuroMorphData << "                    D4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 33 }"			<<std::endl;
	NeuroMorphData << "                    D5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 34 }"			<<std::endl;
	NeuroMorphData << "                    D6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 35 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                    Id 29 }"				<<std::endl;
	NeuroMorphData << "                eRow {"				<<std::endl;
	NeuroMorphData << "                    E1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 37 }"			<<std::endl;
	NeuroMorphData << "                    E2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 38 }"			<<std::endl;
	NeuroMorphData << "                    E3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 39 }"			<<std::endl;
	NeuroMorphData << "                    E4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 40 }"			<<std::endl;
	NeuroMorphData << "                    E5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 41 }"			<<std::endl;
	NeuroMorphData << "                    E6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 42 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                    Id 36 }"				<<std::endl;
	NeuroMorphData << "                greekRow {"				<<std::endl;
	NeuroMorphData << "                    Alpha {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 44 }"			<<std::endl;
	NeuroMorphData << "                    Beta {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 45 }"			<<std::endl;
	NeuroMorphData << "                    Gamma {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 46 }"			<<std::endl;
	NeuroMorphData << "                    Delta {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 47 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                    Id 43 }"				<<std::endl;
	NeuroMorphData << "                Color 0 1 0,"			<<std::endl;
	NeuroMorphData << "                Id 11 }"				<<std::endl;
	NeuroMorphData << "            WhiteMatter {"				<<std::endl;
	NeuroMorphData << "                Color 0.5 1 0.75,"			<<std::endl;
	NeuroMorphData << "                Id 48 }"				<<std::endl;
	NeuroMorphData << "            OtherBarrels {"				<<std::endl;
	NeuroMorphData << "                Color 1 0 1,"			<<std::endl;
	NeuroMorphData << "                Id 49 }"				<<std::endl;
	NeuroMorphData << "            ZAxis {"					<<std::endl;
	NeuroMorphData << "                Color 0 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 50 }"				<<std::endl;
	NeuroMorphData << "            Color 0 1 1,"				<<std::endl;
	NeuroMorphData << "            Id 8 }"					<<std::endl;
	NeuroMorphData << "        Id 0,"					<<std::endl;
	NeuroMorphData << "        Color 0 0 0 }"				<<std::endl;
	NeuroMorphData << "ContentType \"HxSpatialGraph\" }"                 	<<std::endl;
	
	NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 " 		<< std::endl;
	NeuroMorphData << "VERTEX {int GraphLabels } @2 " 			<< std::endl;
	
	NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 " 		<< std::endl;
	NeuroMorphData << "EDGE { int NumEdgePoints } @4 " 			<< std::endl;
	NeuroMorphData << "EDGE { int GraphLabels } @5 " 			<< std::endl;
	
	NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 " 	<< std::endl;
	NeuroMorphData << "POINT { float Radius } @7 " 				<< std::endl;
	
	if(inputSpatialGraph->getNumberOfVertices())
	{
		NeuroMorphData << "\n@1 # Vertices xyz coordinates" 			<< std::endl;
		for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
			NeuroMorphData << (*vertexIt)->coordinates[X_COORD] << " " << (*vertexIt)->coordinates[Y_COORD]  << " " << (*vertexIt)->coordinates[Z_COORD]  << std::endl;
		
		NeuroMorphData << "\n@2 # Vertex Graph Label" << std::endl;
		for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
		{
			NeuroMorphData << (*vertexIt)->label << std::endl;
		}
	}
	
	if(inputSpatialGraph->getNumberOfEdges())
	{
		NeuroMorphData << "\n@3 # Edge Identifiers" << std::endl;
		int last_index = 0;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			NeuroMorphData << (*edgeIt)->edgeConnectivity[0] << " " << (*edgeIt)->edgeConnectivity[1] << std::endl;
		}
		
		NeuroMorphData << "\n@4 # Number of Points per Edge" << std::endl;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			NeuroMorphData <<  (*edgeIt)->numEdgePoints <<std::endl;
		}
		
		NeuroMorphData << "\n@5 # Edge Graph Labels" << std::endl;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			NeuroMorphData << (*edgeIt)->label << std::endl;
		}
	}
	
	if(inputSpatialGraph->getNumberOfPoints())
	{
		NeuroMorphData << "\n@6 # Point xyz coordinates" << std::endl;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			std::list< double * >::iterator contour_it;
			for(contour_it = (*edgeIt)->edgePointCoordinates.begin(); contour_it != (*edgeIt)->edgePointCoordinates.end(); ++contour_it)
			{
				NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
			}
			
		}
		
		NeuroMorphData << "\n@7 # Radius at Point" << std::endl;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->radiusList.size())
			{
				std::list< double >::const_iterator radiusIt;
				for(radiusIt = (*edgeIt)->radiusList.begin(); radiusIt != (*edgeIt)->radiusList.end(); ++radiusIt)
					NeuroMorphData << *radiusIt << std::endl;
			}
			else
				for(int ii = 0; ii < (*edgeIt)->edgePointCoordinates.size(); ++ii)
				{
					NeuroMorphData << (*edgeIt)->radius << std::endl;
				}
		}
	}
	
	NeuroMorphData.close();
};

void Reader::writeSpatialGraphFileFromEdges()
{
	std::vector< Edge * >::iterator edgeIt;
	
	int number_of_edge_points = inputSpatialGraph->getNumberOfPoints();
	
	std::string format = outputFilename;
	format += ".am";
	
	#ifdef DEBUG
	std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
	//std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
	#endif
	
	std::ofstream NeuroMorphData( format.c_str() );
	
	NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
	NeuroMorphData << "# This SpatialGraph file was created by the Neuron Reconstruction Tool NeuroMorph " << std::endl;
	NeuroMorphData << "# NeuroMorph was programmed by Marcel Oberlaender and Philip J. Broser," << std::endl;
	NeuroMorphData << "# Max-Planck-Institute for Medical Research Heidelberg, Germany " << std::endl;
	
	int nrOfVertices = 0;
	for(int ii = 0; ii < inputSpatialGraph->getNumberOfEdges(); ++ii)
	{
		if((*(inputSpatialGraph->edgesPointer()))[ii]->edgeConnectivity[0] == (*(inputSpatialGraph->edgesPointer()))[ii]->edgeConnectivity[1])
			nrOfVertices += 1;
		else
			nrOfVertices += 2;
	}
	
	NeuroMorphData << "define VERTEX " << nrOfVertices << std::endl;
	NeuroMorphData << "define EDGE " << inputSpatialGraph->getNumberOfEdges()  << std::endl;
	NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
	
	NeuroMorphData << "Parameters {GraphLabels {"                          	<<std::endl;
	NeuroMorphData << "        Neuron { "                                	<<std::endl;
	NeuroMorphData << "            Dendrite {"                           	<<std::endl;
	NeuroMorphData << "                ApicalDendrite {"                 	<<std::endl;
	NeuroMorphData << "                    Color 1 0.5 0.5,"          	<<std::endl;
	NeuroMorphData << "                    Id 4 }"                     	<<std::endl;
	NeuroMorphData << "                BasalDendrite {"         		<<std::endl;
	NeuroMorphData << "                    Color 0.8 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                    Id 5 }"				<<std::endl;
	NeuroMorphData << "                Color 1 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 3 }"				<<std::endl;
	NeuroMorphData << "            Axon {"					<<std::endl;
	NeuroMorphData << "                Color 0 0 1,"			<<std::endl;
	NeuroMorphData << "                Id 6 }"				<<std::endl;
	NeuroMorphData << "            Soma {"					<<std::endl;
	NeuroMorphData << "                Color 1 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 7 }"				<<std::endl;
	NeuroMorphData << "            Color 1 0 0,"				<<std::endl;
	NeuroMorphData << "            Id 2 }"					<<std::endl;
	NeuroMorphData << "        Landmark {"					<<std::endl;
	NeuroMorphData << "            Pia {"					<<std::endl;
	NeuroMorphData << "                Color 0 1 0.5,"			<<std::endl;
	NeuroMorphData << "                Id 9 }"				<<std::endl;
	NeuroMorphData << "            Vessel {"				<<std::endl;
	NeuroMorphData << "                Color 1 0.5 0,"			<<std::endl;
	NeuroMorphData << "                Id 10 }"				<<std::endl;
	NeuroMorphData << "            Barrel {"				<<std::endl;
	NeuroMorphData << "                aRow {"				<<std::endl;
	NeuroMorphData << "                    A1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 13 }"			<<std::endl;
	NeuroMorphData << "                    A2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 14 }"			<<std::endl;
	NeuroMorphData << "                    A3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 15 }"			<<std::endl;
	NeuroMorphData << "                    A4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.2 0.2,"		<<std::endl;
	NeuroMorphData << "                        Id 16 }"			<<std::endl;
	NeuroMorphData << "                Color 1 0.2 0.2,"			<<std::endl;
	NeuroMorphData << "                Id 12 }"				<<std::endl;
	NeuroMorphData << "                bRow {"				<<std::endl;
	NeuroMorphData << "                    B1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 18 }"			<<std::endl;
	NeuroMorphData << "                    B2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 19 }"			<<std::endl;
	NeuroMorphData << "                    B3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 20 }"			<<std::endl;
	NeuroMorphData << "                    B4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                        Id 21 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.25 0.25,"		<<std::endl;
	NeuroMorphData << "                    Id 17 }"				<<std::endl;
	NeuroMorphData << "                cRow {"				<<std::endl;
	NeuroMorphData << "                    C1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 23 }"			<<std::endl;
	NeuroMorphData << "                    C2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 24 }"			<<std::endl;
	NeuroMorphData << "                    C3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 25 }"			<<std::endl;
	NeuroMorphData << "                    C4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 26 }"			<<std::endl;
	NeuroMorphData << "                    C5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 27 }"			<<std::endl;
	NeuroMorphData << "                    C6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                        Id 28 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.3 0.3,"		<<std::endl;
	NeuroMorphData << "                    Id 22 }"				<<std::endl;
	NeuroMorphData << "                dRow {"				<<std::endl;
	NeuroMorphData << "                    D1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 30 }"			<<std::endl;
	NeuroMorphData << "                    D2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 31 }"			<<std::endl;
	NeuroMorphData << "                    D3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 32 }"			<<std::endl;
	NeuroMorphData << "                    D4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 33 }"			<<std::endl;
	NeuroMorphData << "                    D5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 34 }"			<<std::endl;
	NeuroMorphData << "                    D6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                        Id 35 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.35 0.35,"		<<std::endl;
	NeuroMorphData << "                    Id 29 }"				<<std::endl;
	NeuroMorphData << "                eRow {"				<<std::endl;
	NeuroMorphData << "                    E1 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 37 }"			<<std::endl;
	NeuroMorphData << "                    E2 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 38 }"			<<std::endl;
	NeuroMorphData << "                    E3 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 39 }"			<<std::endl;
	NeuroMorphData << "                    E4 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 40 }"			<<std::endl;
	NeuroMorphData << "                    E5 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 41 }"			<<std::endl;
	NeuroMorphData << "                    E6 {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                        Id 42 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.4 0.4,"		<<std::endl;
	NeuroMorphData << "                    Id 36 }"				<<std::endl;
	NeuroMorphData << "                greekRow {"				<<std::endl;
	NeuroMorphData << "                    Alpha {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 44 }"			<<std::endl;
	NeuroMorphData << "                    Beta {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 45 }"			<<std::endl;
	NeuroMorphData << "                    Gamma {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 46 }"			<<std::endl;
	NeuroMorphData << "                    Delta {"				<<std::endl;
	NeuroMorphData << "                        Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                        Id 47 }"			<<std::endl;
	NeuroMorphData << "                    Color 1 0.1 0.1,"		<<std::endl;
	NeuroMorphData << "                    Id 43 }"				<<std::endl;
	NeuroMorphData << "                Color 0 1 0,"			<<std::endl;
	NeuroMorphData << "                Id 11 }"				<<std::endl;
	NeuroMorphData << "            WhiteMatter {"				<<std::endl;
	NeuroMorphData << "                Color 0.5 1 0.75,"			<<std::endl;
	NeuroMorphData << "                Id 48 }"				<<std::endl;
	NeuroMorphData << "            OtherBarrels {"				<<std::endl;
	NeuroMorphData << "                Color 1 0 1,"			<<std::endl;
	NeuroMorphData << "                Id 49 }"				<<std::endl;
	NeuroMorphData << "            ZAxis {"					<<std::endl;
	NeuroMorphData << "                Color 0 0 0,"			<<std::endl;
	NeuroMorphData << "                Id 50 }"				<<std::endl;
	NeuroMorphData << "            Color 0 1 1,"				<<std::endl;
	NeuroMorphData << "            Id 8 }"					<<std::endl;
	NeuroMorphData << "        Id 0,"					<<std::endl;
	NeuroMorphData << "        Color 0 0 0 }"				<<std::endl;
	NeuroMorphData << "ContentType \"HxSpatialGraph\" }"                 	<<std::endl;
	
	NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 " 		<< std::endl;
	NeuroMorphData << "VERTEX {int GraphLabels } @2 " 			<< std::endl;
	
	NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 " 		<< std::endl;
	NeuroMorphData << "EDGE { int NumEdgePoints } @4 " 			<< std::endl;
	NeuroMorphData << "EDGE { int GraphLabels } @5 " 			<< std::endl;
	
	NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 " 	<< std::endl;
	NeuroMorphData << "POINT { float Radius } @7 " 				<< std::endl;
	
	if(inputSpatialGraph->getNumberOfEdges())
	{
		NeuroMorphData << "\n@1 # Vertices xyz coordinates" 			<< std::endl;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->edgeConnectivity[0] == (*edgeIt)->edgeConnectivity[1])
				NeuroMorphData << (*edgeIt)->edgePointCoordinates.front()[X_COORD] << " " << (*edgeIt)->edgePointCoordinates.front()[Y_COORD] << " " << (*edgeIt)->edgePointCoordinates.front()[Z_COORD] << std::endl;
			else
			{
				NeuroMorphData << (*edgeIt)->edgePointCoordinates.front()[X_COORD] << " " << (*edgeIt)->edgePointCoordinates.front()[Y_COORD] << " " << (*edgeIt)->edgePointCoordinates.front()[Z_COORD] << std::endl;
				NeuroMorphData << (*edgeIt)->edgePointCoordinates.back()[X_COORD] << " " << (*edgeIt)->edgePointCoordinates.back()[Y_COORD] << " " << (*edgeIt)->edgePointCoordinates.back()[Z_COORD] << std::endl;
			}
		}
		
		NeuroMorphData << "\n@2 # Vertex Graph Label" << std::endl;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->edgeConnectivity[0] == (*edgeIt)->edgeConnectivity[1])
				NeuroMorphData << (*edgeIt)->label << std::endl;
			else
			{
				NeuroMorphData << (*edgeIt)->label << std::endl;
				NeuroMorphData << (*edgeIt)->label << std::endl;
			}
		}
		
		NeuroMorphData << "\n@3 # Edge Identifiers" << std::endl;
		int lastID = 0;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->edgeConnectivity[0] == (*edgeIt)->edgeConnectivity[1])
			{
				NeuroMorphData << lastID << " " << lastID << std::endl;
				++lastID;
			}
			else
			{
				NeuroMorphData << lastID << " " << lastID + 1 << std::endl;
				lastID += 2;
			}
		}
		
		NeuroMorphData << "\n@4 # Number of Points per Edge" << std::endl;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			NeuroMorphData <<  (*edgeIt)->numEdgePoints <<std::endl;
		}
		
		NeuroMorphData << "\n@5 # Edge Graph Labels" << std::endl;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			NeuroMorphData << (*edgeIt)->label << std::endl;
		}
	}
	
	if(inputSpatialGraph->getNumberOfPoints())
	{
		NeuroMorphData << "\n@6 # Point xyz coordinates" << std::endl;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			std::list< double * >::iterator contour_it;
			for(contour_it = (*edgeIt)->edgePointCoordinates.begin(); contour_it != (*edgeIt)->edgePointCoordinates.end(); ++contour_it)
			{
				NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
			}
			
		}
		
		NeuroMorphData << "\n@7 # Radius at Point" << std::endl;
		for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->radiusList.size())
			{
				std::list< double >::const_iterator radiusIt;
				for(radiusIt = (*edgeIt)->radiusList.begin(); radiusIt != (*edgeIt)->radiusList.end(); ++radiusIt)
					NeuroMorphData << *radiusIt << std::endl;
			}
			else
				for(int ii = 0; ii < (*edgeIt)->edgePointCoordinates.size(); ++ii)
				{
					NeuroMorphData << (*edgeIt)->radius << std::endl;
				}
				
		}
	}
	
	NeuroMorphData.close();
};

/* Read in SpatialGraphSetFile
 * originalGraphIndices (1 x 1517439) -> cell ID
 * cellTypeIDs (1 x 1517439) -> Cell Type ID
 * spatialGraphTransforms (1 x 1517439) -> transformation matrix
 * originalGraphFiles (1 x 35482) -> filenames
 * cellTypeIDLabels (1 x 43) -> Map
 */
void Reader::readSpatialGraphSetFile(const char * fname, std::vector< unsigned int >& originalGraphIndices, std::vector< unsigned int >& spatialGraphSetLabels,
									std::vector< double * >& spatialGraphTransforms, std::vector< std::string >& originalGraphFiles,
									std::map< unsigned int, std::string >& cellTypeIDLabels)
{
	std::ifstream inputStream(fname);
	
	if(!inputStream.fail())
	{
		std::cout << "Reading SpatialGraphSet file " << fname << "..." << std::endl;
		unsigned int graphIndicesIndex = -1, graphSetLabelsIndex = -1, transformIndex = -1;
		unsigned int currentIndex = 0;
		unsigned int nGraphs = 0;
		bool parameters = 0, parameters_files = 0, parameters_labels = 0;
		unsigned int brackets = 0, labelBrackets = 0, fileBrackets = 0;
		std::list< unsigned int > cellTypeIDs;
		std::list< std::string > cellTypeLabels;
		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
		{
// 			if(!parameters)
// 			{
// 				std::cout << currentLine << std::endl;
// 				++line;
// 			}
			
			if(currentLine.size())
			{
				if(currentLine.find("#", 0) == 0)
				{
					continue;
				}
				if(currentLine.find("@", 0) == 0)
				{
					char * tmp = new char[currentLine.size() - 1];
					currentLine.copy(tmp, currentLine.size() - 1, 1);
					currentIndex = atoi(tmp);
#ifdef DEBUG
					std::cout << std::endl;
					std::cout << "**********************" << std::endl;
					std::cout << "Reading data section " << currentIndex << std::endl;
					std::cout << "**********************" << std::endl;
#endif
					delete [] tmp;
					continue;
				}
				
				if(currentIndex == 0)
				{
					std::string::size_type loc = currentLine.find("define", 0);
					if(loc != std::string::npos)
					{
						if(currentLine.find("GRAPH", 7) != std::string::npos)
						{
							char * tmp = new char[currentLine.size() - 13];
							currentLine.copy(tmp, currentLine.size() - 13, 13);
							nGraphs = atoi(tmp);
// 							std::cout << "vertex = " << vertex << std::endl;
// 							inputVertices.resize(vertex);
							delete [] tmp;
							continue;
						}
					}
					
					loc = currentLine.find("Parameters", 0);
					if(loc == 0)
					{
						parameters = 1;
						if(currentLine.find("{", 0) != std::string::npos)
							brackets = 1;
						continue;
					}
// 					if(parameters && currentLine.find("{", 0) == std::string::npos && currentLine.find("}", 0) == std::string::npos)
// 						continue;
					if(parameters)
					{
						std::string::size_type startPos = 0;
						for(std::string::size_type bPos = currentLine.find("{", startPos); bPos != std::string::npos; )
						{
							++brackets;
							if(bPos == currentLine.size() - 1)
								break;
							bPos = currentLine.find("{", bPos+1);
						}
						for(std::string::size_type bPos = currentLine.find("}", startPos); bPos != std::string::npos; )
						{
							--brackets;
							if(bPos == currentLine.size() - 1)
								break;
							bPos = currentLine.find("}", bPos+1);
						}
						if(!brackets)
						{
							parameters = 0;
							parameters_files = 0;
							parameters_labels = 0;
							continue;
						}
						
						if(currentLine.find("Files") != std::string::npos)
						{
							parameters_files = 1;
							fileBrackets = brackets - 1;
							continue;
						}
						if(currentLine.find("SpatialGraphSetLabels") != std::string::npos)
						{
							parameters_labels = 1;
							labelBrackets = brackets - 1;
							continue;
						}
						if(brackets <= fileBrackets)
						{
							parameters_files = 0;
							continue;
						}
						if(brackets <= labelBrackets)
						{
							parameters_labels = 0;
							continue;
						}
						
						if(parameters_labels)
						{
							if(currentLine.find("Color") != std::string::npos)
							{
								continue;
							}
							if(currentLine.find("Id") != std::string::npos)
							{
								char idChar[2];
								unsigned int cellTypeID;
								sscanf(currentLine.c_str(), " %2c %d", idChar, &cellTypeID);
								cellTypeIDs.push_back(cellTypeID);
#ifdef DEBUG
								std::flush(std::cout << "currentLine: " << currentLine.c_str() << "\n");
								std::flush(std::cout << "cellTypeID: " << cellTypeID << "\n");
#endif
							}
							if(currentLine.find("{") != std::string::npos)
							{
								size_t delim1 = currentLine.find_first_of(letters_const);
								size_t delim2 = currentLine.find_first_of(whitespace_const, delim1);
								std::string cellTypeLabel = currentLine.substr(delim1, delim2-delim1);
								cellTypeLabels.push_back(cellTypeLabel);
#ifdef DEBUG
								std::flush(std::cout << "currentLine: " << currentLine.c_str() << "\n");
								std::flush(std::cout << "cellTypeLabel: " << cellTypeLabel.c_str() << "| (end of string)\n");
#endif
							}
							continue;
						}
						if(parameters_files)
						{
							if(currentLine.find("Path") != std::string::npos)
							{
								char * pathName = new char[1024];
								sscanf(currentLine.c_str(), " Path %s", pathName);
								originalGraphFiles.push_back(std::string(pathName));
#ifdef DEBUG
// 								std::flush(std::cout << "currentLine: " << currentLine.c_str() << "\n");
								std::flush(std::cout << "Found path: " << pathName << "\r");
#endif
								delete [] pathName;
								continue;
							}
							else if(currentLine.find("File") != std::string::npos)
							{
								// nothing to be done
// 								unsigned int * fileID = new unsigned int;
// 								sscanf(currentLine.c_str(), " File%d { ", fileID);
							}
						}
					}
					
					loc = currentLine.find("GRAPH", 0);
					if(loc == 0)
					{
						if(currentLine.find("OriginalGraphIndices", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							graphIndicesIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("SpatialGraphSetLabels", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							graphSetLabelsIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
						if(currentLine.find("AdditionalTransforms", 0) != std::string::npos)
						{
							loc = currentLine.find("@", 0);
							char * tmp = new char[currentLine.size() - loc - 1];
							currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
							transformIndex = atoi(tmp);
							delete [] tmp;
							continue;
						}
					}
				}
				
				if(currentIndex == graphIndicesIndex)
				{
					unsigned int index;
					sscanf(currentLine.c_str(), " %d ", &index);
					originalGraphIndices.push_back(index);
				}
				
				if(currentIndex == graphSetLabelsIndex)
				{
					unsigned int label;
					sscanf(currentLine.c_str(), " %d ", &label);
					spatialGraphSetLabels.push_back(label);
				}
				
				if(currentIndex == transformIndex)
				{
					double * transform = new double[16];
					sscanf(currentLine.c_str(), " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
						transform, transform+1, transform+2, transform+3, transform+4, transform+5, transform+6, transform+7,
						transform+8, transform+9, transform+10, transform+11, transform+12, transform+13, transform+14, transform+15);
					spatialGraphTransforms.push_back(transform);
#ifdef DEBUG
					std::flush(std::cout << "Found transform: " << *transform << " " << *(transform+1) << " " << *(transform+2) << " " << *(transform+3) << " "
					<< *(transform+4) << " " << *(transform+5) << " " << *(transform+6) << " " << *(transform+7) << " " << *(transform+8) << " " << *(transform+9) << " "
					<< *(transform+10) << " " << *(transform+11) << " " << *(transform+12) << " " << *(transform+13) << " " << *(transform+14) << " " << *(transform+15)
					<< "\r");
#endif
				}
			}
		}
		
#ifdef DEBUG
		std::flush(std::cout << std::endl);
#endif
		std::list< unsigned int >::const_iterator cellTypeIDsIt;
		std::list< std::string >::const_iterator cellTypeLabelsIt;
		for(cellTypeIDsIt = cellTypeIDs.begin(), cellTypeLabelsIt = cellTypeLabels.begin();
			cellTypeIDsIt != cellTypeIDs.end(), cellTypeLabelsIt != cellTypeLabels.end(); ++cellTypeIDsIt, ++cellTypeLabelsIt)
		{
			cellTypeIDLabels.insert(std::pair< unsigned int, std::string >(*cellTypeIDsIt, *cellTypeLabelsIt));
#ifdef DEBUG
			std::flush(std::cout << "Cell type ID/label pair: " << *cellTypeIDsIt << " / " << *cellTypeLabelsIt << std::endl);
#endif
		}
	}
	else
	{
		std::cout << "Error opening SpatialGraphSet file " << fname << "!" << std::endl;
	}
}

void Reader::writeAmiraSurfaceFile(PolyDataPointerType triangleData)
{
	std::string format = outputFilename;
	format += ".surf";
	
	std::ofstream SurfaceData( format.c_str() );
	
	SurfaceData << "# HyperSurface 0.1 ASCII" << std::endl;
	SurfaceData << "" << std::endl;
	SurfaceData << "Parameters {" << std::endl;
	SurfaceData << "\tMaterials {" << std::endl;
	SurfaceData << "\t\tExterior {" << std::endl;
	SurfaceData << "\t\t\tid 0," << std::endl;
	SurfaceData << "\t\t\tColor 1 1 1" << std::endl;
	SurfaceData << "\t\t}" << std::endl;
	SurfaceData << "\t\tUnsortedContours0 {" << std::endl;
	SurfaceData << "\t\t\tid 1," << std::endl;
	SurfaceData << "\t\t\tColor 1 0 0" << std::endl;
	SurfaceData << "\t\t}" << std::endl;
	SurfaceData << "\t}" << std::endl;
	SurfaceData << "\tBoundaryIds {" << std::endl;
	SurfaceData << "\t\tname \"BoundaryConditions\"" << std::endl;
	SurfaceData << "\t}" << std::endl;
	SurfaceData << "}" << std::endl;
	SurfaceData << "" << std::endl;
	
	SurfaceData << "Vertices " << triangleData->GetNumberOfPoints() << std::endl;
	for(int ii = 0; ii < triangleData->GetNumberOfPoints(); ++ii)
	{
// 		double * point = triangleData->GetPoint(ii);
		double point[3];
		triangleData->GetPoint(ii, point);
		SurfaceData << std::fixed << "\t" << point[0] << " " << point[1] << " " << point[2] << std::endl;
	}
	
	SurfaceData << "NBranchingPoints " << 0 << std::endl;
	SurfaceData << "NVerticesOnCurves " << 0 << std::endl;
	SurfaceData << "BoundaryCurves " << 0 << std::endl;
	SurfaceData << "Patches " << 1 << std::endl;
	SurfaceData << "{" << std::endl;
	SurfaceData << "InnerRegion UnsortedContours" << 0 << std::endl;
	SurfaceData << "OuterRegion Exterior" << std::endl;
	SurfaceData << "BoundaryID " << 0 << std::endl;
	SurfaceData << "BranchingPoints " << 0 << std::endl;
	SurfaceData << "" << std::endl;
	SurfaceData << "Triangles " << triangleData->GetNumberOfPolys() << std::endl;
	for(unsigned int ii = 0; ii < triangleData->GetNumberOfPolys(); ++ii)
	{
		vtkIdList * pointIDs = triangleData->GetCell(ii)->GetPointIds();
		SurfaceData << "\t" << pointIDs->GetId(0) + 1 << " " << pointIDs->GetId(1) + 1 << " " << pointIDs->GetId(2) + 1 << std::endl;
	}
	SurfaceData << "}" << std::endl;
	
	SurfaceData.close();
};

void Reader::readHocFileAllPiaContours()
{
	std::ifstream inputStream(inputFilename);

	if(!inputStream.fail())
	{
		std::cout << "Reading hoc file " << inputFilename << std::endl;
		std::list< double * > edgePtList;
		std::list< double > radiusList;
		std::list< unsigned int > edgePtCountList;
		std::list< unsigned int > labelList;
		std::map< std::string, int > segmentInsertOrder;
		std::map< int, std::string > segmentFatherList;
		bool readPts = 0;
		unsigned int edgePtCount = 0;
		int insertCnt = 0;
		int homeBarrel = 0;
		std::list< int >::const_iterator labelIt;
		segmentInsertOrder.insert(std::pair< const char *, int >("soma", -1));

		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof())
			if(currentLine.size())
			{
				// in case neuron is registered
				// to different column:
				// use registered home barrel
				// in NeuroMap header
				if(!homeBarrel && currentLine.find("Registered home barrel") != std::string::npos)
				{
					for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
					{
						int ID = *labelIt;
						if(currentLine.find(int2Labels[ID]) != std::string::npos)
						{
							homeBarrel = ID;
							break;
						}

					}
				}
				// ignore comments
				if(currentLine.find("/*") != std::string::npos && currentLine.find("*/") != std::string::npos)
					continue;
				// read pts belonging to current segment
				if(readPts)
				{
					// ignore spines for now
					if(currentLine.find("Spine") != std::string::npos)
						continue;
					if(currentLine.find("pt3dadd") != std::string::npos)
					{
						size_t nrStart = currentLine.find("(") + 1;
						size_t sep1 = currentLine.find(",", nrStart);
						size_t sep2 = currentLine.find(",", sep1+1);
						size_t sep3 = currentLine.find(",", sep2+1);
						char * xChar = new char[sep1-nrStart+1];
						char * yChar = new char[sep2-sep1];
						char * zChar = new char[sep3-sep2];
						char * rChar = new char[currentLine.size()-sep3];
						currentLine.copy(xChar, sep1-nrStart+1, nrStart);
						currentLine.copy(yChar, sep2-sep1, sep1+1);
						currentLine.copy(zChar, sep3-sep2, sep2+1);
						currentLine.copy(rChar, currentLine.size()-sep3, sep3+1);
						double * tmpCoords = new double[3], radius;
						char ** endptr = new char*;
						tmpCoords[0] = strtod(xChar, endptr);
						tmpCoords[1] = strtod(yChar, endptr);
						tmpCoords[2] = strtod(zChar, endptr);
						radius = strtod(rChar, endptr);
						edgePtList.push_back(tmpCoords);
						radiusList.push_back(radius);
						++edgePtCount;
						delete endptr;
						continue;
					}
					else if(currentLine.find("pt3dadd") == std::string::npos && edgePtCount)
					{
						readPts = 0;
						edgePtCountList.push_back(edgePtCount);
						edgePtCount = 0;
					}
				}
				// case: pia (== alpha)
				/* MODIFIED LINE */
				if((currentLine.find("alpha") != std::string::npos
						|| currentLine.find("Pia") != std::string::npos
						|| currentLine.find("pia") != std::string::npos)
						&& (currentLine.find("create") != std::string::npos))
				{
					labelList.push_back(Pia);
					readPts = 1;
					edgePtCount = 0;
					++insertCnt;
					continue;
				}
				// case: white matter
				if(currentLine.find("WM") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
					labelList.push_back(WhiteMatter);
					readPts = 1;
					edgePtCount = 0;
					++insertCnt;
					continue;
				}
				// case: soma
				if(currentLine.find("soma") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
					labelList.push_back(Soma);
					readPts = 1;
					edgePtCount = 0;
					++insertCnt;
					continue;
				}

				// case: dendrite
				if(currentLine.find("dend") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
					labelList.push_back(Dendrite);
					readPts = 1;
					edgePtCount = 0;
					// insert name
					size_t nameStart = currentLine.find_first_of(letters, currentLine.find_first_of(whitespace));
					size_t nameEnd = currentLine.find_first_of(" \t}", nameStart);
					std::string name(currentLine, nameStart, (nameEnd-nameStart));
					segmentInsertOrder.insert(std::pair< std::string, int >(name, insertCnt));
					++insertCnt;
					continue;
				}

				// case: basal dendrite
				if(currentLine.find("BasalDendrite") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
					labelList.push_back(BasalDendrite);
					readPts = 1;
					edgePtCount = 0;
					// insert name
					size_t nameStart = currentLine.find_first_of(letters, currentLine.find_first_of(whitespace));
					size_t nameEnd = currentLine.find_first_of(" \t}", nameStart);
					std::string name(currentLine, nameStart, (nameEnd-nameStart));
					segmentInsertOrder.insert(std::pair< std::string, int >(name, insertCnt));
					++insertCnt;
					continue;
				}

				// case: apical denrite
				if(currentLine.find("apical") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
					labelList.push_back(ApicalDendrite);
					readPts = 1;
					edgePtCount = 0;
					// insert name
					size_t nameStart = currentLine.find_first_of(letters, currentLine.find_first_of(whitespace));
					size_t nameEnd = currentLine.find_first_of(" \t}", nameStart);
					std::string name(currentLine, nameStart, (nameEnd-nameStart));
					segmentInsertOrder.insert(std::pair< std::string, int >(name, insertCnt));
					++insertCnt;
					continue;
				}
				// case: axon
				if(currentLine.find("axon") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
					labelList.push_back(Axon);
					readPts = 1;
					edgePtCount = 0;
// 					// insert name
					size_t nameStart = currentLine.find_first_of(letters, currentLine.find_first_of(whitespace));
					size_t nameEnd = currentLine.find_first_of(" \t}", nameStart);
					std::string name(currentLine, nameStart, (nameEnd-nameStart));
					segmentInsertOrder.insert(std::pair< const char *, unsigned int >(name.c_str(), insertCnt));
					++insertCnt;
					continue;
				}
				// case: barrel
				for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
				{
					int ID = *labelIt;
					if(currentLine.find(int2Labels[ID]) != std::string::npos && currentLine.find("create") != std::string::npos)
					{
						if(!homeBarrel && currentLine.find("home_barrel") != std::string::npos)
							homeBarrel = ID;
						labelList.push_back(ID);
						readPts = 1;
						edgePtCount = 0;
						++insertCnt;
						break;
					}
				}
				// connectivity
				if(currentLine.find("connect") != std::string::npos)
				{
					// connected to soma => end!
					if(currentLine.find("soma") != std::string::npos)
					{
						segmentFatherList.insert(std::pair< int, const char * >(insertCnt-1, "soma"));
						continue;
					}
					size_t tmp = currentLine.find(")");
					size_t nameStart = currentLine.find_first_of(letters, tmp+1);
					size_t nameEnd = currentLine.find("(", nameStart);
					std::string fatherName(currentLine, nameStart, (nameEnd-nameStart));
					segmentFatherList.insert(std::pair< int, std::string >(insertCnt-1, fatherName));
				}
			}
		// make sure EOF doesn't mess up anything
		if(edgePtCount && edgePtCountList.size() == labelList.size()-1)
			edgePtCountList.push_back(edgePtCount);

		// fill SpatialGraph with input data
		if(!inputSpatialGraph)
			inputSpatialGraph = new AmiraSpatialGraph;
		if(homeBarrel)
			inputSpatialGraph->setHomeBarrel(homeBarrel);
		if(edgePtCountList.size() == labelList.size())
		{
			int segmentCnt = 0;
			std::list< double * >::const_iterator edgePtListIt = edgePtList.begin();
			std::list< double >::const_iterator radiusListIt = radiusList.begin();
			std::list< unsigned int >::const_iterator edgePtCountListIt;
			std::list< unsigned int >::const_iterator labelListIt;
			for(labelListIt = labelList.begin(), edgePtCountListIt = edgePtCountList.begin();
			labelListIt != labelList.end(), edgePtCountListIt != edgePtCountList.end();
			++labelListIt, ++edgePtCountListIt, ++segmentCnt)
			{
				int nrOfEdgePts = *edgePtCountListIt, edgeConnectivity[2], ID = *labelListIt;
				std::list< double * >::const_iterator tmpIt = edgePtListIt;
				std::list< double >::const_iterator tmpIt2 = radiusListIt;
				for(int ii = 0; ii < nrOfEdgePts; ++ii)
					++tmpIt, ++tmpIt2;
				std::list< double * > segmentPtList(edgePtListIt, tmpIt);
				std::list< double > segmentRadiusList(radiusListIt, tmpIt2);
				edgePtListIt = tmpIt, radiusListIt = tmpIt2;
				if(ID == Soma)
				{
					Vertex * newVert1 = new Vertex(segmentPtList.front(), ID);
					Vertex * newVert2 = new Vertex(segmentPtList.back(), ID);
					inputSpatialGraph->addVertex(newVert1);
					inputSpatialGraph->addVertex(newVert2);
					edgeConnectivity[0] = inputSpatialGraph->getNumberOfVertices()-2;
					edgeConnectivity[1] = inputSpatialGraph->getNumberOfVertices()-1;
				} // Connect Neurites To Soma
				else if(ID == Dendrite || ID == ApicalDendrite || ID == BasalDendrite || ID == Axon)
				{
					Vertex * newVert2 = new Vertex(segmentPtList.back(), ID);
					inputSpatialGraph->addVertex(newVert2);
					edgeConnectivity[0] = -1;
					edgeConnectivity[1] = inputSpatialGraph->getNumberOfVertices()-1;
				}
				else // Do NOT connect other structures to Soma
				{
					Vertex * newVert1 = new Vertex(segmentPtList.front(), ID);
					Vertex * newVert2 = new Vertex(segmentPtList.back(), ID);
					inputSpatialGraph->addVertex(newVert1);
					inputSpatialGraph->addVertex(newVert2);
					edgeConnectivity[0] = inputSpatialGraph->getNumberOfVertices()-2;
					edgeConnectivity[1] = inputSpatialGraph->getNumberOfVertices()-1;
				}
				Edge * newSegment  = new Edge(edgeConnectivity, nrOfEdgePts, ID, segmentPtList, segmentRadiusList);
				if(ID == Dendrite || ID == ApicalDendrite || ID == BasalDendrite || ID == Axon)
				{
					newSegment->fatherID = segmentInsertOrder[segmentFatherList[segmentCnt]];
				}
				inputSpatialGraph->addEdge(newSegment);
			}

			for(int ii = 0; ii < inputSpatialGraph->edgesPointer()->size(); ++ii)
			{
				Edge * tmpSegment = inputSpatialGraph->edgesPointer()->at(ii);
				int edgeLabel = tmpSegment->label;

				// Only set up connectivity if neurite
				if(tmpSegment->fatherID > -1 && (edgeLabel == Dendrite || edgeLabel == ApicalDendrite || edgeLabel == BasalDendrite || edgeLabel == Axon))
				{
					double * fatherVertexPt = inputSpatialGraph->edgesPointer()->at(tmpSegment->fatherID)->edgePointCoordinates.back();
					int fatherVertexID = inputSpatialGraph->edgesPointer()->at(tmpSegment->fatherID)->edgeConnectivity[1];
					tmpSegment->edgeConnectivity[0] = fatherVertexID;
				}
				else if(tmpSegment->fatherID == -1 && (edgeLabel == Dendrite || edgeLabel == ApicalDendrite || edgeLabel == BasalDendrite || edgeLabel == Axon))
				{
					if(!inputSpatialGraph->isLabelInSpatialGraph(Soma))
					{
						std::cout << "Error! Cannot set up correct connectivity without soma at root!" << std::endl;
					}
					for(int jj = 0; jj < inputSpatialGraph->edgesPointer()->size(); ++jj)
					{
						Edge * tmpEdge = inputSpatialGraph->edgesPointer()->at(jj);
						if(tmpEdge->label == Soma)
						{
							double * fatherVertexPt = tmpEdge->edgePointCoordinates.back();
							int fatherVertexID = tmpEdge->edgeConnectivity[1];
							double * newEdgePt = new double[3];
							newEdgePt[0] = fatherVertexPt[0];
							newEdgePt[1] = fatherVertexPt[1];
							newEdgePt[2] = fatherVertexPt[2];
							tmpSegment->edgeConnectivity[0] = fatherVertexID;
							tmpSegment->edgePointCoordinates.insert(tmpSegment->edgePointCoordinates.begin(), newEdgePt);
							tmpSegment->radiusList.insert(tmpSegment->radiusList.begin(), tmpSegment->radiusList.front());
							tmpSegment->numEdgePoints += 1;
						}
					}
				}
			}
		}
		else
			std::cout << "Error reading hoc file! Nr of labels is not equal to nr of edges" << std::endl;

		inputStream.close();
	}
};

void Reader::readHocFile()
{
	std::ifstream inputStream(inputFilename);
	
	if(!inputStream.fail())
	{
		std::cout << "Reading hoc file " << inputFilename << std::endl;
//		initializeConstants();
		std::list< double * > edgePtList;
		std::list< double > radiusList;
		std::list< unsigned int > edgePtCountList;
		std::list< unsigned int > labelList;
		std::map< std::string, int > segmentInsertOrder;
		std::map< int, std::string > segmentFatherList;
		bool readPts = 0;
		unsigned int edgePtCount = 0;
		int insertCnt = 0;
		int homeBarrel = 0;
		std::list< int >::const_iterator labelIt;
		segmentInsertOrder.insert(std::pair< const char *, int >("soma", -1));
		
		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof())
			if(currentLine.size())
			{
				// in case neuron is registered
				// to different column:
				// use registered home barrel
				// in NeuroMap header
				if(!homeBarrel && currentLine.find("Registered home barrel") != std::string::npos)
				{
					for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
					{
						int ID = *labelIt;
						if(currentLine.find(int2Labels[ID]) != std::string::npos)
						{
							homeBarrel = ID;
							break;
						}

					}
				}
				// ignore comments
				if(currentLine.find("/*") != std::string::npos && currentLine.find("*/") != std::string::npos)
					continue;
				// read pts belonging to current segment
				if(readPts)
				{
// 					std::cout << "checking pt for label " << labelList.back() << std::endl;
					// ignore spines for now
					if(currentLine.find("Spine") != std::string::npos)
						continue;
					if(currentLine.find("pt3dadd") != std::string::npos)
					{
// 						std::cout << "reading pt for label " << int2Labels[labelList.back()] << std::endl;
						size_t nrStart = currentLine.find("(") + 1;
						size_t sep1 = currentLine.find(",", nrStart);
						size_t sep2 = currentLine.find(",", sep1+1);
						size_t sep3 = currentLine.find(",", sep2+1);
						char * xChar = new char[sep1-nrStart+1];
						char * yChar = new char[sep2-sep1];
						char * zChar = new char[sep3-sep2];
						char * rChar = new char[currentLine.size()-sep3];
						currentLine.copy(xChar, sep1-nrStart+1, nrStart);
						currentLine.copy(yChar, sep2-sep1, sep1+1);
						currentLine.copy(zChar, sep3-sep2, sep2+1);
						currentLine.copy(rChar, currentLine.size()-sep3, sep3+1);
// 						std::cout << xChar << " " << yChar << " " << zChar << " " << rChar << std::endl;
// 						const char * thisLine = currentLine.c_str();
						double * tmpCoords = new double[3], radius;
						char ** endptr = new char*;
						tmpCoords[0] = strtod(xChar, endptr);
						tmpCoords[1] = strtod(yChar, endptr);
						tmpCoords[2] = strtod(zChar, endptr);
						radius = strtod(rChar, endptr);
// 						std::cout << tmpCoords[0] << " " << tmpCoords[1] << " " << tmpCoords[2] << " " << radius << std::endl;
						edgePtList.push_back(tmpCoords);
						radiusList.push_back(radius);
						++edgePtCount;
						delete endptr;
						continue;
					}
					else if(currentLine.find("pt3dadd") == std::string::npos && edgePtCount)
					{
						readPts = 0;
						edgePtCountList.push_back(edgePtCount);
						edgePtCount = 0;
					}
				}
				// case: pia (== alpha)
				if(currentLine.find("alpha") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
					labelList.push_back(Pia);
					readPts = 1;
					edgePtCount = 0;
					++insertCnt;
					continue;
				}
				// case: white matter
				if(currentLine.find("WM") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
					labelList.push_back(WhiteMatter);
					readPts = 1;
					edgePtCount = 0;
					++insertCnt;
					continue;
				}
				// case: soma
				if(currentLine.find("soma") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
#ifdef DEBUG
					std::cout << "Found soma!" << std::endl;
#endif
					labelList.push_back(Soma);
					readPts = 1;
					edgePtCount = 0;
					++insertCnt;
					continue;
				}
#ifdef DEBUG
				if(currentLine.find("create") != std::string::npos)
				{
					std::size_t somaPos = currentLine.find("soma");
					std::size_t createPos = currentLine.find("create");
					std::cout << "*****************************" << std::endl;
					std::cout << "Line:" << std::endl;
					std::cout << currentLine.c_str() << std::endl;
					std::cout << "somaPos:   " << somaPos << std::endl;
					std::cout << "createPos: " << createPos << std::endl;
				}
#endif
				// case: dendrite
				if(currentLine.find("dend") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
#ifdef DEBUG
					std::cout << "Found dendrite!" << std::endl;
#endif
					labelList.push_back(Dendrite);
					readPts = 1;
					edgePtCount = 0;
					// insert name
					size_t nameStart = currentLine.find_first_of(letters, currentLine.find_first_of(whitespace));
					size_t nameEnd = currentLine.find_first_of(" \t}", nameStart);
// 					char * name = new char[nameEnd-nameStart];
// 					currentLine.copy(name, (nameEnd-nameStart), nameStart);
// 					segmentInsertOrder.insert(std::pair< const char *, int >(name, insertCnt));
					std::string name(currentLine, nameStart, (nameEnd-nameStart));
					segmentInsertOrder.insert(std::pair< std::string, int >(name, insertCnt));
					++insertCnt;
					continue;
				}

				// case: basal dendrite
				if(currentLine.find("BasalDendrite") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
#ifdef DEBUG
					std::cout << "Found BasalDendrite!" << std::endl;
#endif
					labelList.push_back(BasalDendrite);
					readPts = 1;
					edgePtCount = 0;
					// insert name
					size_t nameStart = currentLine.find_first_of(letters, currentLine.find_first_of(whitespace));
					size_t nameEnd = currentLine.find_first_of(" \t}", nameStart);
// 					char * name = new char[nameEnd-nameStart];
// 					currentLine.copy(name, (nameEnd-nameStart), nameStart);
// 					segmentInsertOrder.insert(std::pair< const char *, int >(name, insertCnt));
					std::string name(currentLine, nameStart, (nameEnd-nameStart));
					segmentInsertOrder.insert(std::pair< std::string, int >(name, insertCnt));
					++insertCnt;
					continue;
				}

				// case: apical denrite
				if(currentLine.find("apical") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
#ifdef DEBUG
					std::cout << "Found apical!" << std::endl;
#endif
					labelList.push_back(ApicalDendrite);
					readPts = 1;
					edgePtCount = 0;
					// insert name
					size_t nameStart = currentLine.find_first_of(letters, currentLine.find_first_of(whitespace));
					size_t nameEnd = currentLine.find_first_of(" \t}", nameStart);
// 					char * name = new char[nameEnd-nameStart];
// 					currentLine.copy(name, (nameEnd-nameStart), nameStart);
// 					segmentInsertOrder.insert(std::pair< const char *, int >(name, insertCnt));
					std::string name(currentLine, nameStart, (nameEnd-nameStart));
					segmentInsertOrder.insert(std::pair< std::string, int >(name, insertCnt));
					++insertCnt;
					continue;
				}
				// case: axon
				if(currentLine.find("axon") != std::string::npos && currentLine.find("create") != std::string::npos)
				{
#ifdef DEBUG
					std::cout << "Found axon!" << std::endl;
#endif
					labelList.push_back(Axon);
					readPts = 1;
					edgePtCount = 0;
// 					// insert name
					size_t nameStart = currentLine.find_first_of(letters, currentLine.find_first_of(whitespace));
					size_t nameEnd = currentLine.find_first_of(" \t}", nameStart);
					std::string name(currentLine, nameStart, (nameEnd-nameStart));
					segmentInsertOrder.insert(std::pair< const char *, unsigned int >(name.c_str(), insertCnt));
					++insertCnt;
					continue;
				}
				// case: barrel
				for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
				{
					int ID = *labelIt;
					if(currentLine.find(int2Labels[ID]) != std::string::npos && currentLine.find("create") != std::string::npos)
					{
						if(!homeBarrel && currentLine.find("home_barrel") != std::string::npos)
							homeBarrel = ID;
						labelList.push_back(ID);
						readPts = 1;
						edgePtCount = 0;
						++insertCnt;
						break;
					}
				}
				// connectivity
				if(currentLine.find("connect") != std::string::npos)
				{
					// connected to soma => end!
					if(currentLine.find("soma") != std::string::npos)
					{
						segmentFatherList.insert(std::pair< int, const char * >(insertCnt-1, "soma"));
						continue;
					}
					size_t tmp = currentLine.find(")");
					size_t nameStart = currentLine.find_first_of(letters, tmp+1);
					size_t nameEnd = currentLine.find("(", nameStart);
// 					char * name = new char[nameEnd-nameStart];
// 					currentLine.copy(name, (nameEnd-nameStart), nameStart);
// 					segmentFatherList.insert(std::pair< int, const char * >(insertCnt-1, name));
					std::string fatherName(currentLine, nameStart, (nameEnd-nameStart));
					segmentFatherList.insert(std::pair< int, std::string >(insertCnt-1, fatherName));
				}
#ifdef DEBUG
				std::cout << "***********************" << std::endl;
				std::cout << "Line:" << std::endl;
				std::cout << currentLine.c_str() << std::endl;
				std::cout << "Didn't find anything..." << std::endl;
#endif
			}
		// make sure EOF doesn't mess up anything
		if(edgePtCount && edgePtCountList.size() == labelList.size()-1)
			edgePtCountList.push_back(edgePtCount);
		
		//debug
// 		std::cout << "home barrel = " << int2Labels[homeBarrel] << std::endl;
// 		std::cout << "nr of points = " << edgePtList.size() << std::endl;
// 		std::cout << "nr of labels = " << labelList.size() << std::endl;
// 		std::cout << "nr of segments = " << edgePtCountList.size() << std::endl;
		
		// fill SpatialGraph with input data
		if(!inputSpatialGraph)
			inputSpatialGraph = new AmiraSpatialGraph;
		if(homeBarrel)
			inputSpatialGraph->setHomeBarrel(homeBarrel);
		
// 		std::map< int, std::string >::const_iterator fatherListIt;
// 		std::map< std::string, int >::const_iterator segmentOrderIt;
// 		std::cout << "segment father map:" << std::endl;
// 		for(fatherListIt = segmentFatherList.begin(); fatherListIt != segmentFatherList.end(); ++fatherListIt)
// 			std::cout << fatherListIt->first << " -> " << fatherListIt->second << std::endl;
// 		std::cout << "segment order map:" << std::endl;
// 		for(segmentOrderIt = segmentInsertOrder.begin(); segmentOrderIt != segmentInsertOrder.end(); ++segmentOrderIt)
// 			std::cout << segmentOrderIt->first << " -> " << segmentOrderIt->second << std::endl;
		
		if(edgePtCountList.size() == labelList.size())
		{
			int segmentCnt = 0;
			std::list< double * >::const_iterator edgePtListIt = edgePtList.begin();
			std::list< double >::const_iterator radiusListIt = radiusList.begin();
			std::list< unsigned int >::const_iterator edgePtCountListIt;
			std::list< unsigned int >::const_iterator labelListIt;
			for(labelListIt = labelList.begin(), edgePtCountListIt = edgePtCountList.begin();
			labelListIt != labelList.end(), edgePtCountListIt != edgePtCountList.end();
			++labelListIt, ++edgePtCountListIt, ++segmentCnt)
			{
				int nrOfEdgePts = *edgePtCountListIt, edgeConnectivity[2], ID = *labelListIt;
				std::list< double * >::const_iterator tmpIt = edgePtListIt;
				std::list< double >::const_iterator tmpIt2 = radiusListIt;
				for(int ii = 0; ii < nrOfEdgePts; ++ii)
					++tmpIt, ++tmpIt2;
				std::list< double * > segmentPtList(edgePtListIt, tmpIt);
				std::list< double > segmentRadiusList(radiusListIt, tmpIt2);
				edgePtListIt = tmpIt, radiusListIt = tmpIt2;
				if(ID == Soma)
				{
					Vertex * newVert1 = new Vertex(segmentPtList.front(), ID);
					Vertex * newVert2 = new Vertex(segmentPtList.back(), ID);
					inputSpatialGraph->addVertex(newVert1);
					inputSpatialGraph->addVertex(newVert2);
					edgeConnectivity[0] = inputSpatialGraph->getNumberOfVertices()-2;
					edgeConnectivity[1] = inputSpatialGraph->getNumberOfVertices()-1;
				} // Connect Neurites To Soma
				else if(ID == Dendrite || ID == ApicalDendrite || ID == BasalDendrite || ID == Axon)
				{
	// 				Vertex * newVert1 = new Vertex(segmentPtList.front(), ID);
					Vertex * newVert2 = new Vertex(segmentPtList.back(), ID);
					inputSpatialGraph->addVertex(newVert2);
	// 				inputSpatialGraph->addVertex(newVert1), inputSpatialGraph->addVertex(newVert2);
	// 				edgeConnectivity[0] = inputSpatialGraph->getNumberOfVertices()-2, edgeConnectivity[1] = inputSpatialGraph->getNumberOfVertices()-1;
					edgeConnectivity[0] = -1;
					edgeConnectivity[1] = inputSpatialGraph->getNumberOfVertices()-1;
				}
				else // Do NOT connect other structures to Soma
				{
					Vertex * newVert1 = new Vertex(segmentPtList.front(), ID);
					Vertex * newVert2 = new Vertex(segmentPtList.back(), ID);
					inputSpatialGraph->addVertex(newVert1);
					inputSpatialGraph->addVertex(newVert2);
					edgeConnectivity[0] = inputSpatialGraph->getNumberOfVertices()-2;
					edgeConnectivity[1] = inputSpatialGraph->getNumberOfVertices()-1;  
				}
				Edge * newSegment  = new Edge(edgeConnectivity, nrOfEdgePts, ID, segmentPtList, segmentRadiusList);
				if(ID == Dendrite || ID == ApicalDendrite || ID == BasalDendrite || ID == Axon)
				{
// 					std::cout << "label = " << ID << std::endl;
// 					std::cout << "segmentCnt = " << segmentCnt << std::endl;
// 					std::cout << "segmentFatherList[segmentCnt] = " << segmentFatherList[segmentCnt] << std::endl;
// 					std::cout << "fatherID = " << segmentInsertOrder[segmentFatherList[segmentCnt]] << std::endl;
					newSegment->fatherID = segmentInsertOrder[segmentFatherList[segmentCnt]];
				}
				inputSpatialGraph->addEdge(newSegment);
			}
			
			for(int ii = 0; ii < inputSpatialGraph->edgesPointer()->size(); ++ii)
			{
				Edge * tmpSegment = inputSpatialGraph->edgesPointer()->at(ii);
				int edgeLabel = tmpSegment->label;
				
				// Only set up connectivity if neurite
				if(tmpSegment->fatherID > -1 && (edgeLabel == Dendrite || edgeLabel == ApicalDendrite || edgeLabel == BasalDendrite || edgeLabel == Axon))
				{
					double * fatherVertexPt = inputSpatialGraph->edgesPointer()->at(tmpSegment->fatherID)->edgePointCoordinates.back();
					int fatherVertexID = inputSpatialGraph->edgesPointer()->at(tmpSegment->fatherID)->edgeConnectivity[1];
// 					Vertex * newVert = new Vertex(segmentPtList.back(), ID);
// 					inputSpatialGraph->addVertex(newVert);
					tmpSegment->edgeConnectivity[0] = fatherVertexID;
// 					tmpSegment->edgeConnectivity[1] = inputSpatialGraph->getNumberOfVertices()-1;
// 					inputSpatialGraph->addEdge(newSegment);
				}
				else if(tmpSegment->fatherID == -1 && (edgeLabel == Dendrite || edgeLabel == ApicalDendrite || edgeLabel == BasalDendrite || edgeLabel == Axon))
				{
					if(!inputSpatialGraph->isLabelInSpatialGraph(Soma))
					{
						std::cout << "Error! Cannot set up correct connectivity without soma at root!" << std::endl;
					}
					for(int jj = 0; jj < inputSpatialGraph->edgesPointer()->size(); ++jj)
					{
						Edge * tmpEdge = inputSpatialGraph->edgesPointer()->at(jj);
						if(tmpEdge->label == Soma)
						{
							double * fatherVertexPt = tmpEdge->edgePointCoordinates.back();
							int fatherVertexID = tmpEdge->edgeConnectivity[1];
							double * newEdgePt = new double[3];
							newEdgePt[0] = fatherVertexPt[0];
							newEdgePt[1] = fatherVertexPt[1];
							newEdgePt[2] = fatherVertexPt[2];
// 							Vertex * newVert = new Vertex(segmentPtList.back(), ID);
// 							inputSpatialGraph->addVertex(newVert);
							tmpSegment->edgeConnectivity[0] = fatherVertexID;
// 							tmpSegment->edgeConnectivity[1] = inputSpatialGraph->getNumberOfVertices()-1;
							tmpSegment->edgePointCoordinates.insert(tmpSegment->edgePointCoordinates.begin(), newEdgePt);
							tmpSegment->radiusList.insert(tmpSegment->radiusList.begin(), tmpSegment->radiusList.front());
							tmpSegment->numEdgePoints += 1; 
// 							inputSpatialGraph->addEdge(newSegment);
						}
					}
				}
			}
		}
		else
			std::cout << "Error reading hoc file! Nr of labels is not equal to nr of edges" << std::endl;
		
		inputStream.close();
	}
};

// write hoc file in same format as input hoc file
// only point coordinates are different
void Reader::writeHocFile()
{
	int homeBarrel = inputSpatialGraph->getHomeBarrel();
	std::string ofName(outputFilename);
	ofName += "_registered_";
	if(homeBarrel)
	{
		ofName += int2Labels[homeBarrel];
	}
	else
	{
		ofName += "global";
	}
	ofName += ".hoc";
	std::ifstream inputStream(inputFilename);
	std::ofstream outStream(ofName.c_str());
	
	if(!inputStream.fail() && !outStream.fail())
	{
		std::cout << "Writing hoc file " << ofName.c_str() << std::endl;
		std::vector< Edge * >::const_iterator edgeIt = inputSpatialGraph->edgesBegin();
		std::list< double * >::const_iterator edgePtIt = (*edgeIt)->edgePointCoordinates.begin();
		std::list< double >::const_iterator radiusIt = (*edgeIt)->radiusList.begin();
		std::list< int >::const_iterator labelIt;
		
		outStream << "/*---------------------------------------------------------------------------*/" << std::endl;
		outStream << "/* Neuron morphology registered to standard barrel field with NeuroMap       */" << std::endl;
		if(homeBarrel)
		{
			outStream << "/* Registered home barrel: " << int2Labels[homeBarrel] << "                                                */" << std::endl;
		}
		outStream << "/*---------------------------------------------------------------------------*/" << std::endl;
		outStream << std::endl;
		
		bool writeFlag = 1;
		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof())
		{
			// ignore spines for now
			if(currentLine.find("Spine") != std::string::npos)
				continue;
			// ignore daVinci registration
			if(currentLine.find("/* EOF */") != std::string::npos)
				break;
			if(currentLine.find("create") != std::string::npos)
			{
				bool ignoreLabel = 1;
				std::list< const char * >::const_iterator hocLabelIt;
				for(hocLabelIt = hocLabels.begin(); hocLabelIt != hocLabels.end(); ++hocLabelIt)
					if(currentLine.find(*hocLabelIt) != std::string::npos)
						ignoreLabel = 0;
				if(ignoreLabel)
					writeFlag = 0;
				else
					writeFlag = 1;
			}
			if(currentLine.find("pt3dadd") != std::string::npos)
			{
				if(writeFlag && edgeIt != inputSpatialGraph->edgesEnd())
				{
					// readHocFile adds point connecting to soma (skip this point, otherwise all point coordinates are shifted by one, leading to the last point being removed)
					// This is necessary because the output is written like the input hoc file (in order to work it needs the same number of points)
					if( (*edgeIt)->fatherID == -1 && edgePtIt == (*edgeIt)->edgePointCoordinates.begin() && ((*edgeIt)->label == Dendrite || (*edgeIt)->label == BasalDendrite || (*edgeIt)->label == ApicalDendrite || (*edgeIt)->label == Axon ))
					{
						++edgePtIt;
						++radiusIt;
					}
				  
					if(edgePtIt == (*edgeIt)->edgePointCoordinates.end())
					{
						++edgeIt;
						edgePtIt = (*edgeIt)->edgePointCoordinates.begin();
						radiusIt = (*edgeIt)->radiusList.begin();
						
						if((*edgeIt)->fatherID == -1 && ((*edgeIt)->label == Dendrite || (*edgeIt)->label == BasalDendrite || (*edgeIt)->label == ApicalDendrite || (*edgeIt)->label == Axon))
						{
							++radiusIt;
							++edgePtIt;
						}
					}
					
					outStream << "{pt3dadd(" << (*edgePtIt)[0] << "," << (*edgePtIt)[1] << "," << (*edgePtIt)[2] << "," << *radiusIt << ")}" << std::endl;
					
					if(edgePtIt != (*edgeIt)->edgePointCoordinates.end())
					{
						++edgePtIt;
						++radiusIt;
					}
				}
				else if(!writeFlag)
					outStream << currentLine << std::endl;
				continue;
			}
// 			for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
// 			{
// 				int ID = *labelIt;
// 				if(currentLine.find(int2Labels[ID]) != std::string::npos)
// 				{
// 					
// 				}
// 				continue;
// 			}
			//case: nothing special -> just pass through
			outStream << currentLine << std::endl;
		}
		outStream.close();
	}
	else
		std::cout << "Error writing hoc file!" << std::endl;
};

// write hoc file in same format as input hoc file
// only point coordinates are different
void Reader::writeSeparateHocFiles()
{
	int homeBarrel = inputSpatialGraph->getHomeBarrel();
	std::string ofName1(outputFilename);
	std::string ofName2(outputFilename);
	ofName1 += "_neuron_registered_";
	ofName2 += "_landmarks_registered_";
	if(homeBarrel)
	{
		ofName1 += int2Labels[homeBarrel];
		ofName2 += int2Labels[homeBarrel];
	}
	else
	{
		ofName1 += "global";
		ofName2 += "global";
	}
	ofName1 += ".hoc";
	ofName2 += ".hoc";
	std::ifstream inputStream(inputFilename);
	std::ofstream outStream1(ofName1.c_str());
	std::ofstream outStream2(ofName2.c_str());
	
	if(!inputStream.fail() && !outStream1.fail() && !outStream2.fail())
	{
		std::cout << "Writing hoc file " << ofName1.c_str() << std::endl;
		std::vector< Edge * >::const_iterator edgeIt = inputSpatialGraph->edgesBegin();
		std::list< double * >::const_iterator edgePtIt = (*edgeIt)->edgePointCoordinates.begin();
		std::list< double >::const_iterator radiusIt = (*edgeIt)->radiusList.begin();
		std::list< int >::const_iterator labelIt;
		
		outStream1 << "/*---------------------------------------------------------------------------*/" << std::endl;
		outStream1 << "/* Neuron morphology registered to standard barrel field with NeuroMap       */" << std::endl;
		if(homeBarrel)
		{
			outStream1 << "/* Registered home barrel: " << int2Labels[homeBarrel] << "                                                */" << std::endl;
		}
		outStream1 << "/*---------------------------------------------------------------------------*/" << std::endl;
		outStream1 << std::endl;
		outStream2 << "/*---------------------------------------------------------------------------*/" << std::endl;
		outStream2 << "/* Landmarks registered to standard barrel field with NeuroMap               */" << std::endl;
		if(homeBarrel)
		{
			outStream2 << "/* Registered home barrel: " << int2Labels[homeBarrel] << "                                                */" << std::endl;
		}
		outStream2 << "/*---------------------------------------------------------------------------*/" << std::endl;
		outStream2 << std::endl;
		
		bool writeFlag = 1;
		bool writeNeuronFlag = 1;
		bool writeLandmarkFlag = 1;
		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof())
		{
			// ignore spines for now
			if(currentLine.find("Spine") != std::string::npos)
				continue;
			// ignore daVinci registration
			if(currentLine.find("/* EOF */") != std::string::npos)
				break;
			if(currentLine.find("create") != std::string::npos)
			{
				bool ignoreLabel = 1;
				bool neuronLabel = 0;
				bool landmarkLabel = 0;
				std::list< const char * >::const_iterator hocLabelIt;
				for(hocLabelIt = hocNeuronLabels.begin(); hocLabelIt != hocNeuronLabels.end(); ++hocLabelIt)
					if(currentLine.find(*hocLabelIt) != std::string::npos)
					{
						ignoreLabel = 0;
						neuronLabel = 1;
					}
				for(hocLabelIt = hocLandmarkLabels.begin(); hocLabelIt != hocLandmarkLabels.end(); ++hocLabelIt)
					if(currentLine.find(*hocLabelIt) != std::string::npos)
					{
						ignoreLabel = 0;
						landmarkLabel = 1;
					}
				if(ignoreLabel)
				{
					writeFlag = 0;
					writeNeuronFlag = 0;
					writeLandmarkFlag = 0;
				}
				else
				{
					writeFlag = 1;
					if(neuronLabel && !landmarkLabel)
					{
						writeNeuronFlag = 1;
						writeLandmarkFlag = 0;
					}
					else if(!neuronLabel && landmarkLabel)
					{
						writeNeuronFlag = 0;
						writeLandmarkFlag = 1;
					}
					else if(neuronLabel && landmarkLabel)
					{
						std::cout << "Error! Corrupt hoc file format. Abort writing..." << std::endl;
						break;
					}
				}
			}
			if(currentLine.find("pt3dadd") != std::string::npos)
			{
				if(writeFlag && edgeIt != inputSpatialGraph->edgesEnd())
				{
					// readHocFile adds point connecting to soma (skip this point, otherwise all point coordinates are shifted by one, leading to the last point being removed)
					// This is necessary because the output is written like the input hoc file (in order to work it needs the same number of points)
					if( (*edgeIt)->fatherID == -1 && edgePtIt == (*edgeIt)->edgePointCoordinates.begin() && ((*edgeIt)->label == Dendrite || (*edgeIt)->label == ApicalDendrite || (*edgeIt)->label == BasalDendrite || (*edgeIt)->label == Axon))
					{
						++radiusIt;
						++edgePtIt;
					}
					
					if(edgePtIt == (*edgeIt)->edgePointCoordinates.end())
					{
						++edgeIt;
						edgePtIt = (*edgeIt)->edgePointCoordinates.begin();
						radiusIt = (*edgeIt)->radiusList.begin();
						
						if((*edgeIt)->fatherID == -1 && ((*edgeIt)->label == Dendrite || (*edgeIt)->label == ApicalDendrite || (*edgeIt)->label == BasalDendrite || (*edgeIt)->label == Axon))
						{
							++radiusIt;
							++edgePtIt;
						}
					}
					if(writeNeuronFlag)
						outStream1 << "{pt3dadd(" << (*edgePtIt)[0] << "," << (*edgePtIt)[1] << "," << (*edgePtIt)[2] << "," << *radiusIt << ")}" << std::endl;
					else if(writeLandmarkFlag)
						outStream2 << "{pt3dadd(" << (*edgePtIt)[0] << "," << (*edgePtIt)[1] << "," << (*edgePtIt)[2] << "," << *radiusIt << ")}" << std::endl;
					if(edgePtIt != (*edgeIt)->edgePointCoordinates.end())
					{
						++edgePtIt;
						++radiusIt;
					}
				}
//				else if(!writeFlag && writeNeuronFlag && !writeLandmarkFlag)
//					outStream1 << currentLine << std::endl;
//				else if(!writeFlag && !writeNeuronFlag && writeLandmarkFlag)
//					outStream2 << currentLine << std::endl;
//				else if(!writeFlag && !writeNeuronFlag && !writeLandmarkFlag)
//				{
//					outStream1 << currentLine << std::endl;
//					outStream2 << currentLine << std::endl;
//				}
				continue;
			}
// 			for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
// 			{
// 				int ID = *labelIt;
// 				if(currentLine.find(int2Labels[ID]) != std::string::npos)
// 				{
// 					
// 				}
// 				continue;
// 			}
			//case: nothing special -> just pass through
//			outStream << currentLine << std::endl;
			if(writeNeuronFlag && !writeLandmarkFlag)
				outStream1 << currentLine << std::endl;
			else if(!writeNeuronFlag && writeLandmarkFlag)
				outStream2 << currentLine << std::endl;
			else if(writeNeuronFlag && writeLandmarkFlag)
			{
				outStream1 << currentLine << std::endl;
				outStream2 << currentLine << std::endl;
			}
		}
		outStream1.close();
		outStream2.close();
	}
	else
		std::cout << "Error writing hoc file!" << std::endl;
};

PolyDataPointerType Reader::readAmiraSurfaceFile()
{
	std::ifstream inputStream(inputFilename);
	PolyDataPointerType surface = PolyDataPointerType::New();
	PointsPointerType points = PointsPointerType::New();
	TransformPointerType transform = NULL;
	
	if(!inputStream.fail())
	{
		std::string currentLine;
		
		unsigned int currentIndex = 0;
		const unsigned int point = 1, cell = 2;
		unsigned int pointID = 0, cellID = 0;
		
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				std::string::size_type loc1, loc2, loc3;
				
				if(!currentIndex)
				{
					loc1 = currentLine.find("Vertices", 0);
					if(loc1 == 0)
					{
						
						currentIndex = point;
						char * tmp = new char[currentLine.size() - 9];
						currentLine.copy(tmp, currentLine.size() - 9, 9);
						int noOfPoints = atoi(tmp);
						points->SetDataTypeToFloat();
						points->SetNumberOfPoints(noOfPoints);
						delete [] tmp;
					}
					
					loc1 = currentLine.find("Triangles", 0);
					if(loc1 == 0)
					{
						currentIndex = cell;
						char * tmp = new char[currentLine.size() - 10];
						currentLine.copy(tmp, currentLine.size() - 10, 10);
						int noOfCells = atoi(tmp);
						surface->Allocate(1);
						delete [] tmp;
					}
				}
				
				else if(currentIndex == point || currentIndex == cell)
				{
					if(currentLine.find_first_of(letters, 0) != std::string::npos || currentLine.find_first_of(otherChars, 0) != std::string::npos)
					{
						currentIndex = 0;
						continue;
					}
					
					if(currentIndex == point)
					{
						const char * thisLine = currentLine.c_str();
						double * tmpCoords = new double[3];
						char ** endptr = new char*;
						tmpCoords[0] = strtod(thisLine, endptr);
						tmpCoords[1] = strtod(*endptr, endptr);
						tmpCoords[2] = strtod(*endptr, endptr);
						points->SetPoint(pointID, tmpCoords);
						++pointID;
						delete [] tmpCoords, delete endptr;
					}
					if(currentIndex == cell)
					{
						const char * thisLine = currentLine.c_str();
						int * tmpCoords = new int[3];
						char ** endptr = new char*;
						tmpCoords[0] = static_cast< int >(strtol(thisLine, endptr, 10));
						tmpCoords[1] = static_cast< int >(strtol(*endptr, endptr, 10));
						tmpCoords[2] = static_cast< int >(strtol(*endptr, endptr, 10));
						PolygonPointerType poly = PolygonPointerType::New();
						poly->GetPointIds()->SetNumberOfIds(3);
						for(int ii = 0; ii < 3; ++ii)
							poly->GetPointIds()->SetId(ii, tmpCoords[ii]-1);
						surface->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
						delete [] tmpCoords, delete endptr;
					}
					continue;
				}
				
				if(currentLine.find("TransformationMatrix") != std::string::npos)
				{
					
					transform = TransformPointerType::New();
					HomogeneousMatrixPointerType mat = HomogeneousMatrixPointerType::New();
					std::string::size_type loc1, loc2;
					loc1 = currentLine.find_first_of(numbers, 0);
					loc2 = currentLine.find_first_of(signs, 0);
					if(loc2 != std::string::npos)
						if(loc2 < loc1)
							loc1 = loc2;
					char * matChar = new char[currentLine.size()-loc1];
					currentLine.copy(matChar, currentLine.size()-loc1, loc1);
					double tmpMat[16];
					char ** endptr = new char*;
					tmpMat[0] = strtod(matChar, endptr);
					for(int ii = 1; ii < 16; ++ii)
						tmpMat[ii] = strtod(*endptr, endptr);
					for(int ii = 0; ii < 4; ++ii)
						for(int jj = 0; jj < 4; ++jj)
							mat->SetElement(jj, ii, tmpMat[ii*4+jj]);
					transform->SetMatrix(mat);
					transform->Update();
				}
			}
		}
		surface->SetPoints(points);
		surface->Update();
		inputStream.close();
	}
	else
		std::cerr << "Error! Could not read file " << inputFilename << std::endl;
	
	if(!transform)
		return surface;
	else
	{
		TransformFilterType tFilter = TransformFilterType::New();
		tFilter->SetTransform(transform);
		tFilter->SetInput(surface);
		tFilter->Update();
		return tFilter->GetOutput();
	}
};


/****************************************************************************/
/*read Amira landmark file with closed contours                             */
/****************************************************************************/
PointsPointerType Reader::readLandmarkFile()
{
	PointsPointerType landmarkPts = PointsPointerType::New();
	landmarkPts->SetDataTypeToFloat();
	
	std::ifstream inputStream(inputFilename);
	if(!inputStream.fail())
	{
		std::string currentLine;
		
		unsigned int currentIndex = 0;
		const unsigned int point = 1, cell = 2;
		unsigned int pointID = 0, cellID = 0;
		
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				std::string::size_type loc1, loc2, loc3;
				
				if(!currentIndex)
				{
					loc1 = currentLine.find("@1", 0);
					if(loc1 == 0)
					{
						
						currentIndex = point;
						continue;
// 						char * tmp = new char[currentLine.size() - 9];
// 						currentLine.copy(tmp, currentLine.size() - 9, 9);
// 						int noOfPoints = atoi(tmp);
// 						points->SetDataTypeToFloat();
// 						points->SetNumberOfPoints(noOfPoints);
// 						delete [] tmp;
					}
				}
				
				else if(currentIndex == point)
				{
					const char * thisLine = currentLine.c_str();
					double * tmpCoords = new double[3];
					char ** endptr = new char*;
					tmpCoords[0] = strtod(thisLine, endptr);
					tmpCoords[1] = strtod(*endptr, endptr);
					tmpCoords[2] = strtod(*endptr, endptr);
					landmarkPts->InsertNextPoint(tmpCoords);
					++pointID;
					delete endptr;
				}
			}
		}
// 		double * lastPt = new double[3];
// 		landmarkPts->GetPoint(0, lastPt);
// 		landmarkPts->InsertNextPoint(lastPt);
// 		PolygonPointerType poly = PolygonPointerType::New();
// 		poly->GetPointIds()->SetNumberOfIds(landmarkPts->GetNumberOfPoints());
// 		for(int ii = 0; ii < landmarkPts->GetNumberOfPoints(); ++ii)
// 			poly->GetPointIds()->SetId(ii, ii);
// 		PolyDataPointerType polyData = PolyDataPointerType::New();
// 		polyData->Allocate(1);
// 		polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
// 		polyData->SetPoints(landmarkPts);
// 		polyData->Update();
// 		if(!inputSpatialGraph)
// 			inputSpatialGraph = new AmiraSpatialGraph;
// 		inputSpatialGraph->addPolyDataObject(polyData, 0);
	}
	inputStream.close();
	
	return landmarkPts;
};

/****************************************************************************/
/*read Amira landmark file with closed contours                             */
/****************************************************************************/
void Reader::readLandmarkFile(bool applyTransform)
{
	std::ifstream inputStream(inputFilename);
	if(!inputStream.fail())
	{
		PointsPointerType landmarkPts = PointsPointerType::New();
		landmarkPts->SetDataTypeToFloat();
		std::string currentLine;
		
		double ** transformation = new double *[4];
		for(int ii = 0; ii < 4; ++ii)
		{
			transformation[ii] = new double[4];
			for(int jj = 0; jj < 4; ++jj)
			{
				if(ii != jj)
					transformation[ii][jj] = 0;
				else
					transformation[ii][jj] = 1;
			}
		}
		
		unsigned int currentIndex = 0;
		const unsigned int point = 1, cell = 2;
		unsigned int pointID = 0, cellID = 0;
		
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				if(currentLine.find("TransformationMatrix ", 0) != std::string::npos)
				{
// 					std::cout << "found correct section transform parameters!" << std::endl;
					unsigned int count = 0;
					std::string::size_type loc1, loc2, loc3;
					loc1 = currentLine.find_first_of(numbers, 0);
					loc2 = currentLine.find_first_of(signs, 0);
					if(loc2 != std::string::npos)
						if(loc2 < loc1)
							loc1 = loc2;
					loc2 = currentLine.find_first_of(whitespace, loc1 + 1);	//ignores last value: is always 1 anyways
					while(loc2 != std::string::npos && count < 16)
					{
						char * tmp1 = new char[loc2 - loc1];
						currentLine.copy(tmp1, loc2 - loc1, loc1);
						double ftmp1 = atof(tmp1);
						transformation[count%4][count/4]= ftmp1;	// amira files are columns after each other
						loc3 = loc2;
						loc1 = currentLine.find_first_of(numbers, loc3);
						loc2 = currentLine.find_first_of(signs, loc3);
						if(loc2 != std::string::npos)
							if(loc2 < loc1)
								loc1 = loc2;
						loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
						++count;
						delete [] tmp1;
					}
				}
				
				std::string::size_type loc1, loc2, loc3;
				if(!currentIndex)
				{
					loc1 = currentLine.find("@1", 0);
					if(loc1 == 0)
					{
						
						currentIndex = point;
						continue;
// 						char * tmp = new char[currentLine.size() - 9];
// 						currentLine.copy(tmp, currentLine.size() - 9, 9);
// 						int noOfPoints = atoi(tmp);
// 						points->SetDataTypeToFloat();
// 						points->SetNumberOfPoints(noOfPoints);
// 						delete [] tmp;
					}
				}
				
				else if(currentIndex == point)
				{
					const char * thisLine = currentLine.c_str();
					double * tmpCoords = new double[3];
					char ** endptr = new char*;
					tmpCoords[0] = strtod(thisLine, endptr);
					tmpCoords[1] = strtod(*endptr, endptr);
					tmpCoords[2] = strtod(*endptr, endptr);
					landmarkPts->InsertNextPoint(tmpCoords);
					++pointID;
					delete endptr;
				}
			}
		}
// 		double * lastPt = new double[3];
// 		landmarkPts->GetPoint(0, lastPt);
// 		landmarkPts->InsertNextPoint(lastPt);
		PolygonPointerType poly = PolygonPointerType::New();
		poly->GetPointIds()->SetNumberOfIds(landmarkPts->GetNumberOfPoints());
		for(int ii = 0; ii < landmarkPts->GetNumberOfPoints(); ++ii)
			poly->GetPointIds()->SetId(ii, ii);
		PolyDataPointerType polyData = PolyDataPointerType::New();
		polyData->Allocate(1);
		polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
		polyData->SetPoints(landmarkPts);
		polyData->Update();
		if(!inputSpatialGraph)
			inputSpatialGraph = new AmiraSpatialGraph;
		inputSpatialGraph->addPolyDataObject(polyData, 0);
		if(applyTransform)
		{
			inputSpatialGraph->setTransformation(transformation);
			inputSpatialGraph->applyTransformation();
		}
	}
	inputStream.close();
};

void Reader::writeLandmarkFile(PointsPointerType pts)
{
	if(pts)
	{
		std::string ofName(outputFilename);
		ofName += ".landmarkAscii";
		std::ofstream outStream(ofName.c_str());
		
		outStream << "# AmiraMesh 3D ASCII 2.0" << std::endl;
		outStream << std::endl;
		outStream << "define Markers " << pts->GetNumberOfPoints() << std::endl;
		outStream << std::endl;
		outStream << "Parameters {" << std::endl;
		outStream << "\tNumSets 1," << std::endl;
		outStream << "\tContentType \"LandmarkSet\"" << std::endl;
		outStream << "}" << std::endl;
		outStream << std::endl;
		outStream << "Markers { float[3] Coordinates } @1" << std::endl;
		outStream << std::endl;
		outStream << "# Data section follows" << std::endl;
		outStream << "@1" << std::endl;
		
		for(int ii = 0; ii < pts->GetNumberOfPoints(); ++ii)
		{
			double tmp[3];
			pts->GetPoint(ii, tmp);
			outStream << tmp[0] << " " << tmp[1] << " " << tmp[2] << std::endl;
		}
		
		outStream.close();
	}
};

ImageDataPointerType Reader::readScalarField()
{
	ImageDataPointerType field = ImageDataPointerType::New();
	
	std::ifstream inputStream(inputFilename);
	if(!inputStream.fail())
	{
		field->SetNumberOfScalarComponents(1);
		field->SetScalarTypeToDouble();
		
		double spacing[3], bounds[6];
		int extent[6], dims[3], extentDiff[3];
		bool dataSection = 0, hasExtent = 0, hasBounds = 0;
		unsigned long index = 0;
		
		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				std::string::size_type loc1, loc2, loc3;
				
				// set up lattice
				if(!dataSection)
				{
					if(currentLine.find("define") != std::string::npos
						&& currentLine.find("Lattice") != std::string::npos)
					{
						std::size_t startNrs = currentLine.find("Lattice") + 7;
						std::string dimStr = currentLine.substr(startNrs);
						char ** endptr = new char*;
						
						dims[0] = static_cast< int >(strtol(dimStr.c_str(), endptr, 10));
						dims[1] = static_cast< int >(strtol(*endptr, endptr, 10));
						dims[2] = static_cast< int >(strtol(*endptr, endptr, 10));
						
						for(int ii = 0; ii < 3; ++ii)
						{
							extent[2*ii] = 0;
							extent[2*ii+1] = dims[ii] - 1;
							extentDiff[ii] = extent[2*ii+1] - extent[2*ii];
						}
						
						field->SetExtent(extent);
						hasExtent = 1;
						
						delete endptr;
					}
					
					if(currentLine.find("BoundingBox") != std::string::npos)
					{
						std::size_t startNrs = currentLine.find("BoundingBox") + 11;
						std::string bbStr = currentLine.substr(startNrs);
						char ** endptr = new char*;
						
						bounds[0] = strtod(bbStr.c_str(), endptr);
						bounds[1] = strtod(*endptr, endptr);
						bounds[2] = strtod(*endptr, endptr);
						bounds[3] = strtod(*endptr, endptr);
						bounds[4] = strtod(*endptr, endptr);
						bounds[5] = strtod(*endptr, endptr);
						
						field->SetOrigin(bounds[0], bounds[2], bounds[4]);
						hasBounds = 1;
						
						delete endptr;
					}
					
					if(hasExtent && hasBounds)
					{
						for(int ii = 0; ii < 3; ++ii)
							spacing[ii] = (bounds[2*ii+1]-bounds[2*ii])/(extent[2*ii+1]-extent[2*ii]);
						
						field->SetSpacing(spacing);
						field->AllocateScalars();
					}
					
					if(currentLine.find("@1", 0) == 0)
					{
						dataSection = 1;
						continue;
					}
				}
				
				// main data loop
				else
				{
					const char * thisLine = currentLine.c_str();
					double data;
					char ** endptr = new char*;
					data = strtod(thisLine, endptr);
					delete endptr;
					
					int x, y, z;
					z = int(index/(dims[0]*dims[1]));
					y = int(index/dims[0]) - dims[1]*z;
					x = index - (dims[0]*(y + dims[1]*z));
// 					std::flush(std::cout << "x = " << x << " y = " << y << " z = " << z << std::endl);
					
					double * dataPtr = static_cast< double * >(field->GetScalarPointer(x, y, z));
					*dataPtr = data;
					
					++index;
				}
			}
		}
		
		inputStream.close();
		
		field->Update();
		return field;
	}
	else
	{
		std::cout << "Error! Could not open " << inputFilename << std::endl;
		std::cout << "Vector field invalid!" << std::endl;
		return NULL;
	}
}

void Reader::writeScalarField(ImageDataPointerType field)
{
	if(field)
	{
		std::string ofName(outputFilename);
		ofName += ".am";
		std::ofstream outStream(ofName.c_str());
		
		int dims[3], extent[6];
		double bounds[6], spacing[3];
		field->GetDimensions(dims);
		field->GetExtent(extent);
		field->GetBounds(bounds);
		field->GetSpacing(spacing);
		
		outStream << "# AmiraMesh 3D ASCII 2.0" << std::endl;
		outStream << std::endl;
		outStream << "define Lattice ";
		outStream << extent[1]-extent[0]+1 << " ";
		outStream << extent[3]-extent[2]+1 << " ";
		outStream << extent[5]-extent[4]+1 << std::endl;
		outStream << std::endl;
		outStream << "Parameters {" << std::endl;
		outStream << "\tContent \""<< extent[1]-extent[0]+1 << "x" << extent[3]-extent[2]+1 << "x" << extent[5]-extent[4]+1;
		outStream << " float, uniform coordinates\"," << std::endl;
		outStream << "\tBoundingBox ";
		outStream << bounds[0] << " " << bounds[1] << " ";
		outStream << bounds[2] << " " << bounds[3] << " ";
		outStream << bounds[4] << " " << bounds[5] << " ";
		outStream << std::endl;
		outStream << "\tCoordType \"uniform\"" << std::endl;
		outStream << "}" << std::endl;
		outStream << std::endl;
		outStream << "Lattice { float Data } @1" << std::endl;
		outStream << std::endl;
		outStream << "# Data section follows" << std::endl;
		outStream << "@1" << std::endl;
		
//		for(int x = extent[0]; x <= extent[1]; ++x)
//			for(int y = extent[2]; y <= extent[3]; ++y)
//				for(int z = extent[4]; z <= extent[5]; ++z)
		for(int z = extent[4]; z <= extent[5]; ++z)
			for(int y = extent[2]; y <= extent[3]; ++y)
				for(int x = extent[0]; x <= extent[1]; ++x)
				{
					double * vec = static_cast< double * >(field->GetScalarPointer(x, y, z));
					outStream << std::scientific << std::setprecision(15);
					outStream << *vec << std::endl;
					//std::cout << x << "," << y << "," << z << "," << (*vec) << std::endl;
				}
		
		outStream.close();
	}
}

ImageDataPointerType Reader::readVectorField()
{
	ImageDataPointerType field = ImageDataPointerType::New();
	
	std::ifstream inputStream(inputFilename);
	if(!inputStream.fail())
	{
		field->SetNumberOfScalarComponents(3);
		field->SetScalarTypeToDouble();
		
		double spacing[3], bounds[6];
		int extent[6], dims[3], extentDiff[3];
		bool dataSection = 0, hasExtent = 0, hasBounds = 0;
		unsigned long index = 0;
		
		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				std::string::size_type loc1, loc2, loc3;
				
				// set up lattice
				if(!dataSection)
				{
					if(currentLine.find("define") != std::string::npos
						&& currentLine.find("Lattice") != std::string::npos)
					{
						std::size_t startNrs = currentLine.find("Lattice") + 7;
						std::string dimStr = currentLine.substr(startNrs);
						char ** endptr = new char*;
						
						dims[0] = static_cast< int >(strtol(dimStr.c_str(), endptr, 10));
						dims[1] = static_cast< int >(strtol(*endptr, endptr, 10));
						dims[2] = static_cast< int >(strtol(*endptr, endptr, 10));
						
						for(int ii = 0; ii < 3; ++ii)
						{
							extent[2*ii] = 0;
							extent[2*ii+1] = dims[ii] - 1;
							extentDiff[ii] = extent[2*ii+1] - extent[2*ii];
						}
						
						field->SetExtent(extent);
						hasExtent = 1;
						
						delete endptr;
					}
					
					if(currentLine.find("BoundingBox") != std::string::npos)
					{
						std::size_t startNrs = currentLine.find("BoundingBox") + 11;
						std::string bbStr = currentLine.substr(startNrs);
						char ** endptr = new char*;
						
						bounds[0] = strtod(bbStr.c_str(), endptr);
						bounds[1] = strtod(*endptr, endptr);
						bounds[2] = strtod(*endptr, endptr);
						bounds[3] = strtod(*endptr, endptr);
						bounds[4] = strtod(*endptr, endptr);
						bounds[5] = strtod(*endptr, endptr);
						
						field->SetOrigin(bounds[0], bounds[2], bounds[4]);
						hasBounds = 1;
						
						delete endptr;
					}
					
					if(hasExtent && hasBounds)
					{
						for(int ii = 0; ii < 3; ++ii)
							spacing[ii] = (bounds[2*ii+1]-bounds[2*ii])/(extent[2*ii+1]-extent[2*ii]);
						
						field->SetSpacing(spacing);
						field->AllocateScalars();
					}
					
					if(currentLine.find("@1", 0) == 0)
					{
						dataSection = 1;
						continue;
					}
				}
				
				// main data loop
				else
				{
					const char * thisLine = currentLine.c_str();
					double data[3];
					char ** endptr = new char*;
					data[0] = strtod(thisLine, endptr);
					data[1] = strtod(*endptr, endptr);
					data[2] = strtod(*endptr, endptr);
					delete endptr;
					
					int x, y, z;
					z = int(index/(dims[0]*dims[1]));
					y = int(index/dims[0]) - dims[1]*z;
					x = index - (dims[0]*(y + dims[1]*z));
// 					std::flush(std::cout << "x = " << x << " y = " << y << " z = " << z << std::endl);
					
					double * dataPtr = static_cast< double * >(field->GetScalarPointer(x, y, z));
					dataPtr[0] = data[0];
					dataPtr[1] = data[1];
					dataPtr[2] = data[2];
					
					++index;
				}
			}
		}
		
		inputStream.close();
		
		field->Update();
		return field;
	}
	else
	{
		std::cout << "Error! Could not open " << inputFilename << std::endl;
		std::cout << "Vector field invalid!" << std::endl;
		return NULL;
	}
}

ImageDataPointerType Reader::readNVectorField()
{
	ImageDataPointerType field = ImageDataPointerType::New();
	
	std::ifstream inputStream(inputFilename);
	if(!inputStream.fail())
	{
		field->SetScalarTypeToDouble();
		
		double spacing[3], bounds[6];
		int extent[6], dims[3], extentDiff[3], numScalarComponents;
		bool dataSection = 0, hasExtent = 0, hasBounds = 0;
		unsigned long index = 0;
		
		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				std::string::size_type loc1, loc2, loc3;
				
				// set up lattice
				if(!dataSection)
				{
					if(currentLine.find("define") != std::string::npos
						&& currentLine.find("Lattice") != std::string::npos)
					{
						std::size_t startNrs = currentLine.find("Lattice") + 7;
						std::string dimStr = currentLine.substr(startNrs);
						char ** endptr = new char*;
						
						dims[0] = static_cast< int >(strtol(dimStr.c_str(), endptr, 10));
						dims[1] = static_cast< int >(strtol(*endptr, endptr, 10));
						dims[2] = static_cast< int >(strtol(*endptr, endptr, 10));
						
						for(int ii = 0; ii < 3; ++ii)
						{
							extent[2*ii] = 0;
							extent[2*ii+1] = dims[ii] - 1;
							extentDiff[ii] = extent[2*ii+1] - extent[2*ii];
						}
						
						field->SetExtent(extent);
						hasExtent = 1;
						
						delete endptr;
					}
					
					if(currentLine.find("BoundingBox") != std::string::npos)
					{
						std::size_t startNrs = currentLine.find("BoundingBox") + 11;
						std::string bbStr = currentLine.substr(startNrs);
						char ** endptr = new char*;
						
						bounds[0] = strtod(bbStr.c_str(), endptr);
						bounds[1] = strtod(*endptr, endptr);
						bounds[2] = strtod(*endptr, endptr);
						bounds[3] = strtod(*endptr, endptr);
						bounds[4] = strtod(*endptr, endptr);
						bounds[5] = strtod(*endptr, endptr);
						
						field->SetOrigin(bounds[0], bounds[2], bounds[4]);
						hasBounds = 1;
						
						delete endptr;
					}
					
					if(currentLine.find("float") != std::string::npos && currentLine.find("Lattice") != std::string::npos)
					{
						std::size_t startNrs = currentLine.find("float[") + 6;
						std::size_t endNrs = currentLine.find("] Data")-1;
						std::string scalarCompStr = currentLine.substr(startNrs,endNrs);
						
						char ** endptr = new char*;

						numScalarComponents = static_cast< int >(strtod(scalarCompStr.c_str(), endptr));
						field->SetNumberOfScalarComponents(numScalarComponents);
						
						delete endptr;
					}
					
					if(hasExtent && hasBounds)
					{
						for(int ii = 0; ii < 3; ++ii)
							spacing[ii] = (bounds[2*ii+1]-bounds[2*ii])/(extent[2*ii+1]-extent[2*ii]);
						
						field->SetSpacing(spacing);
						field->AllocateScalars();
					}
					
					if(currentLine.find("@1", 0) == 0)
					{
						dataSection = 1;
						continue;
					}
				}
				
				// main data loop
				else
				{
					const char * thisLine = currentLine.c_str();
					double data[numScalarComponents];
					char ** endptr = new char*;
					
					data[0] = strtod(thisLine, endptr);
					for (int i=1; i<numScalarComponents; i++)
					{
						data[i] = strtod(*endptr, endptr);
					}
					delete endptr;
					
					int x, y, z;
					z = int(index/(dims[0]*dims[1]));
					y = int(index/dims[0]) - dims[1]*z;
					x = index - (dims[0]*(y + dims[1]*z));
// 					std::flush(std::cout << "x = " << x << " y = " << y << " z = " << z << std::endl);
					
					double * dataPtr = static_cast< double * >(field->GetScalarPointer(x, y, z));
					for (int i=0; i<numScalarComponents; i++)
					{
						dataPtr[i] = data[i];
					}
					++index;
				}
			}
		}
		
		inputStream.close();
		
		field->Update();
		return field;
	}
	else
	{
		std::cout << "Error! Could not open " << inputFilename << std::endl;
		std::cout << "Vector field invalid!" << std::endl;
		return NULL;
	}
}

void Reader::writeVectorField ( ImageDataPointerType field )
{
	if(field)
	{
		std::string ofName(outputFilename);
		ofName += ".am";
		std::ofstream outStream(ofName.c_str());
		
		int dims[3], extent[6];
		double bounds[6], spacing[3];
		field->GetDimensions(dims);
		field->GetExtent(extent);
		field->GetBounds(bounds);
		field->GetSpacing(spacing);
		
		outStream << "# AmiraMesh 3D ASCII 2.0" << std::endl;
		outStream << std::endl;
		outStream << "define Lattice ";
		outStream << extent[1]-extent[0]+1 << " ";
		outStream << extent[3]-extent[2]+1 << " ";
		outStream << extent[5]-extent[4]+1 << std::endl;
		outStream << std::endl;
		outStream << "Parameters {" << std::endl;
		outStream << "\tContent \""<< extent[1]-extent[0]+1 << "x" << extent[3]-extent[2]+1 << "x" << extent[5]-extent[4]+1;
		outStream << " float[3], uniform coordinates\"," << std::endl;
		outStream << "\tBoundingBox ";
		outStream << bounds[0] << " " << bounds[1] << " ";
		outStream << bounds[2] << " " << bounds[3] << " ";
		outStream << bounds[4] << " " << bounds[5] << " ";
		outStream << std::endl;
		outStream << "\tCoordType \"uniform\"" << std::endl;
		outStream << "}" << std::endl;
		outStream << std::endl;
		outStream << "Lattice { float[3] Data } @1" << std::endl;
		outStream << std::endl;
		outStream << "# Data section follows" << std::endl;
		outStream << "@1" << std::endl;
		
//		for(int x = extent[0]; x <= extent[1]; ++x)
//			for(int y = extent[2]; y <= extent[3]; ++y)
//				for(int z = extent[4]; z <= extent[5]; ++z)
		for(int z = extent[4]; z <= extent[5]; ++z)
			for(int y = extent[2]; y <= extent[3]; ++y)
				for(int x = extent[0]; x <= extent[1]; ++x)
				{
					double * vec = static_cast< double * >(field->GetScalarPointer(x, y, z));
					outStream << std::scientific << std::setprecision(15);
					outStream << vec[0] << " " << vec[1] << " " << vec[2] << " " << std::endl;
				}
		
		outStream.close();
	}
}

void Reader::writeNVectorField ( ImageDataPointerType field )
{
	if(field)
	{
		std::string ofName(outputFilename);
		ofName += ".am";
		std::ofstream outStream(ofName.c_str());
		
		int dims[3], extent[6];
		double bounds[6], spacing[3];
		field->GetDimensions(dims);
		field->GetExtent(extent);
		field->GetBounds(bounds);
		field->GetSpacing(spacing);
		
		int numScalarComponents = field->GetNumberOfScalarComponents(); 
		
		outStream << "# AmiraMesh 3D ASCII 2.0" << std::endl;
		outStream << std::endl;
		outStream << "define Lattice ";
		outStream << extent[1]-extent[0]+1 << " ";
		outStream << extent[3]-extent[2]+1 << " ";
		outStream << extent[5]-extent[4]+1 << std::endl;
		outStream << std::endl;
		outStream << "Parameters {" << std::endl;
		outStream << "\tContent \""<< extent[1]-extent[0]+1 << "x" << extent[3]-extent[2]+1 << "x" << extent[5]-extent[4]+1;
		outStream << " float[" << numScalarComponents << "], uniform coordinates\"," << std::endl;
		outStream << "\tBoundingBox ";
		outStream << bounds[0] << " " << bounds[1] << " ";
		outStream << bounds[2] << " " << bounds[3] << " ";
		outStream << bounds[4] << " " << bounds[5] << " ";
		outStream << std::endl;
		outStream << "\tCoordType \"uniform\"" << std::endl;
		outStream << "}" << std::endl;
		outStream << std::endl;
		outStream << "Lattice { float[" << numScalarComponents << "] Data } @1" << std::endl;
		outStream << std::endl;
		outStream << "# Data section follows" << std::endl;
		outStream << "@1" << std::endl;
		
		for(int z = extent[4]; z <= extent[5]; ++z)
			for(int y = extent[2]; y <= extent[3]; ++y)
				for(int x = extent[0]; x <= extent[1]; ++x)
				{
					double * vec = static_cast< double * >(field->GetScalarPointer(x, y, z));
					outStream << std::scientific << std::setprecision(15);
					for (int s = 0; s<numScalarComponents; s++)
						outStream << vec[s] << " ";
					outStream << "" << std::endl;
				}
		
		outStream.close();
	}
}

void Reader::readConnectionMatrix(ConnectionMatrix* connectome)
{
	std::ifstream inputStream(inputFilename);
	if(!inputStream.fail())
	{
// 		std::cout << "Loading connection matrix from file " << inputFilename << " ..." << std::endl;
		
		/* Time */
		std::time_t start = std::time(0);
		std::time_t prev = start;
		std::cout << "Start: " << std::asctime(std::localtime(&start));
		int counter = 0;
		/* Time */

		int matrixVersion = 0;
		bool readCellIDs = false;
		bool readCellIDsDone = false;
		bool readConnections = false;
// 		HxConnectionMatrix V2 stuff
		bool readPostConnectionNumbers = false;
		unsigned int numberOfPostCells = 0;
		unsigned int postCellIndex = 0;
		std::list< unsigned int > postCellOrder;
		std::list< unsigned int > connectionsPerPostCell;
		std::list< unsigned int >::const_iterator postCellOrderIt;
		std::list< unsigned int >::const_iterator connectionsPerPostCellIt;
		unsigned int connectionsPerPostCellIndex = 0;

#ifdef DEBUG
		unsigned int IDCount = 0;
		unsigned int connectionCount = 0;
#endif
		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				if(!matrixVersion)
				{
					if(currentLine.find("HxConnectionMatrix CSV V1") != std::string::npos)
					{
						matrixVersion = 1;
						std::flush(std::cout << "Loading HxConnectionMatrix V1 from file " << inputFilename << " ..." << std::endl);
						continue;
					}
					if(currentLine.find("HxConnectionMatrix CSV V2") != std::string::npos)
					{
						matrixVersion = 2;
						std::flush(std::cout << "Loading HxConnectionMatrix V2 from file " << inputFilename << " ..." << std::endl);
						continue;
					}
				}
				if(!readCellIDs && !readConnections)
				{
					if(matrixVersion == 1 && currentLine.find("NumberOfCells") != std::string::npos)
					{
						readCellIDs = true;
						continue;
					}
					if(matrixVersion == 2 && currentLine.find("NumberOfMetaDataEntries") != std::string::npos)
					{
						readCellIDs = true;
						continue;
					}
				}
				if(matrixVersion == 2 && currentLine.find("NumberOfPostsynapticCells") != std::string::npos)
				{
					readCellIDs = false;
					readPostConnectionNumbers = true;
// 					char * tmpChar = new char[64];
// 					sscanf(currentLine.currentLine(), " %s , %u ", tmpChar, &numberOfPostCells);
					std::size_t delim = currentLine.find(",");
					std::string tmpNumberPostCellsStr = currentLine.substr(delim+1, currentLine.size()-delim-1);
					numberOfPostCells = atoi(tmpNumberPostCellsStr.c_str());
#ifdef DEBUG
					std::flush(std::cout << currentLine.c_str() << std::endl);
// 					std::flush(std::cout << "tmpChar: " << tmpChar << std::endl);
					std::flush(std::cout << "numberOfPostCells: " << numberOfPostCells << std::endl);
// 					delete [] tmpChar;
#endif
					continue;
				}
				if(matrixVersion == 2 && readPostConnectionNumbers && postCellIndex == numberOfPostCells)
				{
					readPostConnectionNumbers = false;
					readConnections = true;
					postCellOrderIt = postCellOrder.begin();
					connectionsPerPostCellIt = connectionsPerPostCell.begin();
				}
				
				if(readCellIDs)
				{
#ifdef DEBUG
					++IDCount;
#endif
					readCellIDsDone = true;
					unsigned int cellID, cellType, column;
					char * tmpChar = new char[64];
					sscanf(currentLine.c_str(), " %u,%s ", &cellID, tmpChar);
					std::string tmpStr(tmpChar);
					std::size_t delim = tmpStr.find(",");
					std::string tmpCellTypeStr = tmpStr.substr(0, delim);
					std::string tmpColumnStr = tmpStr.substr(delim+1, tmpStr.size()-delim-1);
					cellType = celltypeLabels2Int[tmpCellTypeStr];
					column = labels2Int[tmpColumnStr];
#ifdef DEBUG
					if(!(IDCount%1000))
					{
						std::flush(std::cout << "IDcount: " << IDCount << " cell ID: " << cellID << std::endl);
						std::flush(std::cout << "tmpChar: " << tmpChar << " delim: " << delim << std::endl);
						std::flush(std::cout << "tmp cell type: " << tmpCellTypeStr.c_str() << " tmp column: " << tmpColumnStr.c_str() << std::endl);
						std::flush(std::cout << "cell type: " << cellType << " column: " << column << std::endl);
					}
#endif
					delete [] tmpChar;
					
					connectome->IDColumnCelltypeMap.insert(std::pair< unsigned int, std::pair< unsigned int, unsigned int > >(cellID, std::pair< unsigned int, unsigned int >(column, cellType)));
					
					// check if we already have this cell type/column stored...
					if(std::find(preTypes.begin(), preTypes.end(), cellType) != preTypes.end())
					{
						if(std::find(connectome->preTypes.begin(), connectome->preTypes.end(), cellType) == connectome->preTypes.end())
						{
							connectome->preTypes.push_back(cellType);
							std::vector< unsigned int > IDVec;
							connectome->preTypeIDs.insert(std::pair< unsigned int, std::vector< unsigned int> >(cellType, IDVec));
							connectome->preTypeNumbers.insert(std::pair< unsigned int, unsigned int >(cellType, 0));
						}
						if(std::find(barrelLabels.begin(), barrelLabels.end(), column) != barrelLabels.end())
						{
							if(std::find(connectome->preColumns.begin(), connectome->preColumns.end(), column) ==  connectome->preColumns.end())
							{
								connectome->preColumns.push_back(column);
								std::vector< unsigned int > IDVec;
								connectome->preColumnIDs.insert(std::pair< unsigned int, std::vector< unsigned int> >(column, IDVec));
								connectome->preColumnNumbers.insert(std::pair< unsigned int, unsigned int >(column, 0));
								std::map< unsigned int, unsigned int > cellTypeNrMap;
								connectome->columnCelltypePreNumbers.insert(std::pair< unsigned int, std::map< unsigned int, unsigned int > >(column, cellTypeNrMap));
								connectome->columnCelltypePreNumbers[column].insert(std::pair< unsigned int, unsigned int >(cellType, 0));
							}
							
							connectome->preTypeIDs[cellType].push_back(cellID);
							connectome->preTypeNumbers[cellType]++;
							connectome->preColumnIDs[column].push_back(cellID);
							connectome->preColumnNumbers[column]++;
							connectome->columnCelltypePreNumbers[column][cellType]++;
						}
					}
					if(std::find(postTypes.begin(), postTypes.end(), cellType) != postTypes.end())
					{
						if(std::find(connectome->postTypes.begin(), connectome->postTypes.end(), cellType) == connectome->postTypes.end())
						{
							connectome->postTypes.push_back(cellType);
							std::vector< unsigned int > IDVec;
							connectome->postTypeIDs.insert(std::pair< unsigned int, std::vector< unsigned int> >(cellType, IDVec));
							connectome->postTypeNumbers.insert(std::pair< unsigned int, unsigned int >(cellType, 0));
						}
						if(std::find(barrelLabels.begin(), barrelLabels.end(), column) != barrelLabels.end())
						{
							if(std::find(connectome->postColumns.begin(), connectome->postColumns.end(), column) == connectome->postColumns.end())
							{
								connectome->postColumns.push_back(column);
								std::vector< unsigned int > IDVec;
								connectome->postColumnIDs.insert(std::pair< unsigned int, std::vector< unsigned int> >(column, IDVec));
								connectome->postColumnNumbers.insert(std::pair< unsigned int, unsigned int >(column, 0));
								std::map< unsigned int, unsigned int > cellTypeNrMap;
								connectome->columnCelltypePostNumbers.insert(std::pair< unsigned int, std::map< unsigned int, unsigned int > >(column, cellTypeNrMap));
								connectome->columnCelltypePostNumbers[column].insert(std::pair< unsigned int, unsigned int >(cellType, 0));
							}
							
							connectome->postTypeIDs[cellType].push_back(cellID);
							connectome->postTypeNumbers[cellType]++;
							connectome->postColumnIDs[column].push_back(cellID);
							connectome->postColumnNumbers[column]++;
							connectome->columnCelltypePostNumbers[column][cellType]++;
						}
					}
				}
				if(readPostConnectionNumbers && matrixVersion == 2)
				{
					unsigned int postCellID;
					unsigned int preSynapticCells;
					sscanf(currentLine.c_str(), " %u , %u  ", &postCellID, &preSynapticCells);
					postCellOrder.push_back(postCellID);
					connectionsPerPostCell.push_back(preSynapticCells);
					++postCellIndex;
#ifdef DEBUG
					if(!(postCellIndex%1000))
					{
						std::flush(std::cout << currentLine.c_str() << std::endl);
						std::flush(std::cout << "postCellIndex: " << postCellIndex << " postCellID: " << postCellID << " preSynapticCells: " << preSynapticCells << std::endl);
					}
#endif
				}
				if(readConnections && matrixVersion == 1)
				{
#ifdef DEBUG
					++connectionCount;
#endif
					unsigned int preCellID;
					unsigned int postCellID;
					float innervation;
					sscanf(currentLine.c_str(), " %u , %u , %f ", &preCellID, &postCellID, &innervation);
#ifdef DEBUG
					if(!(connectionCount%100000))
					{
						std::flush(std::cout << "preCellID: " << preCellID << " postCellID: " << postCellID << " innervation: " << innervation << std::endl);
					}
#endif
					connectome->matrix.insert(std::pair< MatrixIndexType, float >(MatrixIndexType(preCellID, postCellID), innervation));
				}
				if(readConnections && matrixVersion == 2)
				{
#ifdef DEBUG
					++connectionCount;
#endif
					if(postCellOrderIt == postCellOrder.end() || connectionsPerPostCellIt == connectionsPerPostCell.end())
					{
						std::flush(std::cout << "Error: trying to read invalid postsynaptic cell connections!" << std::endl);
						connectome = NULL;
						return;
					}
					if(connectionsPerPostCellIndex >= *connectionsPerPostCellIt)
					{
						++postCellOrderIt;
						++connectionsPerPostCellIt;
						connectionsPerPostCellIndex = 0;
					}
					while(*connectionsPerPostCellIt == 0)
					{
						++postCellOrderIt;
						++connectionsPerPostCellIt;
					}
					unsigned int postCellID = *postCellOrderIt;
					unsigned int preCellID;
					float innervation;
					sscanf(currentLine.c_str(), " %u , %f ", &preCellID, &innervation);
					
					unsigned int postCellType = connectome->IDColumnCelltypeMap[postCellID].second;
					std::string postCellTypeString(int2CelltypeLabels[postCellType]);
					if(postCellTypeString.find("axon") != std::string::npos || postCellTypeString.find("VPM") != std::string::npos)
					{
						std::flush(std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl);
						std::flush(std::cout << "ERROR: cell of type " << postCellTypeString.c_str() << " set as postsynaptic!" << std::endl);
						std::flush(std::cout << "preCellID: " << preCellID << " postCellID: " << postCellID << " innervation: " << innervation << std::endl);
						std::flush(std::cout << "connectionsPerPostCellIndex: " << connectionsPerPostCellIndex << std::endl);
						std::flush(std::cout << "*connectionsPerPostCellIt: " << *connectionsPerPostCellIt << std::endl);
						std::flush(std::cout << "*postCellOrderIt: " << *postCellOrderIt << std::endl);
						std::flush(std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl);
						return;
					}
#ifdef DEBUG
					if(!(connectionCount%100000))
					{
						std::flush(std::cout << "preCellID: " << preCellID << " postCellID: " << postCellID << " innervation: " << innervation << std::endl);
						std::flush(std::cout << "connectionsPerPostCellIndex: " << connectionsPerPostCellIndex << std::endl);
						std::flush(std::cout << "*connectionsPerPostCellIt: " << *connectionsPerPostCellIt << std::endl);
						std::flush(std::cout << "*postCellOrderIt: " << *postCellOrderIt << std::endl);
					}
#endif
					connectome->matrix.insert(std::pair< MatrixIndexType, float >(MatrixIndexType(preCellID, postCellID), innervation));
					
					++connectionsPerPostCellIndex;
				}
			}
			else if(!currentLine.size() && readCellIDsDone && !readConnections && matrixVersion == 1)
			{
				readCellIDs = false;
				readConnections = true;
#ifdef DEBUG
				std::cout << "Done reading " << IDCount << " cell IDs..." << std::endl;
#endif
			}

			/* Time */
			std::time_t tmp = std::time(0);
			counter++;

			if (tmp>prev)
			{
				if ((tmp-start) % 60 == 0)
				{
					std::cout << "Time passed: " << tmp-start << " Line: " << counter << std::endl;
					// std::asctime(std::localtime(&start))
				}
				if ((tmp-start) % 3600 == 0)
				{
					std::cout << "Start time: " << std::asctime(std::localtime(&start));
					std::cout << "Current time: " << std::asctime(std::localtime(&tmp));
				}
			}
			prev = tmp;
			/* Time */
		}
#ifdef DEBUG
		std::cout << "Done reading " << connectionCount << " cell connections..." << std::endl;
		std::cout << "Column\tcell type\tnumber" << std::endl;
		std::map< unsigned int, std::map< unsigned int, unsigned int > >::const_iterator colIt;
		for(colIt = connectome->columnCelltypePreNumbers.begin(); colIt != connectome->columnCelltypePreNumbers.end(); ++colIt)
		{
			std::map< unsigned int, unsigned int >::const_iterator typeNrIt;
			for(typeNrIt = colIt->second.begin(); typeNrIt != colIt->second.end(); ++typeNrIt)
			{
				std::cout << int2Labels[colIt->first] << "\t" << int2CelltypeLabels[typeNrIt->first] << "\t" << typeNrIt->second << std::endl;
			}
		}
		for(colIt = connectome->columnCelltypePostNumbers.begin(); colIt != connectome->columnCelltypePostNumbers.end(); ++colIt)
		{
			std::map< unsigned int, unsigned int >::const_iterator typeNrIt;
			for(typeNrIt = colIt->second.begin(); typeNrIt != colIt->second.end(); ++typeNrIt)
			{
				std::cout << int2Labels[colIt->first] << "\t" << int2CelltypeLabels[typeNrIt->first] << "\t" << typeNrIt->second << std::endl;
			}
		}
#endif
	/* Time */
	std::time_t ending = std::time(0);

	std::cout << "Finished Reading ..." << std::endl;
	std::cout << "Start time: " << std::asctime(std::localtime(&start));
	std::cout << "Current time: " << std::asctime(std::localtime(&ending));
	std::cout << "Time passed: " << ending-start << " Total Number of Lines: " << counter << std::endl;
	/* Time */
	}
	else
	{
		std::cout << "Error! Could not read " << inputFilename << std::endl;
	}
}

void Reader::readSynapsesPerCellTable(CellTable* table)
{
	std::ifstream inputStream(inputFilename);
	if(!inputStream.fail())
	{
		std::cout << "Reading AmiraMesh synapse table " << inputFilename << " ..." << std::endl;
		
		std::vector< std::vector< float > > dataTranspose;
		std::vector< float > somaX;
		std::vector< float > somaY;
		std::vector< float > somaZ;
		
		bool foundHeader = false;
		bool foundReferences = false;
		std::map< unsigned int, unsigned int > dataReferences;
		bool totalSynFlag = false;
		const unsigned int synColumnStart = 8;
		int currentColReadFlag = -1;
		int currentRefFlag = -1;
		unsigned int currentRow = 0;
		std::string tmpStr = "";
		std::string currentLine;
		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				if(currentLine.find("define") != std::string::npos)
				{
					continue;
				}
				if(!foundHeader)
				{
					if(currentLine.find("}") != std::string::npos)
					{
						// only checking TOTAL_SYNAPSES table,
						// not AVERAGE_SYNAPSES table
						foundHeader = true;
					}
					if(currentLine.find("NumRows") != std::string::npos)
					{
// 						std::size_t start = currentLine.find("NumRows");
						std::size_t start = currentLine.find_first_of(numbers);
						std::size_t delim = currentLine.find(",");
						std::string tmpNumberPostCellsStr = currentLine.substr(start, delim-start);
						unsigned int numberOfRows = atoi(tmpNumberPostCellsStr.c_str());
#ifdef DEBUG
						std::flush(std::cout << currentLine.c_str() << std::endl);
						std::flush(std::cout << "numberOfRows: " << numberOfRows << std::endl);
#endif
						for(int i = 0; i < numberOfRows; ++i)
						{
							CellTableRow * newRow = new CellTableRow;
							table->rows.push_back(newRow);
						}
					}
					if(currentLine.find("Column") != std::string::npos)
					{
						size_t delim1 = currentLine.find_first_of(numbers);
						size_t delim2 = currentLine.find_first_of("\"");
						size_t delim3 = currentLine.find_last_of("\"");
						std::string columnIDStr = currentLine.substr(delim1, 4);
						unsigned int columnID = atoi(columnIDStr.c_str());
						std::string attributeStr = currentLine.substr(delim2+1, delim3-delim2-1);
#ifdef DEBUG
						std::flush(std::cout << currentLine.c_str() << std::endl);
						std::flush(std::cout << "columnID: " << columnID << std::endl);
						std::flush(std::cout << "attribute: " << attributeStr.c_str() << std::endl);
#endif
						
// 						if(columnID == 7)
// 						{
// 							if(currentLine.find("TOTAL") != std::string::npos)
// 							{
// 								// all good
// 							}
// 						}
						if(columnID >= synColumnStart)
						{
							size_t delim = attributeStr.find_first_of("_");
							std::string columnStr = attributeStr.substr(0, delim);
							std::string tmpStr = attributeStr.substr(delim+1, attributeStr.size()-delim-1);
							// cut off axon at the end
							std::string cellTypeStr;
							if(tmpStr.find("axon") != std::string::npos)
							{
								size_t axonDelim = tmpStr.find("axon");
								cellTypeStr = tmpStr.substr(0, axonDelim);
// #ifdef DEBUG
// 								std::flush(std::cout << "axonDelim: " << axonDelim << std::endl);
// 								std::flush(std::cout << "tmpStr.size(): " << tmpStr.size() << std::endl);
// #endif
							}
							else
							{
								cellTypeStr = tmpStr;
// #ifdef DEBUG
// 								std::flush(std::cout << "axonDelim: " << axonDelim << std::endl);
// 								std::flush(std::cout << "tmpStr.size(): " << tmpStr.size() << std::endl);
// #endif
							}
#ifdef DEBUG
							std::flush(std::cout << "column: " << columnStr.c_str() << std::endl);
							std::flush(std::cout << "tmpStr: " << tmpStr.c_str() << std::endl);
							std::flush(std::cout << "cellType: " << cellTypeStr.c_str() << std::endl);
#endif
							unsigned int column = labels2Int[columnStr];
							unsigned int cellType = celltypeLabels2Int[cellTypeStr];
							table->header.insert(std::pair< ColumnCellTypePair, unsigned int >(ColumnCellTypePair(column, cellType), columnID-synColumnStart));
							std::vector< float > dataColumn;
							dataTranspose.push_back(dataColumn);
						}
					}
					continue;
				}
				if(!foundReferences)
				{
					if(currentLine.find("Table0000Column") != std::string::npos)
					{
						std::string colStr = currentLine.substr(15, 4);
						size_t delim = currentLine.find("@");
						std::string refStr = currentLine.substr(delim+1, currentLine.size()-delim-1);
						unsigned int colNr = atoi(colStr.c_str());
						unsigned int refNr = atoi(refStr.c_str());
#ifdef DEBUG
						std::flush(std::cout << currentLine.c_str() << std::endl);
						std::flush(std::cout << "colNr: " << colNr << std::endl);
						std::flush(std::cout << "refNr: " << refNr << std::endl);
#endif
						dataReferences.insert(std::pair< unsigned int, unsigned int >(refNr, colNr));
					}
					if(dataReferences.size() == table->header.size() + synColumnStart)
					{
						foundReferences = true;
					}
					continue;
				}
				// reading data begins here
				if(currentLine.find("@") == 0)
				{
					currentRow = 0;
					currentRefFlag = atoi(currentLine.substr(1, currentLine.size()-1).c_str());
#ifdef DEBUG
					std::flush(std::cout << currentLine.c_str() << std::endl);
					std::flush(std::cout << "currentRefFlag: " << currentRefFlag << std::endl);
#endif
					if(dataReferences.find(currentRefFlag) != dataReferences.end())
					{
						currentColReadFlag = dataReferences[currentRefFlag];
#ifdef DEBUG
						std::flush(std::cout << "currentColReadFlag: " << currentColReadFlag << std::endl);
#endif
					}
					continue;
				}
				if(dataReferences.find(currentRefFlag) == dataReferences.end())
				{
					continue;
				}
				// CELLID (unsigned int)
				if(currentColReadFlag == 0)
				{
					table->rows[currentRow]->cellID = atoi(currentLine.c_str());
					++currentRow;
				}
				// CELLTYPE (ASCII char 0-terminated)
				if(currentColReadFlag == 1)
				{
					if(atoi(currentLine.c_str()) == 0)
					{
						table->rows[currentRow]->cellType = celltypeLabels2Int[tmpStr];
						tmpStr = "";
						++currentRow;
					}
					else
					{
						tmpStr += (char)atoi(currentLine.c_str());
					}
				}
				// COLUMN (ASCII char 0-terminated)
				if(currentColReadFlag == 2)
				{
					if(atoi(currentLine.c_str()) == 0)
					{
						table->rows[currentRow]->column = labels2Int[tmpStr];
						tmpStr = "";
						++currentRow;
					}
					else
					{
						tmpStr += (char)atoi(currentLine.c_str());
					}
				}
				// SOMA_X (float)
				if(currentColReadFlag == 3)
				{
					somaX.push_back(atof(currentLine.c_str()));
				}
				// SOMA_Y (float)
				if(currentColReadFlag == 4)
				{
					somaY.push_back(atof(currentLine.c_str()));
				}
				// SOMA_Z (float)
				if(currentColReadFlag == 5)
				{
					somaZ.push_back(atof(currentLine.c_str()));
				}
				// INSIDE_COLUMN (bool)
				if(currentColReadFlag == 6)
				{
					table->rows[currentRow]->insideColumn = atoi(currentLine.c_str());
					++currentRow;
				}
				// TOTAL_SYNAPSES (float)
				if(currentColReadFlag == 7)
				{
					table->rows[currentRow]->totalSynapses = atof(currentLine.c_str());
					++currentRow;
				}
// 				if(currentColReadFlag == 7 && totalSynFlag)
// 				{
// 					
// 				}
// 				if(currentColReadFlag == 7 && !totalSynFlag)
// 				{
// 					
// 				}
				// synapses per pre column cell type (float)
				if(currentColReadFlag > 7)
				{
					unsigned int colIndex = currentColReadFlag - synColumnStart;
					dataTranspose[colIndex].push_back(atof(currentLine.c_str()));
				}
			}
		}
		
#ifdef DEBUG
		std::flush(std::cout << "somaX.size(): " << somaX.size() << std::endl);
		std::flush(std::cout << "somaY.size(): " << somaY.size() << std::endl);
		std::flush(std::cout << "somaZ.size(): " << somaZ.size() << std::endl);
		std::flush(std::cout << "dataTranspose[0].size(): " << dataTranspose[0].size() << std::endl);
		std::flush(std::cout << "table->rows.size(): " << table->rows.size() << std::endl);
#endif
		for(int i = 0; i < somaX.size(); ++i)
		{
			std::vector< float > somaLoc;
			somaLoc.push_back(somaX[i]);
			somaLoc.push_back(somaY[i]);
			somaLoc.push_back(somaZ[i]);
			table->rows[i]->somaLocation = somaLoc;
		}
		for(int i = 0; i < table->rows.size(); ++i)
		{
			for(int j = 0; j < dataTranspose.size(); ++j)
			{
				table->rows[i]->synapsesPerPreTypeColumn.push_back(dataTranspose[j][i]);
			}
		}
		
		somaX.clear();
		somaY.clear();
		somaZ.clear();
		dataTranspose.clear();
	}
}

/* read NeuroNet convergence matrix tables
 * Skips vpm and celltypes containing "axon" on postsynaptic site
 * Uses table structe, totalSynapses is set to -1, synapsesPerPreTypeColumn contains convergence values
 */
void Reader::readConvergenceTableCSV(CellTable * table)
{
	std::ifstream inputStream(inputFilename);
	if(!inputStream.fail())
	{
		std::cout << "Reading Convergence csv table " << inputFilename << " ..." << std::endl;

		bool foundHeader = false;
		std::string currentLine;

		while(!std::getline(inputStream, currentLine).eof())
		{
			if(currentLine.size())
			{
				if (currentLine.find("CONVERGENCE") != std::string::npos)
				{
					continue;
				}

				if (!foundHeader)
				{
					// Initialize header
					if (currentLine.find("CELLID") != std::string::npos)
					{
						std::size_t start = currentLine.find("INSIDE_COLUMN");
						std::string tmpStr = currentLine.substr(start+14);
						std::size_t delim = tmpStr.find(",");
						std::string tmp;
						std::string col;
						std::string celltype;
						unsigned int columnID = 0;

						while(delim != std::string::npos)
						{
							tmp = tmpStr.substr(0,delim);
							tmpStr = tmpStr.substr(delim+1);
							delim = tmpStr.find(",");

							col = tmp.substr(0,tmp.find("_"));
							celltype = tmp.substr(tmp.find("_")+1);

							table->header.insert(std::pair< ColumnCellTypePair, unsigned int >(ColumnCellTypePair(labels2Int[col], celltypeLabels2Int[celltype]), columnID));
							columnID++;
							//flush(std::cout << col << " " << celltype << std::endl);
						}

						col = tmpStr.substr(0,tmpStr.find("_"));
						celltype = tmpStr.substr(tmpStr.find("_")+1);
						table->header.insert(std::pair< ColumnCellTypePair, unsigned int >(ColumnCellTypePair(labels2Int[col], celltypeLabels2Int[celltype]), columnID));
						//flush(std::cout << col << " " << celltype << std::endl);

						foundHeader = true;
					}
				}

				// Data
				if (foundHeader)
				{
					// Skip if axon
					if(currentLine.find("axon") != std::string::npos)
					{
						continue;
					}

					// Skip if VPM
					if(currentLine.find("VPM") != std::string::npos)
					{
						continue;
					}

					CellTableRow * newRow = new CellTableRow;
					std::size_t delim = currentLine.find(",");

					// CellID (int)
					std::string tmp = currentLine.substr(0,delim);
					newRow->cellID = atoi(tmp.c_str());
					//std::cout << tmp << ",";

					// CellType (int)
					currentLine= currentLine.substr(delim+1);
					delim = currentLine.find(",");
					tmp = currentLine.substr(0,delim);
					newRow->cellType = celltypeLabels2Int[tmp];
					//std::cout << tmp << ",";

					// Column (int)
					currentLine= currentLine.substr(delim+1);
					delim = currentLine.find(",");
					tmp = currentLine.substr(0,delim);
					newRow->column = labels2Int[tmp];
					//std::cout << tmp << ",";

					// Soma Position (float)
					std::vector< float > somaLocation;
					// X
					currentLine= currentLine.substr(delim+1);
					delim = currentLine.find(",");
					tmp = currentLine.substr(0,delim);
					somaLocation.push_back(atof(tmp.c_str()));
					//std::cout << tmp << ",";
					// Y
					currentLine= currentLine.substr(delim+1);
					delim = currentLine.find(",");
					tmp = currentLine.substr(0,delim);
					somaLocation.push_back(atof(tmp.c_str()));
					//std::cout << tmp << ",";
					// Z
					currentLine= currentLine.substr(delim+1);
					delim = currentLine.find(",");
					tmp = currentLine.substr(0,delim);
					somaLocation.push_back(atof(tmp.c_str()));
					//std::cout << tmp << ",";
					newRow->somaLocation = somaLocation;

					// Inside Column (bool)
					currentLine= currentLine.substr(delim+1);
					delim = currentLine.find(",");
					tmp = currentLine.substr(0,delim);
					newRow->insideColumn = atoi(tmp.c_str());
					//std::cout << tmp << ",";

					// synapsesPerPreTypeColumn // probability (float)
					std::vector< float > convergenceColumn;

					while(delim != std::string::npos)
					{
						currentLine = currentLine.substr(delim+1);
						delim = currentLine.find(",");
						tmp = currentLine.substr(0,delim);
						convergenceColumn.push_back(atof(tmp.c_str()));
						//std::cout << tmp << ",";
					}
					//std::cout << "" << std::endl;

					newRow->synapsesPerPreTypeColumn = convergenceColumn;
					newRow-> totalSynapses = -1;
					table->rows.push_back(newRow);

					if (convergenceColumn.size() != table->header.size())
					{
						std::cout << "ERROR! Something is wrong, header size (" << table->header.size() << ") does not match size of convergence values (" << convergenceColumn.size() << ")" << std::endl;
						break;
					}
				}
			}
		}

		if (!foundHeader)
		{
			std::cout << "No Header containing \"CELLID\" found, file format might be not correct!" << std::endl;
		}

	}
	else
	{
		std::cout << "ERROR! Reading Convergence csv table " << inputFilename << " failed!" << std::endl;
	}
}


void Reader::initializeConstants()
{
	if(this->barrelLabels.size())
		this->barrelLabels.clear();
	this->barrelLabels.push_back(Alpha);
	this->barrelLabels.push_back(A1);
	this->barrelLabels.push_back(A2);
	this->barrelLabels.push_back(A3);
	this->barrelLabels.push_back(A4);
	this->barrelLabels.push_back(Beta);
	this->barrelLabels.push_back(B1);
	this->barrelLabels.push_back(B2);
	this->barrelLabels.push_back(B3);
	this->barrelLabels.push_back(B4);
	this->barrelLabels.push_back(Gamma);
	this->barrelLabels.push_back(C1);
	this->barrelLabels.push_back(C2);
	this->barrelLabels.push_back(C3);
	this->barrelLabels.push_back(C4);
	this->barrelLabels.push_back(C5);
	this->barrelLabels.push_back(C6);
	this->barrelLabels.push_back(Delta);
	this->barrelLabels.push_back(D1);
	this->barrelLabels.push_back(D2);
	this->barrelLabels.push_back(D3);
	this->barrelLabels.push_back(D4);
	this->barrelLabels.push_back(D5);
	this->barrelLabels.push_back(D6);
	this->barrelLabels.push_back(E1);
	this->barrelLabels.push_back(E2);
	this->barrelLabels.push_back(E3);
	this->barrelLabels.push_back(E4);
	this->barrelLabels.push_back(E5);
	this->barrelLabels.push_back(E6);
	if(this->int2Labels.size())
		this->int2Labels.clear();
	this->int2Labels.insert(std::pair< int, const char * >(Alpha, "Alpha"));
	this->int2Labels.insert(std::pair< int, const char * >(A1, "A1"));
	this->int2Labels.insert(std::pair< int, const char * >(A2, "A2"));
	this->int2Labels.insert(std::pair< int, const char * >(A3, "A3"));
	this->int2Labels.insert(std::pair< int, const char * >(A4, "A4"));
	this->int2Labels.insert(std::pair< int, const char * >(Beta, "Beta"));
	this->int2Labels.insert(std::pair< int, const char * >(B1, "B1"));
	this->int2Labels.insert(std::pair< int, const char * >(B2, "B2"));
	this->int2Labels.insert(std::pair< int, const char * >(B3, "B3"));
	this->int2Labels.insert(std::pair< int, const char * >(B4, "B4"));
	this->int2Labels.insert(std::pair< int, const char * >(Gamma, "Gamma"));
	this->int2Labels.insert(std::pair< int, const char * >(C1, "C1"));
	this->int2Labels.insert(std::pair< int, const char * >(C2, "C2"));
	this->int2Labels.insert(std::pair< int, const char * >(C3, "C3"));
	this->int2Labels.insert(std::pair< int, const char * >(C4, "C4"));
	this->int2Labels.insert(std::pair< int, const char * >(C5, "C5"));
	this->int2Labels.insert(std::pair< int, const char * >(C6, "C6"));
	this->int2Labels.insert(std::pair< int, const char * >(Delta, "Delta"));
	this->int2Labels.insert(std::pair< int, const char * >(D1, "D1"));
	this->int2Labels.insert(std::pair< int, const char * >(D2, "D2"));
	this->int2Labels.insert(std::pair< int, const char * >(D3, "D3"));
	this->int2Labels.insert(std::pair< int, const char * >(D4, "D4"));
	this->int2Labels.insert(std::pair< int, const char * >(D5, "D5"));
	this->int2Labels.insert(std::pair< int, const char * >(D6, "D6"));
	this->int2Labels.insert(std::pair< int, const char * >(E1, "E1"));
	this->int2Labels.insert(std::pair< int, const char * >(E2, "E2"));
	this->int2Labels.insert(std::pair< int, const char * >(E3, "E3"));
	this->int2Labels.insert(std::pair< int, const char * >(E4, "E4"));
	this->int2Labels.insert(std::pair< int, const char * >(E5, "E5"));
	this->int2Labels.insert(std::pair< int, const char * >(E6, "E6"));
	if(this->labels2Int.size())
		this->labels2Int.clear();
	this->labels2Int.insert(std::pair< std::string, int >(std::string("Alpha"), Alpha));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("A1"), A1));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("A2"), A2));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("A3"), A3));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("A4"), A4));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("Beta"), Beta));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("B1"), B1));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("B2"), B2));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("B3"), B3));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("B4"), B4));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("Gamma"), Gamma));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("C1"), C1));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("C2"), C2));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("C3"), C3));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("C4"), C4));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("C5"), C5));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("C6"), C6));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("Delta"), Delta));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("D1"), D1));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("D2"), D2));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("D3"), D3));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("D4"), D4));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("D5"), D5));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("D6"), D6));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("E1"), E1));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("E2"), E2));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("E3"), E3));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("E4"), E4));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("E5"), E5));
	this->labels2Int.insert(std::pair< std::string, int >(std::string("E6"), E6));
	if(this->hocLabels.size())
		this->hocLabels.clear();
	this->hocLabels.push_back("alpha");
	this->hocLabels.push_back("WM");
	this->hocLabels.push_back("soma");
	this->hocLabels.push_back("dend");
	this->hocLabels.push_back("apical");
	this->hocLabels.push_back("axon");
	this->hocLabels.push_back("Alpha");
	this->hocLabels.push_back("A1");
	this->hocLabels.push_back("A2");
	this->hocLabels.push_back("A3");
	this->hocLabels.push_back("A4");
	this->hocLabels.push_back("Beta");
	this->hocLabels.push_back("B1");
	this->hocLabels.push_back("B2");
	this->hocLabels.push_back("B3");
	this->hocLabels.push_back("B4");
	this->hocLabels.push_back("Gamma");
	this->hocLabels.push_back("C1");
	this->hocLabels.push_back("C2");
	this->hocLabels.push_back("C3");
	this->hocLabels.push_back("C4");
	this->hocLabels.push_back("C5");
	this->hocLabels.push_back("C6");
	this->hocLabels.push_back("Delta");
	this->hocLabels.push_back("D1");
	this->hocLabels.push_back("D2");
	this->hocLabels.push_back("D3");
	this->hocLabels.push_back("D4");
	this->hocLabels.push_back("D5");
	this->hocLabels.push_back("D6");
	this->hocLabels.push_back("E1");
	this->hocLabels.push_back("E2");
	this->hocLabels.push_back("E3");
	this->hocLabels.push_back("E4");
	this->hocLabels.push_back("E5");
	this->hocLabels.push_back("E6");
	if(this->hocNeuronLabels.size())
		this->hocNeuronLabels.clear();
	this->hocNeuronLabels.push_back("soma");
	this->hocNeuronLabels.push_back("dend");
	this->hocNeuronLabels.push_back("apical");
	this->hocNeuronLabels.push_back("axon");
	if(this->hocLandmarkLabels.size())
		this->hocLandmarkLabels.clear();
	this->hocLandmarkLabels.push_back("alpha");
	this->hocLandmarkLabels.push_back("WM");
	this->hocLandmarkLabels.push_back("Alpha");
	this->hocLandmarkLabels.push_back("A1");
	this->hocLandmarkLabels.push_back("A2");
	this->hocLandmarkLabels.push_back("A3");
	this->hocLandmarkLabels.push_back("A4");
	this->hocLandmarkLabels.push_back("Beta");
	this->hocLandmarkLabels.push_back("B1");
	this->hocLandmarkLabels.push_back("B2");
	this->hocLandmarkLabels.push_back("B3");
	this->hocLandmarkLabels.push_back("B4");
	this->hocLandmarkLabels.push_back("Gamma");
	this->hocLandmarkLabels.push_back("C1");
	this->hocLandmarkLabels.push_back("C2");
	this->hocLandmarkLabels.push_back("C3");
	this->hocLandmarkLabels.push_back("C4");
	this->hocLandmarkLabels.push_back("C5");
	this->hocLandmarkLabels.push_back("C6");
	this->hocLandmarkLabels.push_back("Delta");
	this->hocLandmarkLabels.push_back("D1");
	this->hocLandmarkLabels.push_back("D2");
	this->hocLandmarkLabels.push_back("D3");
	this->hocLandmarkLabels.push_back("D4");
	this->hocLandmarkLabels.push_back("D5");
	this->hocLandmarkLabels.push_back("D6");
	this->hocLandmarkLabels.push_back("E1");
	this->hocLandmarkLabels.push_back("E2");
	this->hocLandmarkLabels.push_back("E3");
	this->hocLandmarkLabels.push_back("E4");
	this->hocLandmarkLabels.push_back("E5");
	this->hocLandmarkLabels.push_back("E6");
	if(this->postTypes.size())
		this->postTypes.clear();
	this->postTypes.push_back(L2);
	this->postTypes.push_back(L34);
	this->postTypes.push_back(L4py);
	this->postTypes.push_back(L4sp);
	this->postTypes.push_back(L4ss);
	this->postTypes.push_back(L5st);
	this->postTypes.push_back(L5tt);
	this->postTypes.push_back(L6cc);
	this->postTypes.push_back(L6ccinv);
	this->postTypes.push_back(L6ct);
	this->postTypes.push_back(SymLocal);
	this->postTypes.push_back(SymLocal1);
	this->postTypes.push_back(SymLocal2);
	this->postTypes.push_back(SymLocal3);
	this->postTypes.push_back(SymLocal4);
	this->postTypes.push_back(SymLocal5);
	this->postTypes.push_back(SymLocal6);
	this->postTypes.push_back(L1);
	this->postTypes.push_back(L23Trans);
	this->postTypes.push_back(L45Sym);
	this->postTypes.push_back(L45Peak);
	this->postTypes.push_back(L56Trans);
	if(this->preTypes.size())
		this->preTypes.clear();
	this->preTypes.push_back(L2axon);
	this->preTypes.push_back(L34axon);
	this->preTypes.push_back(L4pyaxon);
	this->preTypes.push_back(L4spaxon);
	this->preTypes.push_back(L4ssaxon);
	this->preTypes.push_back(L5staxon);
	this->preTypes.push_back(L5ttaxon);
	this->preTypes.push_back(L6ccaxon);
	this->preTypes.push_back(L6ccinvaxon);
	this->preTypes.push_back(L6ctaxon);
	this->preTypes.push_back(SymLocalaxon);
	this->preTypes.push_back(SymLocal1axon);
	this->preTypes.push_back(SymLocal2axon);
	this->preTypes.push_back(SymLocal3axon);
	this->preTypes.push_back(SymLocal4axon);
	this->preTypes.push_back(SymLocal5axon);
	this->preTypes.push_back(SymLocal6axon);
	this->preTypes.push_back(L1axon);
	this->preTypes.push_back(L23Transaxon);
	this->preTypes.push_back(L45Symaxon);
	this->preTypes.push_back(L45Peakaxon);
	this->preTypes.push_back(L56Transaxon);
	this->preTypes.push_back(VPM);
	if(this->celltypeLabels2Int.size())
		this->celltypeLabels2Int.clear();
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L2"),L2));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L34"),L34));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4py"),L4py));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4sp"),L4sp));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4ss"),L4ss));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5st"),L5st));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5tt"),L5tt));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6cc"),L6cc));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccinv"),L6ccinv));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ct"),L6ct));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal"),SymLocal));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal1"),SymLocal1));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal2"),SymLocal2));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal3"),SymLocal3));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal4"),SymLocal4));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal5"),SymLocal5));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal6"),SymLocal6));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L1"),L1));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L23Trans"),L23Trans));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Sym"),L45Sym));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Peak"),L45Peak));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L56Trans"),L56Trans));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L2axon"),L2axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L34axon"),L34axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4pyaxon"),L4pyaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4spaxon"),L4spaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4ssaxon"),L4ssaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5staxon"),L5staxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5ttaxon"),L5ttaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccaxon"),L6ccaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccinvaxon"),L6ccinvaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ctaxon"),L6ctaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocalaxon"),SymLocalaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal1axon"),SymLocal1axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal2axon"),SymLocal2axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal3axon"),SymLocal3axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal4axon"),SymLocal4axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal5axon"),SymLocal5axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal6axon"),SymLocal6axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L1axon"),L1axon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L23Transaxon"),L23Transaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Symaxon"),L45Symaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Peakaxon"),L45Peakaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L56Transaxon"),L56Transaxon));
	this->celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("VPM"),VPM));
	if(this->int2CelltypeLabels.size())
		this->int2CelltypeLabels.clear();
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L2,"L2"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L34,"L34"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4py,"L4py"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4sp,"L4sp"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4ss,"L4ss"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5st,"L5st"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5tt,"L5tt"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6cc,"L6cc"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccinv,"L6ccinv"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ct,"L6ct"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal,"SymLocal"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal1,"SymLocal1"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal2,"SymLocal2"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal3,"SymLocal3"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal4,"SymLocal4"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal5,"SymLocal5"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal6,"SymLocal6"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L1,"L1"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L23Trans,"L23Trans"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Sym,"L45Sym"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Peak,"L45Peak"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L56Trans,"L56Trans"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L2axon,"L2axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L34axon,"L34axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4pyaxon,"L4pyaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4spaxon,"L4spaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4ssaxon,"L4ssaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5staxon,"L5staxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5ttaxon,"L5ttaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccaxon,"L6ccaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccinvaxon,"L6ccinvaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ctaxon,"L6ctaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocalaxon,"SymLocalaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal1axon,"SymLocal1axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal2axon,"SymLocal2axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal3axon,"SymLocal3axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal4axon,"SymLocal4axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal5axon,"SymLocal5axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal6axon,"SymLocal6axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L1axon,"L1axon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L23Transaxon,"L23Transaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Symaxon,"L45Symaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Peakaxon,"L45Peakaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L56Transaxon,"L56Transaxon"));
	this->int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(VPM,"VPM"));
}

void initializeConstants_const()
{
	if(barrelLabels_const.size())
		barrelLabels_const.clear();
	barrelLabels_const.push_back(Alpha);
	barrelLabels_const.push_back(A1);
	barrelLabels_const.push_back(A2);
	barrelLabels_const.push_back(A3);
	barrelLabels_const.push_back(A4);
	barrelLabels_const.push_back(Beta);
	barrelLabels_const.push_back(B1);
	barrelLabels_const.push_back(B2);
	barrelLabels_const.push_back(B3);
	barrelLabels_const.push_back(B4);
	barrelLabels_const.push_back(Gamma);
	barrelLabels_const.push_back(C1);
	barrelLabels_const.push_back(C2);
	barrelLabels_const.push_back(C3);
	barrelLabels_const.push_back(C4);
	barrelLabels_const.push_back(C5);
	barrelLabels_const.push_back(C6);
	barrelLabels_const.push_back(Delta);
	barrelLabels_const.push_back(D1);
	barrelLabels_const.push_back(D2);
	barrelLabels_const.push_back(D3);
	barrelLabels_const.push_back(D4);
	barrelLabels_const.push_back(D5);
	barrelLabels_const.push_back(D6);
	barrelLabels_const.push_back(E1);
	barrelLabels_const.push_back(E2);
	barrelLabels_const.push_back(E3);
	barrelLabels_const.push_back(E4);
	barrelLabels_const.push_back(E5);
	barrelLabels_const.push_back(E6);
	if(int2Labels_const.size())
		int2Labels_const.clear();
	int2Labels_const.insert(std::pair< int, const char * >(Alpha, "Alpha"));
	int2Labels_const.insert(std::pair< int, const char * >(A1, "A1"));
	int2Labels_const.insert(std::pair< int, const char * >(A2, "A2"));
	int2Labels_const.insert(std::pair< int, const char * >(A3, "A3"));
	int2Labels_const.insert(std::pair< int, const char * >(A4, "A4"));
	int2Labels_const.insert(std::pair< int, const char * >(Beta, "Beta"));
	int2Labels_const.insert(std::pair< int, const char * >(B1, "B1"));
	int2Labels_const.insert(std::pair< int, const char * >(B2, "B2"));
	int2Labels_const.insert(std::pair< int, const char * >(B3, "B3"));
	int2Labels_const.insert(std::pair< int, const char * >(B4, "B4"));
	int2Labels_const.insert(std::pair< int, const char * >(Gamma, "Gamma"));
	int2Labels_const.insert(std::pair< int, const char * >(C1, "C1"));
	int2Labels_const.insert(std::pair< int, const char * >(C2, "C2"));
	int2Labels_const.insert(std::pair< int, const char * >(C3, "C3"));
	int2Labels_const.insert(std::pair< int, const char * >(C4, "C4"));
	int2Labels_const.insert(std::pair< int, const char * >(C5, "C5"));
	int2Labels_const.insert(std::pair< int, const char * >(C6, "C6"));
	int2Labels_const.insert(std::pair< int, const char * >(Delta, "Delta"));
	int2Labels_const.insert(std::pair< int, const char * >(D1, "D1"));
	int2Labels_const.insert(std::pair< int, const char * >(D2, "D2"));
	int2Labels_const.insert(std::pair< int, const char * >(D3, "D3"));
	int2Labels_const.insert(std::pair< int, const char * >(D4, "D4"));
	int2Labels_const.insert(std::pair< int, const char * >(D5, "D5"));
	int2Labels_const.insert(std::pair< int, const char * >(D6, "D6"));
	int2Labels_const.insert(std::pair< int, const char * >(E1, "E1"));
	int2Labels_const.insert(std::pair< int, const char * >(E2, "E2"));
	int2Labels_const.insert(std::pair< int, const char * >(E3, "E3"));
	int2Labels_const.insert(std::pair< int, const char * >(E4, "E4"));
	int2Labels_const.insert(std::pair< int, const char * >(E5, "E5"));
	int2Labels_const.insert(std::pair< int, const char * >(E6, "E6"));
	if(hocLabels_const.size())
		hocLabels_const.clear();
	hocLabels_const.push_back("alpha");
	hocLabels_const.push_back("WM");
	hocLabels_const.push_back("soma");
	hocLabels_const.push_back("dend");
	hocLabels_const.push_back("apical");
	hocLabels_const.push_back("axon");
	hocLabels_const.push_back("Alpha");
	hocLabels_const.push_back("A1");
	hocLabels_const.push_back("A2");
	hocLabels_const.push_back("A3");
	hocLabels_const.push_back("A4");
	hocLabels_const.push_back("Beta");
	hocLabels_const.push_back("B1");
	hocLabels_const.push_back("B2");
	hocLabels_const.push_back("B3");
	hocLabels_const.push_back("B4");
	hocLabels_const.push_back("Gamma");
	hocLabels_const.push_back("C1");
	hocLabels_const.push_back("C2");
	hocLabels_const.push_back("C3");
	hocLabels_const.push_back("C4");
	hocLabels_const.push_back("C5");
	hocLabels_const.push_back("C6");
	hocLabels_const.push_back("Delta");
	hocLabels_const.push_back("D1");
	hocLabels_const.push_back("D2");
	hocLabels_const.push_back("D3");
	hocLabels_const.push_back("D4");
	hocLabels_const.push_back("D5");
	hocLabels_const.push_back("D6");
	hocLabels_const.push_back("E1");
	hocLabels_const.push_back("E2");
	hocLabels_const.push_back("E3");
	hocLabels_const.push_back("E4");
	hocLabels_const.push_back("E5");
	hocLabels_const.push_back("E6");
	if(hocNeuronLabels_const.size())
		hocNeuronLabels_const.clear();
	hocNeuronLabels_const.push_back("soma");
	hocNeuronLabels_const.push_back("dend");
	hocNeuronLabels_const.push_back("apical");
	hocNeuronLabels_const.push_back("axon");
	if(hocLandmarkLabels_const.size())
		hocLandmarkLabels_const.clear();
	hocLandmarkLabels_const.push_back("alpha");
	hocLandmarkLabels_const.push_back("WM");
	hocLandmarkLabels_const.push_back("Alpha");
	hocLandmarkLabels_const.push_back("A1");
	hocLandmarkLabels_const.push_back("A2");
	hocLandmarkLabels_const.push_back("A3");
	hocLandmarkLabels_const.push_back("A4");
	hocLandmarkLabels_const.push_back("Beta");
	hocLandmarkLabels_const.push_back("B1");
	hocLandmarkLabels_const.push_back("B2");
	hocLandmarkLabels_const.push_back("B3");
	hocLandmarkLabels_const.push_back("B4");
	hocLandmarkLabels_const.push_back("Gamma");
	hocLandmarkLabels_const.push_back("C1");
	hocLandmarkLabels_const.push_back("C2");
	hocLandmarkLabels_const.push_back("C3");
	hocLandmarkLabels_const.push_back("C4");
	hocLandmarkLabels_const.push_back("C5");
	hocLandmarkLabels_const.push_back("C6");
	hocLandmarkLabels_const.push_back("Delta");
	hocLandmarkLabels_const.push_back("D1");
	hocLandmarkLabels_const.push_back("D2");
	hocLandmarkLabels_const.push_back("D3");
	hocLandmarkLabels_const.push_back("D4");
	hocLandmarkLabels_const.push_back("D5");
	hocLandmarkLabels_const.push_back("D6");
	hocLandmarkLabels_const.push_back("E1");
	hocLandmarkLabels_const.push_back("E2");
	hocLandmarkLabels_const.push_back("E3");
	hocLandmarkLabels_const.push_back("E4");
	hocLandmarkLabels_const.push_back("E5");
	hocLandmarkLabels_const.push_back("E6");
}
