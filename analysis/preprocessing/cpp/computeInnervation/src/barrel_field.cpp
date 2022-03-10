/****************************************************************************/
/*                                                                          */
/* File:      barrel_field.cpp                                              */
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
/* Date: 		06.12.2018                                                  */
/* 																			*/
/* Note: 		Change path to amira files in line 585						*/
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "barrel_field.h"

BarrelField::BarrelField()
{
// 	std::cout << "Loading Standard Barrel Field..." << std::endl;
	initializeConstants();
	readStandardBarrelField();
};

BarrelField::BarrelField(const char * filepath)
{
	initializeConstants();
	readStandardBarrelField(filepath);
};

BarrelField::~BarrelField()
{
	if(avgBarrelField) delete avgBarrelField;
	if(avgPiaSurface) delete avgPiaSurface;
	if(avgWMSurface) delete avgWMSurface;
	if(avgL4UpperSurface) delete avgL4UpperSurface;
	if(avgL4LowerSurface) delete avgL4LowerSurface;
	if(S1ConvexHull) delete S1ConvexHull;
	std::list < int >::const_iterator labelIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(avgColumns.find(ID) != avgColumns.end() && avgColumns[ID]) delete avgColumns[ID];
		if(avgBarrels.find(ID) != avgBarrels.end() && avgBarrels[ID]) delete avgBarrels[ID];
		if(avgAxes.find(ID) != avgAxes.end() && avgAxes[ID]) delete [] avgAxes[ID];
		if(avgCenters.find(ID) != avgCenters.end() && avgCenters[ID]) delete [] avgCenters[ID];
	}
	barrelLabels.clear();
	borderBarrels.clear();
	neighborBarrel.clear();
	avgTopDist.clear();
	avgPiaWMDist.clear();
	avgBarrelArea.clear();
	int2Labels.clear();
	labels2Int.clear();
	avgColumns.clear();
	avgBarrels.clear();
	avgAxes.clear();
	avgCenters.clear();
};

/******************************************************************************/
/*local z axis: nearest neighbor interpolation                                */
/******************************************************************************/
void BarrelField::localZAxis(double x[3], double zAxis[3])
{
	// check if axis vector field can be used
	if(axisVecField)
	{
		double bounds[6], delta[3] = {1E-6, 1E-6, 1E-6};
		axisVecField->GetBounds(bounds);
		if(vtkMath::PointIsWithinBounds(x, bounds, delta))
		{
			double spacing[3], origin[3];
			int xyz[3];
			axisVecField->GetSpacing(spacing);
			axisVecField->GetOrigin(origin);
			for(int ii = 0; ii < 3; ++ii)
				xyz[ii] = int((x[ii] - origin[ii])/spacing[ii]);
			
			double * dataPtr = static_cast< double * >(axisVecField->GetScalarPointer(xyz));
			zAxis[0] = dataPtr[0];
			zAxis[1] = dataPtr[1];
			zAxis[2] = dataPtr[2];
			
			return;
		}
	}
	
	// otherwise use coarser version
	double minDist = 1E06;
	int nearestID;
	for(int ii = 0; ii < avgAxesField->GetNumberOfCells(); ++ii)
	{
		double tmpAxis[3], pt1[3], pt2[3], closestPt[3], t;
		avgAxesField->GetCell(ii)->GetPoints()->GetPoint(0, pt1);
		avgAxesField->GetCell(ii)->GetPoints()->GetPoint(1, pt2);
		for(int jj = 0; jj < 3; ++jj)
			tmpAxis[jj] = pt1[jj] - pt2[jj];
		vtkMath::Normalize(tmpAxis);
		for(int jj = 0; jj < 3; ++jj)
		{
			pt1[jj] += 2000*tmpAxis[jj];
			pt2[jj] -= 2000*tmpAxis[jj];
		}
		double dist = vtkLine::DistanceToLine(x, pt1, pt2, t, closestPt);
		dist = sqrt(dist);
		if(dist < minDist)
		{
			minDist = dist;
			nearestID = ii;
		}
	}
	double axis[3], pt1[3], pt2[3];
	avgAxesField->GetCell(nearestID)->GetPoints()->GetPoint(0, pt1);
	avgAxesField->GetCell(nearestID)->GetPoints()->GetPoint(1, pt2);
	for(int ii = 0; ii < 3; ++ii)
		axis[ii] = pt1[ii] - pt2[ii];
	if(axis[2] < 0)
		axis[0] = -axis[0], axis[1] = -axis[1], axis[2] = -axis[2];
	vtkMath::Normalize(axis);
	for(int ii = 0; ii < 3; ++ii)
		zAxis[ii] = axis[ii];
};

/******************************************************************************/
/*local z axis: nearest neighbor interpolation; also sets nearest neighbor    */
/*found in closestPt                                                          */
/******************************************************************************/
void BarrelField::localZAxis(double x[3], double closestPt[3], double zAxis[3])
{
	// check if axis vector field can be used
	if(axisVecField)
	{
		double bounds[6], delta[3] = {1E-6, 1E-6, 1E-6};
		axisVecField->GetBounds(bounds);
		if(vtkMath::PointIsWithinBounds(x, bounds, delta))
		{
			double spacing[3], origin[3];
			int xyz[3];
			axisVecField->GetSpacing(spacing);
			axisVecField->GetOrigin(origin);
			for(int ii = 0; ii < 3; ++ii)
			{
				xyz[ii] = int((x[ii] - origin[ii])/spacing[ii]);
				closestPt[ii] = xyz[ii]*spacing[ii] + origin[ii];
			}
			
			double * dataPtr = static_cast< double * >(axisVecField->GetScalarPointer(xyz));
			zAxis[0] = dataPtr[0];
			zAxis[1] = dataPtr[1];
			zAxis[2] = dataPtr[2];
			
			return;
		}
	}
	
	// otherwise use coarser version
	double minDist = 1E06;
	int nearestID;
	for(int ii = 0; ii < avgAxesField->GetNumberOfCells(); ++ii)
	{
		double tmpAxis[3], pt1[3], pt2[3], projectedPt[3], t;
		avgAxesField->GetCell(ii)->GetPoints()->GetPoint(0, pt1);
		avgAxesField->GetCell(ii)->GetPoints()->GetPoint(1, pt2);
		for(int jj = 0; jj < 3; ++jj)
			tmpAxis[jj] = pt1[jj] - pt2[jj];
		vtkMath::Normalize(tmpAxis);
		for(int jj = 0; jj < 3; ++jj)
		{
			pt1[jj] += 2000*tmpAxis[jj];
			pt2[jj] -= 2000*tmpAxis[jj];
		}
		double dist = vtkLine::DistanceToLine(x, pt1, pt2, t, projectedPt);
		dist = sqrt(dist);
		if(dist < minDist)
		{
			minDist = dist;
			nearestID = ii;
			closestPt[0] = projectedPt[0], closestPt[1] = projectedPt[1], closestPt[2] = projectedPt[2];
		}
	}
	double axis[3], pt1[3], pt2[3];
	avgAxesField->GetCell(nearestID)->GetPoints()->GetPoint(0, pt1);
	avgAxesField->GetCell(nearestID)->GetPoints()->GetPoint(1, pt2);
	for(int ii = 0; ii < 3; ++ii)
		axis[ii] = pt1[ii] - pt2[ii];
	if(axis[2] < 0)
		axis[0] = -axis[0], axis[1] = -axis[1], axis[2] = -axis[2];
	vtkMath::Normalize(axis);
	for(int ii = 0; ii < 3; ++ii)
		zAxis[ii] = axis[ii];
};

/******************************************************************************/
/*local z axis: nearest neighbor interpolation; also sets nearest neighbor    */
/*found in closestPt                                                          */
/******************************************************************************/
void BarrelField::localZAxis(double x[3], double closestPt[3], bool verbose, double zAxis[3])
{
	// check if axis vector field can be used
	if(axisVecField)
	{
		double bounds[6], delta[3] = {1E-6, 1E-6, 1E-6};
		axisVecField->GetBounds(bounds);
		if(vtkMath::PointIsWithinBounds(x, bounds, delta))
		{
			double spacing[3], origin[3];
			int xyz[3];
			axisVecField->GetSpacing(spacing);
			axisVecField->GetOrigin(origin);
			for(int ii = 0; ii < 3; ++ii)
			{
				xyz[ii] = int((x[ii] - origin[ii])/spacing[ii]);
				closestPt[ii] = xyz[ii]*spacing[ii] + origin[ii];
			}
			
			double * dataPtr = static_cast< double * >(axisVecField->GetScalarPointer(xyz));
			zAxis[0] = dataPtr[0];
			zAxis[1] = dataPtr[1];
			zAxis[2] = dataPtr[2];
			
			return;
		}
	}
	
	// otherwise use coarser version
	double minDist = 1E06;
	int nearestID;
	for(int ii = 0; ii < avgAxesField->GetNumberOfCells(); ++ii)
	{
		double tmpAxis[3], pt1[3], pt2[3], projectedPt[3], t;
		avgAxesField->GetCell(ii)->GetPoints()->GetPoint(0, pt1);
		avgAxesField->GetCell(ii)->GetPoints()->GetPoint(1, pt2);
		for(int jj = 0; jj < 3; ++jj)
			tmpAxis[jj] = pt1[jj] - pt2[jj];
		vtkMath::Normalize(tmpAxis);
		for(int jj = 0; jj < 3; ++jj)
		{
			pt1[jj] += 2000*tmpAxis[jj];
			pt2[jj] -= 2000*tmpAxis[jj];
		}
		double dist = vtkLine::DistanceToLine(x, pt1, pt2, t, projectedPt);
		dist = sqrt(dist);
		if(dist < minDist)
		{
			minDist = dist;
			nearestID = ii;
			closestPt[0] = projectedPt[0], closestPt[1] = projectedPt[1], closestPt[2] = projectedPt[2];
		}
	}
	
	if(verbose)
	{
		std::cout << "minDist =\t" << minDist << std::endl;
		std::cout << "nearestID =\t" << nearestID << std::endl;
		std::cout << "closestPt =\t[" << closestPt[0] << "," << closestPt[1] << "," << closestPt[2] << "]" << std::endl;
	}
	
	double axis[3], pt1[3], pt2[3];
	avgAxesField->GetCell(nearestID)->GetPoints()->GetPoint(0, pt1);
	avgAxesField->GetCell(nearestID)->GetPoints()->GetPoint(1, pt2);
	for(int ii = 0; ii < 3; ++ii)
		axis[ii] = pt1[ii] - pt2[ii];
	if(axis[2] < 0)
		axis[0] = -axis[0], axis[1] = -axis[1], axis[2] = -axis[2];
	vtkMath::Normalize(axis);
	for(int ii = 0; ii < 3; ++ii)
		zAxis[ii] = axis[ii];
};

/******************************************************************************/
/*local z axis: of home barrel                                                */
/******************************************************************************/
void BarrelField::localZAxis(int HBID, double zAxis[3])
{
	double axis[3];
	for(int ii = 0; ii < 3; ++ii)
		axis[ii] = avgBarrels[HBID]->top[ii] - avgBarrels[HBID]->bottom[ii];
	vtkMath::Normalize(axis);
	for(int ii = 0; ii < 3; ++ii)
		zAxis[ii] = axis[ii];
};

void BarrelField::localCoordinateSystem(double x[3], double xAxis[3], double yAxis[3], double zAxis[3])
{
	double dist = -1;
	localZAxis(x, zAxis);
	
	int HBID = closestBarrel(x);
	int NBID = neighborBarrel[HBID];
	
	for(int ii = 0; ii < 3; ++ii)
		xAxis[ii] = avgCenters[NBID][ii] - avgCenters[HBID][ii];
	double correction = vtkMath::Dot(zAxis, xAxis);
	for(int ii = 0; ii < 3; ++ii)
		xAxis[ii] -= correction*zAxis[ii];
	vtkMath::Normalize(xAxis);
	
	vtkMath::Cross(zAxis, xAxis, yAxis);
// 	vtkMath::Normalize(yAxis);
};

void BarrelField::localCylindricalCoordinates(double x[3], double coordinates[3])
{
	double dist = -1;
	double zAxis[3];
	localZAxis(x, zAxis);
	avgPiaSurface->intersectLine(zAxis, x);
	double * piaIntersectPt = avgPiaSurface->getLastIntersectPoint();
	if(piaIntersectPt)
	{
		dist = sqrt(vtkMath::Distance2BetweenPoints(x, piaIntersectPt));
		delete [] piaIntersectPt;
	}
	
	//compute coordinates in cylindrical coordinates
	//measured in home column
	//phi measured relative to neighboring barrel along row
	double rPos, phiPos, zPos = dist;
	double radialSomaPt[3];
	int HBID = closestBarrel(x);
	int NBID = neighborBarrel[HBID];
	
	double t, closestPt[3], colAxis[3], xAxis[3], yAxis[3];
	rPos = vtkLine::DistanceToLine(x, avgColumns[HBID]->top, avgColumns[HBID]->bottom, t, closestPt);
	rPos = sqrt(rPos);
	
	for(int ii = 0; ii < 3; ++ii)
	{
		radialSomaPt[ii] = x[ii] - closestPt[ii];
		colAxis[ii] = avgColumns[HBID]->top[ii] - avgColumns[HBID]->bottom[ii];
		xAxis[ii] = avgCenters[NBID][ii] - avgCenters[HBID][ii];
	}
	vtkMath::Normalize(colAxis);
	double correction = vtkMath::Dot(colAxis, xAxis);
	for(int ii = 0; ii < 3; ++ii)
		xAxis[ii] -= correction*colAxis[ii];
	vtkMath::Normalize(xAxis);
	vtkMath::Cross(colAxis, xAxis, yAxis);
	double xNew = vtkMath::Dot(radialSomaPt, xAxis);
	double yNew = vtkMath::Dot(radialSomaPt, yAxis);
	phiPos = std::atan2(yNew, xNew) + PI;
	
	coordinates[0] = rPos;
	coordinates[1] = phiPos;
	coordinates[2] = zPos;
};

void BarrelField::localCylindricalCoordinates(double x[3], int HBID, double coordinates[3])
{
	double dist = -1;
	double zAxis[3];
	localZAxis(x, zAxis);
	avgPiaSurface->intersectLine(zAxis, x);
	double * piaIntersectPt = avgPiaSurface->getLastIntersectPoint();
	if(piaIntersectPt)
	{
		dist = sqrt(vtkMath::Distance2BetweenPoints(x, piaIntersectPt));
		delete [] piaIntersectPt;
	}
	
	//compute coordinates in cylindrical coordinates
	//measured in home column
	//phi measured relative to neighboring barrel along row
	double rPos, phiPos, zPos = dist;
	double radialSomaPt[3];
	int NBID = neighborBarrel[HBID];
	
	double t, closestPt[3], colAxis[3], xAxis[3], yAxis[3];
	rPos = vtkLine::DistanceToLine(x, avgColumns[HBID]->top, avgColumns[HBID]->bottom, t, closestPt);
	rPos = sqrt(rPos);
	
	for(int ii = 0; ii < 3; ++ii)
	{
		radialSomaPt[ii] = x[ii] - closestPt[ii];
		colAxis[ii] = avgColumns[HBID]->top[ii] - avgColumns[HBID]->bottom[ii];
		xAxis[ii] = avgCenters[NBID][ii] - avgCenters[HBID][ii];
	}
	vtkMath::Normalize(colAxis);
	double correction = vtkMath::Dot(colAxis, xAxis);
	for(int ii = 0; ii < 3; ++ii)
		xAxis[ii] -= correction*colAxis[ii];
	vtkMath::Normalize(xAxis);
	vtkMath::Cross(colAxis, xAxis, yAxis);
	double xNew = vtkMath::Dot(radialSomaPt, xAxis);
	double yNew = vtkMath::Dot(radialSomaPt, yAxis);
	phiPos = std::atan2(yNew, xNew) + PI;
	
	coordinates[0] = rPos;
	coordinates[1] = phiPos;
	coordinates[2] = zPos;
};

int BarrelField::closestBarrel(double x[3])
{
	int HBID = 0;
	double minDist = 1E08;
	std::list< int >::const_iterator labelIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(avgCenters.find(ID) != avgCenters.end() && avgAxes.find(ID) != avgAxes.end())
		{
			double linePt1[3], linePt2[3], closestPt[3], t;
			for(int ii = 0; ii < 3; ++ii)
			{
				linePt1[ii] = avgCenters[ID][ii] + 2000*avgAxes[ID][ii];
				linePt2[ii] = avgCenters[ID][ii] - 2000*avgAxes[ID][ii];
			}
			double dist = vtkLine::DistanceToLine(x, linePt1, linePt2, t, closestPt);
			dist = sqrt(dist);
			if(dist < minDist)
			{
				minDist = dist;
				HBID = ID;
			}
		}
	}
	return HBID;
};

bool BarrelField::isInsideS1 ( double x[3] )
{
	if(S1ConvexHull)
	{
		if(S1ConvexHull->isPointInsideSurface(x))
		{
			// check whether point is also
			// between Pia and WM
			double zAxis[3];
			localZAxis(x, zAxis);
			avgPiaSurface->intersectLine(zAxis, x);
			avgWMSurface->intersectLine(zAxis, x);
			if(avgPiaSurface->isIntersectionFound() && avgWMSurface->isIntersectionFound())
			{
				double piaIntersectPt[3], WMIntersectPt[3], projectedPt[3], t;
				avgPiaSurface->getLastIntersectPoint(piaIntersectPt);
				avgWMSurface->getLastIntersectPoint(WMIntersectPt);
				vtkLine::DistanceToLine(x, piaIntersectPt, WMIntersectPt, t, projectedPt);
				if(t >= 0 && t <= 1)
					return 1;
				return 0;
			}
			else
			{
				std::cout << "Error! Could not determine location relative to Pia/WM!" << std::endl;
				return 0;
			}
		}
		else
			return 0;
	}
	else
	{
		std::cout << "Error! S1 convex hull invalid!" << std::endl;
		return 0;
	}
}

int BarrelField::laminarPosition ( double x[3] )
{
	double zAxis[3];
	localZAxis(x, zAxis);
	// check supragranular
	avgPiaSurface->intersectLine(zAxis, x);
	avgL4UpperSurface->intersectLine(zAxis, x);
	if(avgPiaSurface->isIntersectionFound() && avgL4UpperSurface->isIntersectionFound())
	{
		double piaIntersectPt[3], L4UIntersectPt[3], projectedPt[3], t;
		avgPiaSurface->getLastIntersectPoint(piaIntersectPt);
		avgL4UpperSurface->getLastIntersectPoint(L4UIntersectPt);
		vtkLine::DistanceToLine(x, piaIntersectPt, L4UIntersectPt, t, projectedPt);
		if(t >= 0 && t <= 1)
			return SUPRA;
	}
	// check granular
	avgL4UpperSurface->intersectLine(zAxis, x);
	avgL4LowerSurface->intersectLine(zAxis, x);
	if(avgL4LowerSurface->isIntersectionFound() && avgL4UpperSurface->isIntersectionFound())
	{
		double L4LIntersectPt[3], L4UIntersectPt[3], projectedPt[3], t;
		avgL4LowerSurface->getLastIntersectPoint(L4LIntersectPt);
		avgL4UpperSurface->getLastIntersectPoint(L4UIntersectPt);
		vtkLine::DistanceToLine(x, L4LIntersectPt, L4UIntersectPt, t, projectedPt);
		if(t >= 0 && t <= 1)
			return GRAN;
	}
	// check infragranular
	avgL4LowerSurface->intersectLine(zAxis, x);
	avgWMSurface->intersectLine(zAxis, x);
	if(avgL4LowerSurface->isIntersectionFound() && avgWMSurface->isIntersectionFound())
	{
		double L4LIntersectPt[3], WMIntersectPt[3], projectedPt[3], t;
		avgL4LowerSurface->getLastIntersectPoint(L4LIntersectPt);
		avgWMSurface->getLastIntersectPoint(WMIntersectPt);
		vtkLine::DistanceToLine(x, L4LIntersectPt, WMIntersectPt, t, projectedPt);
		if(t >= 0 && t <= 1)
			return INFRA;
	}
	
	return 0;
}

int BarrelField::insideColumn ( double x[3] )
{
	int closestColumn = 0;
// 	double r, dist;
// 	double closestPt[3], t;
	double minDist = 1E08;
	std::list< int >::const_iterator labelIt;
	for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
	{
		int ID = *labelIt;
		if(avgCenters.find(ID) != avgCenters.end() && avgAxes.find(ID) != avgAxes.end())
		{
			double linePt1[3], linePt2[3], closestPt[3], t;
			for(int ii = 0; ii < 3; ++ii)
			{
				linePt1[ii] = avgCenters[ID][ii] + 2000*avgAxes[ID][ii];
				linePt2[ii] = avgCenters[ID][ii] - 2000*avgAxes[ID][ii];
			}
			double dist, r;
			dist = sqrt(vtkLine::DistanceToLine(x, linePt1, linePt2, t, closestPt));
			if(dist < minDist)
			{
				dist = sqrt(vtkLine::DistanceToLine(x, avgColumns[ID]->top, avgColumns[ID]->bottom, t, closestPt));
				r = sqrt(avgBarrelArea[ID]/PI);
				if(dist <= r && t >= 0  && t <= 1)
				{
					minDist = dist;
					closestColumn = ID;
				}
			}
		}
	}
// 	dist = sqrt(vtkLine::DistanceToLine(x, avgColumns[closestColumn]->top, avgColumns[closestColumn]->bottom, t, closestPt));
// 	r = sqrt(avgBarrelArea[closestColumn]/PI);
// 	if(dist <= r && t >= 0  && t <= 1)
// 		return closestColumn;
// 	return 0;
	return closestColumn;
}

double BarrelField::piaDistance(double x[3])
{
	double zAxis[3];
	localZAxis(x, zAxis);
	avgPiaSurface->intersectLine(zAxis, x);
	if(avgPiaSurface->isIntersectionFound())
	{
		double piaIntersectPt[3];
		avgPiaSurface->getLastIntersectPoint(piaIntersectPt);
		double dist2 = vtkMath::Distance2BetweenPoints(x, piaIntersectPt);
		return sqrt(dist2);
	}
	else
		return -1;
}

void BarrelField::readStandardBarrelField()
{
	// CHANGE PATH HERE!
	const char * filepath = "/home/dudvary/project_src/BarrelField3D/common/";
	readStandardBarrelField(filepath);
}

void BarrelField::readStandardBarrelField(const char * filepath)
{
	std::string filepathStr(filepath);

	std::string bfNameStr = filepathStr + "average_barrel_field_complete_S1.am";
	std::string piaNameStr = filepathStr + "average_barrel_field_pia_complete_S1.surf";
	std::string WMNameStr = filepathStr + "average_barrel_field_WM_complete_S1.surf";
	std::string L4UpperNameStr = filepathStr + "average_barrel_field_L4Upper_complete_S1.surf";
	std::string L4LowerNameStr = filepathStr + "average_barrel_field_L4Lower_complete_S1.surf";
	std::string S1HullNameStr = filepathStr + "average_barrel_field_S1_convex_hull.surf";
	std::string axisFieldNameStr = filepathStr + "axis_field_50mu.am";

	const char * bfName = bfNameStr.c_str();
	const char * piaName = piaNameStr.c_str();
	const char * WMName = WMNameStr.c_str();
	const char * L4UpperName = L4UpperNameStr.c_str();
	const char * L4LowerName = L4LowerNameStr.c_str();
	const char * S1HullName = S1HullNameStr.c_str();
	const char * axisFieldName = axisFieldNameStr.c_str();
	
	Reader * amReader = new Reader(bfName);
	amReader->readSpatialGraphFile(0);
	avgBarrelField = amReader->getSpatialGraph();
	if(avgBarrelField)
	{
		std::list< int >::const_iterator labelIt;
		if(avgColumns.size())
			avgColumns.clear();
		if(avgBarrels.size())
			avgBarrels.clear();
		if(avgAxes.size())
			avgAxes.clear();
		if(avgCenters.size())
			avgCenters.clear();
		for(labelIt = barrelLabels.begin(); labelIt != barrelLabels.end(); ++labelIt)
		{
			int ID = *labelIt;
			if(avgBarrelField->isLabelInSpatialGraph(ID))
			{
				PolyDataPointerType currBarrel = PolyDataPointerType::New();
				if(avgBarrelField->extractLandmark(ID, currBarrel))
				{
					std::map< double, int > cellIDMap;	// sorts cell IDs by z value
					if(currBarrel->GetNumberOfCells() != 4)
					{
						std::cout << "Error! Wrong format of barrel " << int2Labels[ID] << " in Standard Barrel Field! Check SpatialGraph file." << std::endl;
						return;
					}
					for(int ii = 0; ii < 4; ++ii)
						cellIDMap.insert(std::pair< double, int >(currBarrel->GetCell(ii)->GetBounds()[4], ii));
					
					std::map< double, int >::const_iterator cellIDMapIt = cellIDMap.begin();
					std::map< double, int >::const_reverse_iterator cellIDMapRIt = cellIDMap.rbegin();
					
					PolyDataPointerType columnContours = PolyDataPointerType::New();
					PolyDataPointerType barrelContours = PolyDataPointerType::New();
					IdListPointerType columnCellIds = IdListPointerType::New();
					IdListPointerType barrelCellIds = IdListPointerType::New();
					columnContours->Allocate();
					barrelContours->Allocate();
					
					columnCellIds->InsertId(0, cellIDMapRIt->second);
					columnCellIds->InsertId(1, cellIDMapIt->second);
					++cellIDMapRIt, ++cellIDMapIt;
					barrelCellIds->InsertId(0, cellIDMapRIt->second);
					barrelCellIds->InsertId(1, cellIDMapIt->second);
					columnContours->CopyCells(currBarrel, columnCellIds);
					barrelContours->CopyCells(currBarrel, barrelCellIds);
					
					double barrelTop[3], barrelBottom[3], columnTop[3], columnBottom[3], paramCenter[3];
					//parametric center
					int subID;
					double pCoords[3], * weights1, * weights2, * weights3, *weights4;
					weights1 = new double[barrelContours->GetCell(0)->GetNumberOfPoints()];
					weights2 = new double[barrelContours->GetCell(1)->GetNumberOfPoints()];
					weights3 = new double[columnContours->GetCell(0)->GetNumberOfPoints()];
					weights4 = new double[columnContours->GetCell(1)->GetNumberOfPoints()];
					barrelContours->GetCell(0)->GetParametricCenter(pCoords);
					barrelContours->GetCell(0)->EvaluateLocation(subID, pCoords, barrelTop, weights1);
					barrelContours->GetCell(1)->GetParametricCenter(pCoords);
					barrelContours->GetCell(1)->EvaluateLocation(subID, pCoords, barrelBottom, weights2);
					columnContours->GetCell(0)->GetParametricCenter(pCoords);
					columnContours->GetCell(0)->EvaluateLocation(subID, pCoords, columnTop, weights3);
					columnContours->GetCell(1)->GetParametricCenter(pCoords);
					columnContours->GetCell(1)->EvaluateLocation(subID, pCoords, columnBottom, weights4);
					
					Column * newCol = new Column(columnContours, columnTop, columnBottom);
					Column * newBarrel = new Column(barrelContours, barrelTop, barrelBottom);
					avgColumns.insert(std::pair< int, Column * >(ID, newCol));
					avgBarrels.insert(std::pair< int, Column * >(ID, newBarrel));
					double * avgAxis = new double[3], * avgCenter = new double[3];
					for(int ii = 0; ii < 3; ++ii)
					{
						avgAxis[ii] = barrelTop[ii] - barrelBottom[ii];
						avgCenter[ii] = 0.5*(barrelTop[ii] + barrelBottom[ii]);
					}
					vtkMath::Normalize(avgAxis);
					avgAxes.insert(std::pair< int, double * >(ID, avgAxis));
					avgCenters.insert(std::pair< int, double * >(ID, avgCenter));
// 					std::cout << "barrel " << int2Labels[ID] << " top @ [" << newBarrel->top[0] << "," << newBarrel->top[1] << "," << newBarrel->top[2] << "]" << std::endl;
// 					std::cout << "barrel " << int2Labels[ID] << " bottom @ [" << newBarrel->bottom[0] << "," << newBarrel->bottom[1] << "," << newBarrel->bottom[2] << "]" << std::endl;
// 					std::cout << "column " << int2Labels[ID] << " top @ [" << newCol->top[0] << "," << newCol->top[1] << "," << newCol->top[2] << "]" << std::endl;
// 					std::cout << "column " << int2Labels[ID] << " bottom @ [" << newCol->bottom[0] << "," << newCol->bottom[1] << "," << newCol->bottom[2] << "]" << std::endl;
					
					delete [] weights1, delete [] weights2, delete [] weights3, delete [] weights4;
				}
			}
		}
		if(avgBarrelField->isLabelInSpatialGraph(ZAxis))
		{
			avgAxesField = PolyDataPointerType::New();
			if(avgBarrelField->extractLandmark(ZAxis, avgAxesField))
			{
// 				TransformPointerType zScale = TransformPointerType::New();
// 				TransformFilterType zScaleFilter = TransformFilterType::New();
// 				zScale->Scale(3,3,3);
// 				zScaleFilter->SetTransform(zScale);
// 				zScaleFilter->SetInput(avgAxesField);
// 				zScaleFilter->Update();
// 				avgAxesField->DeepCopy(zScaleFilter->GetOutput());
			}
		}
	}
	else
		std::cout << "Error! Could not read Standard Barrel Field file!" << std::endl;
	
	Reader * piaReader = new Reader(piaName);
	avgPia = piaReader->readAmiraSurfaceFile();
	if(!avgPia->GetNumberOfCells())
		std::cout << "Error! Could not read Standard Pia surface file!" << std::endl;
	
	Reader * wmReader = new Reader(WMName);
	avgWM = wmReader->readAmiraSurfaceFile();
	if(!avgWM->GetNumberOfCells())
		std::cout << "Error! Could not read Standard WM surface file!" << std::endl;
	
	Reader * L4UReader = new Reader(L4UpperName);
	avgL4Upper = L4UReader->readAmiraSurfaceFile();
	if(!avgL4Upper->GetNumberOfCells())
		std::cout << "Error! Could not read Standard L4 upper border surface file!" << std::endl;
	
	Reader * L4LReader = new Reader(L4LowerName);
	avgL4Lower = L4LReader->readAmiraSurfaceFile();
	if(!avgL4Lower->GetNumberOfCells())
		std::cout << "Error! Could not read Standard L4 lower border surface file!" << std::endl;
	
	Reader * S1HullReader = new Reader(S1HullName);
	S1HullData = S1HullReader->readAmiraSurfaceFile();
	if(!S1HullData->GetNumberOfCells())
		std::cout << "Error! Could not read Standard S1 convex hull surface file!" << std::endl;
	
	avgPiaSurface = new Surface(avgPia);
	avgWMSurface = new Surface(avgWM);
	avgL4UpperSurface = new Surface(avgL4Upper);
	avgL4LowerSurface = new Surface(avgL4Lower);
	S1ConvexHull = new ClosedSurface(S1HullData);
	
	Reader * axisFieldReader = new Reader(axisFieldName);
	axisVecField = axisFieldReader->readVectorField();
	if(!axisVecField)
		std::cout << "Error! Could not read axis vector field!" << std::endl;
	
	delete amReader;
	delete piaReader;
	delete wmReader;
	delete L4UReader;
	delete L4LReader;
	delete S1HullReader;
	delete axisFieldReader;
};

void BarrelField::initializeConstants()
{
	if(barrelLabels.size())
		barrelLabels.clear();
	barrelLabels.push_back(Alpha);
	barrelLabels.push_back(A1);
	barrelLabels.push_back(A2);
	barrelLabels.push_back(A3);
	barrelLabels.push_back(A4);
	barrelLabels.push_back(Beta);
	barrelLabels.push_back(B1);
	barrelLabels.push_back(B2);
	barrelLabels.push_back(B3);
	barrelLabels.push_back(B4);
	barrelLabels.push_back(Gamma);
	barrelLabels.push_back(C1);
	barrelLabels.push_back(C2);
	barrelLabels.push_back(C3);
	barrelLabels.push_back(C4);
	barrelLabels.push_back(C5);
	barrelLabels.push_back(C6);
	barrelLabels.push_back(Delta);
	barrelLabels.push_back(D1);
	barrelLabels.push_back(D2);
	barrelLabels.push_back(D3);
	barrelLabels.push_back(D4);
	barrelLabels.push_back(D5);
	barrelLabels.push_back(D6);
	barrelLabels.push_back(E1);
	barrelLabels.push_back(E2);
	barrelLabels.push_back(E3);
	barrelLabels.push_back(E4);
	barrelLabels.push_back(E5);
	barrelLabels.push_back(E6);
	if(borderBarrels.size())
		borderBarrels.clear();
	borderBarrels.push_back(Alpha);
	borderBarrels.push_back(A1);
	borderBarrels.push_back(A2);
	borderBarrels.push_back(A3);
	borderBarrels.push_back(A4);
	borderBarrels.push_back(Beta);
	borderBarrels.push_back(B4);
	borderBarrels.push_back(Gamma);
	borderBarrels.push_back(C4);
	borderBarrels.push_back(Delta);
	borderBarrels.push_back(D4);
	borderBarrels.push_back(E1);
	borderBarrels.push_back(E2);
	borderBarrels.push_back(E3);
	borderBarrels.push_back(E4);
	if(neighborBarrel.size())
		neighborBarrel.clear();
	neighborBarrel.insert(std::pair< int, int >(Alpha, A1));
	neighborBarrel.insert(std::pair< int, int >(A1, A2));
	neighborBarrel.insert(std::pair< int, int >(A2, A3));
	neighborBarrel.insert(std::pair< int, int >(A3, A4));
	neighborBarrel.insert(std::pair< int, int >(A4, A3));
	neighborBarrel.insert(std::pair< int, int >(Beta, B1));
	neighborBarrel.insert(std::pair< int, int >(B1, B2));
	neighborBarrel.insert(std::pair< int, int >(B2, B3));
	neighborBarrel.insert(std::pair< int, int >(B3, B4));
	neighborBarrel.insert(std::pair< int, int >(B4, B3));
	neighborBarrel.insert(std::pair< int, int >(Gamma, C1));
	neighborBarrel.insert(std::pair< int, int >(C1, C2));
	neighborBarrel.insert(std::pair< int, int >(C2, C3));
	neighborBarrel.insert(std::pair< int, int >(C3, C4));
	neighborBarrel.insert(std::pair< int, int >(C4, C3));
	neighborBarrel.insert(std::pair< int, int >(Delta, D1));
	neighborBarrel.insert(std::pair< int, int >(D1, D2));
	neighborBarrel.insert(std::pair< int, int >(D2, D3));
	neighborBarrel.insert(std::pair< int, int >(D3, D4));
	neighborBarrel.insert(std::pair< int, int >(D4, D3));
	neighborBarrel.insert(std::pair< int, int >(E1, E2));
	neighborBarrel.insert(std::pair< int, int >(E2, E3));
	neighborBarrel.insert(std::pair< int, int >(E3, E4));
	neighborBarrel.insert(std::pair< int, int >(E4, E3));
	if(int2Labels.size())
		int2Labels.clear();
	int2Labels.insert(std::pair< int, const char * >(Alpha, "Alpha"));
	int2Labels.insert(std::pair< int, const char * >(A1, "A1"));
	int2Labels.insert(std::pair< int, const char * >(A2, "A2"));
	int2Labels.insert(std::pair< int, const char * >(A3, "A3"));
	int2Labels.insert(std::pair< int, const char * >(A4, "A4"));
	int2Labels.insert(std::pair< int, const char * >(Beta, "Beta"));
	int2Labels.insert(std::pair< int, const char * >(B1, "B1"));
	int2Labels.insert(std::pair< int, const char * >(B2, "B2"));
	int2Labels.insert(std::pair< int, const char * >(B3, "B3"));
	int2Labels.insert(std::pair< int, const char * >(B4, "B4"));
	int2Labels.insert(std::pair< int, const char * >(Gamma, "Gamma"));
	int2Labels.insert(std::pair< int, const char * >(C1, "C1"));
	int2Labels.insert(std::pair< int, const char * >(C2, "C2"));
	int2Labels.insert(std::pair< int, const char * >(C3, "C3"));
	int2Labels.insert(std::pair< int, const char * >(C4, "C4"));
	int2Labels.insert(std::pair< int, const char * >(C5, "C5"));
	int2Labels.insert(std::pair< int, const char * >(C6, "C6"));
	int2Labels.insert(std::pair< int, const char * >(Delta, "Delta"));
	int2Labels.insert(std::pair< int, const char * >(D1, "D1"));
	int2Labels.insert(std::pair< int, const char * >(D2, "D2"));
	int2Labels.insert(std::pair< int, const char * >(D3, "D3"));
	int2Labels.insert(std::pair< int, const char * >(D4, "D4"));
	int2Labels.insert(std::pair< int, const char * >(D5, "D5"));
	int2Labels.insert(std::pair< int, const char * >(D6, "D6"));
	int2Labels.insert(std::pair< int, const char * >(E1, "E1"));
	int2Labels.insert(std::pair< int, const char * >(E2, "E2"));
	int2Labels.insert(std::pair< int, const char * >(E3, "E3"));
	int2Labels.insert(std::pair< int, const char * >(E4, "E4"));
	int2Labels.insert(std::pair< int, const char * >(E5, "E5"));
	int2Labels.insert(std::pair< int, const char * >(E6, "E6"));
	if(labels2Int.size())
		labels2Int.clear();
	labels2Int.insert(std::pair< std::string, int >(std::string("Alpha"), Alpha));
	labels2Int.insert(std::pair< std::string, int >(std::string("A1"), A1));
	labels2Int.insert(std::pair< std::string, int >(std::string("A2"), A2));
	labels2Int.insert(std::pair< std::string, int >(std::string("A3"), A3));
	labels2Int.insert(std::pair< std::string, int >(std::string("A4"), A4));
	labels2Int.insert(std::pair< std::string, int >(std::string("Beta"), Beta));
	labels2Int.insert(std::pair< std::string, int >(std::string("B1"), B1));
	labels2Int.insert(std::pair< std::string, int >(std::string("B2"), B2));
	labels2Int.insert(std::pair< std::string, int >(std::string("B3"), B3));
	labels2Int.insert(std::pair< std::string, int >(std::string("B4"), B4));
	labels2Int.insert(std::pair< std::string, int >(std::string("Gamma"), Gamma));
	labels2Int.insert(std::pair< std::string, int >(std::string("C1"), C1));
	labels2Int.insert(std::pair< std::string, int >(std::string("C2"), C2));
	labels2Int.insert(std::pair< std::string, int >(std::string("C3"), C3));
	labels2Int.insert(std::pair< std::string, int >(std::string("C4"), C4));
	labels2Int.insert(std::pair< std::string, int >(std::string("C5"), C5));
	labels2Int.insert(std::pair< std::string, int >(std::string("C6"), C6));
	labels2Int.insert(std::pair< std::string, int >(std::string("Delta"), Delta));
	labels2Int.insert(std::pair< std::string, int >(std::string("D1"), D1));
	labels2Int.insert(std::pair< std::string, int >(std::string("D2"), D2));
	labels2Int.insert(std::pair< std::string, int >(std::string("D3"), D3));
	labels2Int.insert(std::pair< std::string, int >(std::string("D4"), D4));
	labels2Int.insert(std::pair< std::string, int >(std::string("D5"), D5));
	labels2Int.insert(std::pair< std::string, int >(std::string("D6"), D6));
	labels2Int.insert(std::pair< std::string, int >(std::string("E1"), E1));
	labels2Int.insert(std::pair< std::string, int >(std::string("E2"), E2));
	labels2Int.insert(std::pair< std::string, int >(std::string("E3"), E3));
	labels2Int.insert(std::pair< std::string, int >(std::string("E4"), E4));
	labels2Int.insert(std::pair< std::string, int >(std::string("E5"), E5));
	labels2Int.insert(std::pair< std::string, int >(std::string("E6"), E6));
	
	if(avgTopDist.size())
		avgTopDist.clear();
	avgTopDist.insert(std::pair< int, double >(Alpha, 479));
	avgTopDist.insert(std::pair< int, double >(A1, 455));
	avgTopDist.insert(std::pair< int, double >(A2, 467));
	avgTopDist.insert(std::pair< int, double >(A3, 485));
	avgTopDist.insert(std::pair< int, double >(A4, 489));
	avgTopDist.insert(std::pair< int, double >(Beta, 472));
	avgTopDist.insert(std::pair< int, double >(B1, 481));
	avgTopDist.insert(std::pair< int, double >(B2, 490));
	avgTopDist.insert(std::pair< int, double >(B3, 490));
	avgTopDist.insert(std::pair< int, double >(B4, 501));
	avgTopDist.insert(std::pair< int, double >(Gamma, 478));
	avgTopDist.insert(std::pair< int, double >(C1, 478));
	avgTopDist.insert(std::pair< int, double >(C2, 496));
	avgTopDist.insert(std::pair< int, double >(C3, 534));
	avgTopDist.insert(std::pair< int, double >(C4, 557));
	avgTopDist.insert(std::pair< int, double >(C5, 549));
	avgTopDist.insert(std::pair< int, double >(C6, 526));
	avgTopDist.insert(std::pair< int, double >(Delta, 506));
	avgTopDist.insert(std::pair< int, double >(D1, 491));
	avgTopDist.insert(std::pair< int, double >(D2, 526));
	avgTopDist.insert(std::pair< int, double >(D3, 552));
	avgTopDist.insert(std::pair< int, double >(D4, 556));
	avgTopDist.insert(std::pair< int, double >(D5, 567));
	avgTopDist.insert(std::pair< int, double >(D6, 549));
	avgTopDist.insert(std::pair< int, double >(E1, 546));
	avgTopDist.insert(std::pair< int, double >(E2, 557));
	avgTopDist.insert(std::pair< int, double >(E3, 549));
	avgTopDist.insert(std::pair< int, double >(E4, 580));
	avgTopDist.insert(std::pair< int, double >(E5, 566));
	avgTopDist.insert(std::pair< int, double >(E6, 545));
	if(avgPiaWMDist.size())
		avgPiaWMDist.clear();
	avgPiaWMDist.insert(std::pair< int, double >(Alpha, 1600));
	avgPiaWMDist.insert(std::pair< int, double >(A1, 1651));
	avgPiaWMDist.insert(std::pair< int, double >(A2, 1759));
	avgPiaWMDist.insert(std::pair< int, double >(A3, 1825));
	avgPiaWMDist.insert(std::pair< int, double >(A4, 1916));
	avgPiaWMDist.insert(std::pair< int, double >(Beta, 1623));
	avgPiaWMDist.insert(std::pair< int, double >(B1, 1736));
	avgPiaWMDist.insert(std::pair< int, double >(B2, 1815));
	avgPiaWMDist.insert(std::pair< int, double >(B3, 1899));
	avgPiaWMDist.insert(std::pair< int, double >(B4, 1961));
	avgPiaWMDist.insert(std::pair< int, double >(Gamma, 1713));
	avgPiaWMDist.insert(std::pair< int, double >(C1, 1800));
	avgPiaWMDist.insert(std::pair< int, double >(C2, 1892));
	avgPiaWMDist.insert(std::pair< int, double >(C3, 1985));
	avgPiaWMDist.insert(std::pair< int, double >(C4, 2038));
	avgPiaWMDist.insert(std::pair< int, double >(C5, 2036));
	avgPiaWMDist.insert(std::pair< int, double >(C6, 2064));
	avgPiaWMDist.insert(std::pair< int, double >(Delta, 1845));
	avgPiaWMDist.insert(std::pair< int, double >(D1, 1865));
	avgPiaWMDist.insert(std::pair< int, double >(D2, 1957));
	avgPiaWMDist.insert(std::pair< int, double >(D3, 2046));
	avgPiaWMDist.insert(std::pair< int, double >(D4, 2081));
	avgPiaWMDist.insert(std::pair< int, double >(D5, 2071));
	avgPiaWMDist.insert(std::pair< int, double >(D6, 2087));
	avgPiaWMDist.insert(std::pair< int, double >(E1, 1977));
	avgPiaWMDist.insert(std::pair< int, double >(E2, 2096));
	avgPiaWMDist.insert(std::pair< int, double >(E3, 2117));
	avgPiaWMDist.insert(std::pair< int, double >(E4, 2111));
	avgPiaWMDist.insert(std::pair< int, double >(E5, 2087));
	avgPiaWMDist.insert(std::pair< int, double >(E6, 2119));
	if(avgBarrelArea.size())
		avgBarrelArea.clear();
	avgBarrelArea.insert(std::pair< int, double >(Alpha, 8.79E4));
	avgBarrelArea.insert(std::pair< int, double >(A1, 8.55E4));
	avgBarrelArea.insert(std::pair< int, double >(A2, 8.30E4));
	avgBarrelArea.insert(std::pair< int, double >(A3, 6.48E4));
	avgBarrelArea.insert(std::pair< int, double >(A4, 6.50E4));
	avgBarrelArea.insert(std::pair< int, double >(Beta, 10.2E4));
	avgBarrelArea.insert(std::pair< int, double >(B1, 8.74E4));
	avgBarrelArea.insert(std::pair< int, double >(B2, 9.01E4));
	avgBarrelArea.insert(std::pair< int, double >(B3, 7.66E4));
	avgBarrelArea.insert(std::pair< int, double >(B4, 7.89E4));
	avgBarrelArea.insert(std::pair< int, double >(Gamma, 12.7E4));
	avgBarrelArea.insert(std::pair< int, double >(C1, 10.1E4));
	avgBarrelArea.insert(std::pair< int, double >(C2, 10.3E4));
	avgBarrelArea.insert(std::pair< int, double >(C3, 10.4E4));
	avgBarrelArea.insert(std::pair< int, double >(C4, 9.03E4));
	avgBarrelArea.insert(std::pair< int, double >(C5, 7.06E4));
	avgBarrelArea.insert(std::pair< int, double >(C6, 5.63E4));
	avgBarrelArea.insert(std::pair< int, double >(Delta, 14.3E4));
	avgBarrelArea.insert(std::pair< int, double >(D1, 11.2E4));
	avgBarrelArea.insert(std::pair< int, double >(D2, 12.4E4));
	avgBarrelArea.insert(std::pair< int, double >(D3, 11.6E4));
	avgBarrelArea.insert(std::pair< int, double >(D4, 10.3E4));
	avgBarrelArea.insert(std::pair< int, double >(D5, 8.54E4));
	avgBarrelArea.insert(std::pair< int, double >(D6, 6.40E4));
	avgBarrelArea.insert(std::pair< int, double >(E1, 14.6E4));
	avgBarrelArea.insert(std::pair< int, double >(E2, 15.9E4));
	avgBarrelArea.insert(std::pair< int, double >(E3, 15.8E4));
	avgBarrelArea.insert(std::pair< int, double >(E4, 12.3E4));
	avgBarrelArea.insert(std::pair< int, double >(E5, 8.89E4));
	avgBarrelArea.insert(std::pair< int, double >(E6, 7.31E4));
	
	// set up barrel grid
	// i.e., list of neighbor
	// barrel IDs
	std::list< int > alpha;
	alpha.push_back(A1), alpha.push_back(B1), alpha.push_back(Beta);
	std::list< int > beta;
	beta.push_back(Alpha), beta.push_back(B1), beta.push_back(C1), beta.push_back(Gamma);
	std::list< int > gamma;
	gamma.push_back(Beta), gamma.push_back(C1), gamma.push_back(D1), gamma.push_back(Delta);
	std::list< int > delta;
	delta.push_back(Gamma), delta.push_back(D1), delta.push_back(E1);
	std::list< int > a1;
	a1.push_back(Alpha), a1.push_back(B1), a1.push_back(B2), a1.push_back(A2);
	std::list< int > a2;
	a2.push_back(A1), a2.push_back(B1), a2.push_back(B2), a2.push_back(B3), a2.push_back(A3);
	std::list< int > a3;
	a3.push_back(A2), a3.push_back(B2), a3.push_back(B3), a3.push_back(B4), a3.push_back(A4);
	std::list< int > a4;
	a4.push_back(A3), a4.push_back(B3), a4.push_back(B4);
	std::list< int > b1;
	b1.push_back(Beta), b1.push_back(C1), b1.push_back(C2), b1.push_back(B2), b1.push_back(A2), b1.push_back(A1), b1.push_back(Alpha);
	std::list< int > b2;
	b2.push_back(B1), b2.push_back(C1), b2.push_back(C2), b2.push_back(C3), b2.push_back(B3), b2.push_back(A3), b2.push_back(A2), b2.push_back(A1);
	std::list< int > b3;
	b3.push_back(B2), b3.push_back(C2), b3.push_back(C3), b3.push_back(C4), b3.push_back(B4), b3.push_back(A4), b3.push_back(A3), b3.push_back(A2);
	std::list< int > b4;
	b4.push_back(B3), b4.push_back(C3), b4.push_back(C4), b4.push_back(C5), b4.push_back(A4), b4.push_back(A3);
	std::list< int > c1;
	c1.push_back(Gamma), c1.push_back(D1), c1.push_back(D2), c1.push_back(C2), c1.push_back(B2), c1.push_back(B1), c1.push_back(Beta);
	std::list< int > c2;
	c2.push_back(C1), c2.push_back(D1), c2.push_back(D2), c2.push_back(D3), c2.push_back(C3), c2.push_back(B3), c2.push_back(B2), c2.push_back(B1);
	std::list< int > c3;
	c3.push_back(C2), c3.push_back(D2), c3.push_back(D3), c3.push_back(D4), c3.push_back(C4), c3.push_back(B4), c3.push_back(B3), c3.push_back(B2);
	std::list< int > c4;
	c4.push_back(C3), c4.push_back(D3), c4.push_back(D4), c4.push_back(D5), c4.push_back(C5), c4.push_back(B4), c4.push_back(B3);
	std::list< int > c5;
	c5.push_back(C4), c5.push_back(D4), c5.push_back(D5), c5.push_back(D6), c5.push_back(C6), c5.push_back(B4);
	std::list< int > c6;
	c6.push_back(C5), c6.push_back(D5), c6.push_back(D6);
	std::list< int > d1;
	d1.push_back(Delta), d1.push_back(E1), d1.push_back(E2), d1.push_back(D2), d1.push_back(C2), d1.push_back(C1), d1.push_back(Gamma);
	std::list< int > d2;
	d2.push_back(D1), d2.push_back(E1), d2.push_back(E2), d2.push_back(E3), d2.push_back(D3), d2.push_back(C3), d2.push_back(C2), d2.push_back(C1);
	std::list< int > d3;
	d3.push_back(D2), d3.push_back(E2), d3.push_back(E3), d3.push_back(E4), d3.push_back(D4), d3.push_back(C4), d3.push_back(C3), d3.push_back(C2);
	std::list< int > d4;
	d4.push_back(D3), d4.push_back(E3), d4.push_back(E4), d4.push_back(E5), d4.push_back(D5), d4.push_back(C5), d4.push_back(C4), d4.push_back(C3);
	std::list< int > d5;
	d5.push_back(D4), d5.push_back(E4), d5.push_back(E5), d5.push_back(E6), d5.push_back(D6), d5.push_back(C6), d5.push_back(C5), d5.push_back(C4);
	std::list< int > d6;
	d6.push_back(D5), d6.push_back(E5), d6.push_back(E6), d6.push_back(C6), d6.push_back(C5);
	std::list< int > e1;
	e1.push_back(Delta), e1.push_back(E2), e1.push_back(D2), e1.push_back(D1);
	std::list< int > e2;
	e2.push_back(E1), e2.push_back(E3), e2.push_back(D3), e2.push_back(D2), e2.push_back(D1);
	std::list< int > e3;
	e3.push_back(E2), e3.push_back(E4), e3.push_back(D4), e3.push_back(D3), e3.push_back(D2);
	std::list< int > e4;
	e4.push_back(E3), e4.push_back(E5), e4.push_back(D5), e4.push_back(D4), e4.push_back(D3);
	std::list< int > e5;
	e5.push_back(E4), e5.push_back(E6), e5.push_back(D6), e5.push_back(D5), e5.push_back(D4);
	std::list< int > e6;
	e6.push_back(E5), e6.push_back(D6), e6.push_back(D5);
	
	barrelGrid.insert(std::pair< int, std::list< int > >(Alpha, alpha));
	barrelGrid.insert(std::pair< int, std::list< int > >(Beta, beta));
	barrelGrid.insert(std::pair< int, std::list< int > >(Gamma, gamma));
	barrelGrid.insert(std::pair< int, std::list< int > >(Delta, delta));
	barrelGrid.insert(std::pair< int, std::list< int > >(A1, a1));
	barrelGrid.insert(std::pair< int, std::list< int > >(A2, a2));
	barrelGrid.insert(std::pair< int, std::list< int > >(A3, a3));
	barrelGrid.insert(std::pair< int, std::list< int > >(A4, a4));
	barrelGrid.insert(std::pair< int, std::list< int > >(B1, b1));
	barrelGrid.insert(std::pair< int, std::list< int > >(B2, b2));
	barrelGrid.insert(std::pair< int, std::list< int > >(B3, b3));
	barrelGrid.insert(std::pair< int, std::list< int > >(B4, b4));
	barrelGrid.insert(std::pair< int, std::list< int > >(C1, c1));
	barrelGrid.insert(std::pair< int, std::list< int > >(C2, c2));
	barrelGrid.insert(std::pair< int, std::list< int > >(C3, c3));
	barrelGrid.insert(std::pair< int, std::list< int > >(C4, c4));
	barrelGrid.insert(std::pair< int, std::list< int > >(C5, c5));
	barrelGrid.insert(std::pair< int, std::list< int > >(C6, c6));
	barrelGrid.insert(std::pair< int, std::list< int > >(D1, d1));
	barrelGrid.insert(std::pair< int, std::list< int > >(D2, d2));
	barrelGrid.insert(std::pair< int, std::list< int > >(D3, d3));
	barrelGrid.insert(std::pair< int, std::list< int > >(D4, d4));
	barrelGrid.insert(std::pair< int, std::list< int > >(D5, d5));
	barrelGrid.insert(std::pair< int, std::list< int > >(D6, d6));
	barrelGrid.insert(std::pair< int, std::list< int > >(E1, e1));
	barrelGrid.insert(std::pair< int, std::list< int > >(E2, e2));
	barrelGrid.insert(std::pair< int, std::list< int > >(E3, e3));
	barrelGrid.insert(std::pair< int, std::list< int > >(E4, e4));
	barrelGrid.insert(std::pair< int, std::list< int > >(E5, e5));
	barrelGrid.insert(std::pair< int, std::list< int > >(E6, e6));
	
// 	if(avgTopDist.size())
// 		avgTopDist.clear();
// 	avgTopDist.insert(std::pair< int, double >(Alpha, 496));
// 	avgTopDist.insert(std::pair< int, double >(A1, 459));
// 	avgTopDist.insert(std::pair< int, double >(A2, 480));
// 	avgTopDist.insert(std::pair< int, double >(A3, 507));
// 	avgTopDist.insert(std::pair< int, double >(A4, 501));
// 	avgTopDist.insert(std::pair< int, double >(Beta, 492));
// 	avgTopDist.insert(std::pair< int, double >(B1, 489));
// 	avgTopDist.insert(std::pair< int, double >(B2, 491));
// 	avgTopDist.insert(std::pair< int, double >(B3, 469));
// 	avgTopDist.insert(std::pair< int, double >(B4, 491));
// 	avgTopDist.insert(std::pair< int, double >(Gamma, 487));
// 	avgTopDist.insert(std::pair< int, double >(C1, 485));
// 	avgTopDist.insert(std::pair< int, double >(C2, 487));
// 	avgTopDist.insert(std::pair< int, double >(C3, 508));
// 	avgTopDist.insert(std::pair< int, double >(C4, 531));
// 	avgTopDist.insert(std::pair< int, double >(C5, 544));
// 	avgTopDist.insert(std::pair< int, double >(C6, 506));
// 	avgTopDist.insert(std::pair< int, double >(Delta, 525));
// 	avgTopDist.insert(std::pair< int, double >(D1, 491));
// 	avgTopDist.insert(std::pair< int, double >(D2, 544));
// 	avgTopDist.insert(std::pair< int, double >(D3, 547));
// 	avgTopDist.insert(std::pair< int, double >(D4, 549));
// 	avgTopDist.insert(std::pair< int, double >(D5, 541));
// 	avgTopDist.insert(std::pair< int, double >(D6, 559));
// 	avgTopDist.insert(std::pair< int, double >(E1, 552));
// 	avgTopDist.insert(std::pair< int, double >(E2, 568));
// 	avgTopDist.insert(std::pair< int, double >(E3, 530));
// 	avgTopDist.insert(std::pair< int, double >(E4, 571));
// 	avgTopDist.insert(std::pair< int, double >(E5, 533));
// 	avgTopDist.insert(std::pair< int, double >(E6, 538));
// 	if(avgPiaWMDist.size())
// 		avgPiaWMDist.clear();
// 	avgPiaWMDist.insert(std::pair< int, double >(Alpha, 1573));
// 	avgPiaWMDist.insert(std::pair< int, double >(A1, 1630));
// 	avgPiaWMDist.insert(std::pair< int, double >(A2, 1754));
// 	avgPiaWMDist.insert(std::pair< int, double >(A3, 1805));
// 	avgPiaWMDist.insert(std::pair< int, double >(A4, 1856));
// 	avgPiaWMDist.insert(std::pair< int, double >(Beta, 1578));
// 	avgPiaWMDist.insert(std::pair< int, double >(B1, 1726));
// 	avgPiaWMDist.insert(std::pair< int, double >(B2, 1792));
// 	avgPiaWMDist.insert(std::pair< int, double >(B3, 1828));
// 	avgPiaWMDist.insert(std::pair< int, double >(B4, 1911));
// 	avgPiaWMDist.insert(std::pair< int, double >(Gamma, 1681));
// 	avgPiaWMDist.insert(std::pair< int, double >(C1, 1786));
// 	avgPiaWMDist.insert(std::pair< int, double >(C2, 1851));
// 	avgPiaWMDist.insert(std::pair< int, double >(C3, 1929));
// 	avgPiaWMDist.insert(std::pair< int, double >(C4, 1966));
// 	avgPiaWMDist.insert(std::pair< int, double >(C5, 1997));
// 	avgPiaWMDist.insert(std::pair< int, double >(C6, 2017));
// 	avgPiaWMDist.insert(std::pair< int, double >(Delta, 1853));
// 	avgPiaWMDist.insert(std::pair< int, double >(D1, 1825));
// 	avgPiaWMDist.insert(std::pair< int, double >(D2, 1938));
// 	avgPiaWMDist.insert(std::pair< int, double >(D3, 2006));
// 	avgPiaWMDist.insert(std::pair< int, double >(D4, 2054));
// 	avgPiaWMDist.insert(std::pair< int, double >(D5, 2017));
// 	avgPiaWMDist.insert(std::pair< int, double >(D6, 2076));
// 	avgPiaWMDist.insert(std::pair< int, double >(E1, 1957));
// 	avgPiaWMDist.insert(std::pair< int, double >(E2, 2079));
// 	avgPiaWMDist.insert(std::pair< int, double >(E3, 2077));
// 	avgPiaWMDist.insert(std::pair< int, double >(E4, 2086));
// 	avgPiaWMDist.insert(std::pair< int, double >(E5, 2031));
// 	avgPiaWMDist.insert(std::pair< int, double >(E6, 2071));
// 	if(avgBarrelArea.size())
// 		avgBarrelArea.clear();
// 	avgBarrelArea.insert(std::pair< int, double >(Alpha, 7.99E4));
// 	avgBarrelArea.insert(std::pair< int, double >(A1, 7.49E4));
// 	avgBarrelArea.insert(std::pair< int, double >(A2, 8.13E4));
// 	avgBarrelArea.insert(std::pair< int, double >(A3, 6.15E4));
// 	avgBarrelArea.insert(std::pair< int, double >(A4, 6.14E4));
// 	avgBarrelArea.insert(std::pair< int, double >(Beta, 9.52E4));
// 	avgBarrelArea.insert(std::pair< int, double >(B1, 7.96E4));
// 	avgBarrelArea.insert(std::pair< int, double >(B2, 8.16E4));
// 	avgBarrelArea.insert(std::pair< int, double >(B3, 7.60E4));
// 	avgBarrelArea.insert(std::pair< int, double >(B4, 7.84E4));
// 	avgBarrelArea.insert(std::pair< int, double >(Gamma, 11.0E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C1, 9.73E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C2, 10.2E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C3, 9.77E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C4, 8.61E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C5, 7.33E4));
// 	avgBarrelArea.insert(std::pair< int, double >(C6, 5.82E4));
// 	avgBarrelArea.insert(std::pair< int, double >(Delta, 12.0E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D1, 10.8E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D2, 11.7E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D3, 11.0E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D4, 10.1E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D5, 8.10E4));
// 	avgBarrelArea.insert(std::pair< int, double >(D6, 6.26E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E1, 13.3E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E2, 15.1E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E3, 15.0E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E4, 11.9E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E5, 8.76E4));
// 	avgBarrelArea.insert(std::pair< int, double >(E6, 6.70E4));

};
