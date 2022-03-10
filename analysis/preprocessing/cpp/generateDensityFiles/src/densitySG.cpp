/****************************************************************************/
/*                                                                          */
/* File:    densitySG.cpp                                                 	*/
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
/* Date: 	17.10.2018                                                  	*/
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "densitySG.h"

// Constructor with AmiraSpatialGraph
densitySG::densitySG(AmiraSpatialGraph* sg, int label)
{
	this->sg = sg; 
	VoxelSize = 50.0;
	this->label = label; 
	outputpath = ""; 
	epsilon = 1E-5; 
	boolAmirafile = false;
	double sgBounds[6];
	sg->getBoundingBox(label,sgBounds);
	extendBounds(sgBounds);

	if (!sg->isLabelInSpatialGraph(label))
	{
		std::cout << "WARNING! Label " << helper::getNeuriteName(label) << " not found in SpatialGraph!" << std::endl;
	}
}

// Constructor with Path/To/AmiraSpatialGraph
densitySG::densitySG(const char* inputfilename, int label)
{
	Reader * fileReader = new Reader(inputfilename, inputfilename);
	
	std::string fname(inputfilename); 
	if((access(fname.c_str(), F_OK) != -1))
	{
		if (fname.compare(fname.size()-3,3,".am") == 0)
		{
			fileReader->readSpatialGraphFile(0); 
			boolAmirafile = true;
		}
		else if (fname.compare(fname.size()-4,4,".hoc") == 0)
		{
			fileReader->readHocFile();
			boolAmirafile = false;
		}
		else
		{
			std::cout << "Inputfile is neither .hoc nor .am! Cannot read file!" << std::endl; 
		}
		
		sg = fileReader->getSpatialGraph();
	}
	else
	{
		std::cout << "Inputfile not found: " << fname << std::endl; 
	}
  
	VoxelSize = 50.0;
	this->label = label; 
	outputpath = ""; 
	epsilon = 1E-5; 
	double sgBounds[6];
	sg->getBoundingBox(label,sgBounds);
	extendBounds(sgBounds);

	if (!sg->isLabelInSpatialGraph(label))
	{
		std::cout << "WARNING! Label " << helper::getNeuriteName(label) << " not found in SpatialGraph!" << std::endl;
	}
}

ImageDataPointerType densitySG::computeLength()
{
	std::cout << " -- Start Computing Length -- " << std::endl; 
	printProperties();
	
	// For Checking, get total original total length and number of branching points
	double org_length = getLengthInBoxSimpleTotal(bounds); 
	double total_len = 0.0; 

	// Stores Length
	ImageDataPointerType volLength = helper::createImageVolumeNeuroNet(bounds,VoxelSize);
	double origin[3];
	volLength->GetOrigin(origin);

	// Add Intersecting Points to Graph
	addIntersectionPts(volLength);
	
	for(int x = volLength->GetExtent()[0]; x <= volLength->GetExtent()[1]; ++x)
	{  
		for(int y = volLength->GetExtent()[2]; y <= volLength->GetExtent()[3]; ++y)
		{
			for(int z = volLength->GetExtent()[4]; z <= volLength->GetExtent()[5]; ++z)
			{
				// Coordinates x,y,z give center point and not boundingbox limits
				double boundsVoxel[6];
				getVoxelBoundingBox(boundsVoxel,origin,x,y,z,VoxelSize);

				double * len = static_cast< double * >(volLength->GetScalarPointer(x, y, z));
				*len = getLengthInBoxSimple(boundsVoxel);
				total_len += *len; 
			}
		}
	}
	
	volLength->Update();
	
	/* Controlling measurements and computation! If there are deviations, report them!
	 * Get BB of whole spatialGraph (with added intersection pts), check and compare both length measurements */
	double boundsOriginal[6];
	sg->getBoundingBox(label,boundsOriginal);
	double act_len = getLengthInBoxSimpleTotal(boundsOriginal);
		
	if ((fabs(total_len-act_len)>1E-3) || (fabs(org_length-total_len)>1E-3))
	{
		std::cout << "  WARNING! Problem with calculated length!" << std::endl;
		std::cout << "    Actual Length: " << act_len << " - Total Length in Voxels: ";
		std::cout << total_len << " - Original Length: " << org_length << std::endl;
		std::cout << "    Difference: " << fabs(total_len-act_len) << std::endl; 
	}
	
	if (total_len==0.0 || act_len==0.0 || org_length==0.0)
	{
		std::cout << "  WARNING! Length is equal to 0.0! Maybe label does not exist!" << std::endl;
	}
	return volLength;
}

// Calculates the Length of all points within box (Does not take into account intersections!
// Intersecting points have to be computed and added to the SpatialGraph previously)
double densitySG::getLengthInBoxSimpleTotal(double input_bounds[6])
{	
	return getLengthInBoxSimple(input_bounds, true);
}

double densitySG::getLengthInBoxSimple(double input_bounds[6])
{
	return getLengthInBoxSimple(input_bounds, false);
}

// Calculates the Length of all points within box (Does not take into account intersections!
// Intersecting points have to be computed and added to the SpatialGraph previously)
// Prints message if points are outside of bounds
double densitySG::getLengthInBoxSimple(double input_bounds[6], bool showErrorMessage)
{
	// length of edges within clippedBox
	double len = 0; 
	
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = sg->edgesBegin(); edgeIt != sg->edgesEnd(); ++edgeIt)
	{
		// Check only certain label
		if((*edgeIt)->label == label)
		{   
			std::list< double * >::iterator edgeListIt;
			edgeListIt = (*edgeIt)->edgePointCoordinates.begin();
			double * previousPt = *edgeListIt;
			++edgeListIt;
			
			// Go through all points of each Edge starting with the second edge
			for(; edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{  
				double * currentPt = *edgeListIt;
				
				// Check whether Current Coordinates are within clippedBox
				if (inBox(input_bounds,currentPt[X_COORD],currentPt[Y_COORD],currentPt[Z_COORD])
						&& inBox(input_bounds,previousPt[X_COORD],previousPt[Y_COORD],previousPt[Z_COORD]))
				{  
					len += sqrt(pow((currentPt[X_COORD]-previousPt[X_COORD]),2.0) +
								pow((currentPt[Y_COORD]-previousPt[Y_COORD]),2.0) +
								pow((currentPt[Z_COORD]-previousPt[Z_COORD]),2.0));
				}
				else
				{
					if (showErrorMessage)
					{
						std::cout << "  WARNING! Point is not in BoundingBox!" << std::endl;
						std::cout << "   BoundingBox = [" << bounds[0] << " " << bounds[1] << "; ";
						std::cout << bounds[2] << " " << bounds[3] << "; " << bounds[4];
						std::cout << " " << bounds[5] << "]" << std::endl;
						std::cout << "   Current:   " << currentPt[X_COORD] << ", " << currentPt[Y_COORD];
						std::cout << ", " << currentPt[Z_COORD] << std::endl;
						std::cout << "   Previous:   " << previousPt[X_COORD] << ", " << previousPt[Y_COORD];
						std::cout << ", " << previousPt[Z_COORD] << std::endl;
					}
				}
				previousPt = currentPt;
			}
		}
	}
	return len;  
}

// Checks whether Point is within defined bounding box
bool densitySG::inBox(double input_bounds[6], double x, double y, double z)
{
	bool isInBox = ( x >= input_bounds[0]-epsilon &&
					x <= input_bounds[1]+epsilon &&
					y >= input_bounds[2]-epsilon &&
					y <= input_bounds[3]+epsilon &&
					z >= input_bounds[4]-epsilon &&
					z <= input_bounds[5]+epsilon );
	return isInBox; 
}

void densitySG::addIntersectionPts(double input_bounds[6])
{
	ImageDataPointerType volume = helper::createImageVolume(input_bounds, VoxelSize);
	addIntersectionPts(volume);
}

// Add Intersection Points to Spatial Graph in one Voxel defined by bounds!
bool densitySG::addIntersectionPtsInVoxel(double boundsVoxel[6])
{
	vtkBox * BBx = vtkBox::New();
	BBx->SetBounds(boundsVoxel);
	bool sgWithinBounds = false;

	for(std::vector< Edge * >::iterator edgeIt = sg->edgesBegin(); edgeIt != sg->edgesEnd(); ++edgeIt)
	{
		std::list< double * >::iterator edgeListIt;

		// Create Iterator, current one (+1) and previous one (0)
		edgeListIt = (*edgeIt)->edgePointCoordinates.begin();
		double * previousPt = *edgeListIt;
		++edgeListIt;

		// Interpolate Radius
		std::list< double >::iterator edgeRadiusListIt;
		edgeRadiusListIt = (*edgeIt)->radiusList.begin();
		double previousRad = *edgeRadiusListIt;
		++edgeRadiusListIt;

		if ((*edgeIt)->numEdgePoints==1)
		{
			std::cout << "  WARNING! Only One Edge Point Found! " << std::endl;
			std::cout << "    EdgePt: " << previousPt[X_COORD] << ", " << previousPt[Y_COORD];
			std::cout << ", " << previousPt[Z_COORD] << std::endl;
			std::cout << "    EdgeNo: " << (*edgeIt)->fatherID << std::endl;
		}

		bool radiusExist = true;
		if ((*edgeIt)->radiusList.empty() && radiusExist)
		{
			std::cout << "  WARNING! Radius List is empty!" << std::endl;
			radiusExist = false;
		}

		// Go through all points of each Edge starting with the second edge
		for(; edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt, ++edgeRadiusListIt)
		{
			double * currentPt = *edgeListIt;
			double currentRad = *edgeRadiusListIt;

			// Compute Intersection points
			double t1, t2, x1[3], x2[3];
			int plane1,plane2;
			int b = BBx->IntersectWithLine(boundsVoxel,previousPt,currentPt,t1,t2,x1,x2,plane1,plane2);

			// b is zero with line is totally outside of box
			// b is 1 if both pts are inside the box, if the intersect the box or
			// if they are on the edge of the box
			if (b != 0)
			{
				sgWithinBounds = true;

				if (t1>0 && t1<1)
				{
					double * pointerList = new double [3];
					for (int ii=0; ii<3; ii++)
						pointerList[ii] = x1[ii];

					// Add To Edge
					// Increment Number of Edge Points by 1
					(*edgeIt)->numEdgePoints += 1;
					// Insert Coordinates in List
					(*edgeIt)->edgePointCoordinates.insert(edgeListIt,pointerList);

					if (radiusExist)
					{
						double newRadius = currentRad + (previousRad - currentRad) * t1;
						(*edgeIt)->radiusList.insert(edgeRadiusListIt, newRadius);
					}

					//std::cout << "  x1 = [" << x1[0] << " " << x1[1] << " " << x1[2] << "] t1 = " << t1 << " p1 = " << plane1 << std::endl;
				}
				if (t2>0 && t2<1)
				{
					double * pointerList = new double [3];
					for (int ii=0; ii<3; ii++)
						pointerList[ii] = x2[ii];

					// Add To Edge
					// Increment Number of Edge Points by 1
					(*edgeIt)->numEdgePoints += 1;
					// Insert Coordinates in List
					(*edgeIt)->edgePointCoordinates.insert(edgeListIt,pointerList);

					if (radiusExist)
					{
						double newRadius = currentRad + (previousRad - currentRad) * t2;
						(*edgeIt)->radiusList.insert(edgeRadiusListIt, newRadius);
					}
					//std::cout << "  x2 = [" << x2[0] << " " << x2[1] << " " << x2[2] << "] t2 = " << t2 << " p2 = " << plane2 << std::endl;
				}
			}
			previousPt = currentPt;
			previousRad = currentRad;
		}
	}

	return sgWithinBounds;
}

void densitySG::addIntersectionPts(ImageDataPointerType volume)
{
	// Extract Center of BoundingBox Voxels of Volume
	double spacing[3];
	volume->GetSpacing(spacing);
	double origin[3];
	volume->GetOrigin(origin);

	if (spacing[0] != spacing[1] || spacing[1] != spacing [2])
	{
		std::cout << "ERROR! Spacing of Volume is not equal!" << std::endl;
		return;
	}

	if (spacing[0] != VoxelSize)
	{
		VoxelSize = spacing[0];
		std::cout << "WARNING! VoxelSize updated to " << VoxelSize << "!" << std::endl;
	}

	// Iterate over all Boxes
	// range can be negative to positive
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{
		for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
		{
			for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
			{
				double boundsVoxel[6];
				getVoxelBoundingBox(boundsVoxel,origin,x,y,z,VoxelSize);
				bool tmp = addIntersectionPtsInVoxel(boundsVoxel);
			}
		}
	}
}

void densitySG::getVoxelBoundingBox(double boundsVoxel[6],double origin[3],
							int x,int y,int z, double InputVoxelSize)
{
	boundsVoxel[0] = origin[0] + x*InputVoxelSize-0.5*InputVoxelSize;
	boundsVoxel[1] = origin[0] + (x+1)*InputVoxelSize-0.5*InputVoxelSize;
	boundsVoxel[2] = origin[1] + y*InputVoxelSize-0.5*InputVoxelSize;
	boundsVoxel[3] = origin[1] + (y+1)*InputVoxelSize-0.5*InputVoxelSize;
	boundsVoxel[4] = origin[2] + z*InputVoxelSize-0.5*InputVoxelSize;
	boundsVoxel[5] = origin[2] + (z+1)*InputVoxelSize-0.5*InputVoxelSize;
}

// Returns integer in which box x is based on bounding box coordinates
int densitySG::ptInBox(double x, double minX, double maxX)
{
	return std::floor((x - minX)/(maxX - minX));
}

ImageDataPointerType densitySG::computeTP()
{
	std::cout << " -- Start Terminal Points (needs VALIDATION!) -- " << std::endl; 
	printProperties();
	
	if (boolAmirafile)
		std::cout << "      Warning! Assumes TPs are at the end of the spatialGraph edges!" << std::endl; 
	
	double org_TermPts = getTerminalPtInBox(bounds); 
	double total_TermPts = 0.0; 
	
	ImageDataPointerType volTP = helper::createImageVolumeNeuroNet(bounds,VoxelSize);
	double origin[3];
	volTP->GetOrigin(origin);

	for(int x = volTP->GetExtent()[0]; x <= volTP->GetExtent()[1]; ++x)
	{  
		for(int y = volTP->GetExtent()[2]; y <= volTP->GetExtent()[3]; ++y)
		{
			for(int z = volTP->GetExtent()[4]; z <= volTP->GetExtent()[5]; ++z)
			{
				// Coordinates x,y,z give center point and not boundingbox limits
				double boundsVoxel[6];
				getVoxelBoundingBox(boundsVoxel,origin,x,y,z,VoxelSize);

				double * terminalPointer = static_cast< double * >(volTP->GetScalarPointer(x, y, z));
				*terminalPointer = getTerminalPtInBox(boundsVoxel);
				total_TermPts += *terminalPointer; 
			}
		}
	}
	
	volTP->Update();
	
	if (fabs(org_TermPts-total_TermPts)>0)
	{
		std::cout << "  Warning! Problem with Number of Terminal Points!" << std::endl;
		std::cout << "    Total Number of Terminal Points in Voxels: " << total_TermPts;
		std::cout << " - Original Number of Terminal Points: " << org_TermPts << std::endl;
	}
	if (total_TermPts == 0 || org_TermPts == 0)
	{
		std::cout << "  Warning! No Terminal Points were found! Maybe label does not exist!" << std::endl;  
		std::cout << "    Total Number of Terminal Points in Voxels: " << total_TermPts;
		std::cout << " - Original Number of Terminal Points: " << org_TermPts << std::endl;
	}
	
	return volTP; 
}

// Compute theNumber of Terminal Points within specified box
double densitySG::getTerminalPtInBox(double input_bounds[6])
{  
	double NoOfTerminalPts = 0; 
	
	if  (!boolAmirafile)
	{
		std::list< int > branchPtList = getBranchPointIDs(); 
		
		// Go through all EdgeIDs
		for(int ii = 0; ii < sg->getNumberOfEdges(); ++ii)
		{
			bool flagBranch = false; 
			
			// Check whether the current EdgeID is in the list of IDs of Branching Edges
			for(std::list< int >::const_iterator tpListIt = branchPtList.begin();
					tpListIt != branchPtList.end(); ++tpListIt)
			{
				if (ii==(*tpListIt))
				{
					flagBranch = true; 
				}
			}  
			
			// If valid label and not a branching edge, then check whether it's in the specified box. 
			if((*(sg->edgesPointer()))[ii]->label == label && !flagBranch)
			{
				double * tp = (*(sg->edgesPointer()))[ii]->edgePointCoordinates.back();
				if (inBox(input_bounds, tp[X_COORD], tp[Y_COORD], tp[Z_COORD]))
				{
					NoOfTerminalPts++; 
				}
			}
		}
	}
	else 	// WARNING! sometimes TerminalNodes are at the back of edgePointCoordinates and sometimes they are at the front.
		// Thus the algorithm has to be changed accordingly. Just switch from back() to front()
	{  
		for(std::vector< Edge * >::const_iterator edgeIt = sg->edgesBegin();
				edgeIt != sg->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label == label)
			{  
// 				double * pt = (*edgeIt)->edgePointCoordinates.front(); 
				double * pt = (*edgeIt)->edgePointCoordinates.back(); 
				// If that point is in box, check whether this point is the last point of at least two other edges
				if (inBox(input_bounds,pt[X_COORD],pt[Y_COORD],pt[Z_COORD]))
				{  
					int children = 0; 
					for(std::vector< Edge * >::const_iterator edgeIt2 = sg->edgesBegin();
							edgeIt2 != sg->edgesEnd(); ++edgeIt2)
					{
// 						double * pt2 = (*edgeIt2)->edgePointCoordinates.back(); 
						double * pt2 = (*edgeIt2)->edgePointCoordinates.front(); 
						
						if(isSamePoint(pt[X_COORD],pt[Y_COORD],pt[Z_COORD],
								pt2[X_COORD],pt2[Y_COORD],pt2[Z_COORD]))
						{
							children++;
						}
					}
					if (children<1)
					{
						NoOfTerminalPts++;
// 						double terminalLocation[3] = {pt[X_COORD], pt[Y_COORD], pt[Z_COORD]}; 
// 						int oldNrOfPoints = terminalPtMarker->GetNumberOfPoints();
// 						terminalPtMarker->InsertPoint(oldNrOfPoints,terminalLocation);
					}
				}
			}
		}
	}
	
// 	std::string ofNameBorder = "/home/dudvary/Interneurons/L23_IN/validDendrites/artificial/completeCell_density/bipolar_artificial/TPLandmark";
// 	Reader * ScalarWriterBorder = new Reader(ofNameBorder.c_str(), ofNameBorder.c_str());
// 	ScalarWriterBorder->writeLandmarkFile(terminalPtMarker);;
// 	delete ScalarWriterBorder;
	
	return NoOfTerminalPts; 
}

// Check whether two points are identical (given a small epsilon range)
bool densitySG::isSamePoint(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return ((fabs(x1-x2)<epsilon) && (fabs(y1-y2)<epsilon) && (fabs(z1-z2)<epsilon)); 
}

// Computes a list of IDs, specifying which Edge ends with a branching point
std::list< int > densitySG::getBranchPointIDs()
{
	std::map< int, int > nrOfChildBranches;
	std::list< int > branchPtIDs;
	std::vector< Edge * >::const_iterator edgeIt;
		
	for(edgeIt = sg->edgesBegin(); edgeIt != sg->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label == label)
		{  
			int fatherID = (*edgeIt)->fatherID;
			if(nrOfChildBranches.find(fatherID) != nrOfChildBranches.end()) // Check whether fatherID exists, if so, increase count by one
			{
				nrOfChildBranches[fatherID]++;
			}
			else // if not, insert fatherID with corresponding count of 1
			{
				nrOfChildBranches.insert(std::pair< int, int >(fatherID, 1));
			}
		}
	}
	
	std::map< int, int >::const_iterator nrBranchIt;
	for(nrBranchIt = nrOfChildBranches.begin(); nrBranchIt != nrOfChildBranches.end(); ++nrBranchIt)
	{
		if(nrBranchIt->second > 1) // If count of fatherID is larger than 1, then this is a branching point, put fatherID in list
		{
			branchPtIDs.push_back(nrBranchIt->first);
		}
	}
	
	return branchPtIDs;
}

ImageDataPointerType densitySG::computeBP()
{
	std::cout << " -- Start Branch Points (needs VALIDATION!) -- " << std::endl; 
	printProperties();
	
	if (boolAmirafile)
		std::cout << "      Warning! Assumes BPs are at the front of the spatialGraph edges!" << std::endl; 
	
	double org_BranchPts = getBranchingPtInBox(bounds); 
	int total_BranchPts = 0;

	ImageDataPointerType volBP = helper::createImageVolumeNeuroNet(bounds,VoxelSize);
	double origin[3];
	volBP->GetOrigin(origin);

	for(int x = volBP->GetExtent()[0]; x <= volBP->GetExtent()[1]; ++x)
	{  
		for(int y = volBP->GetExtent()[2]; y <= volBP->GetExtent()[3]; ++y)
		{
			for(int z = volBP->GetExtent()[4]; z <= volBP->GetExtent()[5]; ++z)
			{
				// Coordinates x,y,z give center point and not boundingbox limits
				double boundsVoxel[6];
				getVoxelBoundingBox(boundsVoxel,origin,x,y,z,VoxelSize);

				double * branchPointer = static_cast< double * >(volBP->GetScalarPointer(x, y, z));
				*branchPointer = getBranchingPtInBox(boundsVoxel);
				total_BranchPts += *branchPointer; 
			}
		}
	}
	
	volBP->Update();
	
	if (fabs(org_BranchPts-total_BranchPts)>0)
	{
		std::cout << "  Warning! Problem with Number of Branching Points!" << std::endl;
		std::cout << "    Total Number of Branching Points in Voxels: " << total_BranchPts;
		std::cout << " - Original Number of Branching Points: " << org_BranchPts << std::endl;
	}
	if (total_BranchPts == 0 || org_BranchPts == 0)
	{
		std::cout << "  Warning! No Branching Points were found! Maybe label does not exist!" << std::endl;  
		std::cout << "    Total Number of Branching Points in Voxels: " << total_BranchPts;
		std::cout << " - Original Number of Branching Points: " << org_BranchPts << std::endl;
	}
	
	return volBP; 
}

// Compute the Number of BranchingPoints within specified box
int densitySG::getBranchingPtInBox(double input_bounds[6])
{
	int NoOfBranchPts = 0; 
	
	if (!boolAmirafile)
	{
		std::list< int > branchPtList = getBranchPointIDs();
		if (!branchPtList.empty())
		{
			std::list< int >::const_iterator bpListIt;
			for(bpListIt = branchPtList.begin(); bpListIt != branchPtList.end(); ++bpListIt)
			{
				if ((*bpListIt)!=-1)
				{
					double * bp = (*(sg->edgesPointer()))[*bpListIt]->edgePointCoordinates.back();				
					if (inBox(input_bounds, bp[X_COORD], bp[Y_COORD], bp[Z_COORD]))
					{
						NoOfBranchPts++; 
					}
				}
			}
		}
	}
	// WARNING! sometimes TerminalNodes are at the back of edgePointCoordinates and sometimes they are at the front.
	// Thus the algorithm has to be changed accordingly. Just switch from back() to front()
	else // If amira file is loaded, branchingIDs won't work!
	{   
		for(std::vector< Edge * >::const_iterator edgeIt = sg->edgesBegin();
				edgeIt != sg->edgesEnd(); ++edgeIt)
		{
			if((*edgeIt)->label == label)
			{  
				double * pt = (*edgeIt)->edgePointCoordinates.back(); 
				
				// If that point is in box, check whether this point is the last point of at least two other edges
				if (inBox(input_bounds,pt[X_COORD],pt[Y_COORD],pt[Z_COORD]))
				{  
					int children = 0; 
					for(std::vector< Edge * >::const_iterator edgeIt2 = sg->edgesBegin();
							edgeIt2 != sg->edgesEnd(); ++edgeIt2)
					{
						double * pt2 = (*edgeIt2)->edgePointCoordinates.front(); 
						
						if(isSamePoint(pt[X_COORD],pt[Y_COORD],pt[Z_COORD],
								pt2[X_COORD],pt2[Y_COORD],pt2[Z_COORD]))
						{
							children++;
						}
					}
					
					if (children>1)
					{
						NoOfBranchPts++;
					}
				}
			}
		}
	}
	return NoOfBranchPts; 
}

void densitySG::computeParas()
{
	ImageDataPointerType volLen = computeLength(); 
	ImageDataPointerType volTP = computeTP();
	ImageDataPointerType volBP = computeBP(); 
	
	size_t idx = outputpath.find_last_of(".am"); 
	std::string ftmp = outputpath.substr(0,idx-2); 
	
	std::cout << " -- Store Length/BranchPoints/TerminalPoints in " << ftmp << std::endl;

	std::string outfname = ftmp + "_Length.am"; 
	helper::storeImageVolume(outfname.c_str(), volLen);
	outfname = ftmp + "_TP.am"; 
	helper::storeImageVolume(outfname.c_str(), volTP);
	outfname = ftmp + "_BP.am"; 
	helper::storeImageVolume(outfname.c_str(), volBP);
}

// Creates wider Bounding Box
void densitySG::extendBounds(double input_bounds[6])
{
	for (int i = 0; i<3; i++)
	{
		bounds[i*2] = (floor(input_bounds[i*2]/VoxelSize) - 1) * VoxelSize; 	// Min Value
		bounds[i*2+1] = (floor(input_bounds[i*2+1]/VoxelSize) + 2) * VoxelSize; // Max value
	}
}

// Display Information
void densitySG::printProperties()
{
	std::cout << "      Properties: VoxelSize = " << VoxelSize << "; Label = " << helper::getNeuriteName(label);
	std::cout << "; AmiraFile = " << boolAmirafile << std::endl;
	std::cout << "                  BoundingBox(AmiraDensity) = [" << bounds[0] << " " << bounds[1] << "; " << bounds[2] << " ";
	std::cout << bounds[3] << "; " << bounds[4] << " " << bounds[5] << "]" << std::endl;

	double sgBounds[6];
	sg->getBoundingBox(label,sgBounds);
	std::cout << "                  BoundingBox(SG) = [" << sgBounds[0] << " " << sgBounds[1] << "; " << sgBounds[2] << " ";
	std::cout << sgBounds[3] << "; " << sgBounds[4] << " " << sgBounds[5] << "]" << std::endl;
}

/* FUNCTION NOT NEEDED ANYMORE?!
// Are two 3D Points in the same box of a grid defined by minBox[3]
// - pt1[3]: 3D Point Coordinates of point1
// - pt2[3]: 3D Point Coordinates of point2
// - minBox[3]: 3D minimal BoundingBox
// 		-> [50 250 -150] -> checks whether points are within boxes spanned by
//		in x VoxelSize:50:VoxelSize ...
//		in y VoxelSize:250:VoxelSize ...
//		in z VoxelSize:-150:VoxelSize ...
bool densitySG::arePtsInSameBox(double pt1[3], double pt2[3], double minBox[3])
{
	return (arePtsInSameRange(pt1[0],pt2[0],minBox[0],minBox[0]+VoxelSize) &&
			arePtsInSameRange(pt1[1],pt2[1],minBox[1],minBox[1]+VoxelSize) &&
			arePtsInSameRange(pt1[2],pt2[2],minBox[2],minBox[2]+VoxelSize));
}

// Returns boolean whether two points in one dimension are within same bounding box
bool densitySG::arePtsInSameRange(double x1, double x2, double minX, double maxX)
{
	int idx1 = ptInBox(x1,minX,maxX);
	int idx2 = ptInBox(x2,minX,maxX);

	// In Same Box
	if (idx1==idx2)
		return true;

	// Is Poiint on Edge
	bool edge1 = onEdge(x1,minX,maxX);
	bool edge2 = onEdge(x2,minX,maxX);

	// since floor operation have to subtract one
	// e.g. Voxel: [0 50] x1 = 0 -> idx1 = 0
	//		but 0 is also part of Voxel [-50 0] -> idx = -1
	if (!edge1 && !edge2) // Both not on an edge
		return false;
	else if (edge1 && edge2) // Both on an edge
	{
		if (fabs(idx1-idx2)<=1)
			return true;
		else
			return false;
	}
	else if (edge1) // One on edge
	{
		if ((idx1-1)==idx2)
			return true;
		else
			return false;
	}
	else if (edge2) // One on edge
	{
		if ((idx2-1)==idx1)
			return true;
		else
			return false;
	}
	else
	{
		std::cout << "WARNING! densitySG::isPtInSameRange(double x1, double x2, double minX, double maxX)" << std::endl;
		std::cout << "Should never reach here, something is wrong!" << std::endl;
		return false;
	}
}

// Checks whether given point lies on an edge of a boundingbox
bool densitySG::onEdge(double x, double minX, double maxX)
{
	return (std::fmod((x - minX),maxX - minX) == 0);
}

// Computes t values between two points and a box (standing_vector + t*direction_vector)
std::list< double > densitySG::getIntersectAllT(double xBox[2], double yBox[2], double zBox[2],
												double x1, double y1, double z1,
												double x2, double y2, double z2)
{
	std::list < double> t_values;
	double tmp;

	tmp = computeIntersectT(x1-xBox[0],x2-xBox[0]);
	if (tmp != 0.0)
		t_values.push_back(tmp);
	tmp = computeIntersectT(y1-yBox[0],y2-yBox[0]);
	if (tmp != 0.0)
		t_values.push_back(tmp);
	tmp = computeIntersectT(z1-zBox[0],z2-zBox[0]);
	if (tmp != 0.0)
		t_values.push_back(tmp);
	tmp = computeIntersectT(x1-xBox[1],x2-xBox[1]);
	if (tmp != 0.0)
		t_values.push_back(tmp);
	tmp = computeIntersectT(y1-yBox[1],y2-yBox[1]);
	if (tmp != 0.0)
		t_values.push_back(tmp);
	tmp = computeIntersectT(z1-zBox[1],z2-zBox[1]);
	if (tmp != 0.0)
		t_values.push_back(tmp);

	return t_values;
}

// Returns t-value (pt1 + t * (pt2 - pt1))
double densitySG::computeIntersectT(double Dst1, double Dst2)
{
	// Check whether Line is intersecting with given plane
	if ((Dst1*Dst2) >= 0.0)
		return 0.0;
	if (Dst1 ==  Dst2)
		return 0.0;

	// Compute intersecting Point via standing vetor plus weighted direction vector
	double tmp = (-Dst1/(Dst2-Dst1 //+ 1E-8));

	if (tmp>0.0 && tmp<1.0)
		return tmp;
	else
		return 0.0;
}

// Compute intersection point between two points pt1 and pt2 given the relative distance measurement t
void densitySG::getIntersectFromT(double t, double pt1[3], double pt2[3], double intersectPt[3])
{
	for (int i = 0; i<3; i++)
	{
		intersectPt[i] = pt1[i] + (pt2[i]-pt1[i]) * t;
	}
}*/

/* Getter & Setter Methods */
void densitySG::setLabel(int label)
{
	this->label = label; 
}

void densitySG::setVoxelSize(double VoxelSize)
{
	this->VoxelSize = VoxelSize;
	std::cout << "CAUTION! Take care of bounding box!" << std::endl;
	printProperties();
}

void densitySG::setOutputpath(std::string outputpath)
{
	this->outputpath = outputpath; 
}

void densitySG::setBoolAmirafile(bool boolAmirafile)
{
	this->boolAmirafile = boolAmirafile; 
}

void densitySG::setEpsilon(double epsilon)
{
	this->epsilon = epsilon; 
}

void densitySG::setBounds(double inputbounds[6])
{
	for(int i=0; i!=6; i++)
		bounds[i] = inputbounds[i]; 
}

void densitySG::getBounds(double inputbounds[6])
{
	for(int i=0; i!=6; i++)
		inputbounds[i] = bounds[i]; 
}
