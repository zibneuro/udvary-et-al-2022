/****************************************************************************/
/*                                                                          */
/* File:    helper.h                                                      	*/
/*                                                                          */
/* Purpose: Collection of useful functions						        	*/
/*                                                                          */
/* Author: 	Daniel Udvary													*/
/* 			In Silico Brain Sciences										*/
/*			Max Planck Institute for Neurobiology of Behavior – caesar		*/
/*			Ludwig-Erhard-Allee 2											*/
/*			53175 Bonn, Germany												*/
/*                                                                          */
/* e-mail: 	daniel.udvary@mpinb.mpg.de										*/
/*                                                                          */
/* Date: 	02.03.2020                                                  	*/
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "helper.h"

/* Load ImageDataVolume */
ImageDataPointerType helper::loadVolume(const char* inputFilename)
{
	std::string fname(inputFilename);
	if (fname.compare(fname.size()-3,3,".am") != 0)
		fname += ".am";
    
	if((access(fname.c_str(), F_OK) != -1))
	{
		Reader * Reader1 = new Reader(fname.c_str(),fname.c_str()); 
		return Reader1->readScalarField(); 
	}
	else
	{
		std::cout << "Error! Could not load " << fname.c_str() << "!" << std::endl;
	}
	return 0; 
} 

/* Loads ImageDataVolume Histograms (more than one scalar value at each position) */
ImageDataPointerType helper::loadVolumeN(const char* inputFilename)
{
	std::string fname(inputFilename);
	if (fname.compare(fname.size()-3,3,".am") != 0)
		fname += ".am";
  
	if((access(fname.c_str(), F_OK) != -1))
	{
		Reader * Reader1 = new Reader(fname.c_str(),fname.c_str()); 
		return Reader1->readNVectorField(); 
	}
	else
	{
		std::cout << "Error! Could not load " << fname << "!" << std::endl;
	}
	return 0; 
} 

/* Computes Mean of list */
double helper::computeMean(std::list< double > list)
{
	std::list<double>::iterator it = list.begin(); 
	
	double mean = 0.0;
	
	if (list.size()>0)
	{
		for (; it!=list.end(); it++)
		{
			mean += (*it);
		}
		mean /= list.size();
	}
	else
	{
		std::cout << "Error! List has no elements! Mean cannot be computed!" << std::endl; 
	}
	return mean; 
}

/* Computes SD of list */
double helper::computeSTD(std::list< double > list)
{
	std::list<double>::iterator it = list.begin(); 
	
	double mean = computeMean(list); 
	double std = 0.0; 
	
	if (list.size()>0)
	{
		for (; it!=list.end(); it++)
		{
			std += pow((*it)-mean,2.0);
		}
		std = sqrt(std/list.size());
	}
	else
	{
		std::cout << "Error! List has no elements! STD cannot be computed!" << std::endl; 
	}
	return std; 
}

/* Computes Mean and SD of list */
void helper::computeMeanSTD(std::list< double> list, double& mean, double& std)
{
	if (list.size()>0)
	{
		mean = computeMean(list); 
		std = 0.0; 
		std::list<double>::iterator it = list.begin(); 
		for (; it!=list.end(); it++)
		{
			std += pow((*it)-mean,2.0);
		}
		std = sqrt(std/list.size()); 
	}
	else
	{
		std::cout << "Error! List has no elements! Mean and STD cannot be computed!" << std::endl; 
	}
}

void helper::computeMeanSTDMinMax(std::list< double> list, double& mean, double& std, double& minVal, double& maxVal)
{
	if (list.size()>0)
	{
		mean = computeMean(list);
		std = 0.0;
		minVal = std::numeric_limits<double>::infinity();
		maxVal = -std::numeric_limits<double>::infinity();

		std::list<double>::iterator it = list.begin();
		for (; it!=list.end(); it++)
		{
			std += pow((*it)-mean,2.0);

			if (minVal>(*it))
				minVal = (*it);

			if (maxVal<(*it))
				maxVal = (*it);
		}
		std = sqrt(std/list.size());
	}
	else
	{
		std::cout << "Error! List has no elements! Mean and STD cannot be computed!" << std::endl;
	}
}


/* Computes Pearson's Correlation between list1 and list2 */
double helper::computeCorr(std::list< double > list1, std::list< double > list2)
{
	double corrcoef = 0.0; 
	
	if (list1.size()==list2.size() && list1.size()>0)
	{
		double m1, m2, std1, std2;  
		computeMeanSTD(list1,m1,std1); 
		computeMeanSTD(list2,m2,std2); 
		
		std::list<double>::iterator it1 = list1.begin(); 
		std::list<double>::iterator it2 = list2.begin(); 
			
		for (; it1!=list1.end();it1++,it2++)
		{
			corrcoef += ((*it1)-m1)*((*it2)-m2); 
		}
		
		corrcoef /= list1.size(); 
		if (std1<1E-9 || std2<1E-9)
		{
			//std::cout << " Warning! STD is 0.0! Corrcoeff is set to 0.0!" << std::endl;
			corrcoef = 0.0;  
		}
		else
		{
			corrcoef /= (std1*std2);  
		}
	}
	else
	{
		std::cout << " Error! Correlation of two lists with different number of elements cannot be computed!" << std::endl; 
		std::cout << " 		list1.size() = " << list1.size() << " list2.size() = " << list2.size() << std::endl;
	}
	
	return corrcoef; 
}

/* computes Cumulative Sum of inputarray of length */
void helper::computeCumSum(double cumsum[], double inputarray[], int length)
{ 
	for(int i = 0; i<length; i++)
	{
		cumsum[i] = inputarray[i];
		if (i>0)
			cumsum[i] += cumsum[i-1]; 
	}
}

/* Stores ImageDataVolume as .txt file */
void helper::convertImageData2txt(ImageDataPointerType volume, const char* outputFilename)
{
	std::ofstream outStream(outputFilename);
	
	if(!outStream.fail())
	{	  
		int extent[6];
		double CoordOrigin[3]; 
		volume->GetOrigin(CoordOrigin);
		volume->GetExtent(extent);
		
		outStream << extent[0] << "," << extent[1] << "," << extent[2] << "," << extent[3] << "," << extent[4] << "," << extent[5] << std::endl;
		outStream << CoordOrigin[0] << "," << CoordOrigin[1] << "," << CoordOrigin[2] << std::endl;
	  
		for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
		{  
			for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
			{
				for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
				{ 
					double *px = static_cast< double * >(volume->GetScalarPointer(x,y,z));
					outStream << (*px);
					
					if (z==volume->GetExtent()[5])
					{
						outStream << "; ";
					}
					else
					{
						outStream << " "; 
					}
				}
			}
			outStream << "" << std::endl; 
		}
		outStream.close();
	}
	else
	{
		std::cout << "Error! Writing ProfileTxt failed! Path:" <<  outputFilename << std::endl; 
	}
}

/* Corrects ImageDataVolume by setting Origin and extent correctly so that values run from - to + */
ImageDataPointerType helper::correctVolume(ImageDataPointerType volume, double VoxelSize)
{
	double OriginCoord[3];
	int extent[6];
	volume->GetOrigin(OriginCoord);
	volume->GetExtent(extent);
	
	extent[0] += OriginCoord[0]/VoxelSize;
	extent[1] += OriginCoord[0]/VoxelSize;
	extent[2] += OriginCoord[1]/VoxelSize;
	extent[3] += OriginCoord[1]/VoxelSize;
	extent[4] += OriginCoord[2]/VoxelSize;
	extent[5] += OriginCoord[2]/VoxelSize;

	volume->SetOrigin(0,0,0);
	volume->SetExtent(extent);
	
	return volume; 
}

/* Create Image Data Volume (VTK), Input: BoundingBox Coordinates, Output: ImageDataVolume */
ImageDataPointerType helper::createImageVolume(double bounds[6], double VoxelSize)
{
	double spacing[3];
	int dims[6];
	spacing[0] = spacing[1] = spacing[2] = VoxelSize;
	
	calculateExtent(bounds, dims, VoxelSize);
	
	// Create ImageDataVolume, assign initial doubles to 0.0
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	volume->SetExtent(dims);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();
	
	// Initialize each Scalar with 0.0
	for(int z = dims[4]; z <= dims[5]; ++z)
	{
	   for(int y = dims[2]; y <= dims[3]; ++y)
	   {
			for(int x = dims[0]; x <= dims[1]; ++x)
			{
				double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));
				*px = 0.0;
			}
	   }
	}
	volume->Update();
	return volume;
}

/* Create Image Data Volume (VTK), Input: BoundingBox Coordinates, Output: ImageDataVolume
 * NeuroNet: Voxels run from [-50 to 0; 0 to 50], origin is [-50], extends from [0 1] */
ImageDataPointerType helper::createImageVolumeNeuroNet(double bounds[6], double VoxelSize)
{
	double spacing[3];
	spacing[0] = spacing[1] = spacing[2] = VoxelSize;

	int extent[6];
	extent[0] = 0;
	extent[1] = (bounds[1] - bounds[0])/spacing[0];
	extent[2] = 0;
	extent[3] = (bounds[3] - bounds[2])/spacing[1];
	extent[4] = 0;
	extent[5] = (bounds[5] - bounds[4])/spacing[2];

	// Create ImageDataVolume, assign initial doubles to 0.0
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	volume->SetOrigin(bounds[0]+VoxelSize/2,bounds[2]+VoxelSize/2,bounds[4]+VoxelSize/2);
	volume->SetExtent(extent);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();

	// Initialize each Scalar with 0.0
	for(int z = extent[4]; z <= extent[5]; ++z)
	{
	   for(int y = extent[2]; y <= extent[3]; ++y)
	   {
			for(int x = extent[0]; x <= extent[1]; ++x)
			{
				double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));
				*px = 0.0;
			}
	   }
	}
	volume->Update();
	return volume;
}

ImageDataPointerType helper::dimPlusOneVolumne(ImageDataPointerType inputVolume)
{
	// Get Spacing
	double Spacing[3];
	inputVolume->GetSpacing(Spacing[0],Spacing[1],Spacing[2]);

	// Correct Spacing
	double maxSpacing = 0;
	if (Spacing[0]>0 && std::isnan(Spacing[1]) && std::isnan(Spacing[2]))
		maxSpacing = Spacing[0];
	else if (std::isnan(Spacing[0]) && Spacing[1]>0 && std::isnan(Spacing[2]))
		maxSpacing = Spacing[1];
	else if (std::isnan(Spacing[0]) && std::isnan(Spacing[1]) && Spacing[2]>0)
		maxSpacing = Spacing[2];
	else
	{
		std::cout << "ERROR in helper::dimPlusOneVolumne!" << std::endl;
		std::cout << "  Spacing = [" << Spacing[0] << "," << Spacing[1] << "," << Spacing[2] << "]" << std::endl;
		return inputVolume;
	}

	if (maxSpacing==0)
	{
		std::cout << "ERROR in helper::dimPlusOneVolumne! MaxSpacing is zero!" << std::endl;
		std::cout << "  Spacing = [" << Spacing[0] << "," << Spacing[1] << "," << Spacing[2] << "]" << std::endl;
		return inputVolume;
	}

	Spacing[0] = maxSpacing;
	Spacing[1] = maxSpacing;
	Spacing[2] = maxSpacing;

	// Correct Extent by adding one
	int extent[6];
	inputVolume->GetExtent(extent[0],extent[1],extent[2],extent[3],extent[4],extent[5]);
	extent[1] = extent[1] + 1;
	extent[3] = extent[3] + 1;
	extent[5] = extent[5] + 1;

	// Extract Origin
	double origin[3];
	inputVolume->GetOrigin(origin[0],origin[1],origin[2]);

	// Create ImageDataVolume, assign initial doubles to 0.0
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetExtent(extent);
	volume->SetSpacing(Spacing[0],Spacing[1],Spacing[2]);
	volume->SetOrigin(origin[0],origin[1],origin[2]);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();

	// Copy values into new density
	for(int z = extent[4]; z <= extent[5]; ++z)
	{
	   for(int y = extent[2]; y <= extent[3]; ++y)
	   {
			for(int x = extent[0]; x <= extent[1]; ++x)
			{
				double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));

				// Copy values over from inputVolume
				if (z<extent[5] && y<extent[3] && x<extent[1])
				{
					double * px2 = static_cast< double * >(inputVolume->GetScalarPointer(x, y, z));
					*px = (*px2);
				}
				else // set new values in added dimensions to zero
				{
					*px = 0.0;
				}
			}
	   }
	}
	volume->Update();
	return volume;
}

/* Create Image Data Volume (VTK), Input: Dimensions, Output: ImageDataVolume */
ImageDataPointerType helper::createImageVolume(int dims[6], double VoxelSize)
{
	double spacing[3];
	spacing[0] = spacing[1] = spacing[2] = VoxelSize;

	// Create ImageDataVolume, assign initial doubles to 0.0
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	volume->SetExtent(dims);
	volume->SetNumberOfScalarComponents(1);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();
	
	// Initialize each Scalar with 0.0
	for(int z = dims[4]; z <= dims[5]; ++z)
	       for(int y = dims[2]; y <= dims[3]; ++y)
			for(int x = dims[0]; x <= dims[1]; ++x)
			{
				double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));
				*px = 0.0;
			}
	volume->Update();
	return volume;  
}

/* Create Image Data Volume (VTK), Input: Dimensions, Output: ImageDataVolume */
ImageDataPointerType helper::createImageVolume(int dims[6], double VoxelSize, int NumScalarComponents)
{
	double spacing[3];
	spacing[0] = spacing[1] = spacing[2] = VoxelSize;

	// Create ImageDataVolume, assign initial doubles to 0.0
	ImageDataPointerType volume = ImageDataPointerType::New();
	volume->SetSpacing(spacing[0], spacing[1], spacing[2]);
	volume->SetExtent(dims);
	volume->SetNumberOfScalarComponents(NumScalarComponents);
	volume->SetScalarTypeToDouble();
	volume->AllocateScalars();
	
	// Initialize each Scalar with 0.0
	for(int z = dims[4]; z <= dims[5]; ++z)
	       for(int y = dims[2]; y <= dims[3]; ++y)
			for(int x = dims[0]; x <= dims[1]; ++x)
			{
				double * px = static_cast< double * >(volume->GetScalarPointer(x, y, z));
				
				for(int i = 0; i<NumScalarComponents; i++)
					px[i] = 0;
			}
	volume->Update();
	return volume;  
}

/* Draw values from Histogram ImageDataVolume (if binsz set to 0, assume integer values) 
 * returns drawn values as ImageDataVolume
 */
ImageDataPointerType helper::drawFromHist(ImageDataPointerType volume, double VoxelSize, int binSz)
{
	int dim[6]; 
	double OriginCoord[3];
	volume->GetExtent(dim); 
	volume->GetOrigin(OriginCoord);
	// Create new ImageDataVolume
	ImageDataPointerType volumeReturn = createImageVolume(dim,VoxelSize); 
	volumeReturn->SetOrigin(OriginCoord[0],OriginCoord[1],OriginCoord[2]);
	int NumScalarComponents = volume->GetNumberOfScalarComponents(); 
	
	srand(time(NULL));
	int maxValue; 
	int randVal;
	
// 	volume->Print(std::cout); 
// 	volumeReturn->Print(std::cout); 
	
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{  
		for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
		{
			for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
			{    
				double *pxHist = static_cast< double * >(volume->GetScalarPointer(x,y,z)); 	// Histogram
				double *px = static_cast< double * >(volumeReturn->GetScalarPointer(x,y,z)); 
				
				// Compute cumulative sum for histogram
				double cumsum[NumScalarComponents]; 
				computeCumSum(cumsum, pxHist, NumScalarComponents); 
				
				// Draw random number between 1 and maxValue
				maxValue = cumsum[NumScalarComponents-1]; 
				int i = 0; 
				
				if (maxValue>0) // otherwise length is set to zero
				{
					randVal = rand() % maxValue + 1;
					
					while(randVal>cumsum[i] && i<NumScalarComponents)	// find valid index
						i++; 
				}
				
				if (binSz==0)
					*px = double(i); 
				else
				{
					// Get Max value in this voxel (0 25 75 ... )
					*px = double(i*binSz)-binSz/2;
					if ((*px)<0)
						(*px) = 0; 
				}
			}
		}
	}
	
	volumeReturn->Update();
	return volumeReturn; 
}

/* Computes Extend of given Bounds and VoxelSize (Number of Voxels in each Dimensions +/-) */
void helper::calculateExtent(double bounds[6], int extent[6], double VoxelSize)
{
	double spacing[3];
	spacing[0] = spacing[1] = spacing[2] = VoxelSize;
	
	//make sure that maxCoordinates are inside an integer number of cells defined by spacing
	extent[0] = (bounds[0] - spacing[0])/spacing[0]/* - 0.5*/;
	extent[1] = (bounds[1] + spacing[0])/spacing[0]/* + 0.5*/;
	extent[2] = (bounds[2] - spacing[1])/spacing[1]/* - 0.5*/;
	extent[3] = (bounds[3] + spacing[1])/spacing[1]/* + 0.5*/;
	extent[4] = (bounds[4] - spacing[2])/spacing[2]/* - 0.5*/;
	extent[5] = (bounds[5] + spacing[2])/spacing[2]/* + 0.5*/;
}

/* Compute BoundingBox along all dimensions (calls getBoundingBox); 
 * min in [0], max in [1] 
 */
void helper::getBoundingBox(double xBox[2], double yBox[2], double zBox[2],AmiraSpatialGraph * SpatialGraph)
{
	getBoundingBox(X_COORD,xBox,SpatialGraph);
	getBoundingBox(Y_COORD,yBox,SpatialGraph);
	getBoundingBox(Z_COORD,zBox,SpatialGraph);
}

void helper::getBoundingBox(double xBox[2], double yBox[2], double zBox[2],AmiraSpatialGraph * SpatialGraph, int neuriteID)
{
	getBoundingBox(X_COORD,xBox,SpatialGraph,neuriteID);
	getBoundingBox(Y_COORD,yBox,SpatialGraph,neuriteID);
	getBoundingBox(Z_COORD,zBox,SpatialGraph,neuriteID);
}


/* Get BoundingBox along one dimension.int Coordinate between 0 and 2 (X_COORD, Y_COORD, Z_COORD), double Coord returns min and max point */
void helper::getBoundingBox(const int Coordinate, double Coord[2],AmiraSpatialGraph * SpatialGraph)
{
	double maxC = -1E9;
	double minC = 1E9;
	
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = SpatialGraph->edgesBegin(); edgeIt != SpatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label >= Neuron && (*edgeIt)->label <= Soma)
		{
			std::list< double * >::iterator edgeListIt;
			for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{
				double tmpC = (*edgeListIt)[Coordinate];
				if(tmpC > maxC)
					maxC = tmpC;
				if(tmpC < minC)
					minC = tmpC;
			}
		}
	}
	
	Coord[0] = minC;
	Coord[1] = maxC;   
}

/* Get BoundingBox along one dimension.int Coordinate between 0 and 2 (X_COORD, Y_COORD, Z_COORD), double Coord returns min and max point */
void helper::getBoundingBox(const int Coordinate, double Coord[2],AmiraSpatialGraph * SpatialGraph, int neuriteID)
{
	double maxC = -1E9;
	double minC = 1E9;

	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = SpatialGraph->edgesBegin(); edgeIt != SpatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label == neuriteID)
		{
			std::list< double * >::iterator edgeListIt;
			for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{
				double tmpC = (*edgeListIt)[Coordinate];
				if(tmpC > maxC)
					maxC = tmpC;
				if(tmpC < minC)
					minC = tmpC;
			}
		}
	}

	Coord[0] = minC;
	Coord[1] = maxC;
}

/* Stores ImageDataVolume under given name */
void helper::storeImageVolume(const char * outputFilename, ImageDataPointerType volume)
{
	std::string fname(outputFilename);
	if (fname.compare(fname.size()-3,3,".am") == 0)
		fname = fname.substr(0,fname.size()-3);
	
	Reader * Writer = new Reader(fname.c_str(), fname.c_str());
	
	int numScalar = volume->GetNumberOfScalarComponents(); 
	
	if (numScalar==1)
		Writer->writeScalarField(volume);
	else
		Writer->writeNVectorField(volume);
	
	delete Writer; 
}

/* Computes Values of Supragranular, Granular, and Infragranular layer along cortical coord (Y for non-registrated, Z for registrated)*/
void helper::getSGI(ImageDataPointerType volume, double SGI[3], int coord, double VoxelSize)
{
	// Supra:	706 194
	// Granular: 	194 -194
	// Infra:	-194 -1251
	
	SGI[0] = 0;
	SGI[1] = 0;
	SGI[2] = 0;
	
	int coordOthers[2]; 
	if (coord==0)
	{
		coordOthers[0] = 1;
		coordOthers[1] = 2; 
	}
	else if (coord==1)
	{
		coordOthers[0] = 0;
		coordOthers[1] = 2; 
	}
	else if (coord==2)
	{
		coordOthers[0] = 0;
		coordOthers[1] = 1; 
	}
	else
	{
		std::cout << coord << " is not a valid Coordinate [0 2]! 1 (Y) for non-registrated, 2 (Z) for registrated Cells!" << std::endl; 
		return;
	}
	
	std::cout << "Warning! Not perfect, depends on the edges of the voxels within the ImageVolume!" << std::endl; 
	
	for(int c = volume->GetExtent()[coord*2]; c <= volume->GetExtent()[coord*2+1]; ++c)
	{
		double tmp = 0; 
		
		for(int c0 = volume->GetExtent()[coordOthers[0]*2]; c0 <= volume->GetExtent()[coordOthers[0]*2+1]; ++c0)
		{  
			for(int c1 = volume->GetExtent()[coordOthers[1]*2]; c1 <= volume->GetExtent()[coordOthers[1]*2+1]; ++c1)
			{
				double *px;
				if (coord==0)
				{
					px = static_cast< double * >(volume->GetScalarPointer(c,c0,c1));
				}
				else if (coord==1)
				{
					px = static_cast< double * >(volume->GetScalarPointer(c0,c,c1));
				}
				else if (coord==2)
				{
					px = static_cast< double * >(volume->GetScalarPointer(c0,c1,c));
				}
				tmp += (*px); 
			}
		}
		
		if ((c*VoxelSize+VoxelSize/2)<-194) // Infra
		{
			SGI[2] += tmp; 
		}   
		else if ((c*VoxelSize+VoxelSize/2)>=194) // Supra
		{
			SGI[0] += tmp; 
		}  
		else // Granular
		{
			SGI[1] += tmp; 
		}   
	}
	
	std::cout << "SGI: ( " << SGI[0] << " , " << SGI[1] << " , " << SGI[2] << " ) um Length" << std::endl; 
}

/* Checks whether given point coordinates x y z are within ImageDataVolume */
bool helper::isValidPt(int x, int y, int z, ImageDataPointerType volume)
{
	int extent[6];
	volume->GetExtent(extent);	
	return (extent[0]<=x && extent[1]>=x && extent[2]<=y && extent[3]>=y && extent[4]<=z && extent[5]>=z);
}

/* Gets ImageDataPointerType as Input and creates .txt containing the normalized profile along x and z axis */
void helper::getProfileAsTxt(ImageDataPointerType volume, const char * path, double VoxelSize)
{   
	int dim[3]; 
	volume->GetDimensions(dim);

	// Get 2D-Array and Initialize with Zeros
	int x_off = volume->GetExtent()[0]; 
	double x_profile[dim[0]][2];
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{
		x_profile[x-x_off][1] = 0.0; 
		x_profile[x-x_off][0] = x*VoxelSize; 
	}
	
	int y_off = volume->GetExtent()[2]; 
	double y_profile[dim[1]][2];
	for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
	{
		y_profile[y-y_off][1] = 0.0;
		y_profile[y-y_off][0] = y*VoxelSize;
	}
	
	int z_off = volume->GetExtent()[4]; 
	double z_profile[dim[2]][2];
	for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
	{
		z_profile[z-z_off][1] = 0.0;
		z_profile[z-z_off][0] = z*VoxelSize;
	}
	
	// Sum up values in volume
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{  
		for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
		{
			for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
			{ 
				double *px = static_cast< double * >(volume->GetScalarPointer(x,y,z));
				x_profile[x-x_off][1] += (*px); 
				z_profile[z-z_off][1] += (*px); 
				y_profile[y-y_off][1] += (*px); 
			}
		}
	}
	
	// Normalize 
	double x_norm =  1E9/(VoxelSize*VoxelSize*VoxelSize*dim[1]*dim[2]); 	// (50*50*50*y_dim*z_dim)
	double z_norm = 1E9/(VoxelSize*VoxelSize*VoxelSize*dim[0]*dim[1]); 	// (50*50*50*x_dim*y_dim)
	double y_norm = 1E9/(VoxelSize*VoxelSize*VoxelSize*dim[0]*dim[2]); 	// (50*50*50*x_dim*z_dim)
	
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{
		x_profile[x-x_off][1] *= x_norm; 
	}	
	for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
	{
		y_profile[y-y_off][1] *= y_norm;
	}
	for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
	{
		z_profile[z-z_off][1] *= z_norm;
	}
	
	std::ofstream outStream(path);
	
	if(!outStream.fail())
	{	  
		outStream << "x_profile Dim: [" << dim[0] << "," << dim[1] << "," << dim[2] << "]" << std::endl;
	  
		for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
		{
			outStream << x_profile[x-x_off][0] << " " << x_profile[x-x_off][1] << std::endl; 
		}
		
		outStream << "\n" << std::endl; 
		outStream << "y_profile Dim: [" << dim[0] << "," << dim[1] << "," << dim[2] << "]" << std::endl;
		
		for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
		{
			outStream << y_profile[y-y_off][0] << " " << y_profile[y-y_off][1] << std::endl;
 		}

 		outStream << "\n" << std::endl; 
		outStream << "z_profile Dim: [" << dim[0] << "," << dim[1] << "," << dim[2] << "]" << std::endl;
		
		for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
		{
			outStream << z_profile[z-z_off][0] << " " << z_profile[z-z_off][1] << std::endl;
 		}
		outStream.close();
	}
	else
	{
		std::cout << "Error! Writing ProfileTxt failed! Path:" <<  path << std::endl; 
	}
}

/* Computes directory one level up the given directory or file */
std::string helper::getRootFromPath(const char* inputFilenameList)
{
	std::string ofName(inputFilenameList); 
	unsigned found = ofName.find_last_of("/"); 
	std::string rootDir = ofName.substr(0,found);
	rootDir += "/";
	return rootDir; 
}

std::string helper::getFilenameFromPath(const char* inputPathFilename)
{
	// Save truncated Spatial Graph
	std::string tmpStr(inputPathFilename);
	unsigned found = tmpStr.find_last_of("/");
	std::string Filename = tmpStr.substr(found+1);
	return Filename;
}

/* Returns neuriteName of neuriteID */
std::string helper::getNeuriteName(int neuriteID)
{
	switch (neuriteID)
	{  
		case Axon:
			return "axon"; 
			break;
		case Dendrite:
			return "dend";
			break;
		case Soma:
			return "soma";
			break;
		case ApicalDendrite:
			return "apical";
			break;
		case BasalDendrite:
			return "basal";
			break;
		default:
			std::cout<< "neuriteID neither correspondes to axon, soma, dendrite, basal, or apical!" << std::endl;
			break; 
	}
	return ""; 
}

/* Displays array of length len in console / shell */
void helper::printArray(double array[], int len)
{
	for (int i = 0; i<len; i++)
		std::cout << array[i] << " "; 
	std::cout << "" << std::endl; 
}

/* Compute BoundingBox of ImageDataVolume/Density for given VoxelSize (50) */
void helper::getBoundingBox(double xBox[2], double yBox[2], double zBox[2], ImageDataPointerType volume, double VoxelSize)
{
	double box[6];
	volume->GetBounds(box); 
	
	xBox[0] = box[0]-VoxelSize/2; 
	xBox[1] = box[1]+VoxelSize/2; 
	yBox[0] = box[2]-VoxelSize/2; 
	yBox[1] = box[3]+VoxelSize/2; 
	zBox[0] = box[4]-VoxelSize/2; 
	zBox[1] = box[5]+VoxelSize/2; 
}

/* Compute the Z score between a list of artificially created nerons (named in inputFilenameList (.txt) and compares voxelwise with the symmetric image data volumes
 * Compares Length, Terminal and Branch Nodes 
 */
void helper::computeZScore(const char* inputFilenameList, const char* outputPrefix, int neuriteID)
{
	double VoxelSize = 50; 
	std::string outputPrefixString(outputPrefix); 
	std::string rootDir = helper::getRootFromPath(inputFilenameList); 
	std::string rootDirSym = rootDir.substr(0,rootDir.size()-1);
	rootDirSym = helper::getRootFromPath(rootDirSym.c_str()); 
	rootDirSym = rootDirSym.substr(0,rootDirSym.size()-1);
	rootDirSym = helper::getRootFromPath(rootDirSym.c_str()); 
	
	std::ifstream inputStream(inputFilenameList); 
	std::string currentLine;
	int dim[6]; 
	double OriginCoord[3]; 
	std::string neurite = helper::getNeuriteName(neuriteID); 
	double epsilon = 1e-8; 
	
	std::list <ImageDataPointerType> CurrentLengthList; 
	std::list <ImageDataPointerType> CurrentBranchList;
	std::list <ImageDataPointerType> CurrentTerminalList; 
	
	// Read in all the files and store them in four lists (Length Information + Border Information, for dendrite and axon)
	while(!inputStream.eof())
	{
		getline(inputStream,currentLine); 
		if (currentLine.empty())
			break;
		std::string ofName = rootDir + currentLine; 
		ofName.erase(ofName.end()-3,ofName.end()); 
		std::string ofNameLength = ofName + "_length.am";
		std::string ofNameBranch = ofName + "_branch.am"; 
		std::string ofNameTerm = ofName + "_terminal.am"; 
		
		ImageDataPointerType vol1 = helper::loadVolume(ofNameLength.c_str()); 
		CurrentLengthList.push_back(vol1);

		ImageDataPointerType vol2 = helper::loadVolume(ofNameBranch.c_str());
		CurrentBranchList.push_back(vol2);

		ImageDataPointerType vol3 = helper::loadVolume(ofNameTerm.c_str());
		CurrentTerminalList.push_back(vol3);
	}
	inputStream.close(); 
	
	// Load symmetric volumes
	std::string volName= rootDirSym + outputPrefixString + "_" + neurite + "_length_sym.am";
	ImageDataPointerType volumeLengthSym = helper::loadVolume(volName.c_str()); 
	volName= rootDirSym + outputPrefixString + "_" + neurite + "_length_sym_STD.am";
	ImageDataPointerType volumeLengthSymSTD = helper::loadVolume(volName.c_str()); 
	volName= rootDirSym + outputPrefixString + "_" + neurite + "_branch_sym.am";
	ImageDataPointerType volumeBranchSym = helper::loadVolume(volName.c_str()); 
	volName= rootDirSym + outputPrefixString + "_" + neurite + "_branch_sym_STD.am";
	ImageDataPointerType volumeBranchSymSTD = helper::loadVolume(volName.c_str()); 
	volName= rootDirSym + outputPrefixString + "_" + neurite + "_terminal_sym.am";
	ImageDataPointerType volumeTermSym = helper::loadVolume(volName.c_str()); 
	volName= rootDirSym + outputPrefixString + "_" + neurite + "_terminal_sym_STD.am";
	ImageDataPointerType volumeTermSymSTD = helper::loadVolume(volName.c_str()); 
	
	// Get Extent and Coordinates of the Origin 
	volumeLengthSym->GetExtent(dim); 
	volumeLengthSym->GetOrigin(OriginCoord); 

	// New ImageVolumes for the resulting Z scores
	// Length
	ImageDataPointerType volumeLength = helper::createImageVolume(dim,VoxelSize); 
	volumeLength->SetOrigin(OriginCoord[0],OriginCoord[1],OriginCoord[2]);
	
	// Branching Points
	ImageDataPointerType volumeBranch = helper::createImageVolume(dim,VoxelSize); 
	volumeBranch->SetOrigin(OriginCoord[0],OriginCoord[1],OriginCoord[2]);
	
	// Terminal Points
	ImageDataPointerType volumeTerm = helper::createImageVolume(dim,VoxelSize); 
	volumeTerm->SetOrigin(OriginCoord[0],OriginCoord[1],OriginCoord[2]);
	
	for(int x = volumeLength->GetExtent()[0]; x <= volumeLength->GetExtent()[1]; ++x)
	{  
		for(int y = volumeLength->GetExtent()[2]; y <= volumeLength->GetExtent()[3]; ++y)
		{
			for(int z = volumeLength->GetExtent()[4]; z <= volumeLength->GetExtent()[5]; ++z)
			{  
				// Z Store
				double *pxLengthZ = static_cast< double * >(volumeLength->GetScalarPointer(x,y,z)); 
				double *pxBranchZ = static_cast< double * >(volumeBranch->GetScalarPointer(x,y,z)); 
				double *pxTermZ = static_cast< double * >(volumeTerm->GetScalarPointer(x,y,z)); 
				
				// Symmetric Stats
				double *pxLengthMean = static_cast< double * >(volumeLengthSym->GetScalarPointer(x,y,z)); 
				double *pxBranchMean = static_cast< double * >(volumeBranchSym->GetScalarPointer(x,y,z)); 
				double *pxTermMean = static_cast< double * >(volumeTermSym->GetScalarPointer(x,y,z)); 
				double *pxLengthSTD = static_cast< double * >(volumeLengthSymSTD->GetScalarPointer(x,y,z)); 
				double *pxBranchSTD= static_cast< double * >(volumeBranchSymSTD->GetScalarPointer(x,y,z)); 
				double *pxTermSTD = static_cast< double * >(volumeTermSymSTD->GetScalarPointer(x,y,z)); 
				
				std::list<double> lengthList;
				std::list<double> branchList;
				std::list<double> termList;
				
				// Go through the list of ImageDataVolumes
				std::list<ImageDataPointerType>::iterator itBP=CurrentBranchList.begin();
				std::list<ImageDataPointerType>::iterator itTP=CurrentTerminalList.begin();
				
				for (std::list<ImageDataPointerType>::iterator itLen=CurrentLengthList.begin();itLen!=CurrentLengthList.end();++itLen,++itBP,++itTP)
				{
					double tmp = 0; 
					
					// Get current length density, put in list
					double *pxLength = static_cast< double * >((*itLen)->GetScalarPointer(x,y,z));
					if ( fabs((*pxLength)-(*pxLengthMean))<epsilon && fabs(*pxLengthSTD)<epsilon )
						tmp = 0.0; 
					else if (fabs(*pxLengthSTD)<epsilon)
						tmp =  (std::numeric_limits<double>::infinity()); //*((*pxLength)-(*pxLengthMean)); 
					else
						tmp = ((*pxLength)-(*pxLengthMean))/(*pxLengthSTD);
					lengthList.push_back(tmp); 
					
					// Get current number of branching points, put in list
					double *pxBranch = static_cast< double * >((*itBP)->GetScalarPointer(x,y,z));
					if ( fabs((*pxBranch)-(*pxBranchMean))<epsilon && fabs(*pxBranchSTD)<epsilon )
						tmp = 0.0; 
					else if (fabs(*pxBranchSTD)<epsilon)
						tmp =  (std::numeric_limits<double>::infinity()); //*((*pxBranch)-(*pxBranchMean)); 
					else
						tmp = ((*pxBranch)-(*pxBranchMean))/(*pxBranchSTD);
					branchList.push_back(tmp);
					
					// Get current number of terminal points, put in list
					double *pxTerm = static_cast< double * >((*itTP)->GetScalarPointer(x,y,z));
					if ( fabs((*pxTerm)-(*pxTermMean))<epsilon && fabs(*pxTermSTD)<epsilon )
						tmp = 0.0; 
					else if (fabs(*pxTermSTD)<epsilon)
						tmp =  (std::numeric_limits<double>::infinity()); //*((*pxTerm)-(*pxTermMean)); 
					else
						tmp = ((*pxTerm)-(*pxTermMean))/(*pxTermSTD);
					termList.push_back(tmp);
				}
				
				// Mean over all Z scores for this voxel
				*pxLengthZ = helper::computeMean(lengthList); 
				lengthList.clear();
				*pxBranchZ = helper::computeMean(branchList); 
				branchList.clear();
				*pxTermZ = helper::computeMean(termList); 
				termList.clear();
			}
		}
	}
	
	// Store ImageVolumes with z scores
	volumeLength->Update();
	volumeBranch->Update();
	volumeTerm->Update();
	
	// Write ImageDataVolume to Disk
	// Length Information
	std::string outputFilename = rootDir + outputPrefixString + "_" + neurite + "_length_" + "Z";
	helper::storeImageVolume(outputFilename.c_str(), volumeLength);
	// Branching Point Information
	outputFilename = rootDir + outputPrefixString + "_" + neurite + "_branch_" + "Z";
	helper::storeImageVolume(outputFilename.c_str(), volumeBranch); 
	// Terminal Point Information
	outputFilename = rootDir + outputPrefixString + "_" + neurite + "_terminal_" + "Z";
	helper::storeImageVolume(outputFilename.c_str(), volumeTerm); 
}

/* Computes length of spatialGraph neuriteID */
double helper::lengthCell(AmiraSpatialGraph * SpatialGraph, int neuriteID)
{	
	// length of edges within clippedBox
	double len = 0; 
// 	double p = 0; 
	
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = SpatialGraph->edgesBegin(); edgeIt != SpatialGraph->edgesEnd(); ++edgeIt)
	{
		// Check only certain label
		if((*edgeIt)->label == neuriteID)
		{   
			std::list< double * >::iterator edgeListIt;			
			edgeListIt = (*edgeIt)->edgePointCoordinates.begin();
			double * previousPt = *edgeListIt;
			++edgeListIt;

			// Compute Distance between first Vertex and first Point (they should be identical!)
			int tmpVertexID = (*edgeIt)->edgeConnectivity[0];
			Vertex * tmpVertex = SpatialGraph->verticesPointer()->at(tmpVertexID);
			len += sqrt(pow((tmpVertex->coordinates[X_COORD]-previousPt[X_COORD]),2.0)+pow((tmpVertex->coordinates[Y_COORD]-previousPt[Y_COORD]),2.0)+pow((tmpVertex->coordinates[Z_COORD]-previousPt[Z_COORD]),2.0));
			
			// Go through all points of each Edge starting with the second edge
			for(; edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{  
				double * currentPt = *edgeListIt;
				len += sqrt(pow((currentPt[X_COORD]-previousPt[X_COORD]),2.0)+pow((currentPt[Y_COORD]-previousPt[Y_COORD]),2.0)+pow((currentPt[Z_COORD]-previousPt[Z_COORD]),2.0));  
				previousPt = currentPt;
			}
// 			p = p + (*edgeIt)->numEdgePoints;
// 			std::cout << (*edgeIt)->numEdgePoints << std::endl; 
		}
	}
	
// 	std::cout << " >> Edges " << e << " EdgePts " << p << std::endl; 
	return len;  
}

/* Function checks whether Vertex and Edge Points match! (Essential for calculation of Length) */
void helper::checkCell(AmiraSpatialGraph * SpatialGraph, int neuriteID)
{
	double errFront = 0;
	double errBack = 0;

	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = SpatialGraph->edgesBegin(); edgeIt != SpatialGraph->edgesEnd(); ++edgeIt)
	{
		// Check only certain label
		if((*edgeIt)->label == neuriteID)
		{
			std::list< double * >::iterator edgeListIt;

			// Compute Distance between first Vertex and first Point (they should be identical!)
			edgeListIt = (*edgeIt)->edgePointCoordinates.begin();
			double * firstPt = *edgeListIt;
			int tmpVertexID = (*edgeIt)->edgeConnectivity[0];
			Vertex * tmpVertex = SpatialGraph->verticesPointer()->at(tmpVertexID);
			errFront += sqrt(pow((tmpVertex->coordinates[X_COORD]-firstPt[X_COORD]),2.0)+pow((tmpVertex->coordinates[Y_COORD]-firstPt[Y_COORD]),2.0)+pow((tmpVertex->coordinates[Z_COORD]-firstPt[Z_COORD]),2.0));

			// Compute Distance between first Vertex and last Point (they should be identical!)
			edgeListIt = (*edgeIt)->edgePointCoordinates.end();
			edgeListIt--;
			double * lastPt = *edgeListIt;
			tmpVertexID = (*edgeIt)->edgeConnectivity[1];
			tmpVertex = SpatialGraph->verticesPointer()->at(tmpVertexID);
			errBack += sqrt(pow((tmpVertex->coordinates[X_COORD]-lastPt[X_COORD]),2.0)+pow((tmpVertex->coordinates[Y_COORD]-lastPt[Y_COORD]),2.0)+pow((tmpVertex->coordinates[Z_COORD]-lastPt[Z_COORD]),2.0));
		}
	}

	if (errFront > 1)
	{
		std::cout << "-- ISSUE FOUND -- " << std::endl;
		std::cout << "WARNING! The first Edge Point and the first Vertex Point do not match! Introduced length error is " << errFront << "!" << std::endl;
		std::cout << " helper::lengthCell _does_ takes care of that issue when calculating the length!" << std::endl;
	}

	if (errBack > 1)
	{
		std::cout << "-- ISSUE FOUND -- " << std::endl;
		std::cout << "WARNING! The last Edge Point and the last Vertex Point do not match! Introduced length error is " << errBack << "!" << std::endl;
		std::cout << " WARNING! helper::lengthCell _DOES_NOT_ take care of that issue when calculating the length!" << std::endl;
	}

	//std::cout << "Error (front) : " << errFront << " Error (back) : " << errBack << std::endl;
}


/* Returns spatialGraph from given directory */
AmiraSpatialGraph * helper::getSpatialGraph(const char * filename)
{
  	Reader * fileReader = new Reader(filename, filename);
	
	std::string fname(filename); 
	if((access(fname.c_str(), F_OK) != -1))
	{
		if (fname.compare(fname.size()-3,3,".am") == 0)
		{
			fileReader->readSpatialGraphFile(0);
			return fileReader->getSpatialGraph(); 
		}
		else if (fname.compare(fname.size()-4,4,".hoc") == 0)
		{
			fileReader->readHocFile();
			return fileReader->getSpatialGraph(); 
		}
		else
		{
			std::cout << "Inputfile is neither .hoc nor .am! Cannot read file" << std::endl; 
			return 0; 
		}
	}
	else
	{
		std::cout << "Inputfile not found: " << fname << std::endl; 
		return 0; 
	}
}

/* Gets two ImageDataPointerType as Input and creates .txt containing the normalized profile along x and z axis */
void helper::getProfileAsTxtImprov(ImageDataPointerType volume, ImageDataPointerType volumeOld, const char * path, double VoxelSize)
{   
	int dim[3]; 
	volume->GetDimensions(dim);
	
	int x_off = volume->GetExtent()[0]; 
	double profile[dim[0]][5];
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{
		profile[x-x_off][1] = 0.0; 
		profile[x-x_off][2] = 0.0; 
		profile[x-x_off][3] = 0.0; 
		profile[x-x_off][4] = 0.0; 
		profile[x-x_off][0] = x*VoxelSize; 
	}

	// Sum up values in volume
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{  
		for(int y = volume->GetExtent()[2]; y <= volume->GetExtent()[3]; ++y)
		{
			for(int z = volume->GetExtent()[4]; z <= volume->GetExtent()[5]; ++z)
			{ 
				double *px = static_cast< double * >(volume->GetScalarPointer(x,y,z));
				profile[x-x_off][1] += (*px); 
				profile[z-x_off][2] += (*px); 
				
				if (helper::isValidPt(x,y,z,volumeOld))
				{
					double *pxOld = static_cast< double * >(volumeOld->GetScalarPointer(x,y,z));
					profile[x-x_off][3] += (*pxOld); 
					profile[z-x_off+8][4] += (*pxOld);
				}
			}
		}
	}
	
	// Normalize 
	double x_norm =  1E9/(VoxelSize*VoxelSize*VoxelSize*dim[1]*dim[2]); // (50*50*50*y_dim*z_dim)	
	double z_norm = 1E9/(VoxelSize*VoxelSize*VoxelSize*dim[0]*dim[1]); // (50*50*50*y_dim*x_dim)
	for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
	{
		profile[x-x_off][1] *= x_norm; 
		profile[x-x_off][2] *= z_norm;
		profile[x-x_off][3] *= x_norm; 
		profile[x-x_off][4] *= z_norm;
	}
	
	//std::ofstream outStream(name.c_str());
	std::ofstream outStream(path);
	
	if(!outStream.fail())
	{	  
		outStream << "range x z x_old z_old  Dim: [" << dim[0] << "," << dim[1] << "," << dim[2] << "]" << std::endl;
	  
		for(int x = volume->GetExtent()[0]; x <= volume->GetExtent()[1]; ++x)
		{
			outStream << profile[x-x_off][0] << " " << profile[x-x_off][1] << " " << profile[x-x_off][2] << " " << profile[x-x_off][3] << " " << profile[x-x_off][4] << std::endl; 
		}
		outStream.close();
	}
	else
	{
		std::cout << "Error! Writing ProfileTxt failed! Path:" <<  path << std::endl; 
	}
}

/* Scales spatialGraph along cut/shrunk z-dimension. Scales it back to 300um by scaling the BB of the spatialGraph */
void helper::zScaling(AmiraSpatialGraph * spatialGraph)
{
	std::cout << "WARNING! Using default scaling of 300um! L2/3 INs have 350um!" << std::endl;
	zScaling(spatialGraph,300);
};

/* Scales spatialGraph along cut/shrunk z-dimension. Scales it back to original slicing thickness by scaling the BB of the spatialGraph */
void helper::zScaling(AmiraSpatialGraph* spatialGraph, double thickness)
{
	// Get Z coordinates for scaling	
	double r[2]; 
	getBoundingBox(Z_COORD,r,spatialGraph);
	double minZ = r[0];
	double maxZ = r[1];
	
	TransformPointerType scal = TransformPointerType::New();
	
	// Scale z-axis coordinates to Slicing thickness in micrometer
	double zfactor = thickness/(maxZ-minZ);
	scal->Scale(1,1,zfactor);
	
	std::cout << "    z-scaling factor = " << zfactor << " (Slicing thickness = " << thickness << " um)" << std::endl;
	
	spatialGraph->setTransformation(scal);
	spatialGraph->applyTransformation();
}

/* Compute the center of a label in the spatialGraph (most likely center of soma) 
 * can then be aligned with helper::align
 */
void helper::centerOfSpatialGraph(int neuriteID, double centerPt[3], AmiraSpatialGraph * spatialGraph)
{
	PolyDataPointerType structure = PolyDataPointerType::New();
	if(!spatialGraph->extractLandmark(neuriteID, structure))
	{
		std::cout << "Error! Could not find structure with ID " << neuriteID <<" in SpatialGraph!" << std::endl;
		return;
	}
	int subID;
	double pCoords[3], * weights;
	weights = new double[structure->GetCell(0)->GetNumberOfPoints()];
	structure->GetCell(0)->GetParametricCenter(pCoords);
	structure->GetCell(0)->EvaluateLocation(subID, pCoords, centerPt, weights);
	delete [] weights;     
}

/* Aligns spatialGraph to neuriteID (Soma) */
void helper::align(AmiraSpatialGraph * spatialGraph, int neuriteID)
{
	double centerPt[3]; 
	centerOfSpatialGraph(neuriteID, centerPt, spatialGraph); 
		
	// Subtract centerPt from spatialGraph
	for(int i=0; i<3; i++)
	{
	    centerPt[i] *= -1;
	}
  
	TransformPointerType sub = TransformPointerType::New();
	sub->Translate(centerPt);
	spatialGraph->setTransformation(sub);
	spatialGraph->applyTransformation();
}

/* Computes L2 Euclidean Norm */
double helper::norm(double vector[3])
{
	return (sqrt(pow(vector[0],2.0) + pow(vector[1],2.0) + pow(vector[2],2.0))); 
}

/* Draws #numSamples Gaussian Random Numbers given Mean and SD,
 * WARNING! Called in a loop will produce same output because random number generator is set within function
 */
double * helper::randn(double mean, double SD, int numSamples)
{
	gsl_rng_env_setup();	
	//const gsl_rng_type * T = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set (r,time(NULL));
	
	double * randnSample = new double[numSamples];
	
	for (int i = 0; i < numSamples; i++) 
	{
		randnSample[i] = mean + gsl_ran_gaussian_ratio_method(r, SD);
	}
	
	gsl_rng_free(r);
	return randnSample; 
}

/* Draws #numSamples Gaussian Random Numbers given Mean and SD, and random number generator  */
double * helper::randn(double mean, double SD, int numSamples, gsl_rng * r)
{
	double * randnSample = new double[numSamples];
	
	for (int i = 0; i < numSamples; i++) 
	{
		randnSample[i] = mean + gsl_ran_gaussian_ratio_method(r, SD);
	}
	return randnSample; 
}

/* Computes Cross Product between two 3D vectors */
double * helper::cross(double u[3], double v[3])
{
	double * crossProduct = new double[3];
	
	crossProduct[0] = u[1]*v[2] - u[2]*v[1]; 
	crossProduct[1] = u[2]*v[0] - u[0]*v[2]; 
	crossProduct[2] = u[0]*v[1] - u[1]*v[0]; 
	
	return crossProduct; 
}

/* Returns list of IDs to Branch Points in given spatialGraph for given neurite (axon, dendrite, soma)
 * WARNING! Only works for .hoc files. Requires fatherID!
 */
 std::list< int > helper::getBranchPointIDs(AmiraSpatialGraph * spatialGraph, int neuriteID)
{
	std::map< int, int > nrOfChildBranches;
	std::list< int > branchPtIDs;
	std::vector< Edge * >::const_iterator edgeIt;
		
	for(edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label == neuriteID)
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

bool helper::isPtOnEdge(double pt[3], double bounds[6])
{
	bool onEdge = (pt[X_COORD]==bounds[0] || pt[X_COORD]==bounds[1]
					   || pt[Y_COORD]==bounds[2] || pt[Y_COORD]==bounds[3]
					   || pt[Z_COORD]==bounds[4] || pt[Z_COORD]==bounds[5]);
    return onEdge;
}

/* Return whether 3D point is within 3D bounding box
 * [xmin xmax ymin ymax zmin zmax]
 */
bool helper::isPtWithinBounds(double pt[3], double bounds[6])
{
	bool outsideOfBounds = (pt[X_COORD]<bounds[0] || pt[X_COORD]>bounds[1]
					   || pt[Y_COORD]<bounds[2] || pt[Y_COORD]>bounds[3]
					   || pt[Z_COORD]<bounds[4] || pt[Z_COORD]>bounds[5]);

    return !outsideOfBounds;
}

// Computes t values between two points and a box (standing_vector + t*direction_vector)
std::list<double> helper::getIntersectionValue(double bounds[6],double pt1[3], double pt2[3])
{
	std::list < double> t_values;
	for (int i=0; i<2; i++)
	{
		for (int j=0; j<3; j++)
		{
			double tmp = computeIntersectionValue(pt1[j]-bounds[j*2+i],pt2[j]-bounds[j*2+i]);
			if (tmp != 0.0)
				t_values.push_back(tmp);
		}
	}
	return t_values;
}

// Returns t-value (pt1 + t * (pt2 - pt1))
double helper::computeIntersectionValue(double Dst1, double Dst2)
{
	// Check whether Line is intersecting with given plane
	if ((Dst1*Dst2) >= 0.0)
		return 0.0;
	if (Dst1 ==  Dst2)
		return 0.0;

	// Compute intersecting Point via standing vetor plus weighted direction vector
	double tmp = (-Dst1/(Dst2-Dst1 /*+ 1E-8 */));

	if (tmp>0.0 && tmp<1.0)
		return tmp;
	else
		return 0.0;
}

// Compute intersection point between two points pt1 and pt2 given the relative distance measurement t
void helper::getIntersectionPoint(double t, double pt1[3], double pt2[3], double intersectPt[3])
{
	for (int i = 0; i<3; i++)
	{
		intersectPt[i] = pt1[i] + (pt2[i]-pt1[i]) * t;
	}
}

/* Merges Local and NonLocal tree into one tree
 * Returns merged spatialGraph
 * Connects second Point of the NonLocal tree to the first point (SomaLocation) of the Local Tree
 */
AmiraSpatialGraph * helper::mergeLocalNonLocalAxon(AmiraSpatialGraph * sgLocal, AmiraSpatialGraph * sgNonLocal)
{
  	// First - Label the NonLocal Tree
	int currID = 0; 
	int NonLocal = 666; // Dummy
	for(std::vector< Edge * >::iterator edgeIt = sgNonLocal->edgesBegin(); edgeIt != sgNonLocal->edgesEnd(); ++edgeIt, ++currID)
	{
		if((*edgeIt)->label == Axon)
		{
			// Modify Edge Label
			(*(sgNonLocal->edgesPointer()))[currID]->label = NonLocal; 
			
			// Modify corresponding vertices label 
			int tmpVertexID = (*(sgNonLocal->edgesPointer()))[currID]->edgeConnectivity[0]; 
			(*(sgNonLocal->verticesPointer()))[tmpVertexID]->label = NonLocal;  
			
			tmpVertexID = (*(sgNonLocal->edgesPointer()))[currID]->edgeConnectivity[1];
			(*(sgNonLocal->verticesPointer()))[tmpVertexID]->label = NonLocal;
		}
	}
	
	// Second - Merge Local and NonLocal Tree
	AmiraSpatialGraph * finalSG = new AmiraSpatialGraph;  
	finalSG->mergeSpatialGraph(sgLocal);
	finalSG->mergeSpatialGraph(sgNonLocal);
	
	// Get Index of second Vertex (so that the two parts of the axon are connected close to the roots)
	// Somehow idxVertex.front() does not work, maybe because it's not a branching node already
	std::vector<int> idxVertex; 
	int n = 0; 
	for (std::vector< Vertex * >::iterator vertexIter = finalSG->verticesBegin(); vertexIter != finalSG->verticesEnd(); ++vertexIter, ++n)
	{  
		if ((*vertexIter)->label == Axon)
		{
			idxVertex.push_back(n);
		}
		
		if (idxVertex.size()>1) // Delete if clause in case of problems to have previous working version
			break; 
	}
	
	// Third - Connect Local and NonLocal Tree by setting the connections correct
	// Go through SG and modify the first point of the NonLocal SG to be the first point of the Local one and connect it to that point
	// Do this for the edges and vertices
	for(std::vector< Edge * >::iterator edgeIt = finalSG->edgesBegin(); edgeIt != finalSG->edgesEnd(); ++edgeIt)
	{
		if ( (*edgeIt)->label==NonLocal)
		{
// 			(*edgeIt)->edgeConnectivity[0] = n;
// 			(*edgeIt)->edgeConnectivity[0] = idxVertex.front();
			(*edgeIt)->edgeConnectivity[0] = idxVertex.back();
			
			double * tmpCoords = *((*edgeIt)->edgePointCoordinates.begin());
			
// 			Vertex * tmpVertex= finalSG->verticesPointer()->at(n); 
// 			Vertex * tmpVertex= finalSG->verticesPointer()->at(idxVertex.front()); 
			Vertex * tmpVertex= finalSG->verticesPointer()->at(idxVertex.back()); 
			tmpCoords[X_COORD] = tmpVertex->coordinates[X_COORD]; 
			tmpCoords[Y_COORD] = tmpVertex->coordinates[Y_COORD]; 
			tmpCoords[Z_COORD] = tmpVertex->coordinates[Z_COORD]; 
			break; 
		}
	}
	
	for(std::vector< Vertex * >::iterator vertexIter = finalSG->verticesBegin(); vertexIter != finalSG->verticesEnd(); ++vertexIter)
	{	  
		if ( (*vertexIter)->label==NonLocal)
		{
			Vertex * tmpVertex= finalSG->verticesPointer()->at(idxVertex.back()); 
// 			Vertex * tmpVertex= finalSG->verticesPointer()->at(idxVertex.front()); 
// 			Vertex * tmpVertex= finalSG->verticesPointer()->at(n); 
			(*vertexIter)->coordinates[X_COORD] = tmpVertex->coordinates[X_COORD]; 
			(*vertexIter)->coordinates[Y_COORD] = tmpVertex->coordinates[Y_COORD]; 
			(*vertexIter)->coordinates[Z_COORD] = tmpVertex->coordinates[Z_COORD]; 
			break; 
		}
	}
	
	// Fourth - Relabel Edges and Vertices of the NonLocal Tree as Axon
	currID = 0; 
	for(std::vector< Edge * >::iterator edgeIt = finalSG->edgesBegin(); edgeIt != finalSG->edgesEnd(); ++edgeIt, ++currID)
	{
		if((*edgeIt)->label == NonLocal)
		{
			// Modify Edge Label
			(*(finalSG->edgesPointer()))[currID]->label = Axon; 
			
			// Modify corresponding vertices label 
			int tmpVertexID = (*(finalSG->edgesPointer()))[currID]->edgeConnectivity[0]; 
			(*(finalSG->verticesPointer()))[tmpVertexID]->label = Axon;  
			
			tmpVertexID = (*(finalSG->edgesPointer()))[currID]->edgeConnectivity[1];
			(*(finalSG->verticesPointer()))[tmpVertexID]->label = Axon;
		}
	}
	
	for(std::vector< Vertex * >::iterator vertexIter = finalSG->verticesBegin(); vertexIter != finalSG->verticesEnd(); ++vertexIter)
	{	  
		if ( (*vertexIter)->label==NonLocal)
		{
			(*vertexIter)->label=Axon; 
		}
	}
	return finalSG; 
}

/* Merges two spatialGraphs (in silico Axon and in silico Dendrite)
 * Returns merged spatialGraph
 * Connects Axon and Dendrite to Soma!
 * Set Home Barrel to D2
 * Rotates to Cortical Column Axis (z) 
 */
AmiraSpatialGraph * helper::mergeAxonDend(AmiraSpatialGraph * dendSG, AmiraSpatialGraph * axonSG)
{
	AmiraSpatialGraph * finalSG = new AmiraSpatialGraph;  
	
	finalSG->mergeSpatialGraph(dendSG);
	finalSG->mergeSpatialGraph(axonSG);
	
	std::vector<int> idxSomaVertex; 
	int n = 0; 
	for (std::vector< Vertex * >::iterator vertexIter = finalSG->verticesBegin(); vertexIter != finalSG->verticesEnd(); ++vertexIter, ++n)
	{	  
		if ((*vertexIter)->label == Soma)
			idxSomaVertex.push_back(n);
	}
	
	for(std::vector< Edge * >::iterator edgeIt = finalSG->edgesBegin(); edgeIt != finalSG->edgesEnd(); ++edgeIt)
	{	  
		// Remove soma causes former edges connected to soma to be connected with themselves, thus connect them back to the soma and insert one point (equal to soma vertex)
		if ( (*edgeIt)->label==Axon )
		{
			(*edgeIt)->edgeConnectivity[0] = idxSomaVertex.back();
			double * tmpCoords = *((*edgeIt)->edgePointCoordinates.begin());
			
			Vertex * tmpVertex= finalSG->verticesPointer()->at(idxSomaVertex.back()); 
			tmpCoords[X_COORD] = tmpVertex->coordinates[X_COORD]; 
			tmpCoords[Y_COORD] = tmpVertex->coordinates[Y_COORD]; 
			tmpCoords[Z_COORD] = tmpVertex->coordinates[Z_COORD]; 
			break; 
		}
	}
	
	for(std::vector< Vertex * >::iterator vertexIter = finalSG->verticesBegin(); vertexIter != finalSG->verticesEnd(); ++vertexIter)
	{	  
		// Remove soma causes former edges connected to soma to be connected with themselves, thus connect them back to the soma and insert one point (equal to soma vertex)
		if ( (*vertexIter)->label==Axon )
		{
			Vertex * tmpVertex= finalSG->verticesPointer()->at(idxSomaVertex.back()); 
			(*vertexIter)->coordinates[X_COORD] = tmpVertex->coordinates[X_COORD]; 
			(*vertexIter)->coordinates[Y_COORD] = tmpVertex->coordinates[Y_COORD]; 
			(*vertexIter)->coordinates[Z_COORD] = tmpVertex->coordinates[Z_COORD]; 
			break; 
		}
	}
	
	finalSG->setHomeBarrel(D2);
	
	// Rotate in final coordinate system
	switchYZ(finalSG, 1); 
	
	return finalSG; 
}

/* Rotate SpatialGraph so that Y becomes Z (column axis) 
 * Input: 
 * 	-spatialGraph which should be rotated
 * 	- int d: which direction: 
 * 		 1 : y -> z (FINAL, z is column axis)
 * 		-1 : z -> y (ORIGINAL, y is column axis)
 */
void helper::switchYZ(AmiraSpatialGraph * spatialGraph, int d)
{
	if (!(spatialGraph->isLabelInSpatialGraph(Soma)))
	{
		std::cout << "ERROR! Soma was not found in given SpatialGraph! Cannot compute Soma Location for Rotation!" << std::endl;
		return; 
	}
	
	if (d != 1 && d != -1)
	{
		std::cout << "ERROR! d has to be either 1 (make z column axis) or -1 (make y column axis), not " << d << std::endl; 
		return; 
	}
  
	double somaPt[3]; 
	double shiftPt[3]; 
	centerOfSpatialGraph(Soma, somaPt, spatialGraph);
  
	// Set Soma Position to 0 0 0
	TransformPointerType sub1 = TransformPointerType::New();
	shiftPt[0] = -somaPt[0]; 
	shiftPt[1] = -somaPt[1]; 
	shiftPt[2] = -somaPt[2]; 
	sub1->Translate(shiftPt);
	spatialGraph->setTransformation(sub1);
	spatialGraph->applyTransformation();
	
	// Rotate in x-dimension 90°
	TransformPointerType rot = TransformPointerType::New();
	double rotVec[3] = {d,0,0}; 
	rot->RotateWXYZ(90, rotVec);
	spatialGraph->setTransformation(rot);
	spatialGraph->applyTransformation();
	
	// Set back to original Soma Position
	TransformPointerType sub2 = TransformPointerType::New();
	shiftPt[0] = somaPt[0]; 
	shiftPt[1] = somaPt[2]; 
	shiftPt[2] = somaPt[1]; 
	sub2->Translate(shiftPt);
	spatialGraph->setTransformation(sub2);
	spatialGraph->applyTransformation();
}

/* Gets list with all in silico axons and in silico dendrites as input. 
 * Loads the corresponding amira files and then randomly merges the dendrites and axons
 * In silico files should be in same folder as .txt lists. 
 * NumCells tells how many "completeCells" should be generated
 * zPosition specifies position of cell along cortical column (z-axis)
 * outputpath specifies where "completeCells" should be stored
 */
void helper::mergeAxonDendList(const char * axonList, const char * dendList, const char * outputpath, int numCells, double zPosition)
{
	std::string artificialListDend(dendList); 
	std::string artificialListAxon(axonList); 
	std::string artificialPathAxon = getRootFromPath(artificialListAxon.c_str()); 
	std::string artificialPathDend = getRootFromPath(artificialListDend.c_str()); 
	
	// Load dendritic and axonal file names in lists
	std::string currentLine; 
	std::vector<std::string> fnameListDend;
	std::vector<std::string> fnameListAxon;
	
	std::ifstream inputStreamDend(artificialListDend.c_str()); 
	if(!inputStreamDend.fail())
	{
		while(!inputStreamDend.eof())
		{
			getline(inputStreamDend,currentLine); 
			if (currentLine.empty())
				break; 
			currentLine = artificialPathDend + currentLine; 
			fnameListDend.push_back(currentLine);
		}
		inputStreamDend.close();
	}
	else
	{
		std::cout << "Could not find " << artificialListDend << std::endl;
	}
	
	std::ifstream inputStreamAxon(artificialListAxon.c_str()); 
	if(!inputStreamAxon.fail())
	{
		while(!inputStreamAxon.eof())
		{
			getline(inputStreamAxon,currentLine); 
			if (currentLine.empty())
				break; 
			currentLine = artificialPathAxon + currentLine; 
			fnameListAxon.push_back(currentLine);
		}
		inputStreamAxon.close();
	}
	else
	{
		std::cout << "Could not find " << artificialListAxon << std::endl;
	}
	
	helper::mkDir(outputpath); 
	
	srand (time(NULL));
	
	std::vector<std::string>::iterator itDend = fnameListDend.begin();
	std::vector<std::string>::iterator itAxon = fnameListAxon.begin();
	
	int c = 0; 	
	if (numCells > fnameListDend.size() || numCells==0)
		numCells = fnameListDend.size();
	
	if (numCells > fnameListAxon.size() || numCells==0) 
		numCells = fnameListAxon.size(); 
	
	if (fnameListDend.size()==0 || fnameListAxon.size()==0)
		std::cout << "List size is zero! Dendritic List Size: " << fnameListDend.size() << " Axonal List Size: " << fnameListAxon.size() << std::endl;
	
	for (; c != numCells; ++c)
	{
		int idxDend = rand() % fnameListDend.size();
		int idxAxon = rand() % fnameListAxon.size(); 
	  
		// Load Axonal and Dendritic Arbor
		AmiraSpatialGraph * sgDend = helper::getSpatialGraph(fnameListDend.at(idxDend).c_str()); 
		AmiraSpatialGraph * sgAxon = helper::getSpatialGraph(fnameListAxon.at(idxAxon).c_str()); // THIS WAS WRONG? fnameListDend?!?!
		fnameListDend.erase(fnameListDend.begin()+idxDend); 
		fnameListAxon.erase(fnameListAxon.begin()+idxAxon); 
		
		// Merge Dendrite and Axon
		AmiraSpatialGraph * sgFinal = helper::mergeAxonDend(sgDend,sgAxon); 
		
		// L23: 	 325.0 ->  412.0
		// L4: 		 737.5 ->    0.0
		// L5A: 	1012.5 -> -275.5
		// L5B: 	1262.5 -> -525.5
		// L6: 		1650.0 -> -913.0
		double centerPt[3]; 
		centerPt[0] = 0;
		centerPt[1] = 0;
		centerPt[2] = zPosition;
		
		TransformPointerType sub = TransformPointerType::New();
		sub->Translate(centerPt);
	    
		sgFinal->setTransformation(sub);
		sgFinal->applyTransformation();
		
		// Store Complete Cell
		std::stringstream outfnametmp; 
		outfnametmp << outputpath << "_completeCell_" << c;
		std::string outfname = outfnametmp.str();
		
		Reader * amWriter = new Reader((*itAxon).c_str(), outfname.c_str());
		amWriter->setSpatialGraph(sgFinal);
		amWriter->writeSpatialGraphFile();
		
		// Delete SG
		delete sgDend;
		delete sgAxon;
		delete sgFinal; 
		delete amWriter; 
		
		std::cout << outfname << " was created and stored!" << std::endl;
	}
}

/* Creates Directory / Folder in case it does not exist already */
void helper::mkDir(const char * path) 
{
	std::string fname(path); 
	std::string mkdir = "mkdir " + fname; 
	
	if (access(fname.c_str(), F_OK) == -1)
		system(mkdir.c_str());
}

/* executes command line
 * Example:
 * 	char * cmd = ""/bin/bash /home/dudvary/genAxons.sh yeyey""; 
 * 	std::string out = helper:: exec(cmd); 
 * 	std::cout << out << std::endl;
*/
std::string helper::exec(char* cmd)
{
	FILE* pipe = popen(cmd, "r");
	
	if (!pipe) 
		return "ERROR";
	
	char buffer[128];
	std::string result = "";
	
	while(!feof(pipe)) 
	{
		if(fgets(buffer, 128, pipe) != NULL)
			result += buffer;
	}
	
	pclose(pipe);
	return result;
}

/* Returns the output of the ls command to terminal
 * for example the list of all files that match criteria */
std::vector<std::string> helper::getReturnOfLsCmd(std::string path)
{
	std::string cmd =  "ls " + path;
	char * cstr = &cmd[0u];
	std::string txt = exec(cstr);

	std::vector<std::string> filenames;

    std::istringstream iss(txt);
    do
    {
        std::string subs;
        iss >> subs;

        if (subs.size()>0)
        	filenames.push_back(subs);
    } while (iss);
    return filenames;
}

/* Computes BranchOrder for each EdgeID
 * WARNING! Only works for .hoc input! Checks fatherID!
 */
std::map< int, int > helper::getBranchOrder(AmiraSpatialGraph * spatialGraph, int neuriteID)
{
	// Branch Orderfor each Edge/Segment
	std::map< int, int > edgeBO; 
	std::list< int > BranchPointIDs = getBranchPointIDs(spatialGraph, neuriteID); 
	
	int edgeID = 0; 
	for(std::vector< Edge * >::iterator edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt, ++edgeID)
	{ 
		if((*edgeIt)->label == neuriteID)
		{
			// Get Branch Order
			int branchOrder = 0; 
			int fatherID = (*edgeIt)->fatherID;
			while (fatherID != -1)
			{
				branchOrder++; 
				Edge * newEdge = (*(spatialGraph->edgesPointer()))[fatherID];
				fatherID = newEdge->fatherID; 
			}
			edgeBO[edgeID] = branchOrder; 
		}
	}
	return edgeBO; 
}

/* Registration of L5B (Hemberger) cells
 * Available landmarks: Pia, L4 Center, and WM
 * Pia is set at -706um; L4 Center at 0um; WM at 1250um (1957um)
 * Rescale Pia-L4Center and L4Center-WM */
void helper::registrateL5B(const char * hocPath, const char * txtPath, const char * registeredPath)
{
	std::ifstream inputStream(txtPath);
	std::string hocPathStr(hocPath);

	if(!inputStream.fail())
	{
		std::string currentLine;
		std::string registeredPathStr(registeredPath);

		// First Run to Get Global BoundingBox
		while(!inputStream.eof())
		{
			getline(inputStream,currentLine);
			if (currentLine.empty())
				break;

			size_t idx = currentLine.find(".hoc");
			std::string fname = hocPathStr + currentLine.substr(0,idx+4);

			if((access(fname.c_str(), F_OK) != -1))
			{
				// Extract vertical landmark positions of pia, L4 upper border, L4 lower border, WM
				size_t sep1 = currentLine.find(",", idx+5);
				size_t sep2 = currentLine.find(",", sep1+1);
				size_t sep3 = currentLine.find(",", sep2+1);

				std::string val1 = currentLine.substr(idx+5,sep1-idx-5);
				std::string val2 = currentLine.substr(sep1+1,sep2-sep1-1);
				std::string val3 = currentLine.substr(sep2+1,sep3-sep2-1);
				std::string val4 = currentLine.substr(sep3+1);

				bool granular = true;

				if (val1.empty())
				{
					std::cout << "Pia is not defined! Error!" << std::endl;
					continue;
				}
				if (val4.empty())
				{
					std::cout << "WM is not defined! Error!" << std::endl;
					continue;
				}
				if (val2.empty() && val3.empty()) // Granular layer not defined
				{
					granular = false;
				}
				else if (val2.empty() || val3.empty())  // TODO: Only L4 Upper Border or L4 Lower Border defined
				{
					std::cout << "This case is not implemented yet, skipped!" << std::endl;
					continue;
				}

				double wm = NAN;
				double pia = NAN;
				double L4U = NAN;
				double L4L = NAN;
				sscanf(val1.c_str(),"%lf", &pia);
				sscanf(val4.c_str(),"%lf", &wm);

				if (granular)
				{
					sscanf(val2.c_str(),"%lf", &L4U);
					sscanf(val3.c_str(),"%lf", &L4L);
				}

				// Read cells
				std::string currentOut = registeredPathStr + currentLine.substr(0,idx+4);
				Reader * fileReader = new Reader(fname.c_str(), currentOut.c_str());
				fileReader->readHocFile();

				//std::cout << pia << "," << L4U << "," << L4L << "," << wm << std::endl;

				registerIN(fileReader->getSpatialGraph(),pia,L4U,L4L,wm,D2);

//				// Shift Center of Mass of Soma to Origin (0,0,0)
//				align(fileReader->getSpatialGraph(), Soma);
//
//				// Rescale Pia-L4Center and L4Center and WM
//				scaleToD2L4Center(fileReader->getSpatialGraph(),wm,l4center,pia);
//
//				// Rescale along cut dimension to 300um
//				zScaling(fileReader->getSpatialGraph(),300);
//
//				// Check whether the cell penetrates the pia
//				shiftDown(fileReader->getSpatialGraph());
//
//				// Rotation (switch y and z) by rotating 90°
//				switchYZ(fileReader->getSpatialGraph(),1);
//
				//Write File
				fileReader->writeHocFile();
			}
		}
		inputStream.close();
	}
	else
	{
		std::cout << "ERROR! " << txtPath << " not found!" << std::endl;
	}
}

void helper::registerIN(AmiraSpatialGraph * spatialGraph, double pia, double L4U, double L4L, double wm, int colID)
{
	// Align Spatial Graph to Soma
	align(spatialGraph, Soma);

	bool granular = true;

	if (isnan(L4L) && isnan(L4U))
	{
		granular = false;
	}

	// Get average Barrel Field
	BarrelField * BF = new BarrelField();
	std::map< int, Column * > avgColumns = BF->avgColumns;
	Column * refColumn = avgColumns.find(colID)->second;
	double * avgPia = refColumn->top;
	double * avgWM = refColumn->bottom;
	double avgColumnHeigth = refColumn->getHeight();
	double newVerticalSomaPosition = 0;
	double vertScalingFactor = avgColumnHeigth/(pia-wm);

	std::cout << avgPia[2] << "," << avgWM[2] << std::endl;

	// Scale Supragranular, granular, and infragranular differently, but only if the SG, G, and IG ratio is closer to 1 than
	// Column ratio. We want to avoid scaling too much.
	if (granular)
	{
		std::map< int, Column * > avgBarrels = BF->avgBarrels;
		Column * refBarrel = avgBarrels.find(colID)->second;
		double * avgL4U = refBarrel->top;
		double * avgL4L = refBarrel->bottom;

		double vertScalingFactorSG = (avgPia[Z_COORD]-avgL4U[Z_COORD])/(pia-L4U);
		double vertScalingFactorG = (avgL4U[Z_COORD]-avgL4L[Z_COORD])/(L4U-L4L);
		double vertScalingFactorIG = (avgL4L[Z_COORD]-avgWM[Z_COORD])/(L4L-wm);

		double v1 = fabs(vertScalingFactorSG-1);
		double v2 = fabs(vertScalingFactorG-1);
		double v3 = fabs(vertScalingFactorIG-1);

		//if v1<




		//std::cout << "AVG: " << avgPia[Z_COORD] << " " << avgL4U[Z_COORD] << " " << avgL4L[Z_COORD] << " " << avgWM[Z_COORD] << std::endl;
		//std::cout << "AVG: " << avgPia[Z_COORD]-avgL4U[Z_COORD] << " " << avgL4U[Z_COORD]-avgL4L[Z_COORD] << " " << avgL4L[Z_COORD]-avgWM[Z_COORD] << std::endl;

		//std::cout << "CELL: " << pia << " " << L4U << " " << L4L << " " << wm << std::endl;
		//std::cout << "CELL: " << pia-L4U << " " << L4U-L4L << " " << L4L-wm << std::endl;

		//std::cout << vertScalingFactorSG << " " << vertScalingFactorG << " " << vertScalingFactorIG << std::endl;

	}
	else // Scale ColumnHeigth to AverageColumnHeigth
	{
		TransformPointerType scaling = TransformPointerType::New();

		scaling->Scale(1,vertScalingFactor,1);

		std::cout << "Vertical axis scaled by " << vertScalingFactor << " (" << (pia-wm) << "um to " << avgColumnHeigth << "um)" << std::endl;

		spatialGraph->setTransformation(scaling);
		spatialGraph->applyTransformation();

		newVerticalSomaPosition = pia * vertScalingFactor;
		std::cout << "New Pia-Soma-Distance " << newVerticalSomaPosition << " (from " << pia << "um)" << std::endl;
	}

	// Set Soma Position correct
	std::cout << newVerticalSomaPosition << " " << avgPia[2] << std::endl;

	double pt[3] = {0, avgPia[2]-newVerticalSomaPosition, 0};
	TransformPointerType sub = TransformPointerType::New();
	sub->Translate(pt);
	spatialGraph->setTransformation(sub);
	spatialGraph->applyTransformation();

	// Check whether the cell penetrates the pia, if so shift it down
	shiftDown(spatialGraph);

	// Rotation (switch y and z) by rotating 90°
	// z is now vertical column axis, y is slicing axis
	switchYZ(spatialGraph,1);

	spatialGraph->setHomeBarrel(colID);
}

/* DEPRECATED!
 * Registration of amirasave files (no asc!)
 * hocPath: path where hoc files are located 			"/home/dudvary/Interneurons/hoc/"
 * txtPath: path to .txt containing scale factor and offset 	"/home/dudvary/Interneurons/txt/amirasave.txt"
 *  
 * Stores files in subfolder in hocPath named registrated
 * Stored files are properly scaled and rotated ! (y -> z). z = 0 is barrel center, pia at -706um 
 */ 
void helper::registrateAmiraSave(const char * hocPath, const char * txtPath)
{
	std::ifstream inputTxt(txtPath); 
	std::string hocPathStr(hocPath);
	
	if(!inputTxt.fail())
	{
		std::string currentTxt;
		
		// First Run to Get Global BoundingBox
		while(!inputTxt.eof())
		{
			getline(inputTxt,currentTxt); 
			if (currentTxt.empty())
				break; 
			
			// Get offset and scaling from .txt file, convert to double 
			size_t idx1 = currentTxt.find(","); 
			size_t idx2 = currentTxt.find(",",idx1+1); 
			
			std::string fname = currentTxt.substr(0,idx1); 
			std::string tmp = currentTxt.substr(idx1+1,idx2-idx1-1); 
			double scalefactor = atof(tmp.c_str()); 
			
			tmp = currentTxt.substr(idx2+1,currentTxt.size()-idx2-1); 
			double offset = atof(tmp.c_str()); 
			
			// Read File
			std::string currentHoc = hocPathStr + fname; 
			std::string currentOut = hocPathStr + "registrated/" + fname;
			
			Reader * fileReader = new Reader(currentHoc.c_str(), currentOut.c_str());
			fileReader->readHocFile();
			
			// Shift along y
			double somaPt[3]; 
			centerOfSpatialGraph(Soma, somaPt, fileReader->getSpatialGraph()); 
			double pt[3] = {0, offset, 0};
			TransformPointerType sub = TransformPointerType::New();
			sub->Translate(pt);
			fileReader->getSpatialGraph()->setTransformation(sub);
			fileReader->getSpatialGraph()->applyTransformation();
			
			// Scale 
			TransformPointerType scal = TransformPointerType::New();
			scal->Scale(scalefactor,scalefactor,scalefactor);
			fileReader->getSpatialGraph()->setTransformation(scal);
			fileReader->getSpatialGraph()->applyTransformation();
			
			// Rescale along cut dimension to 300um
			zScaling(fileReader->getSpatialGraph(),300);
			
			// Final shift along y -> Barrel Center = 0, Pia = -706
			// Align to soma in x-z
			double ptBC[3] = {-somaPt[0], 706, -somaPt[2]};
			TransformPointerType sub2 = TransformPointerType::New();
			sub2->Translate(ptBC);
			fileReader->getSpatialGraph()->setTransformation(sub2);
			fileReader->getSpatialGraph()->applyTransformation();
			
			// In case the cell penetrates through the pia, shift the cell deeper down
			shiftDown(fileReader->getSpatialGraph());  
			
			std::cout << " >> " << fname << " was shifted by " << offset << " and scaled by " << scalefactor << std::endl; 
			
			// Rotate 
			helper::switchYZ(fileReader->getSpatialGraph(),1); 
			
			//Write File
			fileReader->writeHocFile(); 
		}
		inputTxt.close(); 
	}
	else
		std::cout << "Could not find " << txtPath << std::endl; 
}

/* Shift SpatialGraph into deeper layers if it penetrates through pia (706) */
// see updated function in registerIN.cpp
void helper::shiftDown(AmiraSpatialGraph * spatialGraph)
{
	double bounds[6]; 
	spatialGraph->getBoundingBox(Axon,bounds);
	double DendBounds[6]; 
	spatialGraph->getBoundingBox(Dendrite,DendBounds);
	
	double maxZ = bounds[3]; 
	if (DendBounds[3]>maxZ)
		maxZ = DendBounds[3]; 
	
	if (maxZ>706)
	{
		std::cout << "    OUT OF PIA by " << maxZ-706 << " The Cell is shifted deeper !!!" << std::endl;
		
		double pt[3] = {0, -maxZ+706, 0};
		TransformPointerType sub = TransformPointerType::New();
		sub->Translate(pt);
		spatialGraph->setTransformation(sub);
		spatialGraph->applyTransformation();
	}
}


/* Registration form ASC (i.e. requires .txt file containing center of mass of each layer from ascii (pia, L4/5 border, wm)
 * txtPath = "/home/dudvary/Interneurons/txt/layerASC.txt"; 
 * hocPath = "/home/dudvary/Interneurons/hoc/";
 * rescales layer 1 to 4 to 888um (like in D2) and layer 5 to 6 to 1957.0-888.0
 * shifts coordinates to Barrel Center (706um) 
 * rotates spatialgraphs accordingly (z becomes column axis)
 */
void helper::registrateFromAsc(const char * hocPath, const char * txtPath)
{
	std::ifstream inputStream(txtPath); 
	std::string hocPathStr(hocPath);
	
	if(!inputStream.fail())
	{
		std::string currentLine;
		
		// First Run to Get Global BoundingBox
		while(!inputStream.eof())
		{
			getline(inputStream,currentLine); 
			if (currentLine.empty())
				break; 
			
			size_t idx = currentLine.find(".hoc"); 
			std::string fname = hocPathStr + currentLine.substr(0,idx+4); 
			
			if((access(fname.c_str(), F_OK) != -1))
			{
				size_t sep1 = currentLine.find(",", idx+5);
				size_t sep2 = currentLine.find(",", sep1+1);
			  
				std::string loc1 = currentLine.substr(idx+5,sep1-idx-5); 
				std::string loc2 = currentLine.substr(sep1+1,sep2-sep1-1); 
				std::string loc3 = currentLine.substr(sep2+1); 
				
				double wm; 
				sscanf(loc1.c_str(),"%lf", &wm); 
				double l45; 
				sscanf(loc2.c_str(),"%lf", &l45); 
				double pia; 
				sscanf(loc3.c_str(),"%lf", &pia); 
				
				// Read cells
				std::string currentOut = hocPathStr + "registrated/" + currentLine.substr(0,idx+4);
				Reader * fileReader = new Reader(fname.c_str(), currentOut.c_str());
				fileReader->readHocFile();
				
				// Rescale infragranular (5-6) and supragranular layer (1-3) + layer 4
				scaleToD2(fileReader->getSpatialGraph(),wm,l45,pia); 
				
				// Rescale along cut dimension to 300um
				zScaling(fileReader->getSpatialGraph(),300);
				
				// Check whether the cell penetrates the pia
				shiftDown(fileReader->getSpatialGraph());
				
				// Rotation (switch y and z) by rotating x 90Â°
				switchYZ(fileReader->getSpatialGraph(),1); 
				
				//Write File
				fileReader->writeHocFile(); 
			}
		}
		inputStream.close(); 
	}
}

/* Registration for L4 IN (Koelbl et al., 2013) files where only Pia layer location is known
 * txtPath = "/home/dudvary/Interneurons/txt/layerPiaASC_L4.txt"; 
 * hocPath = "/home/dudvary/Interneurons/hoc/";
 * shifts coordinates to Barrel Center (706um) 
 * rotates spatialgraphs accordingly (z becomes column axis)
 */
void helper::registrateFromPia(const char * hocPath, const char * txtPath)
{
	std::ifstream inputStream(txtPath); 
	std::string hocPathStr(hocPath);
	
	if(!inputStream.fail())
	{
		std::string currentLine;
		
		while(!inputStream.eof())
		{
			getline(inputStream,currentLine); 
			if (currentLine.empty())
				break; 
			
			size_t idx = currentLine.find(".hoc"); 
			std::string fname = hocPathStr + currentLine.substr(0,idx+4); 
			
			if((access(fname.c_str(), F_OK) != -1))
			{
				size_t sep1 = currentLine.find(",", idx+5);			  
				std::string loc1 = currentLine.substr(idx+5,sep1-idx-5); 
				
				double pia; 
				sscanf(loc1.c_str(),"%lf", &pia); 
				
				// Read cells
				std::string currentOut = hocPathStr + "registrated/" + currentLine.substr(0,idx+4);
				Reader * fileReader = new Reader(fname.c_str(), currentOut.c_str());
				fileReader->readHocFile();
				
				// Rescale along cut dimension to 300um
				zScaling(fileReader->getSpatialGraph(),300);
				
				// Final shift along y -> Barrel Center = 0, Pia = -706
				double somaPt[3]; 
				centerOfSpatialGraph(Soma, somaPt, fileReader->getSpatialGraph()); 
				double ptBC[3] = {-somaPt[0], 706-pia, -somaPt[2]};
				TransformPointerType sub2 = TransformPointerType::New();
				sub2->Translate(ptBC);
				fileReader->getSpatialGraph()->setTransformation(sub2);
				fileReader->getSpatialGraph()->applyTransformation();
				
				// Check whether cell penetrates pia
				shiftDown(fileReader->getSpatialGraph());
				
				// Rotation (switch y and z) by rotating x 90Â°
				switchYZ(fileReader->getSpatialGraph(),1); 
				
				//Write File
				fileReader->writeHocFile(); 
			}
		}
		inputStream.close(); 
	}
}

/* Scales SpatialGraph to D2 column proportions
 * Requires spatialGraph, column position of WM, Border of L4 to L5, and Pia
 * Scales area between Pia and L45 to 888.0 and area between L45 and WM to (1957.0-888.0)
 * Aligns spatialGraph so that Barrel Center = 0.0 (Distance between L45 to L4 Center is 194.0
 */
void helper::scaleToD2(AmiraSpatialGraph * spatialGraph, double wm, double l45, double pia)
{
	double dist1 = fabs(pia-l45); 
	double dist2 = fabs(l45-wm); 
	
	double scale1 = 888.0/dist1; 		// Pia to L4/5-Border distance in D2
	double scale2 = (1957.0-888.0)/dist2; 	// L4/5-Border to WM distance in D2
	
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt)
	{
		std::list< double * >::iterator edgeListIt;
		for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
		{
			double tmpY = (*edgeListIt)[Y_COORD];
			if(tmpY > l45) // Above, between L4/5-Border and Pia
			{
				(*edgeListIt)[Y_COORD] = ((*edgeListIt)[Y_COORD] - l45) * scale1 - 194.0; 
// 				std::cout << "Above: " << tmpY << " Pia: " << pia << " L4/5-Border: " << l45 << " WM: " << wm << std::endl; 				
			}
			else // Below, between L4/5-Border and WM
			{
				(*edgeListIt)[Y_COORD] = ((*edgeListIt)[Y_COORD] - l45) * scale2 - 194.0; 
// 				std::cout << "Below: " << tmpY << " Pia: " << pia << " L4/5-Border: " << l45 << " WM: " << wm << std::endl; 
			}
		}
	}
		
	// Shift along y and align soma in x-z (Pia = 0)
	double somaPt[3]; 
	centerOfSpatialGraph(Soma, somaPt, spatialGraph); 
	double pt[3] = {-somaPt[0], 0, -somaPt[2]};
	TransformPointerType sub = TransformPointerType::New();
	sub->Translate(pt);
	spatialGraph->setTransformation(sub);
	spatialGraph->applyTransformation();
}

/* Registration of L23 Cells with no futher Information (Helmstaedter et al., 2009)
   put them in center of L23 (325um from pia, 412.5um from Barrel Center)
   aligned to Soma. 
 */
void helper::registrateL23(const char * hocPath, const char * listPath)
{
 	std::ifstream inputStream(listPath); 
	std::string hocPathStr(hocPath);
	
	if(!inputStream.fail())
	{
		std::string currentLine;
		
		while(!inputStream.eof())
		{
			getline(inputStream,currentLine); 
			if (currentLine.empty())
				break; 
			
			std::string fname = hocPathStr + currentLine; 
			
			// Read cells
			std::string currentOut = hocPathStr + "registrated/" + currentLine;
			Reader * fileReader = new Reader(fname.c_str(), currentOut.c_str());
			fileReader->readHocFile();
			
			// Rescale along cut dimension to 350um
			zScaling(fileReader->getSpatialGraph(),350);
			
			// Final shift along y -> Barrel Center = 0 (737.5), L23 = 325 => 412.5
			double somaPt[3]; 
			centerOfSpatialGraph(Soma, somaPt, fileReader->getSpatialGraph()); 
			double ptBC[3] = {-somaPt[0], 412.5 -somaPt[1], -somaPt[2]};
			TransformPointerType sub2 = TransformPointerType::New();
			sub2->Translate(ptBC);
			fileReader->getSpatialGraph()->setTransformation(sub2);
			fileReader->getSpatialGraph()->applyTransformation();
			
			// If cell exceeds pia boundary, the cell shifted so that no cell parts exceed pia! 
			shiftDown(fileReader->getSpatialGraph()); 
			
			// Rotation (switch y and z) by rotating x 90Â°
			switchYZ(fileReader->getSpatialGraph(),1); 
			
			//Write File
			fileReader->writeHocFile(); 
		}
		inputStream.close(); 
	} 
}


void helper::flipYSpatialGraph(AmiraSpatialGraph * spatialGraph)
{
  	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label >= Neuron && (*edgeIt)->label <= Soma)
		{
			std::list< double * >::iterator edgeListIt;
			for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{
				(*edgeListIt)[Y_COORD] = -(*edgeListIt)[Y_COORD];
			}
		}
	}
}

std::vector<unsigned int> helper::returnCellIDs(const char * inputFilenameList)
{
  	std::vector<unsigned int> CellIDs;
	std::string currentLine;
	std::ifstream inputStream(inputFilenameList);
	if(!inputStream.fail())
	{
		while(!inputStream.eof())
		{
			getline(inputStream,currentLine);
			if (currentLine.empty())
				break;
			unsigned int CellID;
			std::istringstream(currentLine) >> CellID;
			CellIDs.push_back(CellID);
		}
		inputStream.close();
	}
	else
	{
		std::cout << "Could not find " << inputFilenameList << std::endl;
	}

	return CellIDs;
}

std::vector<std::string> helper::returnNames(const char* inputFilenameList, const char* root)
{
  	std::vector<std::string> fnameList;
	std::string currentLine; 
	std::ifstream inputStream(inputFilenameList); 
	if(!inputStream.fail())
	{
		std::string rootStr(root); 
	  
		while(!inputStream.eof())
		{
			getline(inputStream,currentLine); 
			if (currentLine.empty())
				break; 
			currentLine = rootStr + currentLine; 
			fnameList.push_back(currentLine);
		}
		inputStream.close();
	}
	else
	{
		std::cout << "Could not find " << inputFilenameList << std::endl;
	}
	
	return fnameList; 
}

/* Scale in silico Dendrites to compensate for too small dendritic extent compared to in vitro data
 * scaleLat: 0.9736 (1.0649)
 * scaleVert: 1.12 (1.1182)
 */
void helper::scaleDendrites(AmiraSpatialGraph * spatialGraph, double scaleLat, double scaleVert, bool registrated)
{
	TransformPointerType scal = TransformPointerType::New();

	// Scale columnar axis, take care of the rotation (is y columnar axis (non-registrated), or z columnar axis (registrated)?)
	if (registrated)
	{
		scal->Scale(scaleLat,scaleLat,scaleVert);
	}
	else
	{
		scal->Scale(scaleLat,scaleVert,scaleLat);
	}

	spatialGraph->setTransformation(scal);
	spatialGraph->applyTransformation();

	if (spatialGraph->isLabelInSpatialGraph(Axon))
	{
		std::cout << "WARNING! Axon scaled as well!" << std::endl;
	}
}

double helper::computeEuclideanDistance(double pt1[3], double pt2[3])
{
	double dist = sqrt( pow( ( pt1[X_COORD] - pt2[X_COORD] ), 2.0 ) +
						pow( ( pt1[Y_COORD] - pt2[Y_COORD] ), 2.0 ) +
						pow( ( pt1[Z_COORD] - pt2[Z_COORD] ), 2.0 ));
	return dist;
}

// Computes the Surface Area, (specifically for Dendrites)
// Connects primary dendrites with Soma.
// Returns Lateral Surface Area using A = PI * (r1 + r2) * sqrt(height^2 + (r1 - r2)^2)
// 	https://www.quora.com/What-is-the-area-of-a-cylinder-of-different-radius-at-top-and-bottom-respectively
// 			misses Top/Bottom Surface area at sides PI * (r1^2 + r2^2)
//			Total Surface Area = LateralSurfaceArea + Top/BottomSurfaceArea
// with height being the length between two points and r1 and r2 being the radius of the two points
// NOTE: In .hoc files the radius denotes the diameter, thus in amira files, the radius also denotes the diameter!
// If computing the "true" surface area the radius needs to be divided by two.
// If computing the surface area like the NetworkAnalyzer in Amira, the diameter is set to be the radius.
//	Robert: Amira Bug: uses diameter instead of radius; (oesn't matter for end result, but it's affecting
//			the INH PST density -> need to be consistent...)
double helper::computeSurfaceArea(AmiraSpatialGraph * SpatialGraph, int neuriteID)
{
	return helper::computeSurfaceArea(SpatialGraph, neuriteID, false);
}
// If likeAmira true, computes "surface area" as NetworkAnalyzer in Amira does it, by setting the radius = diameter
double helper::computeSurfaceArea(AmiraSpatialGraph * SpatialGraph, int neuriteID, bool likeAmira)
{
	double surfaceArea = 0.0;

	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = SpatialGraph->edgesBegin(); edgeIt != SpatialGraph->edgesEnd(); ++edgeIt)
	{

		if((*edgeIt)->label == neuriteID)
		{
			// Point
			std::list< double * >::iterator edgeListIt;
			edgeListIt = (*edgeIt)->edgePointCoordinates.begin();
			double * previousPt = *edgeListIt;
			++edgeListIt;

			// Radius
			std::list< double >::iterator edgeRadiusListIt;
			edgeRadiusListIt = (*edgeIt)->radiusList.begin();
			double previousRad = *edgeRadiusListIt;
			++edgeRadiusListIt;

			if (likeAmira)
				previousRad = previousRad/2.0;

			// Go through all points of each Edge starting with the second edge
			for(; edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt, ++edgeRadiusListIt)
			{
				double * currentPt = *edgeListIt;
				double currentRad = *edgeRadiusListIt;

				if (likeAmira)
					currentRad = currentRad/2.0;

				double heigth = sqrt(pow((currentPt[X_COORD]-previousPt[X_COORD]),2.0) +
									pow((currentPt[Y_COORD]-previousPt[Y_COORD]),2.0) +
									pow((currentPt[Z_COORD]-previousPt[Z_COORD]),2.0));
				double lateralArea = PI * (currentRad+previousRad) * sqrt(pow(heigth,2.0) +
									pow(currentRad-previousRad,2.0));
				surfaceArea += lateralArea;

				previousPt = currentPt;
				previousRad = currentRad;
			}
		}
	}

	return surfaceArea;
}

bool helper::isInhibitory(int celltype)
{
	return (celltype>=22 && celltype<=45);
}

/* Extract Morphologies from Morphologies.bb/.am file of a certain CellType
 * Input:
 * - spatialGraphSetFilename: path/to/Morphologies.bb file
 * - CellType: CellType ID that should be extracted
 * - outputpath: path where extracted morphologies should be stored to
 * stores morphologies in outputpath
 * stores soma.csv in outputpath (contains soma position of each morphology)
 */
void helper::extractMorphologies(const char * spatialGraphSetFilename, int CellType, const char * outputpath)
{
	double bounds[6];

	bounds[0] = -(std::numeric_limits<double>::infinity());
	bounds[1] = std::numeric_limits<double>::infinity();
	bounds[2] = -(std::numeric_limits<double>::infinity());
	bounds[3] = std::numeric_limits<double>::infinity();
	bounds[4] = -(std::numeric_limits<double>::infinity());
	bounds[5] = std::numeric_limits<double>::infinity();

	extractMorphologies(spatialGraphSetFilename, CellType, outputpath, bounds);
}

/* Extract Morphologies from Morphologies.bb/.am file of a certain CellType
 * Input:
 * - spatialGraphSetFilename: path/to/Morphologies.bb file
 * - CellIDsFilename: path/to/CellIDs.txt file containing all CellIDs to be extracted
 * - outputpath: path where extracted morphologies should be stored to
 * stores morphologies in outputpath
 */
void helper::extractMorphologies(const char * spatialGraphSetFilename, const char * CellIDsFilename, const char * outputpath)
{
	mkDir(outputpath);
	std::string outputpathstr(outputpath);
	std::ifstream inputStream(CellIDsFilename);

	if(!inputStream.fail())
	{
		// Read in SpatialGraphSet
		std::vector< unsigned int > originalGraphIndices;
		std::vector< unsigned int > cellTypeIDs;
		std::vector< double * > spatialGraphTransforms;
		std::vector< std::string > originalGraphFiles;
		std::map< unsigned int, std::string > cellTypeIDLabels;
		Reader::readSpatialGraphSetFile(spatialGraphSetFilename, originalGraphIndices, cellTypeIDs, spatialGraphTransforms, originalGraphFiles, cellTypeIDLabels);

		// Read in CellIDs
		std::string currentLine;

		while(!inputStream.eof()) // go through all CELLIDs
		{
			getline(inputStream,currentLine);
			if (currentLine.empty())
				break;

			int CellID = atoi(currentLine.c_str());

			// Get Filename
			std::string originalName = originalGraphFiles[originalGraphIndices[CellID]];
			originalName = originalName.substr(1, originalName.size()-2);
			std::string loadName = getRootFromPath(spatialGraphSetFilename) + originalName;

			// Outputfilename
			std::ostringstream CellIDStr;   // stream used for the conversion
			CellIDStr << CellID;      // insert the textual representation of 'Number' in the characters in the stream
			std::string saveName = outputpathstr + CellIDStr.str();

			// Load Spatial Graph
			Reader * spatialGraphReader = new Reader(loadName.c_str(), saveName.c_str());
			spatialGraphReader->readSpatialGraphFile(0);
			AmiraSpatialGraph * sg = spatialGraphReader->getSpatialGraph();

			// Apply Transformation
			sg->setTransformation(amiraToVTKTransform(spatialGraphTransforms[CellID]));
			sg->applyTransformation();

			// Store SpatialGraph
			spatialGraphReader->writeSpatialGraphFile();

			delete sg;
			delete spatialGraphReader;
		}
		inputStream.close();
	}
	else
	{
		std::cout << "Error! Reading file with CellIDs failed! Path: " <<  CellIDsFilename << std::endl;
	}
}

/* Extract Morphologies from Morphologies.bb/.am file of a certain CellType
 * Input:
 * - spatialGraphSetFilename: path/to/Morphologies.bb file
 * - CellType: CellType ID that should be extracted
 * - outputpath: path where extracted morphologies should be stored to
 * - bounds[6]: [xmin xmax ymin ymax zmin zmax]: bounding box of to be extracted morphologies.
 * 	 morphologies outside this box will not be extracted/stored
 * stores morphologies in outputpath
 * stores soma.csv in outputpath (contains soma position of each morphology)
 * stores cells.txt in outputpath (contains list of all morphology names)
 */
void helper::extractMorphologies(const char * spatialGraphSetFilename, int CellType, const char * outputpath, double bounds[6])
{
	mkDir(outputpath);
	std::string outputpathstr(outputpath);
	std::string outputFilenameSoma = outputpathstr + "soma.csv";
	std::ofstream outStreamSoma(outputFilenameSoma.c_str());

	std::string outputFilenameCells = outputpathstr + "cells.txt";
	std::ofstream outStreamCells(outputFilenameCells.c_str());

	if(!outStreamSoma.fail() && !outStreamCells.fail())
	{
		std::vector< unsigned int > originalGraphIndices;
		std::vector< unsigned int > cellTypeIDs;
		std::vector< double * > spatialGraphTransforms;
		std::vector< std::string > originalGraphFiles;
		std::map< unsigned int, std::string > cellTypeIDLabels;
		Reader::readSpatialGraphSetFile(spatialGraphSetFilename, originalGraphIndices, cellTypeIDs, spatialGraphTransforms, originalGraphFiles, cellTypeIDLabels);
		int dendCounter = 0;

		outStreamSoma << "OriginalName,StoredName,SomaX,SomaY,SomaZ,CellID," << std::endl;

		for(int i = 0; i < cellTypeIDs.size(); ++i)
		{
			if (cellTypeIDs[i]==CellType)
			{
				// Get Filename
				std::string originalName = originalGraphFiles[originalGraphIndices[i]];
				originalName = originalName.substr(1, originalName.size()-2);

				std::ostringstream dendCounterStr;   // stream used for the conversion
				dendCounterStr << dendCounter;      // insert the textual representation of 'Number' in the characters in the stream

				std::string saveName = outputpathstr + "Dend_" + dendCounterStr.str();
				std::string loadName = getRootFromPath(spatialGraphSetFilename) + originalName;

				// Load Spatial Graph
				Reader * spatialGraphReader = new Reader(loadName.c_str(), saveName.c_str());
				spatialGraphReader->readSpatialGraphFile(0);
				AmiraSpatialGraph * sg = spatialGraphReader->getSpatialGraph();

				// Apply Transformation
				sg->setTransformation(amiraToVTKTransform(spatialGraphTransforms[i]));
				sg->applyTransformation();

				// Get Soma Position
				double centerPt[3];
				centerOfSpatialGraph(Soma, centerPt, sg);

				// Check whether Soma is within bounding box
				if (bounds[0]<=centerPt[0] && bounds[1]>=centerPt[0] && bounds[2]<=centerPt[1] && bounds[3]>=centerPt[1] && bounds[4]<=centerPt[2] && bounds[5]>=centerPt[2])
				{
					// Write out name and soma position
					outStreamSoma << originalName << "," << saveName << "," << centerPt[0] << "," << centerPt[1] << "," << centerPt[2] << "," << i << "," << std::endl;
					outStreamCells << "Dend_" << dendCounterStr.str() << ".am" << std::endl;

					// Store SpatialGraph
					spatialGraphReader->writeSpatialGraphFile();
					dendCounter++;
				}

				delete sg;
				delete spatialGraphReader;
			}
		}
		outStreamSoma.close();
		outStreamCells.close();
	}
	else
	{
		std::cout << "Error! Writing soma.csv failed! Path: " <<  outputFilenameSoma << std::endl;
		std::cout << "Error! Writing cells.txt failed! Path: " <<  outputFilenameCells << std::endl;
	}
}

// Change Label in SpatialGraph (Edge and Vertex)
// for example change BasalDendrite to Dendrite
void helper::changeSpatialGraphLabel(AmiraSpatialGraph * spatialGraph,int oldLabel, int newLabel)
{
	// Change EdgeLabels
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = spatialGraph->edgesBegin(); edgeIt != spatialGraph->edgesEnd(); ++edgeIt)
	{
		if((*edgeIt)->label == oldLabel)
		{
			(*edgeIt)->label = newLabel;
		}
	}

	// Change VertexLabels
	std::vector< Vertex * >::iterator vertexIt;
	for(vertexIt = spatialGraph->verticesBegin(); vertexIt != spatialGraph->verticesEnd(); ++vertexIt)
	{
		if((*vertexIt)->label == oldLabel)
		{
			(*vertexIt)->label = newLabel;
		}
	}
}

TransformPointerType helper::amiraToVTKTransform(double* amiraTransform)
{
	HomogeneousMatrixPointerType mat = HomogeneousMatrixPointerType::New();
	for(int i = 0; i < 4; ++i)
		for(int j = 0; j < 4; ++j)
		{
			mat->SetElement(j, i, amiraTransform[i*4+j]);
		}

	TransformPointerType vtkTransform = TransformPointerType::New();
	vtkTransform->SetMatrix(mat);
	return vtkTransform;
}

void helper::writeListToFile(std::list<double> valList,const char * filename)
{
	std::ofstream TXTWriter;
	TXTWriter.open(filename);
	if(!TXTWriter.fail())
	{
		for (std::list<double>::iterator it=valList.begin(); it != valList.end(); ++it)
			TXTWriter << *it << std::endl;
	}
	else
	{
		std::cout << "ERROR! Writing " << filename << " failed! " << std::endl;
	}
	TXTWriter.close();
}

/* Connection Matrix Helpers */
std::vector <int> helper::getPostCelltypeList()
{
	std::vector <int> PostCelltypeList;
	PostCelltypeList.push_back(L2);
	PostCelltypeList.push_back(L34);
	PostCelltypeList.push_back(L4py);
	PostCelltypeList.push_back(L4sp);
	PostCelltypeList.push_back(L4ss);
	PostCelltypeList.push_back(L5st);
	PostCelltypeList.push_back(L5tt);
	PostCelltypeList.push_back(L6cc);
	PostCelltypeList.push_back(L6ccinv);
	PostCelltypeList.push_back(L6ct);
	PostCelltypeList.push_back(SymLocal1);
	PostCelltypeList.push_back(SymLocal2);
	PostCelltypeList.push_back(SymLocal3);
	PostCelltypeList.push_back(SymLocal4);
	PostCelltypeList.push_back(SymLocal5);
	PostCelltypeList.push_back(SymLocal6);
	PostCelltypeList.push_back(L1);
	PostCelltypeList.push_back(L23Trans);
	PostCelltypeList.push_back(L45Sym);
	PostCelltypeList.push_back(L45Peak);
	PostCelltypeList.push_back(L56Trans);
	return PostCelltypeList;
}

std::vector <int> helper::getPreCelltypeList()
{
	std::vector <int> PreCelltypeList;
	PreCelltypeList.push_back(VPM);
	PreCelltypeList.push_back(L2axon);
	PreCelltypeList.push_back(L34axon);
	PreCelltypeList.push_back(L4pyaxon);
	PreCelltypeList.push_back(L4spaxon);
	PreCelltypeList.push_back(L4ssaxon);
	PreCelltypeList.push_back(L5staxon);
	PreCelltypeList.push_back(L5ttaxon);
	PreCelltypeList.push_back(L6ccaxon);
	PreCelltypeList.push_back(L6ccinvaxon);
	PreCelltypeList.push_back(L6ctaxon);
	PreCelltypeList.push_back(SymLocal1axon);
	PreCelltypeList.push_back(SymLocal2axon);
	PreCelltypeList.push_back(SymLocal3axon);
	PreCelltypeList.push_back(SymLocal4axon);
	PreCelltypeList.push_back(SymLocal5axon);
	PreCelltypeList.push_back(SymLocal6axon);
	PreCelltypeList.push_back(L1axon);
	PreCelltypeList.push_back(L23Transaxon);
	PreCelltypeList.push_back(L45Symaxon);
	PreCelltypeList.push_back(L45Peakaxon);
	PreCelltypeList.push_back(L56Transaxon);
	return PreCelltypeList;
}

std::vector <int> helper::getColumnList()
{
	std::vector <int> ColumnList;
	ColumnList.push_back(A1);
	ColumnList.push_back(A2);
	ColumnList.push_back(A3);
	ColumnList.push_back(A4);
	ColumnList.push_back(B1);
	ColumnList.push_back(B2);
	ColumnList.push_back(B3);
	ColumnList.push_back(B4);
	ColumnList.push_back(C1);
	ColumnList.push_back(C2);
	ColumnList.push_back(C3);
	ColumnList.push_back(C4);
	ColumnList.push_back(D1);
	ColumnList.push_back(D2);
	ColumnList.push_back(D3);
	ColumnList.push_back(D4);
	ColumnList.push_back(E1);
	ColumnList.push_back(E2);
	ColumnList.push_back(E3);
	ColumnList.push_back(E4);
	ColumnList.push_back(Alpha);
	ColumnList.push_back(Beta);
	ColumnList.push_back(Gamma);
	ColumnList.push_back(Delta);
	return ColumnList;
}

std::map< unsigned int, const char * > helper::getInt2CelltypeLabelsMap()
{
	std::map< unsigned int, const char * > int2CelltypeLabels;
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L2,"L2"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L34,"L34"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4py,"L4py"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4sp,"L4sp"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4ss,"L4ss"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5st,"L5st"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5tt,"L5tt"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6cc,"L6cc"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccinv,"L6ccinv"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ct,"L6ct"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal,"SymLocal"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal1,"SymLocal1"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal2,"SymLocal2"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal3,"SymLocal3"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal4,"SymLocal4"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal5,"SymLocal5"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal6,"SymLocal6"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L1,"L1"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L23Trans,"L23Trans"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Sym,"L45Sym"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Peak,"L45Peak"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L56Trans,"L56Trans"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L2axon,"L2axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L34axon,"L34axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4pyaxon,"L4pyaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4spaxon,"L4spaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L4ssaxon,"L4ssaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5staxon,"L5staxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L5ttaxon,"L5ttaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccaxon,"L6ccaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ccinvaxon,"L6ccinvaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L6ctaxon,"L6ctaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocalaxon,"SymLocalaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal1axon,"SymLocal1axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal2axon,"SymLocal2axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal3axon,"SymLocal3axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal4axon,"SymLocal4axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal5axon,"SymLocal5axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(SymLocal6axon,"SymLocal6axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L1axon,"L1axon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L23Transaxon,"L23Transaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Symaxon,"L45Symaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L45Peakaxon,"L45Peakaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(L56Transaxon,"L56Transaxon"));
	int2CelltypeLabels.insert(std::pair< unsigned int, const char * >(VPM,"VPM"));
	return int2CelltypeLabels;
}

std::map< unsigned int, const char * > helper::getInt2ColumnLabelsMap()
{
	std::map< unsigned int, const char * > int2ColumnLabels;
	int2ColumnLabels.insert(std::pair< int, const char * >(Alpha, "Alpha"));
	int2ColumnLabels.insert(std::pair< int, const char * >(A1, "A1"));
	int2ColumnLabels.insert(std::pair< int, const char * >(A2, "A2"));
	int2ColumnLabels.insert(std::pair< int, const char * >(A3, "A3"));
	int2ColumnLabels.insert(std::pair< int, const char * >(A4, "A4"));
	int2ColumnLabels.insert(std::pair< int, const char * >(Beta, "Beta"));
	int2ColumnLabels.insert(std::pair< int, const char * >(B1, "B1"));
	int2ColumnLabels.insert(std::pair< int, const char * >(B2, "B2"));
	int2ColumnLabels.insert(std::pair< int, const char * >(B3, "B3"));
	int2ColumnLabels.insert(std::pair< int, const char * >(B4, "B4"));
	int2ColumnLabels.insert(std::pair< int, const char * >(Gamma, "Gamma"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C1, "C1"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C2, "C2"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C3, "C3"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C4, "C4"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C5, "C5"));
	int2ColumnLabels.insert(std::pair< int, const char * >(C6, "C6"));
	int2ColumnLabels.insert(std::pair< int, const char * >(Delta, "Delta"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D1, "D1"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D2, "D2"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D3, "D3"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D4, "D4"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D5, "D5"));
	int2ColumnLabels.insert(std::pair< int, const char * >(D6, "D6"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E1, "E1"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E2, "E2"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E3, "E3"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E4, "E4"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E5, "E5"));
	int2ColumnLabels.insert(std::pair< int, const char * >(E6, "E6"));
	return int2ColumnLabels;
}

std::map< std::string, unsigned int> helper::getCelltypeLabels2IntMap()
{
	std::map< std::string, unsigned int> celltypeLabels2Int;
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L2"),L2));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L34"),L34));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4py"),L4py));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4sp"),L4sp));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4ss"),L4ss));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5st"),L5st));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5tt"),L5tt));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6cc"),L6cc));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccinv"),L6ccinv));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ct"),L6ct));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal"),SymLocal));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal1"),SymLocal1));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal2"),SymLocal2));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal3"),SymLocal3));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal4"),SymLocal4));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal5"),SymLocal5));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal6"),SymLocal6));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L1"),L1));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L23Trans"),L23Trans));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Sym"),L45Sym));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Peak"),L45Peak));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L56Trans"),L56Trans));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L2axon"),L2axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L34axon"),L34axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4pyaxon"),L4pyaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4spaxon"),L4spaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L4ssaxon"),L4ssaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5staxon"),L5staxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L5ttaxon"),L5ttaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccaxon"),L6ccaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ccinvaxon"),L6ccinvaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L6ctaxon"),L6ctaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocalaxon"),SymLocalaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal1axon"),SymLocal1axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal2axon"),SymLocal2axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal3axon"),SymLocal3axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal4axon"),SymLocal4axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal5axon"),SymLocal5axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("SymLocal6axon"),SymLocal6axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L1axon"),L1axon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L23Transaxon"),L23Transaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Symaxon"),L45Symaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L45Peakaxon"),L45Peakaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("L56Transaxon"),L56Transaxon));
	celltypeLabels2Int.insert(std::pair< std::string, unsigned int >(std::string("VPM"),VPM));
	return celltypeLabels2Int;
}

bool helper::isAxonCellType(int cellTypeID)
{
	if (cellTypeID>=VPM && cellTypeID<=L6ctaxon)
		return true;

	if (cellTypeID>=SymLocalaxon && cellTypeID<=L56Transaxon)
		return true;

	return false;
}

std::map< std::string, unsigned int> helper::getColumnLabels2IntMap()
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
	return ColumnLabels2Int;
}

std::map< unsigned int, unsigned int > helper::getPretype2Posttype()
{
	std::map< unsigned int, unsigned int > pretype2posttype;
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L2axon,L2));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L34axon,L34));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L4pyaxon,L4py));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L4spaxon,L4sp));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L4ssaxon,L4ss));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L5staxon,L5st));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L5ttaxon,L5tt));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L6ccaxon,L6cc));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L6ccinvaxon,L6ccinv));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L6ctaxon,L6ct));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(SymLocalaxon,SymLocal));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(SymLocal1axon,SymLocal1));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(SymLocal2axon,SymLocal2));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(SymLocal3axon,SymLocal3));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(SymLocal4axon,SymLocal4));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(SymLocal5axon,SymLocal5));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(SymLocal6axon,SymLocal6));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L1axon,L1));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L23Transaxon,L23Trans));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L45Symaxon,L45Sym));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L45Peakaxon,L45Peak));
	pretype2posttype.insert(std::pair< unsigned int, unsigned int >(L56Transaxon,L56Trans));
	return pretype2posttype;
}

std::map< unsigned int, unsigned int > helper::getPosttype2Pretype()
{
	std::map< unsigned int, unsigned int > posttype2pretype;
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L2,L2axon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L34,L34axon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L4py,L4pyaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L4sp,L4spaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L4ss,L4ssaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L5st,L5staxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L5tt,L5ttaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L6cc,L6ccaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L6ccinv,L6ccinvaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L6ct,L6ctaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(SymLocal,SymLocalaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(SymLocal1,SymLocal1axon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(SymLocal2,SymLocal2axon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(SymLocal3,SymLocal3axon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(SymLocal4,SymLocal4axon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(SymLocal5,SymLocal5axon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(SymLocal6,SymLocal6axon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L1,L1axon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L23Trans,L23Transaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L45Sym,L45Symaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L45Peak,L45Peakaxon));
	posttype2pretype.insert(std::pair< unsigned int, unsigned int >(L56Trans,L56Transaxon));
	return posttype2pretype;
}

std::map< unsigned int, std::map<int, double> > helper::getCellType2BoutonDensityMap()
{
	std::map<int, double> boutonDensities;
	std::map< unsigned int, std::map<int, double>  > cellType2BoutonDensityMap;

	boutonDensities[INFRA] = 0.28;
	boutonDensities[GRAN] = 0.31;
	boutonDensities[SUPRA] = 0.34;
	cellType2BoutonDensityMap[VPM] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.31;
	boutonDensities[GRAN] = 0.31;
	boutonDensities[SUPRA] = 0.36;
	cellType2BoutonDensityMap[L2axon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.23;
	boutonDensities[GRAN] = 0.25;
	boutonDensities[SUPRA] = 0.25;
	cellType2BoutonDensityMap[L34axon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.18;
	boutonDensities[GRAN] = 0.25;
	boutonDensities[SUPRA] = 0.22;
	cellType2BoutonDensityMap[L4pyaxon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.21;
	boutonDensities[GRAN] = 0.28;
	boutonDensities[SUPRA] = 0.24;
	cellType2BoutonDensityMap[L4spaxon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.24;
	boutonDensities[GRAN] = 0.27;
	boutonDensities[SUPRA] = 0.25;
	cellType2BoutonDensityMap[L4ssaxon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.19;
	boutonDensities[GRAN] = 0.24;
	boutonDensities[SUPRA] = 0.28;
	cellType2BoutonDensityMap[L5staxon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.20;
	boutonDensities[GRAN] = 0.25;
	boutonDensities[SUPRA] = 0.19;
	cellType2BoutonDensityMap[L5ttaxon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.20;
	boutonDensities[GRAN] = 0.27;
	boutonDensities[SUPRA] = 0.26;
	cellType2BoutonDensityMap[L6ccaxon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.29;
	boutonDensities[GRAN] = 0.26;
	boutonDensities[SUPRA] = 0.23;
	cellType2BoutonDensityMap[L6ccinvaxon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.27;
	boutonDensities[GRAN] = 0.27;
	boutonDensities[SUPRA] = 0.27;
	cellType2BoutonDensityMap[L6ctaxon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.20;
	boutonDensities[GRAN] = 0.20;
	boutonDensities[SUPRA] = 0.20;
	cellType2BoutonDensityMap[SymLocal1axon] = boutonDensities;
	cellType2BoutonDensityMap[SymLocal2axon] = boutonDensities;
	cellType2BoutonDensityMap[SymLocal3axon] = boutonDensities;
	cellType2BoutonDensityMap[SymLocal4axon] = boutonDensities;
	cellType2BoutonDensityMap[SymLocal5axon] = boutonDensities;
	cellType2BoutonDensityMap[SymLocal6axon] = boutonDensities;
	cellType2BoutonDensityMap[SymLocalaxon] = boutonDensities;

	boutonDensities.clear();
	boutonDensities[INFRA] = 0.48;
	boutonDensities[GRAN] = 0.48;
	boutonDensities[SUPRA] = 0.48;
	cellType2BoutonDensityMap[L1axon] = boutonDensities;

	return cellType2BoutonDensityMap;
}

std::map< unsigned int, std::map<int, double> > helper::getCellType2SpineDensityMap()
{
	std::map<int, double> spineDensities;
	std::map< unsigned int, std::map<int, double> > cellType2SpineDensityMap;

	// Apical, Basal
	spineDensities[ApicalDendrite] = 1.68;
	spineDensities[BasalDendrite] = 1.68;
	cellType2SpineDensityMap[L2] = spineDensities;

	spineDensities.clear();
	spineDensities[ApicalDendrite] = 1.68;
	spineDensities[BasalDendrite] = 1.68;
	cellType2SpineDensityMap[L34] = spineDensities;

	spineDensities.clear();
	spineDensities[ApicalDendrite] = 1.68;
	spineDensities[BasalDendrite] = 1.17;
	cellType2SpineDensityMap[L4py] = spineDensities;

	spineDensities.clear();
	spineDensities[ApicalDendrite] = 1.17;
	spineDensities[BasalDendrite] = 1.17;
	cellType2SpineDensityMap[L4sp] = spineDensities;

	spineDensities.clear();
	spineDensities[ApicalDendrite] = 1.17;
	spineDensities[BasalDendrite] = 1.17;
	cellType2SpineDensityMap[L4ss] = spineDensities;

	spineDensities.clear();
	spineDensities[ApicalDendrite] = 1.68;
	spineDensities[BasalDendrite] = 1.04;
	cellType2SpineDensityMap[L5st] = spineDensities;

	spineDensities.clear();
	spineDensities[ApicalDendrite] = 1.68;
	spineDensities[BasalDendrite] = 1.04;
	cellType2SpineDensityMap[L5tt] = spineDensities;

	spineDensities.clear();
	spineDensities[ApicalDendrite] = 1.04;
	spineDensities[BasalDendrite] = 1.04;
	cellType2SpineDensityMap[L6cc] = spineDensities;

	spineDensities.clear();
	spineDensities[ApicalDendrite] = 1.04;
	spineDensities[BasalDendrite] = 1.04;
	cellType2SpineDensityMap[L6ct] = spineDensities;

	spineDensities.clear();
	spineDensities[ApicalDendrite] = 1.04;
	spineDensities[BasalDendrite] = 1.04;
	cellType2SpineDensityMap[L6ccinv] = spineDensities;

	return cellType2SpineDensityMap;
}

void helper::getTypeIDs(unsigned int CellType, unsigned int& preCellType, unsigned int& postCellType)
{
	std::map< unsigned int, unsigned int > pretype2posttype = getPretype2Posttype();
	std::map< unsigned int, unsigned int > posttype2pretype = getPosttype2Pretype();

	if(pretype2posttype.count(CellType)>0)
	{
		preCellType = CellType;
		postCellType = pretype2posttype[CellType];
	}
	if(posttype2pretype.count(CellType)>0)
	{
		postCellType = CellType;
		preCellType = posttype2pretype[CellType];
	}
	// VPM?
	if (CellType==VPM)
	{
		preCellType = VPM;
		postCellType = 0;
	}
}
