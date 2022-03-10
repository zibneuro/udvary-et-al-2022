/* alignSpatialGraphToSoma.cpp
 * Aligns SpatialGraph to Soma (0,0,0) in desired dimensions
 * 	saves:
 * 		- aligned Spatial Graph
 *	USAGE:
 *	 ./AlignSpatialGraphToSoma
 *	 	[path/to/SpatialGraphFile.am/.hoc]
 *	 	[outputpath/to/AlignedSpatialGraphFile.am]
 *	 	[AlignInX (0 or 1)] [AlignInY (0 or 1)] [AlignInZ (0 or 1)]
 *
 * --------------------------
 * Author: Daniel Udvary
 * 		   In Silico Brain Sciences										
 *		   Max Planck Institute for Neurobiology of Behavior – caesar		
 *		   Ludwig-Erhard-Allee 2											
 *		   53175 Bonn, Germany
 * e-mail: 	daniel.udvary@mpinb.mpg.de
 * Date: 27 February 2020
 * Version 0
 * --------------------------
 */

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "typedefs.h"
#include "basics.h"
#include "amiraReader.h"
#include "helper.h"
#include <set>
#include <utility>

void printUsage();

int main( int argc , char * argv[])
{
	if (argc==6)
	{
		// Load in Amira Spatial Graph
		const char * inputfilename = argv[1];
		const char * outputfilename = argv[2];

		// Load Spatial Graph
	  	Reader * fileReader = new Reader(inputfilename, outputfilename);

		std::string fname(inputfilename);
		if((access(fname.c_str(), F_OK) != -1))
		{
			if (fname.compare(fname.size()-3,3,".am") == 0)
			{
				fileReader->readSpatialGraphFile(0);
			}
			else if (fname.compare(fname.size()-4,4,".hoc") == 0)
			{
				fileReader->readHocFile();
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

		AmiraSpatialGraph * spatialGraph = fileReader->getSpatialGraph();

		// Initialize boolean
		int xtmp = atoi(argv[3]);
		int ytmp = atoi(argv[4]);
		int ztmp = atoi(argv[5]);
		bool align[3] = {1,1,1};

		// If any of the three inputs is equal to zero, set boolean to false, otherwise its true
		// True (>0): Aligned in that dimension
		// False (==0): Not aligned in that dimension
		if (xtmp==0)
			align[0] = 0;
		if (ytmp==0)
			align[1] = 0;
		if (ztmp==0)
			align[2] = 0;

		//std::cout << align[0] << " " << align[1] << " " << align[2] << std::endl;

		// Get 3D Soma Position
		double centerPt[3];
		helper::centerOfSpatialGraph(Soma, centerPt, spatialGraph);
		std::cout << inputfilename << "," << centerPt[0] << "," << centerPt[1] << "," << centerPt[2] << std::endl;

		// Subtract centerPt from spatialGraph
		for(int i=0; i<3; i++)
		{
		    if (align[i])
			    centerPt[i] *= -1;
		    else
			    centerPt[i] = 0;
		}

		TransformPointerType sub = TransformPointerType::New();
		sub->Translate(centerPt);
		spatialGraph->setTransformation(sub);
		spatialGraph->applyTransformation();

		// Write SpatialGraphfile
		fileReader->setSpatialGraph(spatialGraph);
		fileReader->writeSpatialGraphFile();
	}
	else
	{
		printUsage();
		return 0;
	}
}

void printUsage()
{
	std::cout << "ERROR! Wrong number of input arguments!" << std::endl;
	std::cout << "USAGE:" << std::endl;
	std::cout << " ./AlignSpatialGraphToSoma [path/to/SpatialGraphFile.am/.hoc] [outputpath/to/AlignedSpatialGraphFile.am]";
	std::cout << " [AlignInX (0 or 1)] [AlignInY (0 or 1)] [AlignInZ (0 or 1)]" << std::endl;
}
