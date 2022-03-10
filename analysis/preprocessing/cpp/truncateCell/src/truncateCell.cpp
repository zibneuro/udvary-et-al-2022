/* Truncates Cell Morphology / Emulate Slicing Experiment
 * Takes as an input a spatial graph cell morphology, the bounds of the slice and the orientation of the slice
 * (along x, y, or z axis)
 * Outputs truncated morphology as spatial graph
 * Outputs cell properties such as apical/basal/axon length
 *./TruncateCell [path/to/inputSpatialGraph.am/.hoc]
 *				 [minimum position of slice]
 *				 [maximum position of slice]
 *				 [coordinate of slice (X:0 Y:1 Z:2)]
 *				 [path/to/outputSpatialGraph.am]
 *
 * --------------------------
 * Author: Daniel Udvary
 * 		   In Silico Brain Sciences										
 *		   Max Planck Institute for Neurobiology of Behavior – caesar		
 *		   Ludwig-Erhard-Allee 2											
 *		   53175 Bonn, Germany
 * e-mail: 	daniel.udvary@mpinb.mpg.de
 * Date: 02 August 2018
 * Version 1
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

int main(int argc, const char** argv)
{
	if (argc == 6)
	{
		// Inputs
		const char * inputfilename = argv[1];
		double minLateral = atof(argv[2]);
		double maxLateral = atof(argv[3]);
		int cut_COORD = atoi(argv[4]);
		const char * outputfilename = argv[5]; // "/nas1/Data_daniel/Network/L5/L5_InVitro/L5ttDendrites/"

		// Coordinate Mapper for nicer output
		std::map<int,std::string> coordMapper;
		coordMapper[X_COORD] = "X";
		coordMapper[Y_COORD] = "Y";
		coordMapper[Z_COORD] = "Z";

		// Read in Spatial Graph
		AmiraSpatialGraph * sg = helper::getSpatialGraph(inputfilename);

		// Extract length parameters
		double lenBasalDendrite = helper::lengthCell(sg, BasalDendrite);
		double lenApicalDendrite = helper::lengthCell(sg, ApicalDendrite);
		double lenDendrite = helper::lengthCell(sg, Dendrite);
		double lenAxon = helper::lengthCell(sg, Axon);
		double somaPosition[3];
		helper::centerOfSpatialGraph(Soma, somaPosition, sg);

		// Truncate Spatial Graph
		double bounds[2] = {minLateral, maxLateral};
		sg->truncateSpatialGraph(bounds,cut_COORD);

		// Calculate tissue depth / distance between soma and slicing surface
		double distSomaSlice[2] = {somaPosition[cut_COORD]-minLateral,
									maxLateral-somaPosition[cut_COORD]};

		// Save truncated Spatial Graph
		Reader * amWriterCut = new Reader(outputfilename, outputfilename);
		amWriterCut->setSpatialGraph(sg);
		amWriterCut->writeSpatialGraphFile();

		// Compute Length of Truncated SpatialGraph and Display invivo and truncated length values, and parameters
		double lenBasalDendriteCut = helper::lengthCell(sg, BasalDendrite);
		double lenApicalDendriteCut = helper::lengthCell(sg, ApicalDendrite);
		double lenDendriteCut = helper::lengthCell(sg, Dendrite);
		double lenAxonCut = helper::lengthCell(sg, Axon);

		// Write on terminal
		std::cout << "************************" << std::endl;
		std::cout << "Truncation succeeded" << std::endl;
		std::cout << "Truncated morphology saved as " << outputfilename << std::endl;
		std::cout << "--------------------" << std::endl;
		std::cout << "Parameters of input spatial graph" << std::endl;
		std::cout << " Original length of basal dendrite = " << lenBasalDendrite << std::endl;
		std::cout << " Original length of apical dendrite = " << lenApicalDendrite << std::endl;
		std::cout << " Original length of dendrite = " << lenDendrite << std::endl;
		std::cout << " Original length of axon = " << lenAxon << std::endl;
		std::cout << " Soma position = [" << somaPosition[0] << ", ";
		std::cout << somaPosition[1] << ", " << somaPosition[2] << "]" << std::endl;
		std::cout << "--------------------" << std::endl;
		std::cout << "Parameters of truncated spatial graph" << std::endl;
		std::cout << " Truncated length of basal dendrite = " << lenBasalDendriteCut << std::endl;
		std::cout << " Truncated length of apical dendrite = " << lenApicalDendriteCut << std::endl;
		std::cout << " Truncated length of dendrite = " << lenDendriteCut << std::endl;
		std::cout << " Truncated length of axon = " << lenAxonCut << std::endl;
		std::cout << " Distance between soma and slicing surface (tissue depth) = [";
		std::cout << distSomaSlice[0] << " " << distSomaSlice[1] << "]" << std::endl;
		std::cout << "--------------------" << std::endl;
		std::cout << "Slice configuration" << std::endl;
		std::cout << " Slice position = [" << minLateral << " to " << maxLateral << "]" << std::endl;
		std::cout << " Slice width = " << maxLateral-minLateral << std::endl;
		std::cout << " Slicing axis = " << coordMapper[cut_COORD] << std::endl;
		std::cout << "************************" << std::endl;
	}
	else
	{
		std::cout << "Wrong number of input arguments!" << std::endl;
		printUsage();
	}
	return 0;
}

void printUsage()
{
	std::cout << "USAGE:" << std::endl;
	std::cout << "./TruncateCell";
	std::cout << " [path/to/inputSpatialGraph.am/.hoc]";
	std::cout << " [minimum position of slice]";
	std::cout << " [maximum position of slice]";
	std::cout << " [coordinate of slice (X:0 Y:1 Z:2)]";
	std::cout << " [path/to/outputSpatialGraph.am]" << std::endl;
}
