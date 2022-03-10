/* Generates 3D Amira Density File from SpatialGraphs
 * - Density files of path length per subvolume
 * - Density files of branch points per subvolume
 * - Density files of terminal points per subvolume
 * 
 * USAGE: ./GenerateDensityFiles [Inputfilename] [Outputfilename] 
 * 			[SubcellularCompartment] [Parameter] 
 * 			[SubcellularCompartment: Dendrite, Axon, BasalDendrite, Basal, ApicalDendrite, Apical]
 * 			[Parameter: Length, BP, TP]
 *
 * --------------------------
 * Author: Daniel Udvary
 * 		   In Silico Brain Sciences										
 *		   Max Planck Institute for Neurobiology of Behavior – caesar		
 *		   Ludwig-Erhard-Allee 2											
 *		   53175 Bonn, Germany
 * e-mail: 	daniel.udvary@mpinb.mpg.de
 * Date: 12 May 2021
 * --------------------------
 */
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "typedefs.h"
#include "basics.h"
#include "amiraReader.h"
#include "helper.h"
#include "densitySG.h"
#include <set>
#include <utility>

int main(int argc, const char** argv)
{
	if(argc == 5)
	{
		const char * inputfilename = argv[1];
		const char * outputfilename = argv[2];
		std::string compartment = argv[3];
		std::string parameter = argv[4];

		// Select SubcellularCompartment
		int compartmentID;
		if ((compartment.compare("Dendrite") == 0)
			|| (compartment.compare("dendrite") == 0))
			compartmentID = Dendrite;
		else if ((compartment.compare("Axon") == 0)
				|| (compartment.compare("axon") == 0))
			compartmentID = Axon;
		else if ((compartment.compare("BasalDendrite") == 0)
					|| (compartment.compare("basaldendrite") == 0)
					|| (compartment.compare("Basal") == 0)
					|| (compartment.compare("basal") == 0))
			compartmentID = BasalDendrite;
		else if ((compartment.compare("ApicalDendrite") == 0)
					|| (compartment.compare("apicaldendrite") == 0)
					|| (compartment.compare("Apical") == 0)
					|| (compartment.compare("apical") == 0))
			compartmentID = ApicalDendrite;
		else
		{
			std::cout << "Invalid Subcellular Compartment. Use Dendrite, Axon, BasalDendrite, Basal, ApicalDendrite, or Apical" << std::endl;
			return 0;
		}

		densitySG * dSg1 = new densitySG(inputfilename,compartmentID);
		ImageDataPointerType density1;

		if ((parameter.compare("Length") == 0) || (parameter.compare("length") == 0))
			density1 = dSg1->computeLength();
		else if ((parameter.compare("BP") == 0) || (parameter.compare("bp") == 0)
				|| (parameter.compare("branchpoints") == 0))
			density1 = dSg1->computeBP();
		else if ((parameter.compare("TP") == 0) || (parameter.compare("tp") == 0)
				|| (parameter.compare("terminalpoints") == 0))
			density1 = dSg1->computeTP();
		else
		{
			std::cout << "Invalid Parameter. Use Length, BP [BranchPoints], or TP [TerminalPoints]" << std::endl;
			return 0;
		}

		helper::storeImageVolume(outputfilename,density1);
	}
	else
	{
		std::cout << "ERROR! Wrong number of input arguments!" << std::endl;
		std::cout << "USAGE: ./GenerateDensityFiles [Inputfilename] [Outputfilename] [SubcellularCompartment] [Parameter]" << std::endl;
		std::cout << "       [SubcellularCompartment: Dendrite, Axon, BasalDendrite, Basal, ApicalDendrite, Apical]" << std::endl;
		std::cout << "       [Parameter: Length, BP, TP]" << std::endl;
		return 0;
	}
}
