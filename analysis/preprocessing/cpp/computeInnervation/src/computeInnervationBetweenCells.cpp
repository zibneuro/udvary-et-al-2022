/* ComputeInnervationBetweenCells.cpp
 * Compute Innervation Between Cell Pairs
 * 	saves:
 * 		- bouton density (*_boutonDensity.am)
 * 		- spine density (*_spineDensity.am)
 * 		- innervation density (InnervatonDensity.am)
 * 		- results.txt containing the total number of boutons and spines,
 * 				innervation value, the connection probability
 * 				and mean number of contacts between the cell pair and their distribution
 * USAGE: ./ComputeInnervationBetweenCells [InputfilenamePresynapticCell] [PresynapticCellType]
 * 			[InputfilenamePostsynapticCell] [PostsynapticCellType] [outputpath]
 * 		NOTES:
 * 			Labeling: ApicalDendrite as apical; BasalDendrites as dend; Axon as axon;
 * 			presynaptic cell types: VPM, L2axon, L34axon, ..., L6ccaxon.
 * 			postsynaptic cell types: L2, L34, ..., L6cc.
 *
 * Option 2: Only save Spine or Bouton Density
 * USAGE: ./ComputeInnervationBetweenCells [Option: Spine or Bouton]
 * 			[Inputfilename] [CellType] [outputpath]
 *
 * --------------------------
 * Author: Daniel Udvary
 * 		   In Silico Brain Sciences										
 *		   Max Planck Institute for Neurobiology of Behavior – caesar		
 *		   Ludwig-Erhard-Allee 2											
 *		   53175 Bonn, Germany
 * e-mail: 	daniel.udvary@mpinb.mpg.de
 * Date: 25 March 2021
 * Version 0.2
 * --------------------------
 * 25 March 2021:
 *  - Added option to change path to EXC_POST_ALL.am
 * 	- refined error message when PST is not defined
 * 07 August 2020:
 *  - Added option to only save Spine or Bouton density
 * --------------------------
 */

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "typedefs.h"
#include "basics.h"
#include "amiraReader.h"
#include "barrel_field.h"
#include "helper.h"
#include "densitySG.h"
#include <set>
#include <utility>

ImageDataPointerType computeBoutonDensity(const char * inputfilename, int celltype);
ImageDataPointerType computeSpineDensity(const char * inputfilename, int celltype);
ImageDataPointerType computeInnervationDensity(ImageDataPointerType boutons, ImageDataPointerType spines);
ImageDataPointerType mergeSpineDensities(ImageDataPointerType lenBasal,double basalSpineDensity,
										ImageDataPointerType lenApical,double apicalSpineDensity,
										int foundNeurites);
ImageDataPointerType weigthDensity(ImageDataPointerType density, double scaleFactor);
ImageDataPointerType weigthDensity(ImageDataPointerType density, std::map<int, double> SGIBoutonFactor);
ImageDataPointerType mergeDensities(ImageDataPointerType density1, ImageDataPointerType density2);

double computeSumOfDensity(ImageDataPointerType density);
void getVoxelBoundingBox(double bounds[6],double origin[3],int extent[6],double VoxelSize[3]);
bool isPtInDensity(double voxelCenter[3], ImageDataPointerType density, int pos[3]);
std::vector<double> computeNumberOfContacts(double innervation);
double computeMeanNumberOfContacts(std::vector<double> ContactHist);
void printUsage();

const char * outputfilePSTTotalEXC;

// NeuralNet: Voxels go from [-50 to 0; 0 to 50] and so, BoundingBox is Center [-25 25]
/* 	The origin is BBmin, therefore the extent is positive!
	Bounds:
	  Xmin,Xmax: (-3975, 2925)
	  Ymin,Ymax: (-1725, 4225)
	  Zmin,Zmax: (-2725, 1125)
	Spacing: (50, 50, 50)
	Origin: (-3975, -1725, -2725)
	Dimensions: (139, 120, 78)
	Increments: (1, 0, 0)
	Extent: (0, 138, 0, 119, 0, 77)
*/
int main(int argc, const char** argv)
{
	std::string strVersion = "v0.2 (Date: 25 March 2021)";

	/*  COMPUTE INNERVATION BETWEEN CELL PAIR */
	if(argc == 6 || argc == 7)
	{
		const char * inputfilenameCell1 = argv[1];
		std::string celltypeCell1 = argv[2];
		const char * inputfilenameCell2 = argv[3];
		std::string celltypeCell2 = argv[4];
		const char * outputpath = argv[5];

		if (argc == 7)
			outputfilePSTTotalEXC = argv[6];
		else
			outputfilePSTTotalEXC =
				"/nas1/Data_regger/AXON_SAGA/Axon2/NeuroNet/cache_Vincent_complete_final/Cache/EXToPST.am";

		std::cout << "**********************************************" << std::endl;
		std::cout << "*                                            *" << std::endl;
		std::cout << "*  COMPUTING INNERVATION BETWEEN CELL PAIRS  *" << std::endl;
		std::cout << "*     " << strVersion << "            *" << std::endl;
		std::cout << "*                                            *" << std::endl;
		std::cout << "**********************************************" << std::endl;

		// Get filename only (delete path and file extension)
		std::string inputfilenameCell1Str = helper::getFilenameFromPath(inputfilenameCell1);
		size_t lastindex = inputfilenameCell1Str.find_last_of(".");
		inputfilenameCell1Str = inputfilenameCell1Str.substr(0, lastindex);

		std::string inputfilenameCell2Str = helper::getFilenameFromPath(inputfilenameCell2);
		lastindex = inputfilenameCell2Str.find_last_of(".");
		inputfilenameCell2Str = inputfilenameCell2Str.substr(0, lastindex);

		std::string outputpathStr(outputpath);

		// Get CellTypeIDs
		std::map< std::string, unsigned int> CelltypeLabels2IntMap = helper::getCelltypeLabels2IntMap();
		if (CelltypeLabels2IntMap.count(celltypeCell1)==0)
		{
			std::cout << "ERROR! CellType " << celltypeCell1 << " is unknown!" << std::endl;
			return NULL;
		}
		int celltypeCell1ID = CelltypeLabels2IntMap[celltypeCell1];

		if (CelltypeLabels2IntMap.count(celltypeCell2)==0)
		{
			std::cout << "ERROR! CellType " << celltypeCell2 << " is unknown!" << std::endl;
			return NULL;
		}
		int celltypeCell2ID = CelltypeLabels2IntMap[celltypeCell2];

		/* Compute Bouton Density */
		std::cout << "----------------------------------------------" << std::endl;
		std::cout << "COMPUTING BOUTON DENSITY" << std::endl;
		std::cout << "----------------------------------------------" << std::endl;
		ImageDataPointerType boutonDensity = computeBoutonDensity(inputfilenameCell1,celltypeCell1ID);
		std::string outputFilenameBoutons = outputpathStr + inputfilenameCell1Str + "_boutonDensity.am";
		helper::storeImageVolume(outputFilenameBoutons.c_str(), boutonDensity);
		double numBoutons = computeSumOfDensity(boutonDensity);

		/* Compute Spine Density */
		std::cout << "----------------------------------------------" << std::endl;
		std::cout << "COMPUTING SPINE DENSITY" << std::endl;
		std::cout << "----------------------------------------------" << std::endl;
		ImageDataPointerType spineDensity = computeSpineDensity(inputfilenameCell2,celltypeCell2ID);
		std::string outputFilenameSpines = outputpathStr + inputfilenameCell2Str + "_spineDensity.am";
		helper::storeImageVolume(outputFilenameSpines.c_str(), spineDensity);
		double numSpines = computeSumOfDensity(spineDensity);

		/* Compute Innervation Density */
		std::cout << "----------------------------------------------" << std::endl;
		std::cout << "COMPUTING INNERVATION DENSITY" << std::endl;
		std::cout << "----------------------------------------------" << std::endl;
		ImageDataPointerType InnervationDensity = computeInnervationDensity(boutonDensity, spineDensity);
		std::string outputFilenameInnervation = outputpathStr + "InnervationDensity.am";
		helper::storeImageVolume(outputFilenameInnervation.c_str(), InnervationDensity);
		double I = computeSumOfDensity(InnervationDensity);

		// Results
		std::vector<double> ContactHist = computeNumberOfContacts(I);
		double meanNumContacts = computeMeanNumberOfContacts(ContactHist);

		// Write Results in txt file
		std::ofstream TXTWriter;
		std::string filenameTXT = outputpathStr + "results.txt";
		TXTWriter.open(filenameTXT.c_str());
		if(!TXTWriter.fail())
		{
			// Get Spine Densities
			std::map< unsigned int, std::map<int, double> > SpineDensityMap = helper::getCellType2SpineDensityMap();
			std::map<int, double> spineDensityCellType = SpineDensityMap[celltypeCell2ID];
			// Get Bouton Densities
			std::map< unsigned int, std::map<int, double> > BoutonDensityMap = helper::getCellType2BoutonDensityMap();
			std::map<int, double> SGIBoutonFactor = BoutonDensityMap[celltypeCell1ID];

			TXTWriter << "-----------------------------------------" << std::endl;
			TXTWriter << "computeInnervationBetweenCells.cpp" << std::endl;
			TXTWriter << " " << strVersion << std::endl;
			TXTWriter << "-----------------------------------------" << std::endl;
			TXTWriter << "PresynapticCell: " << inputfilenameCell1 << std::endl;
			TXTWriter << "PresynapticCellType: " << celltypeCell1 << std::endl;
			TXTWriter << "BoutonDensity[Supragranular] = " << SGIBoutonFactor[SUPRA] << std::endl;
			TXTWriter << "BoutonDensity[Granular] = " << SGIBoutonFactor[GRAN] << std::endl;
			TXTWriter << "BoutonDensity[Infragranular] = " << SGIBoutonFactor[INFRA] << std::endl;
			TXTWriter << "-----------------------------------------" << std::endl;
			TXTWriter << "PostsynapticCell: " << inputfilenameCell2 << std::endl;
			TXTWriter << "PostsynapticCellType: " << celltypeCell2 << std::endl;
			TXTWriter << "SpineDensity[Basal] = " << spineDensityCellType[BasalDendrite] << std::endl;
			TXTWriter << "SpineDensity[Apical] = " << spineDensityCellType[ApicalDendrite] << std::endl;
			TXTWriter << "-----------------------------------------" << std::endl;
			TXTWriter << "RESULTS" << std::endl;
			TXTWriter << " Boutons = " << numBoutons << std::endl;
			TXTWriter << " Spines = " << numSpines << std::endl;
			TXTWriter << " Innervation = " << I << std::endl;
			TXTWriter << " Connection Probability = " << 1-exp(-I) << std::endl;
			TXTWriter << " Mean Number of Contacts (if connected) = " << meanNumContacts << std::endl;
			TXTWriter << " Range of Contacts: [numContacts, probability]" << std::endl;
			for(int ii = 0; ii < ContactHist.size(); ++ii)
				TXTWriter << "       [" << ii << ", " << ContactHist.at(ii) << "]" << std::endl;
			TXTWriter << "-----------------------------------------" << std::endl;
		}
		else
		{
			std::cout << "ERROR! Writing " << filenameTXT << " failed! " << std::endl;
		}
		TXTWriter.close();

		std::cout << "----------------------------------------------" << std::endl;
		std::cout << "RESULTS [saved to " << outputpathStr << "]" << std::endl;
		std::cout << "----------------------------------------------" << std::endl;
		std::cout << ">> BOUTONS = " << numBoutons << std::endl;
		std::cout << ">> SPINES = " << numSpines << std::endl;
		std::cout << ">> INNERVATION = " << I << std::endl;
		std::cout << ">> CONNECTION PROBABILITY = " << 1-exp(-I) << std::endl;
		TXTWriter << ">> MEAN NUMBER OF CONTACTS (if connected) = " << meanNumContacts << std::endl;
		std::cout << ">> RANGE OF CONTACTS: [numContacts, probability]" << std::endl;
		for(int ii = 0; ii < ContactHist.size(); ++ii)
			std::cout << "       [" << ii << ", " << ContactHist.at(ii) << "]" << std::endl;
		std::cout << "----------------------------------------------" << std::endl;

	}
	/*  COMPUTE SPINE OR BOUTON DENSITY OF INDIVIDUAL CELL */
	else if (argc == 5)
	{
		std::string option = argv[1];
		const char * inputfilenameCell = argv[2];
		std::string celltypeCell = argv[3];
		const char * outputpath = argv[4];

		int opt;
		if ((option.compare("Spine") == 0)
			|| (option.compare("spine") == 0))
			opt = 1;
		else if ((option.compare("Bouton") == 0)
				|| (option.compare("bouton") == 0))
			opt = 0;
		else
		{
			std::cout << "Invalid Option. Use Bouton or Spine." << std::endl;
			return NULL;
		}

		// Get filename only (delete path and file extension)
		std::string inputfilenameCellStr = helper::getFilenameFromPath(inputfilenameCell);
		size_t lastindex = inputfilenameCellStr.find_last_of(".");
		inputfilenameCellStr = inputfilenameCellStr.substr(0, lastindex);
		std::string outputpathStr(outputpath);

		// Get CellTypeIDs
		std::map< std::string, unsigned int> CelltypeLabels2IntMap = helper::getCelltypeLabels2IntMap();
		if (CelltypeLabels2IntMap.count(celltypeCell)==0)
		{
			std::cout << "ERROR! CellType " << celltypeCell << " is unknown!" << std::endl;
			return NULL;
		}
		int celltypeCellID = CelltypeLabels2IntMap[celltypeCell];

		if (opt==1) // Spine
		{
			std::cout << "**********************************************" << std::endl;
			std::cout << "*                                            *" << std::endl;
			std::cout << "*  COMPUTING SPINE DENSITY OF CELL           *" << std::endl;
			std::cout << "*     " << strVersion << "             *" << std::endl;
			std::cout << "*                                            *" << std::endl;
			std::cout << "**********************************************" << std::endl;

			ImageDataPointerType spineDensity = computeSpineDensity(inputfilenameCell,celltypeCellID);
			std::string outputFilenameSpines = outputpathStr + inputfilenameCellStr + "_spineDensity.am";
			helper::storeImageVolume(outputFilenameSpines.c_str(), spineDensity);
			double numSpines = computeSumOfDensity(spineDensity);
			std::cout << ">> SPINES = " << numSpines << std::endl;

			if (numSpines==0)
			{
				std::cout << "WARNING! NO SPINES FOUND!" << std::endl;
			}
		}
		else // opt==0 // Bouton
		{
			std::cout << "**********************************************" << std::endl;
			std::cout << "*                                            *" << std::endl;
			std::cout << "*  COMPUTING BOUTON DENSITY OF CELL          *" << std::endl;
			std::cout << "*     " << strVersion << "             *" << std::endl;
			std::cout << "*                                            *" << std::endl;
			std::cout << "**********************************************" << std::endl;

			ImageDataPointerType boutonDensity = computeBoutonDensity(inputfilenameCell,celltypeCellID);
			std::string outputFilenameBoutons = outputpathStr + inputfilenameCellStr + "_boutonDensity.am";
			helper::storeImageVolume(outputFilenameBoutons.c_str(), boutonDensity);
			double numBoutons = computeSumOfDensity(boutonDensity);
			std::cout << ">> BOUTONS = " << numBoutons << std::endl;

			if (numBoutons==0)
			{
				std::cout << "WARNING! NO BOUTONS FOUND!" << std::endl;
			}
		}
	}
	else
	{
		printUsage();
	}
}

double computeMeanNumberOfContacts(std::vector<double> ContactHist)
{
	double meanNumContacts = 0;

	// If ContactHistogram has length 1, only NumContact=0 exists
	if (ContactHist.size()==1)
	{
		return meanNumContacts;
	}

	double norm = 0;
	for(int ii = 1; ii < ContactHist.size(); ++ii)
	{
		meanNumContacts += ii*ContactHist.at(ii);
		norm += ContactHist.at(ii);
	}

	return meanNumContacts/norm;
}

/* Compute Probability of Number of Contacts based on Poisson Distribution */
std::vector<double> computeNumberOfContacts(double innervation)
{
	double percentile = 0.99;
	std::vector<double> hist;
	double targetConnectionProb = 1 - exp(-1*innervation);
	double cumulativeConnectionProb = 0;
	unsigned int nSyn = 1;
	hist.push_back(1-targetConnectionProb);
	while(cumulativeConnectionProb < percentile*targetConnectionProb)
	{
		double binProb = pow(innervation, nSyn)/vtkMath::Factorial(nSyn)*exp(-1*innervation);
		hist.push_back(binProb);
		cumulativeConnectionProb += binProb;
		++nSyn;
	}
	return hist;
}

void printUsage()
{
	std::cout << "ERROR! Wrong number of input arguments!" << std::endl;
	std::cout << "USAGE 1: COMPUTE INNERVATION BETWEEN CELL PAIR" << std::endl;
	std::cout << "./ComputeInnervationBetweenCells [InputfilenamePresynapticCell] [PresynapticCellType]";
	std::cout << " [InputfilenamePostsynapticCell] [PostsynapticCellType] [outputpath] [optional: path/to/ExcPostAll.am]" << std::endl;
	std::cout << " NOTES: " << std::endl;
	std::cout << "        Labeling: ApicalDendrite as apical; BasalDendrites as dend; Axon as axon;" << std::endl;
	std::cout << "        presynaptic cell types: VPM, L2axon, L34axon, ..., L6ccaxon." << std::endl;
	std::cout << "        postsynaptic cell types: L2, L34, ..., L6cc." << std::endl;
	std::cout << " " << std::endl;
	std::cout << "USAGE 2: COMPUTE SPINE OR BOUTON DENSITY OF CELL" << std::endl;
	std::cout << "./ComputeInnervationBetweenCells [Option: Spine or Bouton] [Inputfilename] [CellType]";
	std::cout << " [outputpath]" << std::endl;
	std::cout << " NOTES: " << std::endl;
	std::cout << "        Labeling: ApicalDendrite as apical; BasalDendrites as dend; Axon as axon;" << std::endl;
	std::cout << "        presynaptic cell types: VPM, L2axon, L34axon, ..., L6ccaxon." << std::endl;
	std::cout << "        postsynaptic cell types: L2, L34, ..., L6cc." << std::endl;

}

double computeSumOfDensity(ImageDataPointerType density)
{
	double sum = 0;

	for(int x = density->GetExtent()[0]; x <= density->GetExtent()[1]; ++x)
	{
		for(int y = density->GetExtent()[2]; y <= density->GetExtent()[3]; ++y)
		{
			for(int z = density->GetExtent()[4]; z <= density->GetExtent()[5]; ++z)
			{
				double *px = static_cast< double * >(density->GetScalarPointer(x,y,z));
				sum += (*px);
			}
		}
	}
	return sum;
}

ImageDataPointerType computeBoutonDensity(const char * inputfilename, int celltype)
{
	densitySG * dSg = new densitySG(inputfilename,Axon);
	AmiraSpatialGraph * sg = dSg->getSpatialGraph();
	bool axonFound = sg->isLabelInSpatialGraph(Axon);
	ImageDataPointerType lenAxon;
	if (axonFound)
	{
		lenAxon = dSg->computeLength();
	}
	else
	{
		std::cout << "ERROR! No Axon found in SpatialGraph (" << inputfilename << ")" << std::endl;
		return lenAxon;
	}

	// Get Bouton Densities
	std::map< unsigned int, std::map<int, double> > BoutonDensityMap = helper::getCellType2BoutonDensityMap();
	std::map<int, double> SGIBoutonFactor = BoutonDensityMap[celltype];
	std::cout << ">> BoutonDensity used: [SGI] = [" << SGIBoutonFactor[SUPRA] << " ";
	std::cout << SGIBoutonFactor[GRAN] << " " << SGIBoutonFactor[INFRA] << "]" << std::endl;

	ImageDataPointerType boutons = weigthDensity(lenAxon, SGIBoutonFactor);
	return boutons;
}

ImageDataPointerType weigthDensity(ImageDataPointerType density, std::map<int, double> SGIBoutonFactor)
{
	BarrelField * SBF = new BarrelField;

	double VoxelSize[3];
	density->GetSpacing(VoxelSize);
	double origin[3];
	density->GetOrigin(origin);

	for(int x = density->GetExtent()[0]; x <= density->GetExtent()[1]; ++x)
	{
		for(int y = density->GetExtent()[2]; y <= density->GetExtent()[3]; ++y)
		{
			for(int z = density->GetExtent()[4]; z <= density->GetExtent()[5]; ++z)
			{
				// Compute laminar position (supra, granular or infragranular)
				double voxelCenter[3];
				voxelCenter[0] = origin[0] + x*VoxelSize[0];
				voxelCenter[1] = origin[1] + y*VoxelSize[1];
				voxelCenter[2] = origin[2] + z*VoxelSize[2];
				int laminarPos = SBF->laminarPosition(voxelCenter);

				double *px = static_cast< double * >(density->GetScalarPointer(x,y,z));
				*px = (*px) * SGIBoutonFactor[laminarPos];
			}
		}
	}
	return density;
}

ImageDataPointerType computeSpineDensity(const char * inputfilename, int celltype)
{
	densitySG * dSgApical = new densitySG(inputfilename,ApicalDendrite);
	AmiraSpatialGraph * sg = dSgApical->getSpatialGraph();
	bool apicalFound = sg->isLabelInSpatialGraph(ApicalDendrite);
	ImageDataPointerType lenApical;

	// Check basal dendrite based on label "Dendrite"
	densitySG * dSgBasal = new densitySG(inputfilename,Dendrite);
	bool basalFound = sg->isLabelInSpatialGraph(Dendrite);
	ImageDataPointerType lenBasal;

	int foundNeurites = 0; // Apical and Basal exit

	if (apicalFound)
	{
		lenApical = dSgApical->computeLength();
		std::cout << ">> ApicalDendrite found in SpatialGraph!" << std::endl;
	}
	else
	{
		std::cout << ">> No ApicalDendrite found in SpatialGraph (" << inputfilename << ")" << std::endl;
		foundNeurites = BasalDendrite;// Only Basal exists
	}
	if (basalFound)
	{
		lenBasal = dSgBasal->computeLength();
	}
	else
	{
		dSgBasal = new densitySG(inputfilename,BasalDendrite);
		basalFound = sg->isLabelInSpatialGraph(BasalDendrite);
		if (basalFound)
		{
			lenBasal = dSgBasal->computeLength();
		}
	}
	if (!basalFound)
	{
		std::cout << ">> No BasalDendrite found in SpatialGraph (" << inputfilename << ")" << std::endl;
		foundNeurites = ApicalDendrite;// Only Apical exists
	}
	else
	{
		std::cout << ">> BasalDendrite found in SpatialGraph!" << std::endl;
	}
	if (!apicalFound & !basalFound)
	{
		std::cout << "ERROR! No Dendrites found in SpatialGrpah (" << inputfilename << ")" << std::endl;
		return lenBasal;
	}

	// Get Spine Densities
	std::map< unsigned int, std::map<int, double> > SpineDensityMap = helper::getCellType2SpineDensityMap();
	std::map<int, double> spineDensityCellType = SpineDensityMap[celltype];

	double apicalSpineDensity = spineDensityCellType[ApicalDendrite];
	double basalSpineDensity = spineDensityCellType[BasalDendrite];

	std::cout << ">> SpineDensity used: [Apical Basal] = [" << apicalSpineDensity << " " << basalSpineDensity << "]" << std::endl;

	ImageDataPointerType spineDensity = mergeSpineDensities(lenBasal,basalSpineDensity,
											lenApical,apicalSpineDensity,foundNeurites);
	return spineDensity;
}

ImageDataPointerType mergeSpineDensities(ImageDataPointerType lenBasal, double basalSpineDensity,
										ImageDataPointerType lenApical, double apicalSpineDensity,
										int foundNeurites)
{
	if (foundNeurites==BasalDendrite)
	{// Only Basal Dendrite exists
		ImageDataPointerType spinesBasal = weigthDensity(lenBasal,basalSpineDensity);
		return spinesBasal;
	}
	else if (foundNeurites==ApicalDendrite)
	{ 	// Only Apical Dendrite exists
		ImageDataPointerType spinesApical = weigthDensity(lenApical,apicalSpineDensity);
		return spinesApical;
	}

	// Default, apical and basal dendrite exit
	ImageDataPointerType spinesBasal = weigthDensity(lenBasal,basalSpineDensity);
	ImageDataPointerType spinesApical = weigthDensity(lenApical,apicalSpineDensity);
	ImageDataPointerType spineDensity = mergeDensities(spinesApical,spinesBasal);
	return spineDensity;
}

// Compute true bounding box [not centers!]
void getVoxelBoundingBox(double bounds[6],double origin[3],int extent[6],double VoxelSize[3])
{
	bounds[0] = origin[0] - VoxelSize[0]/2;
	bounds[1] = origin[0] + VoxelSize[0]/2 + VoxelSize[0]*extent[1];
	bounds[2] = origin[1] - VoxelSize[1]/2;
	bounds[3] = origin[1] + VoxelSize[1]/2 + VoxelSize[1]*extent[3];
	bounds[4] = origin[2] - VoxelSize[2]/2;
	bounds[5] = origin[2] + VoxelSize[2]/2 + VoxelSize[2]*extent[5];
}

ImageDataPointerType mergeDensities(ImageDataPointerType density1, ImageDataPointerType density2)
{
	double VoxelSize[3];
	density1->GetSpacing(VoxelSize);

	double origin1[3];
	density1->GetOrigin(origin1);
	double origin2[3];
	density2->GetOrigin(origin2);

	int extent1[6];
	density1->GetExtent(extent1);
	int extent2[6];
	density2->GetExtent(extent2);

	double bound1[6];
	double bound2[6];
	getVoxelBoundingBox(bound1,origin1,extent1,VoxelSize);
	getVoxelBoundingBox(bound2,origin2,extent2,VoxelSize);

	// Get Maximum Extent
	double bounds[6];
	for (int i = 0; i<3; i++)
	{
		bounds[i*2] = std::min(bound1[i*2],bound2[i*2]);
		bounds[i*2+1] = std::max(bound1[i*2+1],bound2[i*2+1]);
	}

	ImageDataPointerType mergedDensity = helper::createImageVolumeNeuroNet(bounds,VoxelSize[0]);
	double origin[3];
	mergedDensity->GetOrigin(origin);

	for(int x = mergedDensity->GetExtent()[0]; x <= mergedDensity->GetExtent()[1]; ++x)
	{
	   for(int y = mergedDensity->GetExtent()[2]; y <= mergedDensity->GetExtent()[3]; ++y)
	   {
			for(int z = mergedDensity->GetExtent()[4]; z <= mergedDensity->GetExtent()[5]; ++z)
			{
				double voxelCenter[3];
				voxelCenter[0] = origin[0] + x*VoxelSize[0];
				voxelCenter[1] = origin[1] + y*VoxelSize[1];
				voxelCenter[2] = origin[2] + z*VoxelSize[2];

				double value = 0;

				int pos1[3];
				if (isPtInDensity(voxelCenter,density1,pos1))
				{
					double * px1 = static_cast< double * >(density1->GetScalarPointer(pos1[0],pos1[1],pos1[2]));
					value += *px1;
				}
				int pos2[3];
				if (isPtInDensity(voxelCenter,density2,pos2))
				{
					double * px2 = static_cast< double * >(density2->GetScalarPointer(pos2[0],pos2[1],pos2[2]));
					value += *px2;
				}
				double * px = static_cast< double * >(mergedDensity->GetScalarPointer(x, y, z));
				*px = value;
			}
	   }
	}

	mergedDensity->Update();
	return mergedDensity;
}

bool isPtInDensity(double voxelCenter[3], ImageDataPointerType density, int pos[3])
{
	double origin[3];
	density->GetOrigin(origin);

	double VoxelSize[3];
	density->GetSpacing(VoxelSize);

	for(int x = density->GetExtent()[0]; x <= density->GetExtent()[1]; ++x)
	{
		double xCenter = origin[0] + x*VoxelSize[0];

		if (xCenter!=voxelCenter[0])
			continue;

	   for(int y = density->GetExtent()[2]; y <= density->GetExtent()[3]; ++y)
	   {
			double yCenter = origin[1] + y*VoxelSize[1];

			if (yCenter!=voxelCenter[1])
				continue;

			for(int z = density->GetExtent()[4]; z <= density->GetExtent()[5]; ++z)
			{
				double zCenter = origin[2] + z*VoxelSize[2];

				if (zCenter!=voxelCenter[2])
					continue;

				pos[0] = x;
				pos[1] = y;
				pos[2] = z;
				return true;
			}
	   }
	}
	return false;
}

ImageDataPointerType weigthDensity(ImageDataPointerType density, double scaleFactor)
{
	for(int x = density->GetExtent()[0]; x <= density->GetExtent()[1]; ++x)
	{
		for(int y = density->GetExtent()[2]; y <= density->GetExtent()[3]; ++y)
		{
			for(int z = density->GetExtent()[4]; z <= density->GetExtent()[5]; ++z)
			{
				double *px = static_cast< double * >(density->GetScalarPointer(x,y,z));
				*px = (*px) * scaleFactor;
			}
		}
	}
	return density;
}

ImageDataPointerType computeInnervationDensity(ImageDataPointerType boutonDensity, ImageDataPointerType spineDensity)
{
	ImageDataPointerType PSTEXC = helper::loadVolume(outputfilePSTTotalEXC);

	double origin[3];
	spineDensity->GetOrigin(origin);
	double VoxelSize[3];
	spineDensity->GetSpacing(VoxelSize);

	for(int x = spineDensity->GetExtent()[0]; x <= spineDensity->GetExtent()[1]; ++x)
	{
	   for(int y = spineDensity->GetExtent()[2]; y <= spineDensity->GetExtent()[3]; ++y)
	   {
			for(int z = spineDensity->GetExtent()[4]; z <= spineDensity->GetExtent()[5]; ++z)
			{
				double voxelCenter[3];
				voxelCenter[0] = origin[0] + x*VoxelSize[0];
				voxelCenter[1] = origin[1] + y*VoxelSize[1];
				voxelCenter[2] = origin[2] + z*VoxelSize[2];

				double I = 0;
				double PSTNorm = 0;
				double boutons = 0;
				bool undefinedPSTValue = false;

				int posPST[3];
				if (isPtInDensity(voxelCenter,PSTEXC,posPST))
				{
					double * px1 = static_cast< double * >(PSTEXC->GetScalarPointer(posPST[0],posPST[1],posPST[2]));
					PSTNorm = *px1;
				}
				else
				{
					undefinedPSTValue = true;
				}

				int posBoutons[3];
				if (isPtInDensity(voxelCenter,boutonDensity,posBoutons))
				{
					double * px2 = static_cast< double * >(boutonDensity->GetScalarPointer(posBoutons[0],posBoutons[1],posBoutons[2]));
					boutons = *px2;
				}

				double * px = static_cast< double * >(spineDensity->GetScalarPointer(x, y, z));
				double spines = *px;

				if (PSTNorm>0)
				{
					*px = (spines * boutons)/PSTNorm;
				}
				else if (PSTNorm<=0 && boutons>0 && spines>0)
				{
					*px = (spines * boutons)/spines;

					if (undefinedPSTValue)
					{
						std::cout << "WARNING! Normalized PST Density is not defined (boutons>0 and spine>0). Normalized PST density set to zero at [";
					}
					else
					{
						std::cout << "WARNING! Only normalized PST Density is zero (boutons>0 and spine>0) at [";
					}
					std::cout << voxelCenter[0] << " " << voxelCenter[1] << " " << voxelCenter[2] << "]" << std::endl;
				}
				else
				{
					*px = 0;
				}

				if (PSTNorm<spines)
				{
					std::cout << "WARNING! Normalized PST Density is smaller than the spines of the postsynaptic cell" << std::endl;
					std::cout << "  Normalized PST Density = " << PSTNorm << "; Spines = " << spines <<  "at [";
					std::cout << voxelCenter[0] << " " << voxelCenter[1] << " " << voxelCenter[2] << "]" << std::endl;
				}
			}
	   }
	}
	spineDensity->Update();
	return spineDensity;
}
