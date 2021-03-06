#include <string>
#include <vector>
#include <TMVA/Reader.h>
namespace QQbarAnalysis 
{
	struct MVAInputParameters 
	{
		float btag;
		float   VtxCharge;
		float   nVtxParticles; 
		float VtxMomentum;
		float VtxDistance;
		float   VtxMinHits;
		float VtxMinOffset;
		float VtxChargeBalance;
		
		float JetMaxOffset;
		float JetMaxOffsetMomentum;
		float JetMaxOffsetAngle;
		float   JetMaxOffsetCharge;

		/* data */
	} /* optional variable list */;
	class MVAReader 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			MVAReader (std::string weightpath = "/home/ilc/yokugawa/Bilokin/VertexRecoTest/tmva/");
			virtual ~MVAReader () {};
		//
		//	Methods
		//
			std::vector<float> GetAllTags(MVAInputParameters & parameters);
			float GetTrustTag(MVAInputParameters & parameters);

			float GetTag(MVAInputParameters & parameters, std::string path);
		private:
		//
		//	Data
		//
			std::string myWeightPath;
		//
		//	Private methods
		//
	};
} /* TTBarAnalysis */
