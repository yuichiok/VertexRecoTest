#include <vector>
#include <string>
#include <EVENT/Vertex.h>
#include <EVENT/ReconstructedParticle.h>
#ifndef _VertexTag_hh
#define _VertexTag_hh
namespace TTbarAnalysis
{
		class VertexTag 
		{
			public:
			//
			//	Constants
			//
		
			//
			//	Constructors
			//
				VertexTag (EVENT::Vertex * recoVertex, EVENT::Vertex * mcVertex);
				virtual ~VertexTag () {};
			//
			//	Methods
			//
				EVENT::Vertex * GetVertex();
				EVENT::Vertex * __GetMCVertex();
				int GetGeneration();
				int GetInitialPDG();
				float GetMinimalDistance();
				void SetMinimalDistance(float d);
				int __GetMCTrackNumber();
				int GetTrackNumber();
			private:
			//
			//	Data
			//
				EVENT::Vertex * myRecoVertex;
				EVENT::Vertex * myMCVertex;
				float myMinimalDistanceMC;
				int myInitialPDG;
			//
			//	Private methods
			//
		};
}
#endif
