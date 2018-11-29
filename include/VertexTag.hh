#include <vector>
#include <string>
#include <EVENT/Vertex.h>
#include <EVENT/ReconstructedParticle.h>
#ifndef _VertexTag_hh
#define _VertexTag_hh
namespace QQbarAnalysis
{
		enum TagStatus
		{
			UNKNOWN_TAG,
			EMPTY_TAG,
			MERGED_TAG,
			PRECISE_TAG
		};
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
				VertexTag();
				VertexTag(EVENT::Vertex * mcVertex);
				virtual ~VertexTag () {};
			//
			//	Methods
			//
				EVENT::Vertex * GetVertex();
				EVENT::Vertex * __GetMCVertex();
				void SetStatus(TagStatus status);
				TagStatus GetStatus();
				void SetRecoVertex(EVENT::Vertex * vertex);
				int GetGeneration();
				int GetInitialPDG();
				float GetMinimalDistance();
				void SetMinimalDistance(float d);
				int __GetMCTrackNumber();
				int GetTrackNumber();
				float GetTruthAngle();
				void SetTruthAngle(float theta);
			private:
			//
			//	Data
			//
				EVENT::Vertex * myRecoVertex;
				EVENT::Vertex * myMCVertex;
				TagStatus myStatus;
				float myMinimalDistanceMC;
				int myInitialPDG;
				float myTruthAngle;
			//
			//	Private methods
			//
		};
}
#endif
