#include <vector>
#include <string>
#include <EVENT/Vertex.h>
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
				int GetGeneration();
				int GetInitialPDG();
				float GetMinimalDistance();
				void SetMinimalDistance(float d);
			private:
			//
			//	Data
			//
				EVENT::Vertex * myRecoVertex;
				EVENT::Vertex * myMCVertex;
				float myMinimalDistanceMC;
			//
			//	Private methods
			//
		};
}
#endif
