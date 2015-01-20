#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include "MathOperator.hh"
#include "VertexTag.hh"
#ifndef _VertexRecoOperator_hh
#define _VertexRecoOperator_hh
namespace TTbarAnalysis 
{
	class VertexRecoOperator 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			VertexRecoOperator ();
			VertexRecoOperator (float angle);
			virtual ~VertexRecoOperator () {};
		//
		//	Methods
		//
			std::vector< VertexTag * > * Compare(EVENT::LCCollection * reconstructed, EVENT::LCCollection * mc);
			std::vector< VertexTag * > * CompareDirection(EVENT::LCCollection * reconstructed, EVENT::LCCollection * mc);
			VertexTag * Tag(EVENT::Vertex * vertex, EVENT::LCCollection * mc);
			std::vector< EVENT::Vertex * > & GetUnknownVertexes();
			std::vector< float >  GetDistances();
			bool CheckTaken(EVENT::Vertex * mcVertex);
			int GetID(EVENT::Vertex * mcVertex);
			void SetMinimalDistance(VertexTag * tag, EVENT::LCCollection * mc);
			
		private:
		//
		//	Data
		//
			float myPrecisionCut;
			float myAngleCut;
			std::vector<int> myTakenVertices;
			std::vector<float> myDistances;
			std::vector< EVENT::Vertex * > myUnknownVertexes;

		//
		//	Private methods
		//
			std::vector<float> getAngles(EVENT::Vertex * vertex);
	};
} /* TTbarAnalysis */
#endif
