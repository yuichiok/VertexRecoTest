#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include "MathOperator.hh"
#include "VertexTag.hh"
#include "Particle.hh"
#ifndef _VertexRecoOperator_hh
#define _VertexRecoOperator_hh
namespace QQbarAnalysis 
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
			VertexRecoOperator (float angle, EVENT::Vertex * primary);
			virtual ~VertexRecoOperator () {};
		//
		//	Methods
		//
			std::vector< VertexTag * > * Compare(EVENT::LCCollection * reconstructed, EVENT::LCCollection * mc);
			std::vector< VertexTag * > * CompareDirection(EVENT::LCCollection * reconstructed, EVENT::LCCollection * mc);
			std::vector< Particle > * CompareTracks(EVENT::LCCollection * reconstructed, EVENT::LCCollection * mc, std::vector< EVENT::ReconstructedParticle * > * recomissed = NULL);
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
			EVENT::Vertex * myPrimary;
			double ip[3];
			
			std::vector<int> myTakenVertices;
			std::vector<float> myDistances;
			std::vector< EVENT::Vertex * > myUnknownVertexes;

		//
		//	Private methods
		//
			std::vector<float> getAngles(EVENT::Vertex * vertex);
			std::vector< VertexTag * > * refineTagSelection(std::vector< VertexTag * > * old);
			std::vector< VertexTag * > * refineChain(std::vector< VertexTag * > * refined, std::vector< VertexTag * > & bvtx);
			std::vector< EVENT::ReconstructedParticle * > * getMissedTracks(const std::vector< EVENT::ReconstructedParticle * > & recotracks, const std::vector< EVENT::ReconstructedParticle * > & mctracks);
	};
} /* QQbarAnalysis */
#endif
