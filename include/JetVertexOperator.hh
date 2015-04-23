#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Vertex.h>
#include <UTIL/PIDHandler.h>
#include <UTIL/LCRelationNavigator.h>
#include "MathOperator.hh"
#include "ParticleOperator.hh"
#include "VertexTag.hh"
#include "Particle.hh"
#include "Jet.hh"
#ifndef _JetVertexOperator_hh
#define _JetVertexOperator_hh
namespace TTbarAnalysis 
{
	class JetVertexOperator 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			JetVertexOperator ( );
			virtual ~JetVertexOperator () {};
		//
		//	Methods
		//
			std::vector< Jet * > * TagJets(EVENT::LCCollection * jetcol, EVENT::LCCollection * mccol, EVENT::LCCollection * rel);
			void TagVertices(std::vector< Jet * > * jets, EVENT::LCCollection * mccol, std::vector< EVENT::Vertex * > & lost);
			std::vector< EVENT::ReconstructedParticle * > * GetMissedTracks(std::vector< Jet * > * jets, std::vector< Particle > * converted = NULL);
			std::vector< EVENT::MCParticle * > * GetMissedTracks(EVENT::LCCollection * prongs, EVENT::LCCollection * rel, EVENT::LCCollection * out = NULL);
			std::vector< EVENT::ReconstructedParticle * > * CompareTracks(const std::vector< EVENT::ReconstructedParticle * > & recotracks, const std::vector< EVENT::ReconstructedParticle * > & mctracks);
			void CompareTracks(EVENT::Vertex * mcvertex, const std::vector< EVENT::ReconstructedParticle * > & recotracks, std::vector< EVENT::ReconstructedParticle * > * missedTotal, std::vector< Particle > * convertedTotal = NULL);
			std::vector< EVENT::ReconstructedParticle * > * GetRecoParticles(EVENT::LCCollection * prongs, EVENT::LCCollection * rel);
			
		private:
		//
		//	Data
		//
			std::string myAlgorithmName;
			float myAngleCut;
			double ip[3];
		//
		//	Private methods
		//
			//Jet * composeJet(ReconstructedParticle * particle, std::vector< Vertex * > recovertices, const float tag);
			void produceTags(Jet * jet, std::vector< EVENT::Vertex * > & vertices, std::vector< VertexTag * > * tags = NULL);
			std::vector< Vertex * > * convert(const std::vector< LCObject * > & objs);
			std::vector< VertexTag * > * tagOneVertex(std::vector< EVENT::Vertex * > & mcvertices, Vertex * recovertex);
			bool compareToVertex(EVENT::Vertex * mcvertex,  std::vector< EVENT::Vertex * > * vertices);
	};
}
#endif
