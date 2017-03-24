#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <UTIL/LCRelationNavigator.h>
#include "MathOperator.hh"
#ifndef _ParticleOperator_hh
#define _ParticleOperator_hh
namespace TTbarAnalysis 
{
	class ParticleOperator 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			ParticleOperator ();
			virtual ~ParticleOperator () {};
		//
		//	Methods
		//
			static bool CompareParticles(const EVENT::ReconstructedParticle * particle1, const EVENT::ReconstructedParticle * particle2);
			static bool CompareParticles(EVENT::MCParticle * particle1, EVENT::ReconstructedParticle * particle2);
			static bool CheckForVertex(EVENT::ReconstructedParticle * particle1);
			static bool IsDublicate(EVENT::MCParticle * particle, std::vector< EVENT::MCParticle * > & list);
			static bool IsDublicate(const EVENT::ReconstructedParticle * particle, const std::vector< EVENT::ReconstructedParticle * > & list);
			static std::vector< EVENT::MCParticle * > GetMCParticlesRel(const std::vector< EVENT::ReconstructedParticle * > & secondaries, EVENT::LCCollection * rel, EVENT::LCCollection * trackrel = NULL);
			static EVENT::MCParticle * GetMCParticleTrackRel(const EVENT::ReconstructedParticle * secondaries, EVENT::LCCollection * trackrel);
			static EVENT::ReconstructedParticle * GetRecoParticle(EVENT::MCParticle *  particle, EVENT::LCCollection * rel);
			static EVENT::Track * GetTrackRel(EVENT::MCParticle *  particle, EVENT::LCCollection * rel);
			static EVENT::ReconstructedParticle * ReconstructParticle(EVENT::Track *  track);
			static float GetError(const EVENT::ReconstructedParticle * particle);
		private:
		//
		//	Data
		//
			
		//
		//	Private methods
		//
	};
}
#endif
