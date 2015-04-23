#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
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
			static bool CompareParticles(EVENT::ReconstructedParticle * particle1, EVENT::ReconstructedParticle * particle2);
			static bool CheckForVertex(EVENT::ReconstructedParticle * particle1);
			static bool IsDublicate(EVENT::MCParticle * particle, std::vector< EVENT::MCParticle * > & list);
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
