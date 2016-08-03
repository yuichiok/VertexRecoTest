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
#include "RecoJet.hh"
#ifndef _ComparisonOperator_hh
#define _ComparisonOperator_hh
namespace TTbarAnalysis 
{
	class ComparisonOperator 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			ComparisonOperator (EVENT::LCCollection * pfo , EVENT::LCCollection * egprongs, EVENT::LCCollection * trackrel = NULL);
			virtual ~ComparisonOperator () {};
			void TagJets(std::vector<RecoJet*> * jets, EVENT::LCCollection * mccol, std::vector<float> & angleTags);
			void tagJets(std::vector<RecoJet*> * jets, EVENT::Vertex * vertex, std::vector<float> & angleTags);
		//
		//	Methods
		//
		private:
		//
		//	Data
		//
			std::string myAlgorithmName;
			EVENT::LCCollection * myPFO;
			EVENT::LCCollection * myTrackRel;
			EVENT::LCCollection * myEGProngs;
			float myAngleCut;
			double ip[3];
			
		//
		//	Private methods
		//
	};
}
#endif
