#ifndef _JetOperator_hh
#define _JetOperator_hh

#include "marlin/Processor.h"
#include "lcio.h"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include "Jet.hh"
#include <EVENT/LCCollection.h>
#include <EVENT/LCObject.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/PIDHandler.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/Vertex.h>
#include "MathOperator.hh"

namespace QQbarAnalysis 
{
	class JetOperator 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			JetOperator (float angle, std::string algorithm);
			JetOperator ();
			virtual ~JetOperator () {};
		//
		//	Methods
		//
			std::vector<float> * GetBtags(EVENT::LCCollection * jets);
			std::vector< Jet * > * GetJets(EVENT::LCCollection * jets, EVENT::LCCollection * rel);
			void CompareDirection(std::vector< Jet * > * jets, EVENT::LCCollection * mc);
		private:
		//
		//	Data
		//
			float myAngle;
			std::string myAlgorithmName;
		//
		//	Private methods
		//
			std::vector< EVENT::Vertex * > * convert(const std::vector< EVENT::LCObject * > & objs);
	};
}
#endif
