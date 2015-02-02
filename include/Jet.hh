
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <EVENT/ReconstructedParticle.h>

#ifndef _Jet_hh
#define _Jet_hh
namespace TTbarAnalysis 
{
	class Jet 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			Jet (float btag, float ctag, int number, const double * momentum);
			virtual ~Jet () {};
		//
		//	Methods
		//
			float GetBTag();
			float GetCTag();
			int GetNumberOfVertices();
			const double * GetMomentum();
			int GetMCPDG();
			void SetMCPDG(int pdg);
			const std::vector< EVENT::ReconstructedParticle * > * GetParticles() const;
			void SetParticles(const std::vector< EVENT::ReconstructedParticle * > & particles);
			
		private:
		//
		//	Data
		//
			 float myBTag;
			 float myCTag;
			 int myNumber;
			 int myMCPDG;
			 const double * myMomentum;
			 std::vector<  EVENT::ReconstructedParticle * > * myParticles;
		//
		//	Private methods
		//
	};
}
#endif
