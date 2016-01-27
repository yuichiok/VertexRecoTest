#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#ifndef _Particle_hh
#define _Particle_hh
namespace TTbarAnalysis
{
	class Particle 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			Particle ();
			//Particle (const double * momentum, EVENT::Vertex * vertex, double offset);
			virtual ~Particle () 
			{
				//delete myMomentum;
			};
		//
		//	Methods
		//
			double * GetMomentum();
			void SetMomentum(const double * momentum);
			void SetMass(float momentum);
			EVENT::Vertex * GetVertex();
			void SetVertex(EVENT::Vertex * vertex);
			float GetOffset();
			void SetOffset(float offset);
			float GetTheta();
			float GetMass();
			void SetTheta(float offset);
			int GetGeneration();
			void SetGeneration(int value);
			float GetChi2();
			void SetChi2(float value);
			int GetIsReco();
			void SetIsReco(int value);
			const int * GetHits();
			void SetHits(int * value);
			float GetTruthAngle();
			void SetTruthAngle(float value);
			int GetInteracted();
			void SetInteracted(int value);
			float GetJetBtag();
			void SetJetBtag(float value);
			float GetVertexAngle();
			void SetVertexAngle(float value);
			EVENT::MCParticle * GetMCParticle();
			void SetMCParticle(EVENT::MCParticle * value);

			void Assign(EVENT::ReconstructedParticle * initial);
		private:
		//
		//	Data
		//
			double * myMomentum;
			EVENT::MCParticle * myMCParticle;
			EVENT::Vertex * myVertex;
			float myOffset;
			float myTheta;
			float myMass;
			float myTruthAngle;
			float myVertexAngle;
			float myBtag;
			int myGeneration;
			int myIsReco;
			int myHits[4];
			int myHasDecayed;
			float myChi2;
		//
		//	Private methods
		//
	};
}
#endif
