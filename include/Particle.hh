#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <EVENT/ReconstructedParticle.h>
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
			EVENT::Vertex * GetVertex();
			void SetVertex(EVENT::Vertex * vertex);
			float GetOffset();
			void SetOffset(float offset);
			float GetTheta();
			void SetTheta(float offset);
			void Assign(EVENT::ReconstructedParticle * initial);
		private:
		//
		//	Data
		//
			double * myMomentum;
			EVENT::Vertex * myVertex;
			float myOffset;
			float myTheta;
		//
		//	Private methods
		//
	};
}
#endif
