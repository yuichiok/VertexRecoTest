#include "VertexTag.hh"
#include "MathOperator.hh"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Vertex.h>

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
			Jet ();
			virtual ~Jet () 
			{
				for (unsigned int i = myVertexTags.size(); i > -1; i--) 
				{
					VertexTag * tag = myVertexTags[i];
					delete tag;
				}
			};
		//
		//	Methods
		//
			float GetBTag();
			float GetCTag();
			void SetBTag(float value);
			void SetCTag(float value);
			void AddVertexTag(VertexTag * tag);
			const std::vector< VertexTag * > & GetVertexTags() const;
			void SetRecoVertices(std::vector<  EVENT::Vertex * > * vertices);
			std::vector<  EVENT::Vertex * > * GetRecoVertices();
			int GetNumberOfVertices();
			int GetNumberOfVertexParticles();
			float GetHadronCharge();
			float GetHadronMomentum();
			float GetHadronMass();
			float GetHadronCostheta();
			float GetHadronDistance();
			float GetMinVtxProbability();
			float GetMaxVtxChi2();
			int __GetGenNumberOfVertexParticles();
			int __GetGenCharge();
			const double * GetMomentum();
			int GetMCPDG();
			void SetMCPDG(int pdg);
			const std::vector< EVENT::ReconstructedParticle * > * GetParticles() const;
			void SetParticles(const std::vector< EVENT::ReconstructedParticle * > & particles);
			float GetTagAngle();
			void SetTagAngle(float a);
		private:
		//
		//	Data
		//
			 float myBTag;
			 float myCTag;
			 float myTagAngle;
			 int myNumber;
			 int myMCPDG;
			 std::vector< VertexTag * > myVertexTags;
			 const double * myMomentum;
			 std::vector<  EVENT::ReconstructedParticle * > * myParticles;
			 std::vector<  EVENT::Vertex * > * myRecoVertices;
		//
		//	Private methods
		//
	};
}
#endif
