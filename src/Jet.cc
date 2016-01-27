
#include "Jet.hh"
using EVENT::ReconstructedParticle;
using std::vector;
namespace TTbarAnalysis
{
	Jet::Jet (float btag,float ctag, int number,const double * momentum)
	{
		myBTag = btag;
		myCTag = ctag;
		myNumber = number;
		myMomentum = momentum;
		myMCPDG = 0;
		myParticles = NULL;
		myRecoVertices = NULL;
	}
	Jet::Jet ()
	{
		myRecoVertices = NULL;
		myMCPDG = 0;
		myParticles = NULL;
		
	}
	void Jet::SetRecoVertices(std::vector<  EVENT::Vertex * > * vertices)
	{
		myRecoVertices = vertices;
	}
	std::vector<  EVENT::Vertex * > * Jet::GetRecoVertices()
	{
		return myRecoVertices;
	}
	void Jet::AddVertexTag(VertexTag * tag)
	{
		myVertexTags.push_back(tag);
	}
	const vector< VertexTag * > & Jet::GetVertexTags() const
	{
		return myVertexTags;
	}
	void Jet::SetBTag(float value)
	{
		myBTag = value;
	}
	void Jet::SetCTag(float value)
	{
		myCTag = value;
	}
	float Jet::GetBTag()
	{
		return myBTag;
	}
	float Jet::GetCTag()
	{
		return myCTag;
	}
	int Jet::GetNumberOfVertices()
	{
		if (myRecoVertices) 
		{
			return myRecoVertices->size();
		}
		return myNumber;
	}
	int Jet::GetNumberOfVertexParticles()
	{
		int sum = -1;
		if (myRecoVertices) 
		{
			sum = 0;
			for (unsigned int i = 0; i < myRecoVertices->size(); i++) 
			{
				sum += myRecoVertices->at(i)->getAssociatedParticle()->getParticles().size();
			}
		}
		return sum;
	}
	float Jet::GetHadronCharge()
	{
		float charge = -5.0;
		if (myRecoVertices) 
		{
			charge = 0.0;
			for (unsigned int i = 0; i < myRecoVertices->size(); i++) 
			{
				charge += myRecoVertices->at(i)->getAssociatedParticle()->getCharge();
			}
		}
		return charge;
	}
	float Jet::GetHadronMomentum()
	{
		float momentum = -1.0;
		if (myRecoVertices) 
		{
			momentum = 0.0;
			for (unsigned int i = 0; i < myRecoVertices->size(); i++)
			{
				momentum += MathOperator::getModule(myRecoVertices->at(i)->getAssociatedParticle()->getMomentum()); // CRUNCH!!!
			}
		}
		return momentum;
	}
	float Jet::GetHadronMass()
	{
		float mass = -1.0;
		if (myRecoVertices) 
		{
			 mass = 0.0;
			 for (int i = 0; i < myRecoVertices->size(); i++) 
			 {
			 	mass += myRecoVertices->at(i)->getAssociatedParticle()->getMass();
			 }
		}
		return mass;
	}
	float Jet::GetHadronDistance()
	{
		float distance = 1000.0;
		if (myRecoVertices && myRecoVertices->size() > 0) 
		{
			for (int i = 0; i < myRecoVertices->size(); i++) 
			{
				float current = MathOperator::getModule(myRecoVertices->at(i)->getPosition());
				if (current < distance) 
				{
					distance = current;
				}
			}
		}
		else 
		{
			distance = -1.0;
		}
		return distance;
	}
	const double * Jet::GetMomentum()
	{
		return myMomentum;
	}
	const vector< ReconstructedParticle * > * Jet::GetParticles() const
	{
		return myParticles;
	}
	void Jet::SetParticles(const vector< ReconstructedParticle * > & particles)
	{
		myParticles = new vector< ReconstructedParticle * >( particles);
	}
	
	int Jet::GetMCPDG()
	{
		return myMCPDG;
	}
	void Jet::SetMCPDG(int pdg)
	{
		myMCPDG = pdg;
	}

}
